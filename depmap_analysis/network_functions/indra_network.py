import re
import inspect
import logging
from itertools import product
from collections import defaultdict
from time import time, gmtime, strftime
from numpy import trunc

import requests
import networkx as nx
from networkx import NodeNotFound, NetworkXNoPath, NetworkXError

from indra_db.client.readonly.mesh_ref_counts import get_mesh_ref_counts

from indra.config import CONFIG_DICT
from indra.databases import get_identifiers_url
from indra.assemblers.indranet.net import default_sign_dict
from indra.explanation.pathfinding.util import signed_nodes_to_signed_edge, \
    path_sign_to_signed_nodes
from indra.explanation.model_checker.model_checker import \
    signed_edges_to_signed_nodes
from indra.explanation.pathfinding.pathfinding import shortest_simple_paths, \
    bfs_search, open_dijkstra_search
from depmap_analysis.network_functions import famplex_functions as ff
from depmap_analysis.network_functions import net_functions as nf
from depmap_analysis.network_functions.net_functions import \
    INT_PLUS, INT_MINUS

bfs_kwargs = inspect.signature(bfs_search).parameters.keys()

logger = logging.getLogger('indra network')

GRND_URI = None
try:
    GRND_URI = CONFIG_DICT['GILDA_URL']
except KeyError:
    logger.warning('Indra Grounding service not available. Add '
                   'GILDA_URL to `indra/config.ini`')

MAX_PATHS = 50
TIMEOUT = 60  # Timeout in seconds
MIN_TIMEOUT = 2
MAX_TIMEOUT = 120
MAX_SIGNED_PATH_LEN = 7
EMPTY_RESULT = {'paths_by_node_count': {'forward': {}, 'backward': {},
                                        'path_hashes': []},
                'common_targets': [],
                'shared_regulators': [],
                'common_parents': {},
                'timeout': False,
                'node_not_found': False}
MANDATORY = ['stmt_filter', 'node_filter',
             'path_length', 'weighted', 'bsco', 'fplx_expand',
             'k_shortest', 'curated_db_only', 'two_way']
USER_OVERRIDE = False


def _truncate(n):
    return float(trunc(n * 100) / 100)


class MissingParametersError(Exception):
    """Raise for missing query parameters"""
    pass


class IndraNetwork:
    """Handle searches and graph output of the INDRA DB network"""
    def __init__(self, indra_dir_graph=nx.DiGraph(),
                 indra_multi_dir_graph=nx.MultiDiGraph(),
                 indra_sign_edge_graph=nx.MultiDiGraph(),
                 indra_sign_node_graph=None):
        self.nx_dir_graph_repr = indra_dir_graph
        self.nx_md_graph_repr = indra_multi_dir_graph
        self.sign_node_graph_repr = signed_edges_to_signed_nodes(
            indra_sign_edge_graph) if indra_sign_node_graph\
            is None else indra_sign_node_graph
        self.sign_edge_graph_repr = indra_sign_edge_graph
        self.nodes = indra_dir_graph.nodes \
            if not nx.is_empty(indra_dir_graph) \
            else (indra_sign_edge_graph.nodes
                  if not nx.is_empty(indra_sign_edge_graph)
                  else dict())
        self.signed_nodes = self.sign_node_graph_repr.nodes
        self.dir_edges = indra_dir_graph.edges
        self.mdg_edges = indra_multi_dir_graph.edges
        self.signed_node_edges = self.sign_node_graph_repr.edges
        self.signed_edges = indra_sign_edge_graph.edges
        self.node_by_uri = indra_dir_graph.graph.get('node_by_uri', None)
        self.node_by_ns_id = indra_dir_graph.graph.get('node_by_ns_id', None)
        self.MAX_PATHS = MAX_PATHS
        self.TIMEOUT = TIMEOUT
        self.MANDATORY = MANDATORY
        self.small = False
        self.verbose = 0
        self.query_receive_time = 0.0
        self.query_timed_out = False

    def handle_query(self, **kwargs):
        """Handles path query from client. Returns query result.

        Notes:
            - Parameters for yet to be implemented functionalities are not
              mandatory and have no effect on the path search if provided
            - Parameter combinations that will trigger some default responses:
                * 'path_length' AND 'weighted':
                    Defaults to unrestricted weighted search, i.e.
                    path_length'=False and 'weighted'=True. Placing
                    the path length constrain on top of the weigthed search
                    can take a substantial amount of time

        The query is a json-friendly key-value structure contained in kwargs
        with the following parameters:

        Parameters
        ----------
        source: str
            the source node for the path
        target: str
            the target for the path
        stmt_filter: [str]
            a list of valid indra statement types or FamPlex child-parent
            connections (as 'fplx') *to exclude* in the path
        node_filter: [str]
            a list of node namespaces *to include* in the path
        node_blacklist: [str]
            a list of node names to ignore. If a path contains a node in this
            list, the path will be discarded.
        edge_hash_blacklist: [str/int]
            a list of statement hashes (as strings or ints) to ignore. If an
            edge statement hash is found in this list, it will be discarded
            from the assembled edge list.
        cull_best_node: [int]
            a positive integer. Every x valid paths, cull the node with the
            highest (weighted) degree from the network. This increases the
            variety of paths found and reduces the impact of nodes with
            extremely high connectivity in the network.
        path_length: int|False
            a positive integer stating the number of edges that should be in
            the returned path. If False, return paths with any number of edges.
        sign: None|str [None|'plus'|'minus']
            If 'no_sign' or None, do regular unsigned graph search over the
            directed graph. If 'plus'/'minus', only paths with overall
            up/down regulation will be returned.
        weighted: Bool
            If True, do a weighted path search. Weights in the network are
            assigned as -log(belief score).
        bsco: 0 <= float <= 1.0
            Belief Score Cut-Off, a positive decimal number < 1.0 indicating
            at what belief score an edge statement should be ignored
        direct_only: Bool **currently not implemented**
            Placeholder for future implementation of allowing to filter edges
            on the annotation 'direct' in indra statements
        curated_db_only: Bool
            Filter results to only allow edges that are sourced from curated
            databases
        fplx_expand: Bool
            If True, when no path is found in the initial search, look for
            paths between the parents of the source and target
        k_shortest: Bool|int
            An integer stating the maximum number of directed paths to return
            in the result. The maximum allowed value is 50. If False,
            the maximum number of paths returned will be set to the maximum
            allowed value.
        user_timeout : float
            A decimal specifying the number of seconds to use for timeout. If
            not provided, the default of 30 seconds is used.
        two_way: Bool
            If True, search path both ways, i.e. search A->B and B->A
        mesh_ids : list
            List of MeSH IDs related to the hashes used for filtering edges
        strict_mesh_id_filtering : Bool
            If true, the provided MeSH IDs specify the edges used exclusively
            in path finding; otherwise, they are used in weight calculation
        const_c : int
            Constant used in MeSH IDs-based weight calculation
        const_tk : int
            Constant used in MeSH IDs-based weight calculation

        Returns
        -------
        result : dict('paths_by_node_count'={ksp_forward, ksp_backward},
                      'common_targets'=ct,
                      'common_parents'=cp)
            A dict containing the results from each path search and a flag
            for timeout:
                ksp_forward : dict(int)
                    Dict keyed by node count with the results of directed path
                    search from source to target
                ksp_backward : dict(int)
                    Dict keyed by node count with the results of directed path
                    search from target to source
                ct : list[dict('target')]
                    List of dicts keyed by common target name, sorted on
                    highest lowest belief score
                sr : list[dict('regulator')]
                    List of dicts keyed by shared regulator name, sorted on
                    highest lowest belief score
                cp : dict
                    Dict with result of common parents search together with
                    the ns:id pairs used to resolve the query
                timeout : Bool
                    True if the query timed out
                node_not_found: Bool | str
                    An error message if a queried node is not present in the
                    network, otherwise False
        """
        self.query_receive_time = time()
        self.query_timed_out = False
        logger.info('Query received at %s' %
                    strftime('%Y-%m-%d %H:%M:%S (UTC)',
                             gmtime(self.query_receive_time)))

        if not self.sanity_check(**kwargs):
            # todo Add detailed test with info of network stats like number
            #  of nodes, edges, last updated, available network types
            return EMPTY_RESULT
        if not (kwargs.get('source', False) or kwargs.get('target', False)):
            raise KeyError('At least one of "source" or "target" has to be '
                           'provided')
        if not all([key in kwargs for key in self.MANDATORY]):
            miss = [key in kwargs for key in self.MANDATORY].index(False)
            raise KeyError('Missing mandatory parameter "%s"' %
                           self.MANDATORY[miss])
        options = translate_query(kwargs)

        # If open ended search, skip common parents, common targets,
        # common regulators (for now).

        k_shortest = kwargs.pop('k_shortest', None)
        self.MAX_PATHS = k_shortest if k_shortest else MAX_PATHS
        user_timeout = kwargs.pop('user_timeout', None)
        if user_timeout:
            if user_timeout < MIN_TIMEOUT:
                logger.warning('Resetting timeout to minimum value (%d)' %
                               MIN_TIMEOUT)
                self.TIMEOUT = MIN_TIMEOUT
            elif user_timeout > MAX_TIMEOUT:
                logger.warning('Resetting timeout to maximum value (%d)' %
                               MAX_TIMEOUT)
                self.TIMEOUT = MAX_TIMEOUT
            else:
                self.TIMEOUT = user_timeout
        else:
            self.TIMEOUT = TIMEOUT
        logger.info('Query translated to: %s' % repr(options))
        logger.info('Looking for no more than %d paths' % self.MAX_PATHS)

        ksp_backward = {}
        boptions = options.copy()
        boptions['source'] = options.get('target')
        boptions['target'] = options.get('source')
        node_not_found = False
        try:
            # Special case: 1 or 2 unweighted, unsigned edges only, non-open
            # search
            if not options['weight'] and \
                    not (bool(options.get('target')) ^
                         bool(options.get('source'))) and \
                    options['path_length'] in [1, 2]:
                ksp_forward = self._unweighted_direct(**options)
                if options['two_way']:
                    ksp_backward = self._unweighted_direct(**boptions)
            else:
                ksp_forward = self.find_shortest_paths(**options)
                if options['two_way']:
                    ksp_backward = self.find_shortest_paths(**boptions)
        except (NodeNotFound, NetworkXError) as e:
            # Make sure we're catching the correct NetworkXError
            if isinstance(e, NetworkXError):
                patt = r'The node (.*) is not in the digraph.'
                m = re.search(patt, str(e))
                if m is None:
                    raise e
            logger.info('No paths found, trying to ground source and target')
            ksp_forward = self.grounding_fallback(**options)
            if options['two_way']:
                ksp_backward = self.grounding_fallback(**boptions)
            node_not_found = str(e)

        if options.get('source') and options.get('target'):
            ct = self.find_common_targets(**options)
            sr = self.find_shared_regulators(**options) if\
                options.get('shared_regulators', False) else []
            cp = self.get_common_parents(**options)
        else:
            ct = EMPTY_RESULT['common_targets']
            cp = EMPTY_RESULT['common_parents']
            sr = EMPTY_RESULT['shared_regulators']

        if not ksp_forward and not ksp_backward and not ct and not sr and\
                not cp.get('common_parents', []):
            ckwargs = options.copy()
            bckwargs = boptions.copy()
            if kwargs['fplx_expand'] and options.get('source') and \
                    options.get('target'):

                logger.info('No directed path found, looking for paths '
                            'connected by common parents of source and/or '
                            'target')
                ksp_forward = self.try_parents(**ckwargs)
                if options['two_way']:
                    ksp_backward = self.try_parents(**bckwargs)
                if self.verbose > 2:
                    logger.info('Parents search result: %s' %
                                repr(ksp_forward))

        if not ksp_forward and not ksp_backward:
            logger.info('No directed path found')
        if not options['weight'] and not options['mesh_ids']:
            if ksp_forward:
                # Sort the results in ksp_forward if non-weighted search
                ksp_forward = self._sort_stmts(ksp_forward)
            if ksp_backward:
                # Sort the results in ksp_forward if non-weighted search
                ksp_backward = self._sort_stmts(ksp_backward)
        fwd_hshs = list_all_hashes(ksp_forward) if ksp_forward else []
        bwd_hshs = list_all_hashes(ksp_backward) if ksp_backward else []
        all_path_hashes = fwd_hshs + bwd_hshs
        return {'paths_by_node_count': {'forward': ksp_forward,
                                        'backward': ksp_backward,
                                        'path_hashes': all_path_hashes},
                'common_targets': ct,
                'shared_regulators': sr,
                'common_parents': cp,
                'timeout': self.query_timed_out,
                'node_not_found': node_not_found}

    @staticmethod
    def sanity_check(**options):
        """Checks for some possible gotchas in query"""
        # Check for test
        if options.get('test', False):
            logger.info('Query handling test passed')
            return False
        # Check if source or target is in blacklist
        if options.get('source') in options.get('node_blacklist') or\
                options.get('target') in options.get('node_blacklist'):
            logger.warning('Source and/or target is blacklisted!')
            # Return True so path is returned anyway
            return True

        # Check non-resolving query
        # sns, sid, snn = nf.gilda_normalization(options['source'],
        #                                    gilda_retry=True)
        # tns, tid, tnn = nf.gilda_normalization(options['target'])
        # if (sns and sns.lower() not in options['node_filter']) or \
        #         (tns and tns.lower() not in options['node_filter']):
        #     if sns.lower() not in options['node_filter']:
        #         logger.warning('%s not among accepted nodes' % sns)
        #     if tns.lower() not in options['node_filter']:
        #         logger.warning('%s not among accepted nodes' % tns)
        #     return False

        return True

    def grounding_fallback(self, **ckwargs):
        """Retry search with alternative names found by grounding service"""
        if ckwargs.get("mesh_ids"):
            return None
        logger.info('Expanding search using grounding service')
        org_source = ckwargs.get('source')
        org_target = ckwargs.get('target')
        # ToDo:
        #  -establish grounding priority when scores are equal between
        #   groundings
        #  -Add ignore list for new ns added for UP:PRO

        # Get groundings
        if org_source:
            src_groundings = requests.post(GRND_URI,
                                           json={'text': org_source}).json()
        else:
            src_groundings = {}
        if org_target:
            trgt_groundings = requests.post(GRND_URI,
                                            json={'text': org_target}).json()
        else:
            trgt_groundings = {}

        # Loop combinations of source and target groundings, break if
        # anything found

        # org target with sources (ckwargs['target'] is unaltered here)
        if src_groundings and not trgt_groundings:
            for src in src_groundings:
                if src['term']['entry_name'] == org_source:
                    continue
                ckwargs['source'] = src['term']['entry_name']
                ksp = self.find_shortest_paths(**ckwargs)
                if ksp:
                    return ksp

        # org source with targets
        if not src_groundings and trgt_groundings:
            ckwargs['source'] = org_source
            for trgt in trgt_groundings:
                if trgt['term']['entry_name'] == org_target:
                    continue
                ckwargs['target'] = trgt['term']['entry_name']
                ksp = self.find_shortest_paths(**ckwargs)
                if ksp:
                    return ksp

        # all source groundings with all target groundings
        if src_groundings and trgt_groundings:
            for src, trgt in product(src_groundings, trgt_groundings):
                if trgt['term']['entry_name'] == org_target and \
                        src['term']['entry_name'] == org_source:
                    continue
                ckwargs['source'] = src['term']['entry_name']
                ckwargs['target'] = trgt['term']['entry_name']
                ksp = self.find_shortest_paths(**ckwargs)
                if ksp:
                    return ksp

        if self.verbose:
            if not src_groundings and not trgt_groundings:
                logger.info('No groundings for source or target')
            else:
                logger.info('No paths found between grounding alternatives')
        return {}

    def try_parents(self, **ckwargs):
        """Retry search with sources' and targets' parents

        Search for paths between combinations of the parents of source and
        target.
        """
        source = ckwargs['source']
        target = ckwargs['target']

        if self.verbose > 1:
            logger.info('Parents search: source=%s, target=%s' %
                        (ckwargs['source'], ckwargs['target']))

        # Get closures for source and target
        source_parents = self._get_parents(source)
        target_parents = self._get_parents(target)
        if self.verbose > 3:
            logger.info('Got source_parents: %s' %
                        repr(source_parents))
            logger.info('Got target_parents: %s' %
                        repr(target_parents))

        # First try current source with all target parents
        if target_parents and not source_parents:
            for tp_uri in target_parents:
                ckwargs['target'] = self.node_by_uri[tp_uri]
                if self.verbose > 4:
                    logger.info('Parents search: source=%s, target=%s' %
                                (ckwargs['source'], ckwargs['target']))
                ksp = self.find_shortest_paths(**ckwargs)
                if ksp:
                    return ksp

        # Then, try current target with all source parents
        if source_parents and not target_parents:
            for sp_uri in source_parents:
                ckwargs['source'] = self.node_by_uri[sp_uri]
                if self.verbose > 4:
                    logger.info('Parents search: source=%s, target=%s' %
                                (ckwargs['source'], ckwargs['target']))
                ksp = self.find_shortest_paths(**ckwargs)
                if ksp:
                    return ksp

        # Lastly try all possible pairs of source and target parents
        if source_parents and target_parents:
            for sp_uri, tp_uri in product(source_parents,
                                          target_parents):
                ckwargs['source'] = self.node_by_uri[sp_uri]
                ckwargs['target'] = self.node_by_uri[tp_uri]
                if self.verbose > 4:
                    logger.info('Parents search: source=%s, target=%s' %
                                (ckwargs['source'], ckwargs['target']))
                ksp = self.find_shortest_paths(**ckwargs)
                if ksp:
                    return ksp

        # If we get this far, no path was found
        return {}

    def find_shortest_path(self, source, target, **options):
        """Returns a list of nodes representing a shortest path"""
        try:
            return self._loop_paths(
                source=source,
                target=target,
                paths_gen=nx.shortest_path(
                    G=self.nx_dir_graph_repr,
                    source=source,
                    target=target,
                    weight=options['weight']
                ),
                **options)
        except NetworkXNoPath as e:
            logger.warning(repr(e))
            return {}

    def _unweighted_direct(self, **options):
        logger.info('Doing unweighted path search for %d-edge paths' %
                    options['path_length'])
        if options['path_length'] == 1:
            return self._one_edge_path(**options)
        elif options['path_length'] == 2:
            return self._two_edge_path(**options)
        return {}

    def _one_edge_path(self, source, target, **options):
        res = {}
        path = [source, target]
        hash_path = []
        if options['sign'] is None:
            if self.dir_edges.get((source, target), None):
                if self.verbose > 1:
                    logger.info('Found direct path from %s to %s' %
                                (source, target))
                hash_path = self._get_hash_path(path=path, source=source,
                                                target=target, **options)
        elif options['sign'] is not None:
            int_sign = nf.SIGNS_TO_INT_SIGN[options['sign']]
            if self.signed_edges.get((source, target, int_sign)):
                if self.verbose > 1:
                    logger.info('Found direct signed path from %s to %s' %
                                (source, target))
                hash_path = self._get_hash_path(path=path, source=source,
                                                target=target,
                                                edge_signs=[int_sign],
                                                graph_type='signed',
                                                **options)
        if hash_path and all(hash_path):
            pd = {'stmts': hash_path,
                  'path': path,
                  'cost': str(self._get_cost(path)),
                  'sort_key': str(self._get_sort_key(path, hash_path))}
            res = {2: [pd]}
        return res

    def _two_edge_path(self, source, target, **options):

        def _paths_genr(s, t, imts, ign_nodes, ign_hashes):
            for i in imts:
                if isinstance(i, tuple) and i[0] in ign_nodes or \
                        i in ign_nodes:
                    continue
                else:
                    yield [s, i, t]

        # Loop the set of all intermediate nodes
        ign_options = {'ign_nodes': options.get('node_blacklist', []),
                       'ign_hashes': options.get('edge_hash_blacklist', [])}
        if options['sign'] is not None:
            sign_source, signed_target =\
                path_sign_to_signed_nodes(source, target, options['sign'])
            sign_interm = set(self.sign_node_graph_repr.succ[sign_source]) & \
                set(self.sign_node_graph_repr.pred[signed_target])
            paths_gen = _paths_genr(sign_source, signed_target, sign_interm,
                                    **ign_options)
        else:
            intermediates = set(self.nx_dir_graph_repr.succ[source]) & \
                set(self.nx_dir_graph_repr.pred[target])
            paths_gen = _paths_genr(source, target, intermediates,
                                    **ign_options)
        res = defaultdict(list)
        added_paths = 0
        graph_type = 'digraph' if options['sign'] is None else 'signed'
        for _path in paths_gen:
            if added_paths >= self.MAX_PATHS:
                logger.info('Found all %d shortest paths, returning results.' %
                            self.MAX_PATHS)
                return res
            if time() - self.query_receive_time > self.TIMEOUT:
                logger.info('Reached timeout (%d s) before finding all %d '
                            'paths. Returning search.' % (self.TIMEOUT,
                                                          MAX_PATHS))
                self.query_timed_out = True
                return res
            if options['sign'] is None:
                path = _path
                edge_signs = None
            else:
                path = [n[0] for n in _path]
                edge_signs = [signed_nodes_to_signed_edge(s, t)[2]
                              for s, t in zip(_path[:-1], _path[1:])]
            hash_path = self._get_hash_path(path=path, source=source,
                                            target=target,
                                            edge_signs=edge_signs,
                                            graph_type=graph_type,
                                            **options)
            if hash_path and all(hash_path):
                if self.verbose > 1:
                    logger.info('Adding stmts and path from %s to path list' %
                                repr(hash_path))
                pd = {'stmts': hash_path,
                      'path': path,
                      'cost': str(self._get_cost(path)),
                      'sort_key': str(self._get_sort_key(path, hash_path))}
                res[3].append(pd)
                added_paths += 1
        return res

    def find_shortest_paths(self, source=None, target=None, **options):
        """Return a list of shortest paths with their support in ascending
        order.

        The paths are ranked by cost minimization (weighted search) or
        number of edges (unweighted search).

        If weighted, use
        indra.explanation.pathfinding.shortest_simple_paths

        If open ended, i.e. only source or target is provided, use
        self.open_bfs
        """
        try:
            blacklist_options =\
                {'ignore_nodes': options.get('node_blacklist', None)}
            if bool(source) ^ bool(target):
                logger.info('Doing open ended %sbreadth first search' %
                            ('signed ' if options.get('sign') is not None
                             else ''))
                if source:
                    # Open downstream search
                    start_node = source
                    reverse = False
                else:
                    # Open upstream search
                    start_node = target
                    reverse = True

                if options['strict_mesh_id_filtering'] or \
                        (not options['mesh_ids'] and not options['weight']):
                    return self.open_bfs(start_node=start_node,
                                         reverse=reverse, **options,
                                         **blacklist_options)
                else:
                    return self.open_dijkstra(start_node=start_node,
                                              reverse=reverse, **options,
                                              **blacklist_options)
            else:
                logger.info('Doing simple %spath search' %
                            ('weighted ' if options['weight'] else ''))
            if options['mesh_ids']:
                hash_mesh_dict = get_mesh_ref_counts(options['mesh_ids'],
                                                     require_all=False)
                related_hashes = hash_mesh_dict.keys()
                ref_counts_from_hashes = _get_ref_counts_func(hash_mesh_dict)
            else:
                related_hashes = None
                ref_counts_from_hashes = None
            strict = options['strict_mesh_id_filtering']
            if options['sign'] is None:
                # Do unsigned path search
                search_graph = self.nx_dir_graph_repr
                paths = shortest_simple_paths(search_graph,
                                              source, target,
                                              options['weight'],
                                              hashes=related_hashes,
                                              ref_counts_function=
                                              ref_counts_from_hashes,
                                              strict_mesh_id_filtering=strict,
                                              const_c=options['const_c'],
                                              const_tk=options['const_tk'],
                                              **blacklist_options)
                subj = source
                obj = target
            else:
                # Generate signed nodes from query's overall sign
                (src, src_sign), (trgt, trgt_sign) = \
                    path_sign_to_signed_nodes(
                        source, target, options['sign']
                    )
                # Get signed nodes for source and target
                subj = (src, nf.SIGNS_TO_INT_SIGN[src_sign])
                obj = (trgt, nf.SIGNS_TO_INT_SIGN[trgt_sign])
                # Generate signed blacklisted nodes
                signed_blacklisted_nodes = []
                for n in options.get('node_blacklist', []):
                    signed_blacklisted_nodes += [(n, INT_PLUS), (n, INT_MINUS)]
                search_graph = self.sign_node_graph_repr
                paths = shortest_simple_paths(
                    search_graph, subj, obj, options['weight'],
                    ignore_nodes=signed_blacklisted_nodes,
                    hashes=related_hashes,
                    ref_counts_function=ref_counts_from_hashes,
                    strict_mesh_id_filtering=options[
                        'strict_mesh_id_filtering'],
                    const_c=options['const_c'],
                    const_tk=options['const_tk'])

            return self._loop_paths(graph=search_graph, source=subj,
                                    target=obj, paths_gen=paths, **options)

        except NetworkXNoPath as err:
            logger.warning(repr(err))
            return {}

    def open_bfs(self, start_node, reverse=False, depth_limit=2,
                 path_limit=None, terminal_ns=None, max_per_node=5, **options):
        """Return paths and their data starting from source

        Parameters
        ----------
        start_node : str
            Node to start search from
        reverse : bool
            If True, let source be the start of an upstream search.
            Default: False
        depth_limit : int
            The maximum allowed depth (number of edges). Default: 2
        path_limit : int
            The maximum number of paths the generator should yield. Note
            that this different from max_results which determines how many
            paths to return from the path generator. Default: None (no limit).
        terminal_ns : list[str]
            Force a path to terminate when any of the namespaces in this
            list are encountered. Default: ['chebi', 'pubchem'].
        max_per_node : int
            The maximum number of times a node can be a parent to leaf
            nodes. Default: 5
        options : **kwargs
            For a full list of options see
            indra.explanation.pathfinding::bfs_search
            Notable options:
                -max_results : int
                    The maximum number of results to return. Default: 50.
                -sign : int
                    If sign is present as a kwarg, it specifies the sign of
                    leaf node in the path, i.e. whether the leaf node is up-
                    or down regulated.
                -mesh_ids : list[str]
                    List of MeSH IDs relevance to which is considered if strict
                    filtering by statement hashes is required
                -strict_mesh_id_filtering : bool
                    If True, only consider edges relevant to provided hashes

        Returns
        -------
        list
            List of dicts with results. Each dict has the format:
                {'path': <tuple of nodes>,
                 'stmts': <data supporting the paths>}
        """
        # Signed search
        if options.get('sign') is not None:
            graph = self.sign_node_graph_repr
            signed_node_blacklist = []
            for node in options.get('node_blacklist', []):
                signed_node_blacklist.extend([(node, INT_MINUS),
                                              (node, INT_PLUS)])
            options['node_blacklist'] = signed_node_blacklist

            # Assign the correct sign to source:
            # If search is downstream, source is the first node and the
            # search must always start with + as node sign. The leaf node
            # sign (i.e. the end of the path) in this case will then be
            # determined by the requested sign.
            # If reversed search, the source is the last node and can have
            # + or - as node sign depending on the requested sign.
            starting_node = get_signed_node(start_node, options['sign'] or
                                            None, reverse)

        # Normal search
        else:
            graph = self.nx_dir_graph_repr
            starting_node = start_node

        # Set default terminal_ns
        if terminal_ns is None:
            terminal_ns = ['chebi', 'pubchem']

        # Limit search scope
        #  -Ensure finite search: e.g. Can't have unlimited depth together
        #   with no path limit and no max results limit
        if depth_limit is None and path_limit is None and options.get(
                'max_results', 0) > 10000:
            raise ValueError('Limitless search detected: depth_limit is '
                             'None, path_limit is None and max_results > '
                             '10000, aborting')

        depth_limit = options['path_length'] - 1 if options.get('path_length')\
            else depth_limit
        if self.verbose > 1:
            logger.info(f'Depth limit set to {depth_limit}')

        # Get the bfs options from options
        bfs_options = {k: v for k, v in options.items() if k in bfs_kwargs}

        if options['mesh_ids']:
            hash_mesh_dict = get_mesh_ref_counts(options['mesh_ids'],
                                                 require_all=False)
            related_hashes = hash_mesh_dict.keys()
            allowed_edges = \
                {graph.graph['edge_by_hash'][h] for h in
                 related_hashes if h in graph.graph['edge_by_hash']}
            def allow_edge(u, v):
                return (u, v) in allowed_edges
        else:
            related_hashes = []
            def allow_edge(u, v): return True

        bfs_gen = bfs_search(g=graph, source_node=starting_node,
                             reverse=reverse, depth_limit=depth_limit,
                             path_limit=path_limit, max_per_node=max_per_node,
                             terminal_ns=terminal_ns, hashes=related_hashes,
                             allow_edge=allow_edge, **bfs_options)
        return self._loop_open_paths(graph, bfs_gen, source_node=start_node,
                                     reverse=reverse, **options)

    def open_dijkstra(self, start_node, reverse=False,
                      terminal_ns=None, ignore_nodes=None, **options):
        """Do Dijkstra search from a given node and yield paths

        Parameters
        ----------
        start_node : node
            Node in the graph to start from.
        reverse : bool
            If True go upstream from source, otherwise go downstream. Default:
            False.
        terminal_ns : list[str]
            Force a path to terminate when any of the namespaces in this list
            are encountered and only yield paths that terminate at these
            namepsaces
        ignore_nodes : list[str]
            Paths containing nodes from this list are excluded
        **options: **kwargs
            Options to pass along to open_dijkstra_search. Notable options:
                depth_limit : int
                    Stop when all paths with this many edges have been found.
                    Default: 2.
                path_limit : int
                    The maximum number of paths to return. Default: no limit.

        Yields
        ------
        path : tuple(node)
            Paths in the bfs search starting from `source`.
        """
        # Signed search
        if options.get('sign') is not None:
            graph = self.sign_node_graph_repr
            signed_node_blacklist = []
            for node in options.get('node_blacklist', []):
                signed_node_blacklist.extend([(node, INT_MINUS),
                                              (node, INT_PLUS)])
            options['node_blacklist'] = signed_node_blacklist
            starting_node = get_signed_node(start_node, options['sign'],
                                            reverse)

        # Normal search
        else:
            graph = self.nx_dir_graph_repr
            starting_node = start_node

        if options['mesh_ids']:
            hash_mesh_dict = get_mesh_ref_counts(options['mesh_ids'],
                                                 require_all=False)
            related_hashes = hash_mesh_dict.keys()
            ref_counts_from_hashes = _get_ref_counts_func(hash_mesh_dict)
        else:
            related_hashes = None
            ref_counts_from_hashes = None

        # Set weight style: regular or context
        weight = 'context_weight' if _is_context_weighted(
            mesh_id_list=options['mesh_ids'],
            strict_filtering=options['strict_mesh_id_filtering']
        ) else options['weight']

        dijkstra_gen = open_dijkstra_search(graph, starting_node,
                                            reverse=reverse,
                                            hashes=related_hashes,
                                            terminal_ns=terminal_ns,
                                            weight=weight,
                                            ref_counts_function=
                                            ref_counts_from_hashes,
                                            ignore_nodes=ignore_nodes,
                                            const_c=options['const_c'],
                                            const_tk=options['const_tk'])
        return self._loop_open_paths(graph, dijkstra_gen,
                                     source_node=starting_node,
                                     reverse=reverse, **options)

    def _loop_open_paths(self, graph, open_path_gen, source_node, reverse,
                         **options):
        result = defaultdict(list)
        max_results = int(options['max_results']) \
            if options.get('max_results') is not None else self.MAX_PATHS
        added_paths = 0
        _ = options.pop('source', None)
        _ = options.pop('target', None)

        collect_weights = _get_collect_weights_func(graph, **options)

        # Loop paths
        while True:
            try:
                path = next(open_path_gen)
                # Skip to specific length if requested; does not apply if
                # doing regular weighted or context weighted search
                if options['path_length'] and \
                        not _is_weighted(weight=options['weight'],
                                         mesh_ids=options['mesh_ids'],
                                         strict_mesh_id_filtering=options[
                                             'strict_mesh_id_filtering']):
                    while len(path) != options['path_length']:
                        path = next(open_path_gen)

                # Reverse path if reverse search
                path = path[::-1] if reverse else path

                # Collect weights
                weights = collect_weights(path)

            except StopIteration:
                logger.info('Reached StopIteration, all BFS paths found, '
                            'breaking')
                break

            # Handle signed path
            if options.get('sign') is not None:
                edge_signs = [signed_nodes_to_signed_edge(s, t)[2]
                              for s, t in zip(path[:-1], path[1:])]
                path = [n[0] for n in path]
                graph_type = 'signed'
            else:
                edge_signs = None
                graph_type = 'digraph'
            hash_path = self._get_hash_path(path=path, source=source_node,
                                            weights=weights,
                                            edge_signs=edge_signs,
                                            graph_type=graph_type,
                                            **options)

            # Assemble results
            if hash_path and all(hash_path):
                result[len(path)].append({
                    'path': path,
                    'stmts': hash_path,
                    'sort_key': str(self._get_sort_key(path, hash_path,
                                                       edge_signs)),
                    'cost': str(self._get_cost(path, edge_signs))
                })
                added_paths += 1

                if added_paths >= max_results:
                    logger.info('Max bfs paths found, returning')
                    return result

        return result

    def multi_regulators_targets(self, list_of_regulators=None,
                                 list_of_targets=None, **options):
        if not (bool(list_of_regulators) ^ bool(list_of_targets)):
            logger.warning('Must provide either of targets OR regulators. '
                           'None or both are not allowed.')
            return {}
        if list_of_targets:
            return self.direct_interactors_multi(
                list_of_targets=list_of_targets, **options)
        else:
            return self.direct_interactors_multi(
                list_of_regulators=list_of_regulators, **options)

    def direct_interactors_multi(self, list_of_regulators=None,
                                 list_of_targets=None, **options):
        """Returns a list of statement data that connect list_of_targets for
        upstream regulators

        Parameters
        ----------
        list_of_regulators : list(str)
            A list of nodes to look for shared targets for
        list_of_targets : list(str)
            A list of nodes to look for shared regulators for
        options : kwargs
            Options have to include (see self.handle_query for explanations):
                *node_filter
                *bsco
                *stmt_filter
                *curated_db_only

        Returns
        -------
        dict
            Dictionary containing the results:
                {'targets': <list of targets>,
                 'regulators': <list of regulators>,
                 'stmt_data': <dict of statements data for the edges>,
                 'stmt_hashes': <dict of statement hashes>}
            If there are no results, an empty dictionary is returned.
        """
        def grounding_filter(intrcts):
            # Make sure we have grounded node names
            if not all([n in self.nodes for n in intrcts]):
                not_in_network = [n for n in intrcts
                                  if n not in self.nodes]
                in_network = [n for n in intrcts if n in self.nodes]
                grounded = []
                for trgt in not_in_network:
                    _, _, name = get_top_ranked_name(trgt)
                    if name is None:
                        # abort if one of the targets is ungroundable
                        logger.warning('Target %s is ungroundable' % trgt)
                        return []
                    if name not in self.nodes:
                        logger.warning(
                            'Target %s (grounded to %s) is not a node '
                            'in the graph' % (trgt, name))
                    grounded.append(name)
                assert len(grounded) + len(in_network) == \
                    len(intrcts)
                return in_network + grounded
            else:
                return intrcts

        allowed_ns = options['node_filter']

        input_interactors = list_of_targets if list_of_targets else \
            list_of_regulators

        input_interactors = grounding_filter(input_interactors)

        # Get the intersection of all direct interactors
        first = input_interactors[0]
        other_interactors = set(self.nx_dir_graph_repr.pred[first]) if \
            list_of_targets else set(self.nx_dir_graph_repr.succ[first])
        for interactor in input_interactors[1:]:
            if list_of_targets:
                other_interactors.intersection_update(
                    set(self.nx_dir_graph_repr.pred[interactor])
                )
            else:
                other_interactors.intersection_update(
                    set(self.nx_dir_graph_repr.succ[interactor])
                )

        # Filter to allowed ns
        allowed_other_interactors = {n for n in other_interactors
                                     if self.nodes[n]['ns'].lower() in
                                     allowed_ns}

        if allowed_other_interactors:
            # If targets were input interactors
            if list_of_targets:
                targets = input_interactors
                regulators = allowed_other_interactors
                ign_nodes = 'targets'
            # If regulators were input interactors
            else:
                targets = allowed_other_interactors
                regulators = input_interactors
                ign_nodes = 'regulators'
            return self._loop_direct_regulators_multi(
                targets=targets,
                regulators=regulators,
                ign=ign_nodes,
                **options
            )
        return {}

    def _loop_direct_regulators_multi(self, targets, regulators,
                                      ign, **options):
        result = defaultdict(list)
        all_hashes = []
        for reg in regulators:
            data = {}
            hashes = []
            for target in targets:
                # get hash path for each target-regulator pair
                ign_node = reg if ign == 'regulators' else target
                hash_path = self._get_hash_path(path=[reg, target],
                                                source=ign_node, **options)
                if hash_path and hash_path[0]:
                    result['2'].append({
                        'path': [reg, target],
                        'stmts': hash_path,
                        'sort_key': str(self._get_sort_key([reg, target],
                                                           hash_path,
                                                           None)),
                        'cost': str(self._get_cost([reg, target]))
                    })
                    data[target] = hash_path[0]
                    # The hash path will be a list of len 1 since we only
                    # have one direct edge
                    for key, dl in hash_path[0].items():
                        if key not in {'subj', 'obj'}:
                            hashes.extend([d['stmt_hash'] for d in dl])
            if hashes:
                all_hashes.extend(hashes)

        return {
            'result': {
                'paths_by_node_count': {
                    'forward': result,
                    'path_hashes': all_hashes
                },
            },
            'targets': list(targets),
            'regulators': list(regulators)
        }

    def find_shared_regulators(self, source, target, **options):
        """Returns a list of statement data that explain shared regulators
        for source and target"""
        if source in self.nodes and target in self.nodes:
            source_pred = set(self.nx_dir_graph_repr.pred[source].keys())
            target_pred = set(self.nx_dir_graph_repr.pred[target].keys())
            common = source_pred & target_pred
            if common:
                try:
                    return self._loop_shared_regulators(shared_regs=common,
                                                        source=source,
                                                        target=target,
                                                        **options)
                except NetworkXNoPath as e:
                    logger.warning(repr(e))

        return []

    def _loop_shared_regulators(self, shared_regs, source, target,
                                **options):
        """Order shared regulators by lowest highest belief score"""
        ordered_regulators = []
        added_regulators = 0
        for sr in shared_regs:
            paths1 = self._get_hash_path(target=source, path=[sr, source],
                                         **options)
            paths2 = self._get_hash_path(target=target, path=[sr, target],
                                         **options)
            if paths1 and paths2 and paths1[0] and paths2[0]:
                paths1_stmts = []
                for k, v in paths1[0].items():
                    if k not in {'subj', 'obj', 'weight_to_show'}:
                        paths1_stmts.extend(v)
                paths2_stmts = []
                for k, v in paths2[0].items():
                    if k not in {'subj', 'obj', 'weight_to_show'}:
                        paths2_stmts.extend(v)
                max_belief1 = max([st['belief'] for st in paths1_stmts])
                max_belief2 = max([st['belief'] for st in paths2_stmts])
                ordered_regulators.append({
                    sr: [paths1, paths2],
                    'lowest_highest_belief': min(max_belief1, max_belief2)
                })
                added_regulators += 1
                if added_regulators >= self.MAX_PATHS:
                    if self.verbose:
                        logger.info('Max number of shared regulators '
                                    'reached. Breaking loop')
                    break
        if ordered_regulators:
            return sorted(ordered_regulators,
                          key=lambda m: m['lowest_highest_belief'],
                          reverse=True)
        else:
            return []

    def find_common_targets(self, source, target, **options):
        """Returns a list of statement data that explain common targets for
        source and target"""
        if source in self.nodes and target in self.nodes:
            source_succ = set(self.nx_dir_graph_repr.succ[source].keys())
            target_succ = set(self.nx_dir_graph_repr.succ[target].keys())
            common = source_succ & target_succ
            if common:
                try:
                    return self._loop_common_targets(common_targets=common,
                                                     source=source,
                                                     target=target,
                                                     **options)
                except NetworkXNoPath as e:
                    logger.warning(repr(e))

        return []

    def _loop_common_targets(self, common_targets, source, target, **options):
        """Order common_targets targets by lowest belief in pair."""
        ordered_commons = []
        added_targets = 0
        for ct in common_targets:
            paths1 = self._get_hash_path(source=source, path=[source, ct],
                                         **options)
            paths2 = self._get_hash_path(source=target, path=[target, ct],
                                         **options)
            if paths1 and paths2 and paths1[0] and paths2[0]:
                paths1_stmts = []
                for k, v in paths1[0].items():
                    if k not in {'subj', 'obj', 'weight_to_show'}:
                        paths1_stmts.extend(v)
                paths2_stmts = []
                for k, v in paths2[0].items():
                    if k not in {'subj', 'obj', 'weight_to_show'}:
                        paths2_stmts.extend(v)
                max_belief1 = max([st['belief'] for st in paths1_stmts])
                max_belief2 = max([st['belief'] for st in paths2_stmts])
                ordered_commons.append({
                    ct: [paths1, paths2],
                    'lowest_highest_belief': min(max_belief1, max_belief2)
                })
                added_targets += 1
                if added_targets >= self.MAX_PATHS:
                    if self.verbose:
                        logger.info('Max number of common targets reached. '
                                    'Breaking loop')
                    break
        if ordered_commons:
            return sorted(ordered_commons,
                          key=lambda k: k['lowest_highest_belief'],
                          reverse=True)
        else:
            return []

    def _loop_paths(self, graph, source, target, paths_gen, **options):
        # len(path) = edge count + 1
        sign = options['sign']
        graph_type = 'digraph' if sign is None else 'signed'
        if graph_type == 'signed':
            logger.info('Looping signed paths')
        path_len = options['path_length'] + 1 if \
            options['path_length'] and not options['weight'] else False
        result = defaultdict(list)
        prev_path = None
        added_paths = 0
        skipped_paths = 0
        culled_nodes = set()
        culled_edges = set()  # Currently unused, only operate on node level

        collect_weights = _get_collect_weights_func(graph, **options)

        while True:
            # Check if we found k paths
            if added_paths >= self.MAX_PATHS:
                logger.info('Found all %d shortest paths, returning results.'
                            % self.MAX_PATHS)
                return result
            if time() - self.query_receive_time > self.TIMEOUT:
                logger.info('Reached timeout (%d s) before finding all %d '
                            'shortest paths. Returning search.' %
                            (self.TIMEOUT, MAX_PATHS))
                self.query_timed_out = True
                return result
            # Check if we have to cull the best node, this is the case
            # if the modulo is 1, meaning that in the *following* path we
            # want another node culled
            send_values = None
            if sign is None and (added_paths % options.get(
                    'cull_best_node', float('NaN')) == 1 and
                    prev_path is not None and len(prev_path['path']) >= 3):
                degrees = self.nx_dir_graph_repr.degree(
                    prev_path['path'][1:-1], options.get('weight', None))
                node_highest_degree = max(degrees, key=lambda x: x[1])[0]
                culled_nodes.add(node_highest_degree)
                send_values = (culled_nodes, culled_edges)
                if self.verbose > 1:
                    logger.info('Culled nodes: %s' % repr(culled_nodes))
            try:
                if sign is not None:
                    signed_path_nodes = next(paths_gen)
                    weights = collect_weights(signed_path_nodes)
                    path = [n[0] for n in signed_path_nodes]
                    edge_signs = [signed_nodes_to_signed_edge(s, t)[2]
                                  for s, t in zip(signed_path_nodes[:-1],
                                                  signed_path_nodes[1:])]
                else:
                    # Get next path and send culled nodes and edges info for
                    # the path in the following iteration
                    path = paths_gen.send(send_values)
                    weights = collect_weights(path)
                    edge_signs = None
            except StopIteration:
                logger.info('Reached StopIteration: all paths found. '
                            'breaking.')
                break
            # Todo: skip to correct length here already
            hash_path = self._get_hash_path(source=source, target=target,
                                            path=path, weights=weights,
                                            edge_signs=edge_signs,
                                            graph_type=graph_type,
                                            **options)

            if hash_path and all(hash_path):
                if self.verbose > 1:
                    logger.info('Adding stmts and path from %s to path list' %
                                repr(hash_path))
                pd = {'stmts': hash_path,
                      'path': path,
                      'weight_to_show': weights,
                      'cost': str(self._get_cost(path, edge_signs)),
                      'sort_key': str(self._get_sort_key(path, hash_path,
                                                         edge_signs))}
                if not path_len or (path_len and path_len == len(path)):
                    result[len(path)].append(pd)
                    prev_path = pd
                    added_paths += 1
                elif path_len and len(path) < path_len:
                    continue
                elif path_len and len(path) > path_len:
                    if self.verbose > 1:
                        logger.info('Max path length reached, returning '
                                    'results.')
                    return result
                elif options['weight'] and graph_type == 'signed' and\
                        len(path) > 7:
                    logger.warning('Extremely long signed paths detected, '
                                   'aborting for to avoid long wait for '
                                   'Networkx to return further paths.')
                    return result
                else:
                    logger.warning('This option should not happen')
            else:
                skipped_paths += 1
        if self.verbose > 2:
            logger.info('Done looping paths. Returning result: %s' %
                        repr(result))
        return result

    def has_path(self, source, target):
        """Return true if there is a path from source to target"""
        return nx.has_path(self.nx_dir_graph_repr, source, target)

    def get_common_parents(self, **options):
        """Find common parents between source and target"""
        # Try, in order:
        #   1. ns:id from node dict
        #   2. ns:id from grounding service
        #   3. go with original node name and try HGNC and FPLX

        source_ns, source_id, target_ns, target_id = None, None, None, None

        # Source
        if options['source'] in self.nodes:
            source_id = self.nodes[options['source']]['id']
            source_ns = self.nodes[options['source']]['ns']
        else:
            source_ns, source_id, source_norm_name = \
                nf.gilda_normalization(options['source'])
            if not source_id:
                source_id = options['source']

        # Target
        if options['target'] in self.nodes:
            target_id = self.nodes[options['target']]['id']
            target_ns = self.nodes[options['target']]['ns']
        else:
            target_ns, target_id, target_norm_name = \
                nf.gilda_normalization(options['target'])
            if not target_id:
                target_id = options['target']

        # Initialize result dict
        cp_results = {'source_ns': source_ns, 'source_id': source_id,
                      'target_ns': target_ns, 'target_id': target_id,
                      'common_parents': []}
        cp = set()

        # Try different combinations of ns combinations

        # If both source and target are given
        if source_ns and target_ns:
            if source_ns.lower() in options['node_filter'] and \
                    target_ns.lower() in options['node_filter']:
                if self.verbose > 1:
                    logger.info('Looking for common parents using namespaces '
                                'found in network')
                cp = ff.common_parent(ns1=source_ns, id1=source_id,
                                      ns2=target_ns, id2=target_id)
            else:
                logger.info(f'The namespaces for {source_ns} and/or '
                            f'{target_ns} are not in node filter: '
                            f'({", ".join(options["node_filter"])})'
                            f'Aborting common parent search.')
                cp_results['common_parents'] = []
                return cp_results

        # If only target ns is given
        if not source_ns and target_ns:
            if target_ns.lower() in options['node_filter']:
                if self.verbose > 1:
                    logger.info('No namespace found for %s, trying HGNC and '
                                'FPLX.' % source_id)
                for sns in ['HGNC', 'FPLX']:
                    if sns.lower() not in options['node_filter']:
                        continue
                    else:
                        cp = ff.common_parent(ns1=sns, id1=source_id,
                                              ns2=target_ns, id2=target_id)
                        if cp:
                            if self.verbose:
                                logger.info('Found common parents with source '
                                            'ns %s' % sns)
                            break
            else:
                logger.info('The namespaces for %s is not in node filter. '
                            'Aborting common parent search.' % target_id)
                cp_results['common_parents'] = []
                return cp_results

        # If only source ns is given
        if not target_ns and source_ns:
            if source_ns.lower() in options['node_filter']:
                if self.verbose > 1:
                    logger.info('No namespace found for %s, trying HGNC and '
                                'FPLX.' % target_id)
                for tns in ['HGNC', 'FPLX']:
                    if tns.lower() not in options['node_filter']:
                        continue
                    else:
                        cp = ff.common_parent(ns1=source_ns, id1=source_id,
                                              ns2=tns, id2=target_id)
                        if cp:
                            if self.verbose:
                                logger.info('Found common parents with source '
                                            'ns %s' % tns)
                            break
            else:
                logger.info('The namespaces for %s is not in node filter. '
                            'Aborting common parent search.' % source_id)
                cp_results['common_parents'] = []
                return cp_results

        # If no namespaces exist
        if not source_ns and not target_ns:
            if self.verbose > 1:
                logger.info('No namespaces found for %s and %s, trying HGNC '
                            'and FPLX' % (source_id, target_id))
            for source_ns in ['HGNC', 'FPLX']:
                if source_ns.lower() not in options['node_filter']:
                    continue
                for target_ns in ['HGNC', 'FPLX']:
                    if target_ns.lower() not in options['node_filter']:
                        continue
                    cp = ff.common_parent(ns1=source_ns, id1=source_id,
                                          ns2=target_ns, id2=target_id)
                    if cp:
                        break

        if not cp:
            logger.info('No common parents found')
            cp_results['common_parents'] = []
            return cp_results
        else:
            cp_list = [(ns, _id, ff.get_identifiers_url(ns, _id))
                       for ns, _id in cp]
            cp_results['common_parents'] = sorted(cp_list,
                                                  key=lambda t: (t[0], t[1]))
            return cp_results

    def _get_edges(self, s, o, edge_sign=None, graph='digraph'):
        """Return edges from one of the loaded graphs in a uniform format"""
        if graph == 'multi':
            if not self.mdg_edges:
                raise nx.NetworkXException('MultiDiGraph not loaded')
            i = 0
            edge = self.mdg_edges.get((s, o, i))
            while edge:
                yield edge
                i += 1
                edge = self.mdg_edges.get((s, o, i))

        elif graph == 'signed':
            if edge_sign is None:
                raise nx.NetworkXException('Argument edge_sign needs to be '
                                           'specified to get signed edge.')
            sign = nf.SIGNS_TO_INT_SIGN[edge_sign]
            i = 0
            try:
                stmt_edge = self.signed_edges[(s, o, sign)][
                    'statements'][i]
                while stmt_edge:
                    yield stmt_edge
                    i += 1
                    stmt_edge = self.signed_edges[(s, o, sign)][
                        'statements'][i]
            except IndexError:
                return

        else:
            i = 0
            try:
                stmt_edge = self.dir_edges[(s, o)]['statements'][i]
                while stmt_edge:
                    yield stmt_edge
                    i += 1
                    stmt_edge = self.dir_edges[(s, o)]['statements'][i]
            except IndexError:
                return

    def _get_hash_path(self, path, source=None, target=None,
                       edge_signs=None, graph_type='digraph', weights=None,
                       **options):
        """Return a list of n-1 lists of dicts containing of stmts connecting
        the n nodes in path"""
        hash_path = []
        es = edge_signs if edge_signs else [None]*(len(path)-1)
        weights = weights if weights else ['N/A']*(len(path)-1)
        if self.verbose:
            logger.info('Building evidence for path %s' % str(path))
        for subj, obj, edge_sign, w in zip(path[:-1], path[1:], es, weights):
            # Check node filter, but ignore source or target nodes
            # e.g., check node_filter IFF source != subj AND target != obj
            if (source != subj and target != subj and
                self.nodes[subj]['ns'].lower() not in options['node_filter'])\
                or \
                (source != obj and target != obj and
                 self.nodes[obj]['ns'].lower() not in options['node_filter']):
                if self.verbose:
                    logger.info('Node namespace %s or %s not part of '
                                'acceptable namespaces %s' %
                                (self.nodes[subj]['ns'],
                                 self.nodes[obj]['ns'],
                                 options['node_filter']))
                return []

            # Initialize edges dict
            edges = {}

            # Get first edge statement
            edge_stmts = self._get_edges(subj, obj, edge_sign, graph_type)
            if self.verbose > 3:
                logger.info('First edge stmt %s' % repr(next(edge_stmts)))

            # Exhaustively loop through all edge statements
            for edge_stmt in edge_stmts:
                # If edge statement passes, append to edges list
                if self._pass_stmt(edge_stmt, **options):
                    # convert hash to string for javascript compatibility
                    edge_stmt['stmt_hash'] = str(edge_stmt['stmt_hash'])
                    # ToDo english assemble statements per type
                    try:
                        edges[edge_stmt['stmt_type']].append(edge_stmt)
                    except KeyError:
                        edges['subj'] = subj
                        edges['obj'] = obj
                        edges['weight_to_show'] = str(w)
                        edges[edge_stmt['stmt_type']] = [edge_stmt]
                    if self.verbose > 3:
                        logger.info('edge stmt passed filter, appending to '
                                    'edge list.')
                        logger.info('Next edge stmt %s' % repr(edge_stmt))

            # If edges list contains anything, append to hash_path list
            if edges:
                if self.verbose > 4:
                    logger.info('Appending %s to hash path list' % repr(edges))
                hash_path.append(edges)
            else:
                return []
        if self.verbose > 1 and len(hash_path) > 0:
            logger.info('Returning hash path: %s' % repr(hash_path))
        return hash_path

    def _pass_stmt(self, edge_stmt, **options):
        """Returns True if edge_stmt passes the filters below"""
        # Failsafe for empty statements
        if not edge_stmt:
            logger.warning('No edge statement')
            return False

        # Filter belief score
        if edge_stmt['belief'] < options['bsco']:
            if self.verbose:
                logger.info('Did not pass belief score')
            return False

        # Filter statement type
        if edge_stmt['stmt_type'].lower() in options['stmt_filter']:
            if self.verbose > 4:
                logger.info('statement type %s found in filter %s'
                            % (edge_stmt['stmt_type'],
                               str(options['stmt_filter'])))
            return False

        # Filter out statements with only readers as sources
        if options['curated_db_only'] and not edge_stmt['curated']:
            return False

        # Filter stmt hash
        if options.get('edge_hash_blacklist', None) and \
                edge_stmt['stmt_hash'] in options['edge_hash_blacklist']:
            if self.verbose > 3:
                logger.info('hash %s is blacklisted, skipping' %
                            edge_stmt['stmt_hash'])
            return False

        # Return True is all filters were passed
        return True

    def _get_cost(self, path, edge_signs=None, direct=True, key='weight'):
        if direct:
            # Return sum of aggregated weights per edge
            if edge_signs:
                return sum(self.signed_edges[(s, o, es)][key] for
                           s, o, es in zip(path[:-1], path[1:], edge_signs))
            else:
                return sum(self.dir_edges[(s, o)][key] for s, o in
                           zip(path[:-1], path[1:]))
        else:
            # Return sum of averaged weights per stmts
            cost = 0
            for s, o in zip(path[:-1], path[1:]):
                ew = [e[key] for e in self._get_edges(s, o, direct)]
                cost += sum(ew)/len(ew)
            return cost

    def _aggregated_path_belief(self, path, edge_signs):
        if edge_signs:
            belief_list = [self.signed_edges[e]['belief'] for e in
                           zip(path[:-1], path[1:], edge_signs)]
        else:
            belief_list = [self.dir_edges[e]['belief'] for e in
                           zip(path[:-1], path[1:])]
        return nf.ag_belief_score(belief_list)

    def _get_sort_key(self, path, hash_path, edge_signs=None, method=None):
        """Calculate a number to sort the path on

        `Method` allows to specify the calculation"""

        # Default: aggregated path belief score
        sort_key = self._aggregated_path_belief(path, edge_signs)
        return sort_key

    @staticmethod
    def _sort_stmts(ksp):
        for pl in ksp:
            res_list = ksp[pl]
            ksp[pl] = sorted(res_list, key=lambda pd: pd['sort_key'],
                             reverse=True)
        return ksp

    def _uri_by_node(self, node):
        """Return the fplx URI for the provided node"""
        # Check existence of node outside function
        node_id = self.nodes[node]['id']
        node_ns = self.nodes[node]['ns']
        return get_identifiers_url(node_ns, node_id)

    def _get_parents(self, node):
        if self.nodes.get(node):
            db_id = node
            ns = self.nodes[node]['ns']

            true_ns, true_id, norm_name = nf.gilda_normalization(db_id)
            if true_ns and true_id:
                ns, db_id = true_ns, true_id
            parents = nf.bio_ontology.get_parents(ns, db_id)
            return {get_identifiers_url(*p) for p in parents}
        else:
            return set()


def _get_collect_weights_func(graph, **options):
    if options['mesh_ids']:
        if options['strict_mesh_id_filtering']:
            def _f(path):
                return ['N/A'] * (len(path) - 1)
        else:
            def _f(path):
                return [_truncate(graph[u][v]['context_weight'])
                        for u, v in zip(path[:-1], path[1:])]
    else:
        if options.get('weight', None):
            def _f(path):
                return [_truncate(graph[u][v][options['weight']])
                        for u, v in zip(path[:-1], path[1:])]
        else:
            def _f(path):
                return ['N/A'] * (len(path) - 1)
    return _f


def _get_ref_counts_func(hash_mesh_dict):
    def _func(graph, u, v):
        # Get hashes for edge
        hashes = [d['stmt_hash'] for d in graph[u][v]['statements']]

        # Get all relevant mesh counts
        dicts = [hash_mesh_dict.get(h, {'': 0, 'total': 1})
                 for h in hashes]

        # Count references
        ref_counts = sum(sum(v for k, v in d.items() if k != 'total')
                         for d in dicts)
        total = sum(d['total'] for d in dicts)
        return ref_counts, total if total else 1
    return _func


def get_top_ranked_name(name, context=None):
    """get top ranked result from gilda

    Parameters
    ----------
    name : str
        Provide a name to be grounded
    context : str
        Optional. Provide context for the name provided

    Returns
    -------
    tuple(str)
        A tuple of ns, id, grounded name
    """
    none_triple = (None, )*3
    if not GRND_URI:
        logger.warning('Grounding service URL not set')
        return none_triple
    req_json = {'text': name}
    if context:
        req_json['context'] = context
    res = requests.post(GRND_URI, json=req_json)
    if res.status_code == 200 and res.json():
        top_res = res.json()[0]
        topname = top_res['term']['entry_name']
        topns = top_res['term']['db']
        topid = top_res['term']['id']
        return topns, topid, topname
    elif res.status_code != 200:
        logger.warning('Got status code %d from gilda, no result')
    else:
        logger.warning('No result from gilda for %s' % name)
    return none_triple


def _open_ended_common_search(G, node_set_queue, allowed_ns=None,
                              max_depth=4):
    # Search common upstreams of all nodes in initial_nodes

    # Get first node and its upstream nodes as a set
    first_node = node_set_queue[-1][0]
    upstreams = set(G.pred[first_node])

    # Get the intersection of all upstreams
    for node in node_set_queue[-1][1:]:
        upstreams.intersection_update(set(G.pred[node]))

    # Filter out the upstream set to allowed_ns if present
    if allowed_ns:
        upstreams = {n for n in upstreams if G.nodes[n]['ns'].lower() in
                     allowed_ns}

    if upstreams:
        node_set_queue.append(list(upstreams))

    if not upstreams or len(node_set_queue) >= max_depth:
        return node_set_queue

    else:
        return _open_ended_common_search(G, node_set_queue, allowed_ns,
                                         max_depth)


def translate_query(query_json):
    """Translate query json"""
    options = {k: v for k, v in query_json.items()  # Handled below
               if k not in ['sign', 'weighted']}
    if 'sign_dict' not in options:
        options['sign_dict'] = default_sign_dict
    for k, v in query_json.items():
        if k == 'weighted':
            logger.info('Doing %sweighted path search' % 'un' if not v
                        else '')
            options['weight'] = 'weight' if v else None
        if k == 'sign':
            # Positive regulation: 0; Negative regulation: 1;
            options[k] = INT_PLUS if nf.SIGN_TO_STANDARD.get(v) == '+' else \
                (INT_MINUS if nf.SIGN_TO_STANDARD.get(v) == '-' else None)
        if k == 'edge_hash_blacklist' and options.get(k) and \
                isinstance(options[k][0], int):
            options[k] = [str(i) for i in options[k]]
        if k in ['node_filter', 'stmt_filter']:
            options[k] = [s.lower() for s in options[k]]
        if k == "cull_best_node":
            options[k] = int(v) if v >= 1 else float('NaN')
    return options


def _is_context_weighted(mesh_id_list, strict_filtering):
    if mesh_id_list and not strict_filtering:
        return True
    return False


def _is_weighted(weight, mesh_ids, strict_mesh_id_filtering, **kwargs):
    if mesh_ids:
        ctx_w = _is_context_weighted(mesh_id_list=mesh_ids,
                                     strict_filtering=strict_mesh_id_filtering)
        return bool(weight) or ctx_w
    else:
        return bool(weight)


def list_all_hashes(ksp_results):
    hash_set = set()
    for path_length in ksp_results:
        for res in ksp_results[path_length]:
            stmt_dict_list = res['stmts']
            for stmt_dict in stmt_dict_list:
                for stmt_type in stmt_dict:
                    if stmt_type in ('subj', 'obj', 'weight_to_show'):
                        continue
                    meta_list = stmt_dict[stmt_type]
                    hash_set.update([item['stmt_hash'] for item in
                                     meta_list])
    return list(hash_set)


def get_signed_node(node, sign, reverse):
    """Given sign and direction, return a node

    Assign the correct sign to the source node:
    If search is downstream, source is the first node and the search must
    always start with + as node sign. The leaf node sign (i.e. the end of
    the path) in this case will then be determined by the requested sign.

    If reversed search, the source is the last node and can have
    + or - as node sign depending on the requested sign.
    """
    if sign is None:
        return node
    else:
        # Upstream: return asked sign
        if reverse:
            return node, sign
        # Downstream: return positive node
        else:
            return node, INT_PLUS
