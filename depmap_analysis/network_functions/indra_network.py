import logging
import requests
import networkx as nx
from itertools import product
from collections import defaultdict
from networkx import NodeNotFound, NetworkXNoPath
from time import time, gmtime, strftime

from indra.config import CONFIG_DICT

import depmap_network_functions as dnf
from depmap_analysis.network_functions.network_functions import ag_belief_score, shortest_simple_paths

logger = logging.getLogger('indra network')

GRND_URI = None
try:
    GRND_URI = CONFIG_DICT['INDRA_GROUNDING_SERVICE_URL']
except KeyError:
    logger.warning('Indra Grounding service not available. Add '
                   'INDRA_GROUNDING_SERVICE_URL to `indra/config.ini`')

MAX_PATHS = 50
TIMEOUT = 30  # Timeout in seconds


class IndraNetwork:
    """Handle searches and graph output of the INDRA DB network"""
    def __init__(self, indra_dir_graph=nx.DiGraph(),
                 indra_multi_dir_graph=nx.MultiDiGraph()):
        self.nx_dir_graph_repr = indra_dir_graph
        self.nx_md_graph_repr = indra_multi_dir_graph
        self.nodes = self.nx_dir_graph_repr.nodes
        self.dir_edges = indra_dir_graph.edges
        self.mdg_edges = indra_multi_dir_graph.edges
        self.ehm = indra_dir_graph.graph.get('entity_hierarchy_manager', None)
        self.node_by_uri = indra_dir_graph.graph.get('node_by_uri', None)
        self.MAX_PATHS = MAX_PATHS
        self.small = False
        self.verbose = 0
        self.query_recieve_time = 0.0
        self.query_timed_out = False

    def handle_query(self, **kwargs):
        """Handles path query from client. Returns query result.

        The query is a json-friendly key-value structure contained in kwargs
        with the following parameters:

        (Note that parameters that are not yet implemented are not mandatory
        and have no effect on the path search if provided)

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
        path_length: int|False
            a positive integer stating the number of edges that should be in
            the returned path. If False, return paths with any number of edges.
        sign: str ['no_sign'|'plus'|'minus'] **currently not implemented**
            Placeholder for future implementation of path searches in signed
            graphs
        weighted: Bool
            If True, do a weighted path search. Weights in the network are
            assigned as -log(belief score)
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

        Returns
        -------
        result : dict('paths_by_node_count'=ksp,
                      'common_targets'=ct,
                      'common_parents'=cp)
            A dict containing the results from each path search and a flag
            for timeout:
            ksp : dict(int)
                Dict keyed by node count with the results of directed paths
            ct : dict('target')
                List of dicts keyed by common target name, sorted on highest
                lowest belief score
            cp : dict
                Dict with result of common parents search together which
                ns:id pairs were used to resolve the query
            timeout : Bool
                True if the query timed out
        """
        self.query_recieve_time = time()
        self.query_timed_out = False
        logger.info('Query received at %s' %
                    strftime('%Y-%m-%d %H:%M:%S (UTC)',
                             gmtime(self.query_recieve_time)))
        if not self.sanity_check(**kwargs):
            return {'paths_by_node_count': {},
                    'common_targets': [],
                    'common_parents': {},
                    'timeout': False}
        mandatory = ['source', 'target', 'stmt_filter', 'node_filter',
                     'path_length', 'weighted', 'bsco', 'fplx_expand',
                     'k_shortest', 'curated_db_only']
        if not all([key in kwargs for key in mandatory]):
            miss = [key in kwargs for key in mandatory].index(False)
            raise KeyError('Missing mandatory parameter %s' % mandatory[miss])
        options = {k: v for k, v in kwargs.items()  # Handled below
                   if k not in ['sign', 'weighted']}
        for k, v in kwargs.items():
            if k == 'weighted':
                logger.info('Doing %sweighted path search' % 'un' if not v
                            else '')
                options['weight'] = 'weight' if v else None
            if k == 'sign':
                options[k] = 1 if v == 'plus' \
                    else (-1 if v == 'minus' else 0)
            if k == 'edge_hash_blacklist' and options.get(k) and \
                    isinstance(options[k][0], int):
                options[k] = [str(i) for i in options[k]]
            if k in ['node_filter', 'stmt_filter']:
                options[k] = [s.lower() for s in options[k]]
        k_shortest = kwargs.pop('k_shortest', None)
        self.MAX_PATHS = k_shortest if k_shortest else MAX_PATHS
        logger.info('Query translated to: %s' % repr(options))
        logger.info('Looking for no more than %d paths' % self.MAX_PATHS)
        ksp = self.find_shortest_paths(**options)
        if not ksp:
            ckwargs = options.copy()
            if kwargs['fplx_expand']:

                logger.info('No directed path found, looking for paths '
                            'connected by common parents of source and/or '
                            'target')
                ksp = self.try_parents(**ckwargs)
                if self.verbose > 2:
                    logger.info('Got parents search result: %s' % repr(ksp))

            if not ksp and GRND_URI:
                ksp = self.grounding_fallback(**ckwargs)
        if not ksp:
            logger.info('No directed path found')
        if ksp and not options['weight']:
            # Sort the results in ksp if non-weighted search
            ksp = self._sort_stmts(ksp)
        ct = self.find_common_targets(**options)
        cp = self.get_common_parents(**options)
        return {'paths_by_node_count': ksp,
                'common_targets': ct,
                'common_parents': cp,
                'timeout': self.query_timed_out}

    @staticmethod
    def sanity_check(**options):
        """Checks for some possible gotchas in query"""
        # Check non-resolving query
        sns, sid = dnf._ns_id_from_name(options['source'])
        tns, tid = dnf._ns_id_from_name(options['target'])
        if (sns and sns.lower() not in options['node_filter']) or \
                (tns and tns.lower() not in options['node_filter']):
            if sns.lower() not in options['node_filter']:
                logger.warning('%s not among accepted nodes' % sns)
            if tns.lower() not in options['node_filter']:
                logger.warning('%s not among accepted nodes' % tns)
            return False

        return True

    def grounding_fallback(self, **ckwargs):
        """Retry search with alternative names found by grounding service"""
        if self.verbose:
            logger.info('Expanding search using grounding service')
        org_source = ckwargs['source']
        org_target = ckwargs['target']

        # ToDo establish grounding priority when scores are equal between
        #  groundings

        # Get groundings
        src_groundings = requests.post(GRND_URI,
                                       json={'text': org_source}).json()
        trgt_groundings = requests.post(GRND_URI,
                                        json={'text': org_target}).json()

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
            return self._loop_paths(nx.shortest_path(
                self.nx_dir_graph_repr, source, target, options['weight']),
                **options)
        except NodeNotFound or NetworkXNoPath:
            return {}

    def find_shortest_paths(self, source, target, **options):
        """Returns a list of shortest paths in ascending order"""
        try:
            logger.info('Doing simple %s path search' % 'weigthed'
                        if options['weight'] else '')
            blacklist_options = {}
            blacklist_options['ignore_nodes'] = options.get('node_blacklist', None)
            paths = shortest_simple_paths(self.nx_dir_graph_repr,
                                     source, target, options['weight'],
                                     **blacklist_options)
            # paths = nx.all_shortest_paths(self.nx_md_graph_repr,
            #                               source, target, options['weight'])
            return self._loop_paths(paths, **options)
        except NodeNotFound as e:
            logger.warning(repr(e))
            return {}
        except NetworkXNoPath as e:
            logger.warning(repr(e))
            return {}

    def find_common_targets(self, source, target, **options):
        """Returns a list of statement(?) pairs that explain common targets
        for source and target"""
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
                except NodeNotFound as e:
                    logger.warning(repr(e))
                except NetworkXNoPath as e:
                    logger.warning(repr(e))

        return []

    def _loop_common_targets(self, common_targets, source, target, **options):
        """Order common_targets targets by lowest bs in pair."""
        ordered_commons = []
        added_targets = 0
        for ct in common_targets:
            paths1 = self._get_hash_path(path=[source, ct], **options)
            paths2 = self._get_hash_path(path=[target, ct], **options)
            if paths1 and paths2 and paths1[0] and paths2[0]:
                max_bs1 = max([st['bs'] for st in paths1[0]])
                max_bs2 = max([st['bs'] for st in paths2[0]])
                ordered_commons.append({
                    ct: [sorted(paths1[0],
                                key=lambda k: k['bs'],
                                reverse=True),
                         sorted(paths2[0],
                                key=lambda k: k['bs'],
                                reverse=True)],
                    'lowest_highest_belief': min(max_bs1, max_bs2)
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

    def _loop_paths(self, paths_gen, **options):
        # len(path) = edge count + 1
        path_len = options['path_length'] + 1 if \
            options['path_length'] and not options['weight'] else False
        result = defaultdict(list)
        added_paths = 0
        skipped_paths = 0
        for path in paths_gen:
            # Check if we found k paths
            if added_paths >= self.MAX_PATHS:
                logger.info('Found all %d shortest paths, returning results.' %
                            self.MAX_PATHS)
                return result
            if time() - self.query_recieve_time > TIMEOUT:
                logger.warning('Reached timeout (%d s) before finding all %d '
                               'shortest paths. Returning search.' %
                               (TIMEOUT, MAX_PATHS))
                self.query_timed_out = True
                return result

            hash_path = self._get_hash_path(path, **options)
            if hash_path and all(hash_path):
                if self.verbose > 1:
                    logger.info('Adding stmts and path from %s to path list' %
                                repr(hash_path))
                pd = {'stmts': hash_path,
                      'path': path,
                      'cost': str(self._get_cost(path)),
                      'sort_key': str(self._get_sort_key(path, hash_path))}
                if not path_len or (path_len and path_len == len(path)):
                    result[len(path)].append(pd)
                    added_paths += 1
                elif path_len and len(path) < path_len:
                    continue
                elif path_len and len(path) > path_len:
                    if self.verbose > 1:
                        logger.info('Max path length reached, returning '
                                    'results.')
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
            source_ns, source_id = dnf._ns_id_from_name(options['source'])
            if not source_id:
                source_id = options['source']

        # Target
        if options['target'] in self.nodes:
            target_id = self.nodes[options['target']]['id']
            target_ns = self.nodes[options['target']]['ns']
        else:
            target_ns, target_id = dnf._ns_id_from_name(options['target'])
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
                cp = dnf.common_parent(ns1=source_ns, id1=source_id,
                                       ns2=target_ns, id2=target_id)
            else:
                logger.info('The namespaces for %s and/or %s are not in node '
                            'filter. Aborting common parent search.' %
                            (source_id, target_id))
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
                        cp = dnf.common_parent(ns1=sns, id1=source_id,
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
                        cp = dnf.common_parent(ns1=source_ns, id1=source_id,
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
                    cp = dnf.common_parent(ns1=source_ns, id1=source_id,
                                           ns2=target_ns, id2=target_id)
                    if cp:
                        break

        if not cp:
            logger.info('No common parents found')
            cp_results['common_parents'] = []
            return cp_results
        else:
            cp_results['common_parents'] = sorted(list(cp))
            return cp_results

    def _get_edge(self, s, o, index, simple_graph):
        """Return edges from DiGraph or MultiDigraph in a uniform format"""
        if simple_graph:
            try:
                stmt_edge = self.dir_edges.get((s, o))['stmt_list'][index]
            except IndexError:
                # To keep it consistent with below Multi DiGraph implementation
                stmt_edge = None
            return stmt_edge
        else:
            return self.mdg_edges.get((s, o, index))

    def _get_hash_path(self, path, simple_graph=True, **options):
        """Return a list of n-1 lists of dicts containing of stmts connecting
        the n nodes in path. If simple_graph is True, query edges from DiGraph
        and not from MultiDiGraph representation"""
        hash_path = []
        if self.verbose:
            logger.info('Building evidence for path %s' % str(path))
        for subj, obj in zip(path[:-1], path[1:]):
            # Check node filter
            if self.nodes[subj]['ns'].lower() not in \
                    options['node_filter'] or self.nodes[obj]['ns'].lower() \
                    not in options['node_filter']:
                if self.verbose:
                    logger.info('Node namespace %s or %s not part of '
                                'acceptable namespaces %s' %
                                (self.nodes[subj]['ns'],
                                 self.nodes[obj]['ns'],
                                 options['node_filter']))
                return []

            # Initialize edges list, statement index
            edges = []
            e = 0

            # Get first edge statement
            edge_stmt = self._get_edge(subj, obj, e, simple_graph)
            if self.verbose > 3:
                logger.info('First edge stmt %s' % repr(edge_stmt))

            # Exhaustively loop through all edge statments
            while edge_stmt:

                # If edge statement passes, append to edges list
                if self._pass_stmt(edge_stmt, **options):
                    # convert hash to string for javascript compatability
                    edge_stmt['stmt_hash'] = str(edge_stmt['stmt_hash'])
                    edges.append({**edge_stmt,
                                  'subj': subj,
                                  'obj': obj})
                    if self.verbose > 3:
                        logger.info('edge stmt passed filter, appending to '
                                    'edge list.')
                        logger.info('Next edge stmt %s' % repr(edge_stmt))

                # Incr statement index, get next edge statement
                e += 1
                edge_stmt = self._get_edge(subj, obj, e, simple_graph)

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
        """Returns True if edge_stmt passes the below filters"""
        # Failsafe for empty statements
        if not edge_stmt:
            logger.warning('No edge statement')
            return False

        # Filter belief score
        if edge_stmt['bs'] < options['bsco']:
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

    def _get_cost(self, path, direct=True):
        if direct:
            # Return sum of aggregated weights per edge
            return sum(self.dir_edges[(s, o)]['weight'] for s, o in
                       zip(path[:-1], path[1:]))
        else:
            # Return sum of averaged weights per stmts
            cost = 0
            for s, o in zip(path[:-1], path[1:]):
                ew = []
                e = self._get_edge(s, o, len(ew), direct)
                while e:
                    ew.append(e['weight'])
                    e = self._get_edge(s, o, len(ew), direct)
                cost += sum(ew)/len(ew)
            return cost

    def _aggregated_path_belief(self, path):
        bs_list = [self.dir_edges[e]['bs'] for e in zip(path[:-1], path[1:])]
        return ag_belief_score(bs_list)

    def _get_sort_key(self, path, hash_path, method=None):
        """Calculate a number to sort the path on

        `Method` allows to specify the calculation"""

        # Default: aggregated path belief score
        sort_key = self._aggregated_path_belief(path)
        return sort_key

    @staticmethod
    def _sort_stmts(ksp):
        for l in ksp:
            res_list = ksp[l]
            ksp[l] = sorted(res_list,
                            key=lambda pd: pd['sort_key'],
                            reverse=True)
        return ksp

    def _uri_by_node(self, node):
        """Return the fplx URI for the provided node"""
        # Check existence of node outside function
        node_id = self.nodes[node]['id']
        node_ns = self.nodes[node]['ns']
        return self.ehm.get_uri(id=node_id, ns=node_ns)

    def _get_parents(self, node):
        if self.nodes.get(node):
            id = node
            ns = self.nodes[node]['ns']

            true_ns, true_id = dnf._ns_id_from_name(id)
            if true_ns and true_id:
                return self.ehm.get_parents(uri=self.ehm.get_uri(true_ns,
                                                                 true_id))
            return self.ehm.get_parents(uri=self.ehm.get_uri(ns, id))
        else:
            return set()
