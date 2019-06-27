import json
import logging
import argparse
import requests
import numpy as np
import networkx as nx
from os import path
from jinja2 import Template
from subprocess import call
from datetime import datetime
from itertools import product
from networkx import NodeNotFound
from flask import Flask, request, abort, Response

from indra_db.util import dump_sif
from indra.config import CONFIG_DICT

import depmap_network_functions as dnf
from util.io_functions import _pickle_open, _dump_it_to_pickle

app = Flask(__name__)

logger = logging.getLogger('INDRA GDE API')

HERE = path.dirname(path.abspath(__file__))
CACHE = path.join(HERE, '_cache')

TEST_MDG_CACHE = path.join(CACHE, 'test_mdg_network.pkl')
INDRA_MDG_CACHE = path.join(CACHE,
                            'nx_bs_fam_multi_digraph_db_dump_20190417.pkl')
TEST_DG_CACHE = path.join(CACHE, 'test_dir_network.pkl')
INDRA_DG_CACHE = path.join(CACHE, 'nx_bs_fam_dir_graph_db_dump_20190417.pkl')

GRND_URI = None
try:
    GRND_URI = CONFIG_DICT['INDRA_GROUNDING_SERVICE_URL']
except KeyError:
    logger.warning('Indra Grounding service not available. Add '
                   'INDRA_GROUNDING_SERVICE_URL to `indra/config.ini`')

MAX_PATHS = 50
MAX_SKIP = 1000


def _todays_date():
    return datetime.now().strftime('%Y%m%d')


# Create a template object from the template file, load once
def _load_template(fname):
    template_path = path.join(HERE, fname)
    with open(template_path, 'rt') as f:
        template_str = f.read()
        template = Template(template_str)
    return template


QUERY = _load_template('query.html')


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
            a list of valid indra statement types or FamPlex (as 'isa' or
            'part_of') child-parent connections *to exclude* in the path
        node_filter: [str]
            a list of node namespaces *to include* in the path
        node_blacklist: [str]
            a list of node names to ignore. If a path contains a node in this
            list the path will be discarded.
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
        curated_db_only: Bool **currently not implemented**
            Placeholder for future implementation allowing filtering to edges
            that are sourced from curated databases only
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
        """
        mandatory = ['source', 'target', 'stmt_filter', 'node_filter',
                     'path_length', 'weighted', 'bsco', 'fplx_expand',
                     'k_shortest']
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
        if not ksp and GRND_URI:
            ckwargs = options.copy()
            ksp = self.grounding_fallback(**ckwargs)
        if not ksp and kwargs['fplx_expand']:
            logger.info('No directed path found, looking for paths '
                        'connected by common parents of source and/or '
                        'target')
            ksp = self.try_parents(**ckwargs)
            if self.verbose > 2:
                logger.info('Got parents search result: %s' % repr(ksp))
        elif not ksp:
            logger.info('No directed path found')
        if ksp and not options['weight']:
            # Sort the results in ksp if non-weighted search
            ksp = self._sort_stmts(ksp)
        ct = self.find_common_targets(**options)
        cp = list(self.get_common_parents(**options))
        return {'paths_by_node_count': ksp,
                'common_targets': ct,
                'common_parents': cp}

    def grounding_fallback(self, **ckwargs):
        """Retry search with alternative names found by grounding service"""
        if self.verbose:
            logger.info('Expanding search using grounding service')
        org_source = ckwargs['source']
        org_target = ckwargs['target']

        # Get groundings
        src_groundings = requests.post(GRND_URI,
                                       json={'text': org_source}).json()
        trgt_groundings = requests.post(GRND_URI,
                                        json={'text': org_target}).json()

        # Loop combinations of source and target groundings, break if
        # anything found

        # org target with sources
        if src_groundings and not trgt_groundings:
            for src in src_groundings:
                ckwargs['source'] = src['entry']['entry_name']
                ksp = self.find_shortest_paths(**ckwargs)
                if ksp:
                    return ksp

        # org source with targets
        if not src_groundings and trgt_groundings:
            ckwargs['source'] = org_source
            for trgt in trgt_groundings:
                ckwargs['target'] = trgt['entry']['entry_name']
                ksp = self.find_shortest_paths(**ckwargs)
                if ksp:
                    return ksp

        # all source groundings with all target groundings
        if src_groundings and trgt_groundings:
            for src, trgt in product(src_groundings, trgt_groundings):
                ckwargs['source'] = src['entry']['entry_name']
                ckwargs['target'] = trgt['entry']['entry_name']
                ksp = self.find_shortest_paths(**ckwargs)
                if ksp:
                    return ksp

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

    def find_shortest_path(self, source, target, weight, **options):
        """Returns a list of nodes representing a shortest path"""
        try:
            return self._loop_paths(nx.shortest_path(self.nx_dir_graph_repr,
                                                     source, target, weight),
                                    **options)
        except NodeNotFound or nx.NetworkXNoPath:
            return {}

    def find_shortest_paths(self, source, target, weight, **options):
        """Returns a list of shortest paths in ascending order"""
        try:
            logger.info('Doing simple %s path search' % 'weigthed' if weight
                        else '')
            paths = nx.shortest_simple_paths(self.nx_dir_graph_repr,
                                     source, target, weight)
            # paths = nx.all_shortest_paths(self.nx_md_graph_repr,
            #                               source, target, weight)
            return self._loop_paths(paths, **options)
        except nx.NodeNotFound as e:
            logger.warning(repr(e))
            return {}
        except nx.NetworkXNoPath as e:
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
                except nx.NodeNotFound as e:
                    logger.warning(repr(e))
                except nx.NetworkXNoPath as e:
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
        path_len = options['path_length'] + 1 if options['path_length'] else \
            False
        result = {}
        added_paths = 0
        skipped_paths = 0
        for path in paths_gen:
            # Check if we found k paths
            if added_paths >= self.MAX_PATHS:
                logger.info('Found all %d shortest paths, returning results.' %
                            self.MAX_PATHS)
                return result
            if skipped_paths > MAX_SKIP:
                logger.warning('Reached MAX_SKIP (%d) before finding all %d '
                               'shortest paths. Returning search.' %
                               (MAX_SKIP, MAX_PATHS))
                return result

            hash_path = self._get_hash_path(path, **options)
            if hash_path and all(hash_path):
                if self.verbose > 1:
                    logger.info('Adding stmts and path from %s to path list' %
                                repr(hash_path))
                pd = {'stmts': hash_path, 'path': path,
                      'cost': str(self._get_cost(path))}
                try:
                    if not path_len:
                        result[len(path)].append(pd)
                        added_paths += 1
                    elif path_len and len(path) < path_len:
                        continue
                    elif path_len and len(path) == path_len:
                        result[len(path)].append(pd)
                        added_paths += 1
                    elif path_len and len(path) > path_len:
                        if self.verbose > 1:
                            logger.info('Max path length reached, returning '
                                        'results.')
                        return result
                    else:
                        logger.warning('This option should not happen')
                except KeyError:
                    try:
                        if path_len and len(path) == path_len:
                            result[len(path)] = [pd]
                            added_paths += 1
                        elif not path_len:
                            result[len(path)] = [pd]
                            added_paths += 1
                    except KeyError as ke:
                        logger.warning('Unexpected KeyError: ' + repr(ke))
                        raise ke
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
        source_id = options['source']
        source_ns = None
        target_id = options['target']
        target_ns = None
        cp = {}

        # Get ns
        if options['source'] in self.nodes:
            source_ns = self.nodes[options['source']]['ns']
        if options['target'] in self.nodes:
            target_ns = self.nodes[options['target']]['ns']

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
                return {}

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
                            break
            else:
                logger.info('The namespaces for %s is not in node filter. '
                            'Aborting common parent search.' % target_id)
                return {}

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
                            break
            else:
                logger.info('The namespaces for %s is not in node filter. '
                            'Aborting common parent search.' % source_id)
                return {}

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
            return {}
        else:
            return {'source_ns': source_ns, 'source_id': source_id,
                    'target_ns': target_ns, 'target_id': target_id,
                    'common_parents': sorted(list(cp))}

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
            elif options.get('node_blacklist', None):
                if subj in options['node_blacklist'] or obj in \
                        options['node_blacklist']:
                    if self.verbose:
                        logger.info('%s or %s part of node blacklist, '
                                    'skipping path' % (subj, obj))
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
        # Failsafe for empty statements are sent
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

    @staticmethod
    def _sort_stmts(ksp):
        for l in ksp:
            res_list = ksp[l]
            ksp[l] = sorted(res_list,
                            key=lambda pd: np.longfloat(pd['cost']),
                            reverse=False)
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
            return self.ehm.get_parents(uri=self.ehm.get_uri(ns, id))
        else:
            return set()


def dump_indra_db(path='.'):
    base_name = 'db_dump_' + _todays_date()
    if path is not '.':
        path_base_name = path + base_name
    else:
        path_base_name = base_name
    stmts_file = path_base_name + '.pkl'
    dataframe_file = path_base_name + '_dataframe.pkl'
    csv_file = path_base_name + '.csv'
    files = ' '.join((stmts_file, dataframe_file, csv_file))
    cmd = 'python ' + dump_sif.__file__ + ' ' + files
    logger.info('Executing subprocess: %s' % cmd)
    try:
        retcode = call(cmd, shell=True)
        if retcode < 0:
            logger.warning('Script was terminated by signal: ' + str(-retcode))
        else:
            logger.info('Script finished ' + str(retcode))
    except OSError as e:
        logger.error('Script failed: ' + repr(e))

    return stmts_file, dataframe_file, csv_file


def load_indra_graph(dir_graph_path, multi_digraph_path, update=False):
    global INDRA_DG_CACHE, INDRA_MDG_CACHE
    if update:
        stmts_file, dataframe_file, csv_file = dump_indra_db()
        indra_dir_graph = dnf.nx_digraph_from_sif_dataframe(dataframe_file)
        indra_multi_digraph = dnf.nx_digraph_from_sif_dataframe(dataframe_file,
                                                                multi=True)
        logging.info('Dumping latest indra db snapshot to pickle')
        _dump_it_to_pickle(dir_graph_path, indra_dir_graph)
        INDRA_DG_CACHE = path.join(CACHE, dir_graph_path)
        _dump_it_to_pickle(multi_digraph_path, indra_multi_digraph)
        INDRA_MDG_CACHE = path.join(CACHE, multi_digraph_path)
    else:
        logger.info('Loading indra networks %s and %s' %
                    (dir_graph_path, multi_digraph_path))
        indra_dir_graph = _pickle_open(dir_graph_path)
        indra_multi_digraph = _pickle_open(multi_digraph_path)
        logger.info('Finished loading indra networks.')
    return indra_dir_graph, indra_multi_digraph


if path.isfile(INDRA_DG_CACHE) and path.isfile(
        INDRA_MDG_CACHE):
    indra_network = IndraNetwork()
else:
    # Here should dump new cache instead, but raise error for now

    raise FileExistsError(
        'Could not find one or both of %s and %s' %
        (INDRA_DG_CACHE, INDRA_MDG_CACHE))


@app.route('/')
@app.route('/query')
def get_query_page():
    """Loads the query page"""
    return QUERY.render()


@app.route('/query/submit', methods=['POST'])
def process_query():
    """Processing queries to the indra network"""
    # Print inputs.
    logger.info('Got query')
    logger.info('Incoming Args -----------')
    logger.info(repr(request.args))
    logger.info('Incoming Json ----------------------')
    logger.info(str(request.json))
    logger.info('------------------------------------')

    try:
        result = indra_network.handle_query(**request.json.copy())
        if not result:
            logger.info('Query returned with no path found')
        res = {'result': result}
        if indra_network.verbose > 5:
            logger.info('Result: %s' % str(res))
        return Response(json.dumps(res), mimetype='application/json')

    except KeyError as e:
        # return 400 badly formatted
        if indra_network.verbose:
            logger.exception(e)
        else:
            logger.warning('Missing parameters in query')
        abort(Response('Missing parameters', 400))

    except ValueError as e:
        # Bad values in json, but entry existed
        if indra_network.verbose:
            logger.exception(e)
        else:
            logger.warning('Badly formatted json')
        abort(Response('Badly formatted json', 400))

    except Exception as e:
        logger.exception(e)
        logger.warning('Unhandled internal error, see above error messages')
        abort(Response('Server error during handling of query', 500))


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Run the Indra DepMap Service.')
    parser.add_argument('--host', default='127.0.0.1')
    parser.add_argument('--port', default=5000, type=int)
    parser.add_argument('--test', action='store_true')
    parser.add_argument('--cache', nargs=2, help='Provide files for the '
        'networks instead of reading the default parameters. '
        'Usage: --cache <DiGraph pickle> <MultiDiGraph pickle>')
    parser.add_argument('-v', '--verbose', action='count', default=0)
    args = parser.parse_args()

    if args.cache:
        logger.info('Loading provided network files')
        INDRA_DG_CACHE = args.cache[0]
        INDRA_MDG_CACHE = args.cache[1]
    elif args.test:
        logger.info('Running test network')
        INDRA_DG_CACHE = TEST_DG_CACHE
        INDRA_MDG_CACHE = TEST_MDG_CACHE

    indra_network = \
        IndraNetwork(*load_indra_graph(INDRA_DG_CACHE,
                                       INDRA_MDG_CACHE))
    if args.test:
        indra_network.small = True
        indra_network.verbose = args.verbose if args.verbose else 1
    if args.verbose:
        logger.info('Verbose level %d' % args.verbose)
        indra_network.verbose = args.verbose
    app.run(host=args.host, port=args.port)
