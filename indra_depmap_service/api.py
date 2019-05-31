import json
import logging
import argparse
import networkx as nx
from os import path
from jinja2 import Template
from subprocess import call
from datetime import datetime
from itertools import product
from networkx import NodeNotFound
from flask import Flask, request, abort, Response

from indra_db.util import dump_sif

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

MAX_PATHS = 50
MAX_PATH_LEN = 4


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
        self.dir_edges = self.nx_dir_graph_repr.edges
        self.mdg_edges = self.nx_md_graph_repr.edges
        self.ehm = indra_dir_graph.graph.get('entity_hierarchy_manager', None)
        self.node_by_uri = indra_dir_graph.graph.get('node_by_uri', None)
        self.MAX_PATHS = MAX_PATHS
        self.MAX_PATH_LEN = MAX_PATH_LEN
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
            a list of valid indra statement or FamPlex types *to include* in
            the path
        node_filter: [str]
            a list of node namespaces *to include* in the path
        path_length: int <=4
            a positive integer <= 4 stating the maximum number of edges in
            the path
        spec_len_only: Bool
            If True, only search for paths with number of edges given by
            path_lenth
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
        simple: Bool
            If True, do a simple path search
        k_shortest: Bool|int
            An integer stating the maximum number of directed paths to return
            in the result. The maximum allowed value is 50. If False,
            the maximum number of paths returned will be set to the maximum
            allowed value.

        Returns
        -------
        """
        logger.info('Handling query: %s' % repr(kwargs))
        mandatory = ['source', 'target', 'stmt_filter', 'node_filter',
                     'path_length', 'spec_len_only', 'weighted',
                     'bsco', 'direct_only', 'fplx_expand',
                     'simple', 'k_shortest']
        if not all([key in kwargs for key in mandatory]):
            miss = [key in kwargs for key in mandatory].index(False)
            raise KeyError('Missing mandatory parameter %s' % mandatory[miss])
        options = {k: v for k, v in kwargs.items()
                   if k not in ['path_length', 'sign']}  # Handled below
        for k, v in kwargs.items():
            if k == 'weighted':
                logger.info('Doing weighted path search') if v \
                    else logger.info('Doing unweighted path search')
            if k == 'path_length':
                options['path_length'] = -1 if v == 'no_limit' else int(v)
            if k == 'sign':
                options['sign'] = 1 if v == 'plus' \
                    else (-1 if v == 'minus' else 0)
        k_shortest = kwargs.pop('k_shortest', None)
        self.MAX_PATHS = k_shortest if k_shortest else MAX_PATHS
        logger.info('Looking for no more than %d paths' % self.MAX_PATHS)

        # Todo MultiDiGrap can't do simple graphs: resolve by loading
        #  both a MultiDiGraph and a simple DiGraph - find the simple
        #  paths in the DiGraph and check them in the Multi-DiGraph.
        ksp = self.find_shortest_paths(**options)
        if not ksp:
            if kwargs['fplx_expand']:
                logger.info('No directed path found, looking for paths '
                            'connected by common parents of source and/or '
                            'target')
                ckwargs = options.copy()
                ksp = self.try_parents(**ckwargs)
                if self.verbose > 2:
                    logger.info('Got parents search result: %s' % repr(ksp))
            else:
                logger.info('No directed path found')
        ct = self.find_common_targets(**options)
        cp = self.get_common_parents(**options)
        return {**ksp, 'common_targets': ct, 'common_parents': cp}

    def try_parents(self, **ckwargs):
        """Retry search with sources' and targets' parents

        Search for paths between combinations of the parents of source and
        target.
        """
        source = ckwargs['source']
        target = ckwargs['target']

        if self.verbose > 1:
            logger.info('Parents search: source=%s, target=%s' % \
                        (ckwargs['source'], ckwargs['target']))

        # Get closures for source and target
        source_parent_closure = self._get_closure(source)
        target_parent_closure = self._get_closure(target)
        if self.verbose > 3:
            logger.info('Got source_parent_closure: %s' %
                        repr(source_parent_closure))
            logger.info('Got target_parent_closure: %s' %
                        repr(target_parent_closure))

        # Base case: no further closures found, return empty dict
        if not source_parent_closure and not target_parent_closure:
            return {}

        # First try current source with all target parents
        if target_parent_closure:
            for tp_uri in target_parent_closure:
                ckwargs['target'] = self.node_by_uri[tp_uri]
                if self.verbose > 4:
                    logger.info('Parents search: source=%s, target=%s' % \
                                (ckwargs['source'], ckwargs['target']))
                ksp = self.find_shortest_paths(**ckwargs)
                if ksp:
                    return ksp

        # Then, try current target with all source parents
        if source_parent_closure:
            for sp_uri in source_parent_closure:
                ckwargs['source'] = self.node_by_uri[sp_uri]
                if self.verbose > 4:
                    logger.info('Parents search: source=%s, target=%s' % \
                                (ckwargs['source'], ckwargs['target']))
                ksp = self.find_shortest_paths(**ckwargs)
                if ksp:
                    return ksp

        # Lastly try all possible pairs of source and target parents
        if source_parent_closure and target_parent_closure:
            for sp_uri, tp_uri in product(source_parent_closure,
                                          target_parent_closure):
                    ckwargs['source'] = self.node_by_uri[sp_uri]
                    ckwargs['target'] = self.node_by_uri[tp_uri]
                    if self.verbose > 4:
                        logger.info('Parents search: source=%s, target=%s' % \
                                    (ckwargs['source'], ckwargs['target']))
                    ksp = self.find_shortest_paths(**ckwargs)
                    if ksp:
                        return ksp

        # If we get this far, no path was found
        return {}

    def find_shortest_path(self, source, target, weight=None, simple=False,
                           **kwargs):
        """Returns a list of nodes representing a shortest path"""
        try:
            if not simple:
                path = nx.shortest_path(self.nx_dir_graph_repr, source, target,
                                        weight)
                return {len(path): [{'path': path,
                    'stmts': self._get_hash_path(path, **kwargs)}]}
            else:
                return self._find_shortest_simple_paths(source, target,
                                                        weight, **kwargs)
        except NodeNotFound or nx.NetworkXNoPath:
            return {}

    def find_shortest_paths(self, source, target, weight=None, simple=True,
                            **kwargs):
        """Returns a list of shortest paths in ascending order"""
        path_len = kwargs['path_length']
        if path_len > self.MAX_PATH_LEN:
            logger.warning('path_len > MAX_PATH_LEN, resetting path_len to '
                           'MAX_PATH_LEN (%d).' % self.MAX_PATH_LEN)
            path_len = self.MAX_PATH_LEN
        try:
            if not simple:
                logger.info('Doing non-simple path search')
                # paths = nx.all_shortest_paths(self.nx_dir_graph_repr,
                #                               source, target, weight)
                paths = nx.all_shortest_paths(self.nx_md_graph_repr,
                                              source, target, weight)
                return self._loop_paths(paths, path_len, **kwargs)
            else:
                logger.info('Doing simple path search')
                return self._find_shortest_simple_paths(source, target,
                                                        weight, **kwargs)
        except nx.NodeNotFound as e:
            logger.warning(repr(e))
            return {}
        except nx.NetworkXNoPath as e:
            logger.warning(repr(e))
            return {}

    def _find_shortest_simple_paths(self, source, target, weight=None,
                                    **kwargs):
        """Returns a list of shortest simple paths in ascending order"""
        path_len = kwargs['path_length']
        if path_len > self.MAX_PATH_LEN:
            logger.warning('path_len > MAX_PATH_LEN, resetting path_len to '
                           'MAX_PATH_LEN (%d).' % self.MAX_PATH_LEN)
            path_len = self.MAX_PATH_LEN
        simple_paths = nx.shortest_simple_paths(self.nx_dir_graph_repr,
                                                source, target, weight)
        return self._loop_paths(simple_paths, path_len, **kwargs)

    def find_common_targets(self,**kwargs):
        """Returns a list of statement(?) pairs that explain common targets
        for source and target"""
        if kwargs['source'] in self.nodes and kwargs['target'] in self.nodes:
            source_succ = set(self.nx_dir_graph_repr.succ[
                                  kwargs['source']].keys())
            target_succ = set(self.nx_dir_graph_repr.succ[
                                  kwargs['target']].keys())
            common = source_succ & target_succ
            if common:
                try:
                    return self._loop_common_targets(common_targets=common,
                                                     **kwargs)
                except nx.NodeNotFound as e:
                    logger.warning(repr(e))
                except nx.NetworkXNoPath as e:
                    logger.warning(repr(e))

        return []

    def _loop_common_targets(self, common_targets, **kwargs):
        """Order common_targets targets by lowest bs in pair."""
        ordered_commons = []
        source = kwargs['source']
        target = kwargs['target']
        added_targets = 0
        for ct in common_targets:
            paths1 = self._get_hash_path(path=[source, ct], **kwargs)
            paths2 = self._get_hash_path(path=[target, ct], **kwargs)
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

    def _loop_paths(self, paths_gen, path_len, **kwargs):
        len_only = kwargs['spec_len_only']
        result = {'paths_by_node_count': {}}
        added_paths = 0
        for path in paths_gen:
            # Check if we found k paths
            if added_paths >= self.MAX_PATHS:
                logger.info('Found k shortest paths')
                return result
            # Check if path exceeds MAX_PATH_LEN
            if len(path) >= self.MAX_PATH_LEN:
                if not result['paths_by_node_count']:
                    logger.info('No paths shorther than %d found.' %
                                self.MAX_PATH_LEN)
                    return {}
                logger.info('Reached longest allowed path length. Returning '
                            'results')
                return result
            hash_path = self._get_hash_path(path, **kwargs)
            if hash_path and all(hash_path):
                if self.verbose > 1:
                    logger.info('Adding stmts and path from %s to path list' %
                                repr(hash_path))
                pd = {'stmts': hash_path, 'path': path}
                try:
                    if not len_only:
                        result['paths_by_node_count'][len(path)].append(pd)
                        added_paths += 1
                    elif len_only and len(path) == path_len:
                        result['paths_by_node_count'][len(path)].append(pd)
                        added_paths += 1
                    elif len_only and len(path) != path_len:
                        continue
                    else:
                        logger.warning('This option should not happen')
                except KeyError:
                    try:
                        if len_only and len(path) == path_len:
                            result['paths_by_node_count'][len(path)] = [pd]
                            added_paths += 1
                        elif not len_only:
                            result['paths_by_node_count'][len(path)] = [pd]
                            added_paths += 1
                    except KeyError as ke:
                        logger.warning('Unexpected KeyError: ' + repr(ke))
                        raise ke
        if self.verbose > 2:
            logger.info('Done looping paths. Returning result: %s' %
                        repr(result))
        return result

    def has_path(self, source, target):
        """Return true if there is a path from source to target"""
        return nx.has_path(self.nx_dir_graph_repr, source, target)

    def get_common_parents(self, **kwargs):
        """Find common parents between source and target"""
        source_id = kwargs['source']
        source_ns = None
        target_id = kwargs['target']
        target_ns = None
        cp = {}

        # Get ns
        if kwargs['source'] in self.nodes:
            source_ns = self.nodes[kwargs['source']]['ns']
        if kwargs['target'] in self.nodes:
            target_ns = self.nodes[kwargs['target']]['ns']

        # Try different combinations of ns combinations

        # If both source and target are given
        if source_ns and target_ns:
            if source_ns in kwargs['node_filter'] and \
                    target_ns in kwargs['node_filter']:
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
            if target_ns in kwargs['node_filter']:
                if self.verbose > 1:
                    logger.info('No namespace found for %s, trying HGNC and '
                                'FPLX.' % source_id)
                for sns in ['HGNC', 'FPLX']:
                    if sns not in kwargs['node_filter']:
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
            if source_ns in kwargs['node_filter']:
                if self.verbose > 1:
                    logger.info('No namespace found for %s, trying HGNC and '
                                'FPLX.' % target_id)
                for tns in ['HGNC', 'FPLX']:
                    if tns not in kwargs['node_filter']:
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
                if source_ns not in kwargs['node_filter']:
                    continue
                for target_ns in ['HGNC', 'FPLX']:
                    if target_ns not in kwargs['node_filter']:
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

    def _get_edge(self, s, o, index, directed):
        """Return edges from DiGraph or MultiDigraph in a uniform format"""
        if directed:
            try:
                stmt_edge = self.dir_edges.get((s, o))['stmt_list'][index]
            except IndexError:
                # To keep it consistent with below Multi DiGraph implementation
                stmt_edge = None
            return stmt_edge
        else:
            return self.mdg_edges.get((s, o, index))

    def _get_hash_path(self, path, simple_dir=True, **kwargs):
        """Return a list of n-1 lists of dicts containing of stmts connecting
        the n nodes in path. If simple_dir is True, query edges from DiGraph
        and not from MultiDiGraph representation"""
        hash_path = []
        if self.verbose:
            logger.info('Building evidence for path %s' % str(path))
        for subj, obj in zip(path[:-1], path[1:]):
            # Check node filter
            if self.nodes[subj]['ns'] not in kwargs['node_filter'] \
                    or self.nodes[obj]['ns'] not in kwargs['node_filter']:
                if self.verbose:
                    logger.info('Node namespace %s or %s not part of '
                                'acceptable namespaces %s' %
                                (self.nodes[subj]['ns'],
                                 self.nodes[obj]['ns'],
                                 kwargs['node_filter']))
                return []

            # Initialize edges list, statement index
            edges = []
            e = 0

            # Get first edge statement
            edge_stmt = self._get_edge(subj, obj, e, simple_dir)
            if self.verbose > 3:
                logger.info('First edge stmt %s' % repr(edge_stmt))

            # Exhaustively loop through all edge statments
            while edge_stmt:

                # If edge statement passes, append to edges list
                if self._pass_stmt(subj, obj, edge_stmt, **kwargs):
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
                edge_stmt = self._get_edge(subj, obj, e, simple_dir)

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

    def _pass_stmt(self, subj, obj, edge_stmt, **kwargs):
        """Returns True if edge_stmt passes the below filters"""
        # Failsafe for empty statements are sent
        if not edge_stmt:
            if self.verbose:
                logger.info('No edge statement')
            return False

        # Filter belief score
        if edge_stmt['bs'] < kwargs['bsco']:
            if self.verbose:
                logger.info('Did not pass belief score')
            return False

        # Filter statement type
        if edge_stmt['stmt_type'].lower() not in kwargs['stmt_filter']:
            if self.verbose > 4:
                logger.info('statement type %s not found in filter %s'
                            % (edge_stmt['stmt_type'],
                               str(kwargs['stmt_filter'])))
            return False

        # Return True is all filters were passed
        return True

    def _uri_by_node(self, node):
        """Check existence of node outside function"""
        node_id = self.nodes[node]['id']
        node_ns = self.nodes[node]['ns']
        return self.ehm.get_uri(id=node_id, ns=node_ns)

    def _get_closure(self, node):
        if self.nodes.get(node):
            return set(self.ehm.isa_or_partof_closure.get(self._uri_by_node(
                node), []))
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
    parser.add_argument('-v', '--verbose', action='count', default=0)
    args = parser.parse_args()

    if args.test:
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
