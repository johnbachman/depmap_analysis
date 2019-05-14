import json
import logging
import argparse
import networkx as nx
from os import path
from jinja2 import Template
from subprocess import call
from datetime import datetime
from networkx import NodeNotFound
from flask import Flask, request, abort, Response, redirect, url_for

from indra_db.util import dump_sif

import depmap_network_functions as dnf
from util.io_functions import _pickle_open, _dump_it_to_pickle

app = Flask(__name__)

logger = logging.getLogger('INDRA GDE API')

HERE = path.dirname(path.abspath(__file__))
CACHE = path.join(HERE, '_cache')

TEST_MDG_CACHE = path.join(CACHE, 'test_mdg_network.pkl')
INDRA_MDG_CACHE = path.join(CACHE, 'nx_bs_multi_digraph_db_dump_20190417.pkl')
TEST_DG_CACHE = path.join(CACHE, 'test_dir_network.pkl')
INDRA_DG_CACHE = path.join(CACHE, 'nx_bs_dir_graph_db_dump_20190417.pkl')

MAX_PATHS = 100
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
        self.MAX_PATHS = MAX_PATHS
        self.MAX_PATH_LEN = MAX_PATH_LEN
        self.small = False
        self.verbose = 0

    def handle_query(self, **kwargs):
        """Handles path query from client. Returns query result."""
        logger.info('Handling query: %s' % repr(kwargs))
        # possible keys (* = mandatory):
        # *'source', *'target', 'path_length', 'spec_len_only', 'sign',
        # 'weighted', 'direct_only', 'curated_db_only'
        keys = kwargs.keys()
        # 'source' & 'target' are mandatory
        if 'source' not in keys or 'target' not in keys:
            raise KeyError('Missing mandatory parameters "source" or '
                           '"target"')
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
        logger.info('Lookng for no more than %d paths' % self.MAX_PATHS)

        # Todo MultiDiGrap can't do simple graphs: resolve by loading
        #  both a MultiDiGraph and a simple DiGraph - find the simple
        #  paths in the DiGraph and check them in the Multi-DiGraph.

        return self.find_shortest_paths(**options)

    def find_shortest_path(self, source, target, weight=None, simple=False,
                           **kwargs):
        """Returns a list of nodes representing a shortest path"""
        try:
            if not simple:
                path = nx.shortest_path(self.nx_dir_graph_repr, source, target,
                                        weight)
                return {len(path): [{

                    'stmts': self._get_hash_path(path, **kwargs),
                    'path': path
                }]}
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
        source_succ = set(self.nx_dir_graph_repr.succ[kwargs['source']].keys())
        target_succ = set(self.nx_dir_graph_repr.succ[kwargs['target']].keys())
        common = source_succ & target_succ

        if common:
            return self._loop_common_targets(common_targets=common, **kwargs)

        return []

    def _loop_common_targets(self, common_targets, **kwargs):
        """Order common_targets targets by lowest bs in pair."""
        ordered_commons = []
        source = kwargs['source']
        target = kwargs['target']
        for ct in common_targets:
            paths1 = self._get_hash_path(path=[source, ct], **kwargs)[0]
            paths2 = self._get_hash_path(path=[target, ct], **kwargs)[0]
            if paths1 and paths2:
                max_bs1 = max([st['bs'] for st in paths1])
                max_bs2 = max([st['bs'] for st in paths2])
                ordered_commons.append({
                    ct: [sorted(paths1, key=lambda k: k['bs'], reverse=True),
                         sorted(paths2, key=lambda k: k['bs'], reverse=True)],
                    'lowest_highest_belief': min(max_bs1, max_bs2)
                })
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
            # Check if we path exceeds MAX_PATH_LEN
            if len(path) >= self.MAX_PATH_LEN:
                if not result['paths_by_node_count']:
                    logger.info('No paths shorther than %d found.' %
                                self.MAX_PATH_LEN)
                    return {}
                logger.info('Reached longest allowed path length. Returning '
                            'results')
                return result
            hash_path = self._get_hash_path(path, **kwargs)
            if self.verbose > 1:
                logger.info('Got hash path: %s' % repr(hash_path))
            if hash_path:
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
        if self.verbose:
            logger.info('Done loopngi paths. Returning result: %s' %
                        repr(result))
        return result

    def has_path(self, source, target):
        """Return true if there is a path from source to target"""
        return nx.has_path(self.nx_dir_graph_repr, source, target)

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
        """Return a list of n-1 lists of dicts containing of stmts connected
        by the n nodes in the input path. if simple_dir is True, query edges
        from directed graph and not from MultiDiGraph representation"""
        hash_path = []
        if self.verbose:
            logger.info('Building evidence for path %s' % str(path))
        for n in range(len(path) - 1):
            edges = []
            subj = path[n]
            obj = path[n+1]
            if self.nodes[subj]['ns'] in kwargs['node_filter'] \
                    or self.nodes[obj]['ns'] in kwargs['node_filter']:
                if self.verbose:
                    logger.info('Node namespace %s or %s filtered out using '
                                '%s' % (self.nodes[subj]['ns'],
                                        self.nodes[obj]['ns'],
                                        kwargs['node_filter']))
                return []
            e = 0
            edge_stmt = self._get_edge(subj, obj, e, simple_dir)
            if self.verbose > 2:
                logger.info('edge stmt %s' % repr(edge_stmt))
            while edge_stmt:
                if self._pass_stmt(subj, obj, edge_stmt, **kwargs):
                    if self.verbose > 3:
                        logger.info('edge stmt passed filter, appending to '
                                    'edge list.')
                    # convert hash to string for javascript compatability
                    edge_stmt['stmt_hash'] = str(edge_stmt['stmt_hash'])
                    edges.append({**edge_stmt,
                                  'subj': subj,
                                  'obj': obj})
                e += 1
                edge_stmt = self._get_edge(subj, obj, e, simple_dir)
            if self.verbose:
                logger.info('Appending %s to hash path list' % repr(edges))
            hash_path.append(edges)
        if self.verbose and len(hash_path) > 0:
            logger.info('Returning hash path: %s' % repr(hash_path))
        return hash_path

    def _pass_stmt(self, subj, obj, edge_stmt, **kwargs):
        if not edge_stmt:
            if self.verbose:
                logger.info('No edge statement')
            return False
        if edge_stmt['bs'] < kwargs['bsco']:
            if self.verbose:
                logger.info('Did not pass belief score')
            return False
        if edge_stmt['stmt_type'].lower() in kwargs['stmt_filter']:
            if self.verbose:
                logger.info('statement type %s filtered out as part filter %s'
                            % (edge_stmt['stmt_type'],
                               str(kwargs['stmt_filter'])))
            return False
        return True


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
        logger.info('Result: %s' % str(res))
        return Response(json.dumps(res), mimetype='application/json')

    except KeyError:
        # return 400 badly formatted
        logger.warning('Missing parameters in query')
        abort(Response('Missing parameters', 400))

    except ValueError:
        # Bad values in json, but entry existed
        logger.warning('Badly formatted json')
        abort(Response('Badly formatted json', 400))
    except Exception as e:
        # Anything else: probably bug or networkx error, not the user's fault
        logger.warning('Unhandled internal error: ' + repr(e))
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
    if args.verbose:
        logger.info('Verbose level %d' % args.verbose)
        indra_network.verbose = args.verbose
    app.run(host=args.host, port=args.port)
