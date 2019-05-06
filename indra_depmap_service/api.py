import json
import logging
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

INDRA_MDG_NETWORK_CACHE = path.join(CACHE,
                                    'nx_bs_multi_digraph_db_dump_20190417.pkl')
INDRA_DG_NETWORK_CACHE = path.join(CACHE,
                                   'nx_bs_dir_graph_db_dump_20190417.pkl')
MAX_NUM_PATH = 10
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
    def __init__(self, indra_dir_graph, indra_multi_dir_graph):
        self.nx_dir_graph_repr = indra_dir_graph
        self.nx_md_graph_repr = indra_multi_dir_graph
        self.nodes = self.nx_dir_graph_repr.nodes
        self.dir_edges = self.nx_dir_graph_repr.edges
        self.mdg_edges = self.nx_md_graph_repr.edges
        self.MAX_NUM_PATH = MAX_NUM_PATH
        self.MAX_PATH_LEN = MAX_PATH_LEN

    def handle_query(self, **kwargs):
        """Handles path query from client. Returns query result."""
        logger.info('Handling query: %s' % repr(kwargs))
        try:
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

            # Todo MultiDiGrap can't do simple graphs: resolve by loading
            #  both a MultiDiGraph and a simple DiGraph - find the simple
            #  paths in the DiGraph and check them in the Multi-DiGraph.

            return self.find_shortest_paths(**options)
        except Exception as e:
            logger.warning('Exception: ', repr(e))

    def find_shortest_path(self, source, target, weight=None, simple=False,
                           **kwargs):
        """Returns a list of nodes representing a shortest path"""
        try:
            if not simple:
                path = nx.shortest_path(self.nx_dir_graph_repr, source, target,
                                        weight)
                return {len(path): [{
                    'stmts': self._get_hash_path(path),
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

    def _loop_paths(self, paths_gen, path_len, **kwargs):
        len_only = kwargs['spec_len_only']
        belief_cutoff = kwargs['bsco']
        result = {'paths_by_node_count': {}}
        for n, path in enumerate(paths_gen):
            hash_path = self._get_hash_path(path, belief_cutoff)
            if not hash_path:
                return {}
            pd = {'stmts': hash_path, 'path': path}
            try:
                if not len_only and \
                        len(result['paths_by_node_count'][len(path)]) \
                        < self.MAX_NUM_PATH:
                    result['paths_by_node_count'][len(path)].append(pd)
                elif len_only and \
                        len(result['paths_by_node_count'][len(path)]) \
                        < self.MAX_NUM_PATH \
                        and len(path) == path_len:
                    result['paths_by_node_count'][len(path)].append(pd)
                elif len(result['paths_by_node_count'][len(path)]) \
                        >= self.MAX_NUM_PATH:
                    continue
            except KeyError:
                try:
                    if len_only and len(path) == path_len:
                        result['paths_by_node_count'][len(path)] = [pd]
                    elif not len_only:
                        result['paths_by_node_count'][len(path)] = [pd]
                except KeyError as ke:
                    logger.warning('Unexpected KeyError: ' + repr(ke))
                    raise ke
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
                stmt_edge = None  # To keep it consistent with Multi DiGraph
            return stmt_edge
        else:
            return self.mdg_edges.get((s, o, index))

    def _get_hash_path(self, path, belief_cutoff, simple_dir=True):
        """Return a list of n-1 lists of dicts containing of stmts connected
        by the n nodes in the input path. if simple_dir is True, query edges
        from directed graph and not from MultiDiGraph representation"""

        hash_path = []
        for n in range(len(path) - 1):
            edges = []
            e = 0
            edge_stmt = self._get_edge(path[n], path[n + 1], e, simple_dir)
            if edge_stmt['bs'] < belief_cutoff:
                return []
            while edge_stmt:
                # convert hash to string for javascript compatability
                edge_stmt['stmt_hash'] = str(edge_stmt['stmt_hash'])
                edges.append({**edge_stmt, 'subj': path[n], 'obj': path[n + 1]})
                e += 1
                edge_stmt = self._get_edge(path[n], path[n + 1], e, simple_dir)
                if edge_stmt and edge_stmt['bs'] < belief_cutoff:
                    return []
            hash_path.append(edges)
        return hash_path


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
    global INDRA_DG_NETWORK_CACHE, INDRA_MDG_NETWORK_CACHE
    if update:
        stmts_file, dataframe_file, csv_file = dump_indra_db()
        indra_dir_graph = dnf.nx_digraph_from_sif_dataframe(dataframe_file)
        indra_multi_digraph = dnf.nx_digraph_from_sif_dataframe(dataframe_file,
                                                                multi=True)
        logging.info('Dumping latest indra db snapshot to pickle')
        _dump_it_to_pickle(dir_graph_path, indra_dir_graph)
        INDRA_DG_NETWORK_CACHE = path.join(CACHE, dir_graph_path)
        _dump_it_to_pickle(multi_digraph_path, indra_multi_digraph)
        INDRA_MDG_NETWORK_CACHE = path.join(CACHE, multi_digraph_path)
    else:
        logger.info('Loading indra networks %s and %s' %
                    (dir_graph_path, multi_digraph_path))
        indra_dir_graph = _pickle_open(dir_graph_path)
        indra_multi_digraph = _pickle_open(multi_digraph_path)
        logger.info('Finished loading indra networks.')
    return indra_dir_graph, indra_multi_digraph


if path.isfile(INDRA_DG_NETWORK_CACHE) and path.isfile(
        INDRA_MDG_NETWORK_CACHE):
    indra_network = \
        IndraNetwork(*load_indra_graph(INDRA_DG_NETWORK_CACHE,
                                       INDRA_MDG_NETWORK_CACHE))
else:
    # Here should dump new cache instead, but raise error for now

    raise FileExistsError(
        'Could not find one or both of %s and %s' %
        (INDRA_DG_NETWORK_CACHE, INDRA_MDG_NETWORK_CACHE))


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
        abort(Response('Missing parameters', 400))

    except ValueError:
        # Bad values in json, but entry existed
        abort(Response('Badly formatted json', 400))
    finally:
        # Anything else: bug or networkx error, not the user's fault
        abort(Response('Error handling query', 500))


if __name__ == '__main__':
    app.run()
