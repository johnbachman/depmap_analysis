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

INDRA_NETWORK_CACHE = path.join(CACHE, 'nx_multi_digraph_db_dump_20190417.pkl')
MAX_NUM_PATH = 10
MAX_PATH_LEN = 6


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
    def __init__(self, indra_graph):
        self.nx_graph_repr = indra_graph
        self.nodes = self.nx_graph_repr.nodes
        self.edges = self.nx_graph_repr.edges
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

            # Todo MultiDiGrap can't do simple graphs

            return self.find_shortest_paths(**options)
        except Exception as e:
            logger.warning('Exception: ', repr(e))

    def find_shortest_path(self, source, target, weight=None, simple=False):
        """Returns a list of nodes representing a shortest path"""
        result = []
        try:
            if not simple:
                path = nx.shortest_path(self.nx_graph_repr, source, target,
                                        weight)
                result.append({'stmts': self._get_hash_path(path),
                               'path': path})
                return result
            else:
                return self._find_shortest_simple_paths(source, target,
                                                        weight, 1)
        except NodeNotFound or nx.NetworkXNoPath:
            return []

    def find_shortest_paths(self, source, target, weight=False, simple=True,
                            **kwargs):
        """Returns a list of len <= self.MAX_NUM_PATH of shortest paths"""
        try:
            if not simple:
                try:
                    result = []
                    for n, path in enumerate(nx.all_shortest_paths(
                            self.nx_graph_repr, source, target, weight)):
                        if n > self.MAX_NUM_PATH:
                            logger.info('Max number of paths exceeded, '
                                        'breaking.')
                            return result
                        result.append({'stmts': self._get_hash_path(path),
                                       'path': path})
                except nx.NetworkXNoPath:
                    return []

            else:
                return self._find_shortest_simple_paths(source, target, weight)

        except nx.NodeNotFound:
            return []

        return result

    def _find_shortest_simple_paths(self, source, target, weight=None,
                                    max_paths_found=0):
        if max_paths_found == 0 or max_paths_found > self.MAX_NUM_PATH:
            max_paths_found = self.MAX_NUM_PATH
        result = []
        if source in self.nodes and target in self.nodes:
            try:
                paths = nx.shortest_simple_paths(self.nx_graph_repr, source,
                                                 target, weight)
                path_len = 0
                for c, path in enumerate(paths):
                    try:
                        if path_len == 0:
                            path_len = len(path)
                        elif path_len != 0 and len(path) > path_len:
                            break

                        if c > max_paths_found:
                            logger.info('Max number of paths exceeded, '
                                        'breaking.')
                        result.append({
                            'stmts': self._get_hash_path(path),
                            'path': path
                        })
                    except KeyError as e:
                        logger.warning(repr(e))
                        continue
            except IndexError as e:
                logger.warning(repr(e))
            except nx.NetworkXNoPath as e:
                logger.warning(repr(e))

        return result

    def has_path(self, source, target):
        return nx.has_path(self.nx_graph_repr, source, target)

    def _get_hash_path(self, path):
        hash_path = []
        for n in range(len(path) - 1):
            edges = []
            e = 0
            edge = self.edges.get((path[n], path[n + 1], e))
            while edge:
                edges.append({**edge, 'subj': path[n], 'obj': path[n + 1]})
                e += 1
                edge = self.edges.get((path[n], path[n + 1], e))
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
        logger.error('Script failed: ' + e.strerror)

    return stmts_file, dataframe_file, csv_file


def load_indra_graph(graph_path, update=False):
    if update:
        stmts_file, dataframe_file, csv_file = dump_indra_db()
        indra_graph = dnf.nx_directed_graph_from_sif_dataframe(dataframe_file)
        logging.info('Dumping latest indra db snapshot to pickle')
        _dump_it_to_pickle(graph_path, indra_graph)
    else:
        logger.info('Loading indra network...')
        indra_graph = _pickle_open(graph_path)
        logger.info('Finished loading indra network...')
    return indra_graph


if path.isfile(INDRA_NETWORK_CACHE):
    indra_network = IndraNetwork(load_indra_graph(INDRA_NETWORK_CACHE))
else:
    # Here should dump new cache instead, but raise error for now

    raise FileExistsError('Could not find file: ' + INDRA_NETWORK_CACHE)


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
    logger.info('Args -----------')
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


if __name__ == '__main__':
    app.run()
