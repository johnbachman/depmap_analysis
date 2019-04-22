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

INDRA_NETWORK_CACHE = path.join(CACHE, 'nx_dir_graph_db_dump_20190417.pkl')
MAX_PATH = 10


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
        self.MAX_PATH = MAX_PATH

    def find_shortest_path(self, source, target, weight=None, simple=True):
        """Returns a list of nodes representing a shortest path"""
        try:
            if not simple:
                return nx.shortest_path(self.nx_graph_repr,
                                        source, target, weight)
            else:
                return [] ###################
        except NodeNotFound or nx.NetworkXNoPath:
            return [] ###################

    def _find_shortest_simple_path(self, source, target, weight=None):
        return [] ###################

    def find_shortest_paths(self, source, target, weight=None, simple=True):
        """Returns a list of len <= self.MAX_PATH of shortest paths"""
        try:
            if not simple:
                return nx.all_shortest_paths(self.nx_graph_repr, source,
                                             target, weight)
            else:
                return self._find_shortest_simple_paths(source, target, weight)
        except NodeNotFound or nx.NetworkXNoPath:
            return [] ###################

    def _find_shortest_paths(self, source, target, weight=None,
                             max_len=0):
        if max_len == 0:
            max_len = self.MAX_PATH
        return [] ###################

    def _find_shortest_simple_paths(self, source, target, weight=None,
                                    max_paths_found=0):
        if max_paths_found == 0:
            max_paths_found = self.MAX_PATH
        elif max_paths_found > self.MAX_PATH:
            max_paths_found = self.MAX_PATH
        result = []
        if source in self.nodes and target in self.nodes:
            try:
                paths = nx.shortest_simple_paths(self.nx_graph_repr, source,
                                                 target, weight)
                for c, path in enumerate(paths):
                    path_len = 0
                    try:
                        if path_len == 0:
                            path_len = len(path)
                        elif path_len != 0 and len(path) > path_len:
                            break

                        if c > max_paths_found:
                            logger.info('Max number of paths exceeded, '
                                        'breaking.')
                        hash_path = []
                        for n in range(len(path)-1):
                            hashes = {}
                            for t in self.edges[(path[n],
                                                 path[n+1])]['stmt_list']:
                                hashes[t[0]] = t[2]
                            hash_path.append(hashes)
                        result.append({
                            'stmts': hash_path,
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
    logger.info(request.args)
    logger.info('Json -----------')
    logger.info(str(request.json))
    logger.info('------------------')

    try:
        source = request.json['source']
        target = request.json['target']
        logger.info('Querying indra network for path between %s and %s' %
                    (source, target))
        result = indra_network.find_shortest_path(source, target)
        if not result:
            logger.info('Query returned with no path found')
        res = {'result': result}
        logger.info('Result: %s' % str(res))
        return Response(json.dumps(res), mimetype='application/json')

    except KeyError:
        # return 400 badly formatted
        abort(Response('Badly formatted json/missing parameters', 400))


if __name__ == '__main__':
    app.run()
