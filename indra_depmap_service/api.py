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
from networkx import NodeNotFound, NetworkXNoPath
from time import time, gmtime, strftime
from flask import Flask, request, abort, Response

from indra_db.util import dump_sif
from indra.config import CONFIG_DICT

import depmap_network_functions as dnf
from depmap_analysis.util.io_functions import _pickle_open, _dump_it_to_pickle
from depmap_analysis.network_functions.indra_network import IndraNetwork

app = Flask(__name__)

logger = logging.getLogger('INDRA GDE API')

HERE = path.dirname(path.abspath(__file__))
CACHE = path.join(HERE, '_cache')

TEST_MDG_CACHE = path.join(CACHE, 'test_mdg_network.pkl')
INDRA_MDG_CACHE = path.join(CACHE,
                            'nx_bs_fam_multi_digraph_db_refresh_20190702.pkl')
TEST_DG_CACHE = path.join(CACHE, 'test_dir_network.pkl')
INDRA_DG_CACHE = path.join(CACHE, 'nx_bs_fam_dir_graph_db_refresh_20190702.pkl')

GRND_URI = None
try:
    GRND_URI = CONFIG_DICT['INDRA_GROUNDING_SERVICE_URL']
except KeyError:
    logger.warning('Indra Grounding service not available. Add '
                   'INDRA_GROUNDING_SERVICE_URL to `indra/config.ini`')

MAX_PATHS = 50
TIMEOUT = 30  # Timeout in seconds


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
        if 'test' in request.json and request.json.get('test', False):
            return Response(json.dumps({'paths_by_node_count': {}, 
                                        'common_targets': [],
                                        'common_parents': [] 
                                        }), mimetype='application/json')
        result = indra_network.handle_query(**request.json.copy())
        logger.info('Query resolved at %s' %
                    strftime('%Y-%m-%d %H:%M:%S (UTC)', gmtime(time())))
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
