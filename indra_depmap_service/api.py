"""INDRA Causal Network Search API"""
import os
import json
import pickle
import logging
import argparse
from sys import argv
from os import path, makedirs
from time import time, gmtime, strftime

from jinja2 import Template
from flask import Flask, request, abort, Response, render_template

from indra.config import CONFIG_DICT
from depmap_analysis.network_functions.indra_network import IndraNetwork
from depmap_analysis.util.io_functions import pickle_open, dump_it_to_pickle

from .util import load_indra_graph, _get_s3_client,\
    get_queryable_stmt_types, API_PATH as HERE, CACHE, INDRA_DG, \
    INDRA_DG_CACHE, INDRA_SEG, INDRA_SEG_CACHE, INDRA_SNG_CACHE, TEST_DG_CACHE

app = Flask(__name__)
app.config['DEBUG'] = True

logger = logging.getLogger('INDRA Network Search API')

S3_BUCKET = 'depmap-analysis'


FILES = {
    'dir_graph_path': INDRA_DG_CACHE if path.isfile(INDRA_DG_CACHE)
    else None,
    # 'multi_digraph_path': INDRA_MDG_CACHE if path.isfile(INDRA_MDG_CACHE)
    # else None,
    'multi_digraph_path': None,
    'sign_edge_graph_path': INDRA_SEG_CACHE if path.isfile(INDRA_SEG_CACHE)
    else None,
    'sign_node_graph_path': INDRA_SNG_CACHE if path.isfile(INDRA_SNG_CACHE)
    else None
}

VERBOSITY = int(os.environ.get('VERBOSITY', 0))
API_DEBUG = int(os.environ.get('API_DEBUG', 0))
if API_DEBUG:
    logger.info('API_DEBUG set to %d' % API_DEBUG)

GRND_URI = None
try:
    GRND_URI = CONFIG_DICT['GILDA_URL']
except KeyError:
    logger.warning('Indra Grounding service not available. Add '
                   'GILDA_URL to `indra/config.ini`')

MAX_PATHS = 50
TIMEOUT = 30  # Timeout in seconds


# Create a template object from the template file, load once
def _load_template(fname):
    template_path = path.join(HERE, fname)
    with open(template_path, 'rt') as f:
        template_str = f.read()
        template = Template(template_str)
    return template


# Load network
indra_network = IndraNetwork()
if path.isfile(INDRA_DG_CACHE):
    if API_DEBUG:
        logger.info('Debugging API, no network will be loaded...')
    elif argv[0].split('/')[-1].lower() != 'api.py':
        indra_network = IndraNetwork(*load_indra_graph(**FILES))
else:
    # Try to find file(s) on s3
    try:
        logger.info('%s not found locally, trying to get file from s3...' %
                    INDRA_DG)
        makedirs(CACHE, exist_ok=True)
        s3 = _get_s3_client(unsigned=True)
        logger.info('Caching network to %s' % CACHE)
        dg_key = 'indra_db_files/' + INDRA_DG
        dg_obj = s3.get_object(Bucket=S3_BUCKET, Key=dg_key)
        dg_net = pickle.loads(dg_obj['Body'].read())
        dump_it_to_pickle(INDRA_DG_CACHE, dg_net)

        if FILES['sign_edge_graph_path'] is None:
            seg_key = 'indra_db_file/' + INDRA_SEG
            seg_obj = s3.get_object(Bucket=S3_BUCKET, Key=seg_key)
            seg_net = pickle.loads(seg_obj['Body'].read())
            dump_it_to_pickle(INDRA_SEG_CACHE, seg_net)
        else:
            seg_net = pickle_open(FILES['sign_edge_graph_path'])
        if argv[0].split('/')[-1].lower() != 'api.py':
            indra_network = IndraNetwork(indra_dir_graph=dg_net,
                                         indra_sign_edge_graph=seg_net)
    except Exception as e:
        logger.error('Could not find %s or %s locally or on s3' %
                     (INDRA_DG, INDRA_SEG))
        raise e

# Set verbosity
if VERBOSITY > 0 and\
        argv[0].split('/')[-1].lower() != 'api.py':
    logger.info('Setting verbosity to %d' % VERBOSITY)
    indra_network.verbose = VERBOSITY


@app.route('/')
@app.route('/query')
def get_query_page():
    """Loads the query page"""
    node_name_spaces = ['CHEBI', 'FPLX', 'GO', 'HGNC', 'HMDB', 'MESH',
                        'PUBCHEM']
    stmt_types = get_queryable_stmt_types()
    has_signed_graph = bool(len(indra_network.signed_nodes))
    return render_template('query_template.html',
                           stmt_types=stmt_types,
                           node_name_spaces=node_name_spaces,
                           has_signed_graph=has_signed_graph)


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
        # Test api by POSTing {'test': 'api'} to '/query/submit'
        if 'test' in request.json and request.json.get('test', False) and \
                request.json['test'] == 'api':
            logger.info('api test successful')
            return Response(json.dumps({'result': 'api test passed'}),
                            mimetype='application/json')
        result = indra_network.handle_query(**request.json.copy())
        logger.info('Query resolved at %s' %
                    strftime('%Y-%m-%d %H:%M:%S (UTC)', gmtime(time())))
        if not result or not all(result.values()):
            if API_DEBUG:
                logger.info('API_DEBUG is set to "True" so no network is '
                            'loaded, perhaps you meant to turn it off? '
                            'Run "export API_DEBUG=0" in your terminal to do '
                            'so and then restart the flask service')
            else:
                logger.info('Query returned with no path found')
        res = {'result': result}
        if indra_network.verbose > 5:
            logger.info('Result: %s' % str(res))
        return Response(json.dumps(res), mimetype='application/json')

    except KeyError as e:
        # return 400 badly formatted
        if indra_network.verbose:
            logger.exception(e)
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


# For those running with "python api.py"
if __name__ == '__main__':
    parser = argparse.ArgumentParser('Run the Indra DepMap Service.')
    parser.add_argument('--host', default='127.0.0.1')
    parser.add_argument('--port', default=5000, type=int)
    parser.add_argument('--test', action='store_true')
    parser.add_argument('--cache', nargs='+', help='Provide the file(s) for '
        'the network instead of reading the default parameters. Use "none" as '
        'a placeholder for a file that should be ignored. The DiGraph file '
        'path is mandatory.'
        'Usage: --cache <DiGraph pickle> [MultiDiGraph pickle|None] '
        '[SignedGraphModelChecker pickle|None]')
    parser.add_argument('-v', '--verbose', action='count', default=0)
    args = parser.parse_args()
    if API_DEBUG:
        pass
    elif args.cache:
        logger.info('Loading provided network files')
        dg_file = args.cache[0]
        mdg_file = args.cache[1] if len(args.cache) in (2, 3, 4) and\
            args.cache[1].lower() != 'none' else None
        seg_file = args.cache[2] if len(args.cache) >= 3 and\
            args.cache[2].lower() != 'none' else None
        sng_file = args.cache[3] if len(args.cache) > 3 and\
            args.cache[3].lower() != 'none' else None
        try:
            indra_network = \
                IndraNetwork(*load_indra_graph(dg_file, mdg_file))
        except Exception as e:
            logger.warning('Could not load the provided files. Reverting to '
                           'default network...')

    elif args.test:
        logger.info('Running test network')
        dg_file = TEST_DG_CACHE
        try:
            indra_network = \
                IndraNetwork(*load_indra_graph(dg_file))
        except Exception as e:
            logger.warning('Could not load the provided files. Reverting to '
                           'default network...')
    else:
        indra_network = IndraNetwork(
            *load_indra_graph(**FILES))

    if args.test:
        indra_network.small = True
        indra_network.verbose = args.verbose if args.verbose else 1
    if args.verbose:
        logger.info('Verbose level %d' % args.verbose)
        indra_network.verbose = args.verbose
    logger.info('Running app from "if __name__ == \'__main__\'"')
    app.run(host=args.host, port=args.port)
