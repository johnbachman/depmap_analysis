"""INDRA Causal Network Search API"""
import os
import re
import json
import pickle
import logging
import argparse
import networkx as nx
from fnvhash import fnv1a_32
from os import path, makedirs, environ
import requests
import platform
from sys import argv
from fnvhash import fnv1a_32
from os import path, makedirs
from os import path, makedirs, environ, stat
from datetime import datetime
from time import time, gmtime, strftime

from jinja2 import Template
from indra_db.util.dump_sif import load_db_content, make_dataframe, \
    NS_PRIORITY_LIST as NS_LIST
from flask import Flask, request, abort, Response, render_template, jsonify,\
    session

from indra.config import CONFIG_DICT
from indra.util.aws import get_s3_client
from indralab_web_templates.path_templates import path_temps

from depmap_analysis.network_functions.indra_network import IndraNetwork,\
    EMPTY_RESULT
from depmap_analysis.util.io_functions import pickle_open, \
    dump_it_to_pickle, json_open, dump_it_to_json

from .util import load_indra_graph, get_queryable_stmt_types, API_PATH as \
    HERE, CACHE, INDRA_DG, INDRA_DG_CACHE, INDRA_SEG, INDRA_SEG_CACHE, \
    INDRA_SNG_CACHE, TEST_DG_CACHE, dump_query_result_to_s3

app = Flask(__name__)
app.register_blueprint(path_temps)
app.config['SECRET_KEY'] = environ.get('NETWORK_SEARCH_SESSION_KEY', '')
app.config['DEBUG'] = False
app.config['SECRET_KEY'] = os.environ['SESSION_KEY']

logger = logging.getLogger('INDRA Network Search API')

JSON_CACHE = path.join(HERE, '_json_res')
STATIC = path.join(HERE, 'static')
S3_BUCKET = 'depmap-analysis'
STMT_HASH_CACHE = {}


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

STMTS_FROM_HSH_URL = os.environ.get('INDRA_DB_HASHES_URL')
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
DT_YmdHMS_ = '%Y-%m-%d-%H-%M-%S'
DT_YmdHMS = '%Y%m%d%H%M%S'
DT_Ymd = '%Y%m%d'
RE_YmdHMS_ = r'\d{4}\-\d{2}\-\d{2}\-\d{2}\-\d{2}\-\d{2}'
RE_YYYYMMDD = r'\d{8}'
EMPTY_RESULT = {'paths_by_node_count': {'forward': {}, 'backward': {}},
                'common_targets': [],
                'common_parents': {},
                'timeout': False}


# Load network
indra_network = IndraNetwork()


def _todays_date():
    return datetime.now().strftime(DT_Ymd)


def _get_earliest_date(file):
    """Returns creation or modification timestamp of file"""
    # https://stackoverflow.com/questions/237079/
    # how-to-get-file-creation-modification-date-times-in-python
    if platform.system().lower() == 'windows':
        return path.getctime(file)
    else:
        st = stat(file)
        try:
            return st.st_birthtime
        except AttributeError:
            return st.st_mtime


def get_date_from_str(date_str, dt_format):
    """Returns a datetime object from a datestring of format FORMAT"""
    return datetime.strptime(date_str, dt_format)


def strip_out_date(keystring, re_format):
    """Strips out datestring of format re_format from a keystring"""
    try:
        return re.search(re_format, keystring).group()
    except AttributeError:
        logger.warning('Can\'t parse string %s for date using regex pattern '
                       '%s' % (keystring, re_format))
        return None


def _is_empty_result(res):
    for k, v in res.items():
        if k is not 'timeout' and EMPTY_RESULT[k] != v:
            return False
    return True


def _get_query_resp_fstr(query_hash):
    qf = path.join(JSON_CACHE, 'query_%s.json' % query_hash)
    rf = path.join(JSON_CACHE, 'result_%s.json' % query_hash)
    return qf, rf


def _list_chunk_gen(lst, size=1000):
    """Given list, generate chunks <= size"""
    n = max(1, size)
    return (lst[k:k+n] for k in range(0, len(lst), n))


def sorted_json_string(json_thing):
    """Produce a string that is unique to a json's contents."""
    if isinstance(json_thing, str):
        return json_thing
    elif isinstance(json_thing, (tuple, list)):
        return '[%s]' % (','.join(sorted(sorted_json_string(s)
                                         for s in json_thing)))
    elif isinstance(json_thing, dict):
        return '{%s}' % (','.join(sorted(k + sorted_json_string(v)
                                         for k, v in json_thing.items())))
    elif isinstance(json_thing, (int, float)):
        return str(json_thing)
    else:
        raise TypeError('Invalid type: %s' % type(json_thing))


def _get_query_hash(query_json):
    """Create an FNV-1a 32-bit hash from the query json and model_id."""
    return fnv1a_32(sorted_json_string(query_json).encode('utf-8'))


def _check_existence_and_date(fname):
    pass


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
        s3 = get_s3_client(unsigned=True)
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


def handle_query(**json_query):
    """Handle queries to the indra network"""
    result = indra_network.handle_query(**json_query.copy())
    logger.info('Query resolved at %s' %
                strftime('%Y-%m-%d %H:%M:%S (UTC)', gmtime(time())))
    if _is_empty_result(result):
        if API_DEBUG:
            logger.info('API_DEBUG is set to "True" so no network is '
                        'loaded, perhaps you meant to turn it off? '
                        'Run "export API_DEBUG=0" in your terminal to do '
                        'so and then restart the flask service')
        if result['timeout']:
            logger.info('Query timed out with no path found')
        else:
            logger.info('Query returned with no path found')
    res = {'result': result}
    if indra_network.verbose > 5:
        logger.info('Result: %s' % str(res))
    return res


@app.route('/')
@app.route('/query')
def get_query_page():
    """Loads or responds to queries submitted on the query page"""
    logger.info('Got query')
    logger.info('Incoming Args -----------')
    logger.info(repr(request.args))
    logger.info('Incoming Json ----------------------')
    logger.info(str(request.json))
    logger.info('------------------------------------')

    qh = session.get('query')
    rf = os.path.join(JSON_CACHE, 'result_%s.json' % qh) if qh else False
    qf = os.path.join(JSON_CACHE, 'query_%s.json' % qh) if qh else False
    stmt_types = get_queryable_stmt_types()
    has_signed_graph = bool(len(indra_network.signed_nodes))

    if all([qh, rf, qf, os.path.isfile(rf)]):
        with open(rf) as fr:
            search_results = json.load(fr)
        with open(qf) as fq:
            query = json.load(fq)
            source = query['source']
            target = query['target']
    qh = session.get('query_hash')
    rf = path.join(JSON_CACHE, 'result_%s.json' % qh) if qh else False
    qf = path.join(JSON_CACHE, 'query_%s.json' % qh) if qh else False
    if all([qh, rf, qf, path.isfile(rf)]):
        search_results = json_open(rf)
        query = json_open(qf)
        source = query['source']
        target = query['target']
    else:
        search_results = {'result': EMPTY_RESULT}
        query = {}
        source = ''
        target = ''
    node_name_spaces = ['CHEBI', 'FPLX', 'GO', 'HGNC', 'HMDB', 'MESH',
                        'PUBCHEM']
    return render_template('query_template.html',
                           stmt_types=stmt_types,
                           node_name_spaces=list(node_name_spaces),
                           has_signed_graph=has_signed_graph,
                           search_results=json.dumps(search_results),
                           old_query=json.dumps(query),
                           source=source,
                           target=target)


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
        query_json = request.json.copy()
        result = handle_query(**request.json)
        logger.info('Query resolved at %s' %
                    strftime('%Y-%m-%d %H:%M:%S (UTC)', gmtime(time())))
        if indra_network.verbose > 5:
            logger.info('Result: %s' % str(result))
        if request.json.get('format') and \
                request.json['format'].lower() == 'html':
            ignore_keys = ['format']
            qc = {k: v for k, v in query_json.items() if k not in ignore_keys}
            qh = _get_query_hash(qc)
            if _is_empty_result(result):
                if API_DEBUG:
                    logger.info('API_DEBUG is set to "True" so no network is '
                                'loaded, perhaps you meant to turn it off? '
                                'Run "export API_DEBUG=0" in your terminal to do '
                                'so and then restart the flask service')
                else:
                    logger.info('Query returned with no path found')
                download_link = ''
            else:
                # Upload to public S3 and get the download link
                json_fname = '%s_result.json' % qh
                download_link = dump_query_result_to_s3(json_fname, result)

            all_path_hashes = result['paths_by_node_count'].get('path_hashes', [])
            session['query_hash'] = qh
            STMT_HASH_CACHE[qh] = all_path_hashes
            res = {'result': result,
                   'download_link': download_link,
                   'query_hash': qh}
            # Return response with redirect URL
            return Response(json.dumps({'redirURL': '/query'}),
                            mimetype='application/json')
        else:
            return Response(json.dumps(result), mimetype='application/json')

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


@app.route('/multi_interactors', methods=['POST'])
def multi_interactors():
    logger.info('Got request for multi interactors')
    logger.info('Incoming Json ----------------------')
    logger.info(str(request.json))
    logger.info('------------------------------------')
    if not request.json:
        return abort(415, '')
    query_json = request.json
    if 'help' in query_json:
        return jsonify({
            'required': ['targets XOR regulators'],
            'optional': ['allowed_ns (all)', 'belief_cutoff (0)',
                         'skip_stmt_types ([])',
                         'db_only (False)']
        })
    if not (bool(query_json.get('targets')) ^
            bool(query_json.get('regulators'))):
        abort(Response('Missing required parameter "targets"', 415))

    # Requires the following options:
    # *node_filter
    # *bsco
    # *stmt_filter
    # *curated_db_only
    allowed_ns = [ns.lower() for ns in query_json.get('allowed_ns', [])]
    default_ns = list(map(lambda s: s.lower(), NS_LIST))

    if not allowed_ns:
        allowed_ns = default_ns

    if not set(allowed_ns).issubset(set(default_ns)):
        abort(Response('One or more of the provided ns in "allowed_ns" is '
                       'not part of the standard ns. Provided ns list: %s. '
                       'Allowed ns list: %s' %
                       (str(allowed_ns), str(default_ns)), 415))

    options = {
        'node_filter': allowed_ns,
        'bsco': float(query_json.get('belief_cutoff', 0)),
        'stmt_filter': query_json.get('skip_stmt_types', []),
        'curated_db_only': bool(query_json.get('db_only', False))
    }

    if query_json.get('targets'):
        options['list_of_targets'] = query_json['targets']
    else:
        options['list_of_regulators'] = query_json['regulators']

    try:
        result = indra_network.multi_regulators_targets(**options)
        return jsonify(result)
    except Exception as err:
        logger.warning('Error handling multi interactors query')
        logger.exception(err)
        abort(Response('Internal server error handling multi interactors '
                       'query', 500))


@app.route('/node', methods=['POST'])
def node_check():
    logger.info('Got request for node check')
    logger.info('Incoming Args -----------')
    logger.info(repr(request.args))
    logger.info('Incoming Json ----------------------')
    logger.info(str(request.json))
    logger.info('------------------------------------')

    node = request.json.get('node')

    if node:
        in_network = node in indra_network.nodes
        return jsonify({'node': node, 'in_network': in_network})
    return jsonify({'node': None, 'in_network': False})


@app.route('/nodes', methods=['POST'])
def nodes_check():
    logger.info('Got request to check multiple nodes')
    logger.info('Incoming Args -----------')
    logger.info(repr(request.args))
    logger.info('Incoming Json ----------------------')
    logger.info(str(request.json))
    logger.info('------------------------------------')

    nodes = request.json.get('nodes')

    if nodes:
        result = {node: node in indra_network.nodes for node in nodes}
        return jsonify(result)
    return jsonify({})


@app.route('/stmts_download/stmts.json')
def stmts_download():
    """Getting statement jsons from a list of hashes"""
    # Print inputs.
    logger.info('Got request for statements')
    logger.info('Incoming Args -----------')
    logger.info(repr(request.args))
    logger.info('Incoming Json ----------------------')
    logger.info(str(request.json))
    logger.info('------------------------------------')

    query_hash = session.get('query_hash', '')
    stmt_hashes = STMT_HASH_CACHE.pop(query_hash, [])

    if not STMTS_FROM_HSH_URL:
        logger.error('No URL for statement download set')
        return abort(500)

    if not stmt_hashes:
        logger.info('No hashes provided, returning')
        return jsonify([])

    logger.info('Got %d hashes' % len(stmt_hashes))

    stmt_list = []
    hash_list_iter = _list_chunk_gen(stmt_hashes)
    for hash_list in hash_list_iter:
        res = requests.post(STMTS_FROM_HSH_URL, json={'hashes': hash_list})
        if res.status_code == 200:
            stmts = res.json().get('statements')
            if stmts:
                stmt_list += [stmts]
    logger.info('Returning %d json statements' % len(stmt_list))
    return Response(json.dumps(stmt_list), mimetype='application/json',
                    content_type='attachment')


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
        '[SignedEdgeGraph pickle|None] [SignedNodeGraph pickle|None]')
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
                IndraNetwork(*load_indra_graph(dg_file, mdg_file, sng_file,
                                               seg_file))
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
