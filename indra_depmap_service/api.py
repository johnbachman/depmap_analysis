"""INDRA Causal Network Search API"""
import requests
from sys import argv
from os import makedirs, environ
from time import time, gmtime, strftime

from indra_db.util.dump_sif import NS_PRIORITY_LIST as NS_LIST
from flask import Flask, request, abort, Response, render_template, jsonify,\
    session, url_for

from indra.config import CONFIG_DICT
from indralab_web_templates.path_templates import path_temps

from depmap_analysis.network_functions.indra_network import IndraNetwork, \
    EMPTY_RESULT
from depmap_analysis.network_functions.net_functions import SIGNS_TO_INT_SIGN

from .util import *

app = Flask(__name__)
app.register_blueprint(path_temps)
app.config['SECRET_KEY'] = environ.get('NETWORK_SEARCH_SESSION_KEY', '')
app.config['DEBUG'] = False
app.config['SECRET_KEY'] = environ['SESSION_KEY']
STMT_HASH_CACHE = {}

logger = logging.getLogger('INDRA Network Search API')

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

STMTS_FROM_HSH_URL = environ.get('INDRA_DB_HASHES_URL')
VERBOSITY = int(environ.get('VERBOSITY', 0))
API_DEBUG = int(environ.get('API_DEBUG', 0))
if API_DEBUG:
    logger.info('API_DEBUG set to %d' % API_DEBUG)

if not STMTS_FROM_HSH_URL:
    if API_DEBUG:
        logger.error('No URL for statement download set')
    else:
        raise ValueError('No URL for statement download set. Set it '
                         '"INDRA_DB_HASHES_URL" ')

GRND_URI = None
try:
    GRND_URI = CONFIG_DICT['GILDA_URL']
except KeyError:
    logger.warning('Indra Grounding service not available. Add '
                   'GILDA_URL to `indra/config.ini`')

MAX_PATHS = 50
TIMEOUT = 30  # Timeout in seconds
EMPTY_RESULT = {'paths_by_node_count': {'forward': {}, 'backward': {}},
                'common_targets': [],
                'common_parents': {},
                'timeout': False}


# Load network
indra_network = IndraNetwork()


def _is_empty_result(res):
    for k, v in res.items():
        if k is not 'timeout' and EMPTY_RESULT[k] != v:
            return False
    return True


if path.isfile(INDRA_DG_CACHE):
    INDRANET_DATE = datetime.utcfromtimestamp(get_earliest_date(
        INDRA_DG_CACHE))
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
        dg_obj = s3.get_object(Bucket=NET_BUCKET, Key=dg_key)
        dg_net = pickle.loads(dg_obj['Body'].read())
        dump_it_to_pickle(INDRA_DG_CACHE, dg_net)

        if FILES['sign_edge_graph_path'] is None:
            seg_key = 'indra_db_files/' + INDRA_SEG
            seg_obj = s3.get_object(Bucket=NET_BUCKET, Key=seg_key)
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


@app.route('/health')
def health():
    return jsonify({'status': 'pass'})


@app.route('/')
@app.route('/query')
def get_query_page():
    """Loads or responds to queries submitted on the query page"""
    logger.info('Got query')
    logger.info('Incoming Args -----------')
    logger.info(repr(request.args))

    stmt_types = get_queryable_stmt_types()
    has_signed_graph = bool(len(indra_network.signed_nodes))

    # Get query hash from parameters or session
    qh = request.args.get('query') or session.get('query_hash') or ''
    if qh:
        # Get query hash
        logger.info('Got query hash %s' % str(qh))
        old_results = check_existence_and_date_s3(qh)

        # Get result json
        res_json_key = old_results.get('result_json_key')
        results_json = read_query_json_from_s3(res_json_key) if res_json_key\
            else {}

        # Get query json
        query_json_key = old_results.get('query_json_key')
        query_json = read_query_json_from_s3(query_json_key) if \
            query_json_key else {}

        source = query_json.get('source', '')
        target = query_json.get('target', '')
    else:
        results_json = {'result': EMPTY_RESULT}
        query_json = {}
        source = ''
        target = ''
    return render_template('query_template.html',
                           query_hash=qh,
                           stmt_types=stmt_types,
                           node_name_spaces=list(NS_LIST),
                           has_signed_graph=has_signed_graph,
                           old_result=json.dumps(results_json),
                           old_query=json.dumps(query_json),
                           source=source,
                           target=target)


@app.route('/query/submit', methods=['POST'])
def process_query():
    """Processing queries to the indra network"""
    # Print inputs.
    logger.info('Got network search query')
    logger.info('Incoming Json ----------------------')
    logger.info(str(request.json))
    logger.info('------------------------------------')
    # Separate the requests that want JSON vs HTML, HTML is for POST
    # requests coming from a webpage and they should get back a redirect
    # with a query hash in the parameters

    try:
        # Test api by POSTing {'test': 'api'} to '/query/submit'
        if 'test' in request.json and request.json.get('test', False) and \
                request.json['test'] == 'api':
            logger.info('api test successful')
            return Response(json.dumps({'result': 'api test passed'}),
                            mimetype='application/json')

        # Get query json and query hash
        query_json = request.json.copy()
        ignore_keys = ['format']
        qc = {k: v for k, v in query_json.items() if k not in ignore_keys}
        qh = get_query_hash(qc)

        cached_files = check_existence_and_date_s3(query_hash=qh)
        if cached_files.get('result_json_key'):
            qjs3_key = cached_files['result_json_key']
            logger.info('Result found on s3: %s' % qjs3_key)
            result = read_query_json_from_s3(qjs3_key)
        # Files not cached on s3, run new query
        else:
            # Do new query
            # JS expects the following json for result:
            # {'results': '<indra_network.handle_query()>'
            #  'query_hash': '<query hash>'}
            result = handle_query(**request.json)

            # Empty result
            if _is_empty_result(result['result']):
                if API_DEBUG:
                    logger.info('API_DEBUG is set to "True" so no network is '
                                'loaded, perhaps you meant to turn it off? '
                                'Run "export API_DEBUG=0" in your terminal '
                                'to do so and then restart the flask service')
                else:
                    logger.info('Query returned with no path found')
                s3_query = ''
                result['query_hash'] = qh
                result['path_hashes'] = []
                query_json_fname = '%s_query.json' % qh
                s3_query_res = dump_query_result_to_s3(
                    filename=query_json_fname,
                    json_obj=query_json
                )

            # Non empty new result
            else:
                # Upload to public S3 and get the download link
                all_path_hashes = \
                    result['result']['paths_by_node_count'].get(
                        'path_hashes', [])
                result['query_hash'] = qh
                result['path_hashes'] = all_path_hashes

                # Upload the query itself
                query_json_fname = '%s_query.json' % qh
                s3_query = dump_query_result_to_s3(
                    filename=query_json_fname,
                    json_obj=query_json
                )
                # Upload query result
                res_json_fname = '%s_result.json' % qh
                s3_query_res = dump_query_result_to_s3(
                    filename=res_json_fname,
                    json_obj=result
                )
                logger.info('Uploaded query and results to %s and %s' %
                            (s3_query, s3_query_res))
        result['redirURL'] = url_for('get_query_page', query=qh)
        if request.json.get('format') and \
                request.json['format'].lower() == 'html':
            logger.info('HTML requested, sending redirect url')
            return url_for('get_query_page', query=qh)
        else:
            logger.info('Regular POST detected, sedning json back')
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


@app.route('/bfs_search', methods=['POST'])
def breadth_search():
    logger.info('Got request for breadth first search')
    logger.info('Incoming Json ----------------------')
    logger.info(str(request.json))
    logger.info('------------------------------------')
    if not request.json:
        return abort(415, '')
    query_json = request.json

    # Make lowercase
    allowed_ns = [ns.lower() for ns in query_json.get('node_filter', [])]
    default_ns = list(map(lambda s: s.lower(), NS_LIST))

    if not allowed_ns:
        allowed_ns = default_ns

    if not set(allowed_ns).issubset(set(default_ns)):
        abort(Response('One or more of the provided ns in "node_filter" is '
                       'not part of the standard ns. Provided ns list: %s. '
                       'Allowed ns list: %s' %
                       (str(allowed_ns), str(default_ns)), 415))

    if not query_json.get('source'):
        abort(Response('Missing required parameter "source"', 415))

    sign = SIGNS_TO_INT_SIGN[query_json.get('sign')]

    # If reversed, search upstream instead of downstream from source
    options = {
        'source': query_json['source'],
        'reverse': query_json.get('reverse', False),
        'depth_limit': int(query_json.get('depth_limit', 2)),
        'path_limit': int(query_json.get('path_limit', 100)),
        'node_filter': allowed_ns,
        'node_blacklist': query_json.get('node_blacklist', []),
        'bsco': float(query_json.get('belief_cutoff', 0)),
        'stmt_filter': query_json.get('skip_stmt_types', []),
        'curated_db_only': bool(query_json.get('db_only', False)),
        'terminal_ns': query_json.get('terminal_ns',  ['chebi', 'pubchem']),
        'max_results': query_json.get('max_results', 50),
        'max_per_node': query_json.get('max_per_node', 5),
        'sign': sign
    }
    try:
        results = indra_network.open_bfs_search(**options)
        return jsonify(results)
    except Exception as err:
        logger.warning('Exception handling open bfs search')
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

    query_hash = request.args.get('query', '')
    cached_files = check_existence_and_date_s3(query_hash=query_hash)
    if cached_files.get('result_json_key'):
        results_json = read_query_json_from_s3(
            s3_key=cached_files['result_json_key'])
        stmt_hashes = results_json.get('path_hashes', [])
    else:
        stmt_hashes = []

    if not STMTS_FROM_HSH_URL:
        logger.error('No URL for statement download set')
        return abort(500)

    if not stmt_hashes:
        logger.info('No hashes provided, returning')
        return jsonify({'result': 'no statement hashes found for query "%s"'
                                  % query_hash})

    logger.info('Got %d hashes' % len(stmt_hashes))

    stmt_list = []
    hash_list_iter = list_chunk_gen(stmt_hashes, size=1000)
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
