import os
import json
import boto3
import pickle
import logging
import argparse
import networkx as nx
from os import path, makedirs
from datetime import datetime
from botocore import UNSIGNED
from botocore.client import Config
from time import time, gmtime, strftime

from jinja2 import Template
from flask import Flask, request, abort, Response, render_template, url_for
from indra_db.util.dump_sif import load_db_content, make_dataframe, NS_LIST
from indra.config import CONFIG_DICT

from depmap_analysis.network_functions import net_functions as nf
from depmap_analysis.network_functions.indra_network import IndraNetwork
from depmap_analysis.util.io_functions import pickle_open, dump_it_to_pickle

app = Flask(__name__)
app.config['DEBUG'] = True

logger = logging.getLogger('INDRA Network Search API')

HERE = path.dirname(path.abspath(__file__))
CACHE = path.join(HERE, '_cache')
STATIC = path.join(HERE, 'static')
S3_BUCKET = 'depmap-analysis'
INDRA_MDG = 'indranet_multi_digraph_latest.pkl'
INDRA_DG = 'indranet_dir_graph_latest.pkl'
INDRA_SG_MC = 'indranet_sign_graph_mc_latest.pkl'

TEST_MDG_CACHE = path.join(CACHE, 'test_mdg_network.pkl')
INDRA_MDG_CACHE = path.join(CACHE, INDRA_MDG)
TEST_DG_CACHE = path.join(CACHE, 'test_dir_network.pkl')
INDRA_DG_CACHE = path.join(CACHE, INDRA_DG)
INDRA_SG_MC_CACHE = path.join(CACHE, INDRA_SG_MC)

VERBOSITY = int(os.environ.get('VERBOSITY', 0))

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


def _get_s3_client(unsigned=True):
    if unsigned:
        return boto3.client('s3',
                            config=Config(signature_version=UNSIGNED))
    else:
        return boto3.client('s3')


# Create a template object from the template file, load once
def _load_template(fname):
    template_path = path.join(HERE, fname)
    with open(template_path, 'rt') as f:
        template_str = f.read()
        template = Template(template_str)
    return template


def load_indra_graph(dir_graph_path, multi_digraph_path=None,
                     sign_graph_mc_path=None, update=False,
                     belief_dict=None, strat_ev_dict=None,
                     include_entity_hierarchies=True, verbosity=0):
    """Return a nx.DiGraph and nx.MultiDiGraph representation an INDRA DB dump

    If update is True, make a fresh snapshot from the INDRA DB.
    WARNING: this typically requires a lot of RAM and might slow down your
    computer significantly.
    """
    global INDRA_DG_CACHE, INDRA_MDG_CACHE, INDRA_SG_MC_CACHE
    indra_multi_digraph = nx.MultiDiGraph()
    signed_graph_mc = None
    if update:
        df = make_dataframe(True, load_db_content(True, NS_LIST))
        options = {'df': df,
                   'belief_dict': belief_dict,
                   'strat_ev_dict': strat_ev_dict,
                   'include_entity_hierarchies': include_entity_hierarchies,
                   'verbosity': verbosity}
        indra_dir_graph = nf.sif_dump_df_to_nx_digraph(**options)
        dump_it_to_pickle(dir_graph_path, indra_dir_graph)
        INDRA_DG_CACHE = path.join(CACHE, dir_graph_path)
        if multi_digraph_path:
            indra_multi_digraph = nf.sif_dump_df_to_nx_digraph(
                graph_type='multidigraph', **options)
            dump_it_to_pickle(multi_digraph_path, indra_multi_digraph)
            INDRA_MDG_CACHE = path.join(CACHE, multi_digraph_path)
        if sign_graph_mc_path:
            signed_graph_mc = nf.sif_dump_df_to_nx_digraph(
                graph_type='signed', **options)
    else:
        logger.info('Loading indra network%s %s %s' %
                    ('s' if multi_digraph_path else '',
                     dir_graph_path.split()[-1],
                     'and ' + multi_digraph_path.split()[-1] if
                     multi_digraph_path else ''))
        indra_dir_graph = pickle_open(dir_graph_path)
        if multi_digraph_path:
            indra_multi_digraph = pickle_open(multi_digraph_path)
        if sign_graph_mc_path:
            signed_graph_mc = pickle_open(sign_graph_mc_path)
        logger.info('Finished loading indra networks.')
    return indra_dir_graph, indra_multi_digraph, signed_graph_mc


if path.isfile(INDRA_DG_CACHE):
    if path.isfile(INDRA_SG_MC_CACHE):
        indra_network = IndraNetwork(
            *load_indra_graph(dir_graph_path=INDRA_DG_CACHE,
                              sign_graph_mc_path=INDRA_SG_MC_CACHE))
    else:
        indra_network = IndraNetwork(*load_indra_graph(INDRA_DG_CACHE))
else:
    # Try to find file(s) on s3
    try:
        logger.info('%s not found locally, trying to get file from s3...' %
                    INDRA_DG_CACHE)
        s3 = _get_s3_client(unsigned=True)

        dg_key = 'indra_db_files/' + INDRA_DG
        dg_obj = s3.get_object(Bucket=S3_BUCKET, Key=dg_key)
        dg_net = pickle.loads(dg_obj['Body'].read())

        sg_key = 'indra_db_file/' + INDRA_SG_MC
        sg_obj = s3.get_object(Bucket=S3_BUCKET, Key=sg_key)
        sg_net = pickle.loads(sg_obj['Body'].read())

        logger.info('Caching network to %s' % CACHE)
        try:
            makedirs(CACHE, exist_ok=True)
            dump_it_to_pickle(path.join(CACHE, INDRA_DG), dg_net)
            dump_it_to_pickle(path.join(CACHE, INDRA_SG_MC), sg_net)
        except Exception as e:
            logger.warning('Could not dump file to pickle')
        indra_network = IndraNetwork(indra_dir_graph=dg_net,
                                     indra_sign_graph_mc=sg_net)
    except Exception as e:
        logger.error('Could not find %s or %s locally or on s3' %
                     (INDRA_DG, INDRA_SG_MC))
        raise e

# Set verbosity
if VERBOSITY > 0:
    logger.info('Setting verbosity to %d' % VERBOSITY)
indra_network.verbose = VERBOSITY


@app.route('/')
@app.route('/query')
def get_query_page():
    """Loads the query page"""
    return render_template('query_template.html')


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

    if args.cache:
        logger.info('Loading provided network files')
        dg_file = args.cache[0]
        mdg_file = args.cache[1] if len(args.cache) in (2, 3) and\
            args.cache[1].lower() != 'none' else None
        sgmc_file = args.cache[2] if len(args.cache) >= 3 and\
            args.cache[2].lower() != 'none' else None
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

    if args.test:
        indra_network.small = True
        indra_network.verbose = args.verbose if args.verbose else 1
    if args.verbose:
        logger.info('Verbose level %d' % args.verbose)
        indra_network.verbose = args.verbose
    logger.info('Running app from "if __name__ == \'__main__\'"')
    app.run(host=args.host, port=args.port)
