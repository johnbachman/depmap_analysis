"""Utility functions for the INDRA Causal Network Search API in api.py"""
import boto3
import pickle
import logging
import argparse
import networkx as nx
from os import path
from botocore import UNSIGNED
from datetime import datetime
from botocore.client import Config

from indra.util.aws import get_s3_client, get_s3_file_tree
from indra_db.util.dump_sif import load_db_content, make_dataframe, NS_LIST
from indra.statements import get_all_descendants, Activation, Inhibition, \
    IncreaseAmount, DecreaseAmount, AddModification, RemoveModification

from depmap_analysis.network_functions import net_functions as nf
from depmap_analysis.util.io_functions import pickle_open, dump_it_to_pickle


logger = logging.getLogger('INDRA Network Search util')

API_PATH = path.dirname(path.abspath('./api.py'))
CACHE = path.join(API_PATH, '_cache')
STATIC = path.join(API_PATH, 'static')

INDRA_MDG = 'indranet_multi_digraph_latest.pkl'
INDRA_DG = 'indranet_dir_graph_latest.pkl'
INDRA_SNG = 'indranet_sign_node_graph_latest.pkl'
INDRA_SEG = 'indranet_sign_edge_graph_latest.pkl'

TEST_MDG_CACHE = path.join(CACHE, 'test_mdg_network.pkl')
INDRA_MDG_CACHE = path.join(CACHE, INDRA_MDG)
TEST_DG_CACHE = path.join(CACHE, 'test_dir_network.pkl')
INDRA_DG_CACHE = path.join(CACHE, INDRA_DG)
INDRA_SNG_CACHE = path.join(CACHE, INDRA_SNG)
INDRA_SEG_CACHE = path.join(CACHE, INDRA_SEG)

SIF_BUCKET = 'bigmech'
NET_BUCKET = 'depmap-analysis'


def _todays_date():
    return datetime.now().strftime('%Y%m%d')


# Copied from emmaa_service/api.py
def get_queryable_stmt_types():
    """Return Statement class names that can be used for querying."""
    def _get_sorted_descendants(cls):
        return sorted(_get_names(get_all_descendants(cls)))

    def _get_names(classes):
        return [s.__name__ for s in classes]

    stmt_types = \
        _get_names([
            Activation, Inhibition, IncreaseAmount, DecreaseAmount
        ]) + \
        _get_sorted_descendants(AddModification) + \
        _get_sorted_descendants(RemoveModification)
    return stmt_types


def load_indra_graph(dir_graph_path, multi_digraph_path=None,
                     sign_node_graph_path=None, sign_edge_graph_path=None,
                     update=False, belief_dict=None, strat_ev_dict=None,
                     include_entity_hierarchies=True, verbosity=0):
    """Return a nx.DiGraph and nx.MultiDiGraph representation an INDRA DB dump

    If update is True, make a fresh snapshot from the INDRA DB.
    WARNING: this typically requires a lot of RAM and might slow down your
    system significantly.
    """
    global INDRA_DG_CACHE, INDRA_MDG_CACHE, INDRA_SNG_CACHE, INDRA_SEG_CACHE
    indra_multi_digraph = nx.MultiDiGraph()
    indra_signed_edge_graph = nx.MultiDiGraph()
    indra_signed_node_graph = nx.DiGraph()

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
        if sign_node_graph_path or sign_edge_graph_path:
            indra_signed_edge_graph, indra_signed_node_graph = \
                nf.sif_dump_df_to_nx_digraph(graph_type='signed', **options)
    else:
        logger.info('Loading indra network representations from pickles')
        indra_dir_graph = pickle_open(dir_graph_path)
        if multi_digraph_path:
            indra_multi_digraph = pickle_open(multi_digraph_path)
        if sign_edge_graph_path:
            indra_signed_edge_graph = pickle_open(sign_edge_graph_path)
        if sign_node_graph_path:
            indra_signed_node_graph = pickle_open(sign_node_graph_path)
        logger.info('Finished loading indra networks.')
    return indra_dir_graph, indra_multi_digraph, indra_signed_edge_graph,\
        indra_signed_node_graph


def _get_latest_files():
    necc_files = ['belief_dict', 'db_dump_df', 'strat_ev']
    s3 = get_s3_client(unsigned=False)
    tree = get_s3_file_tree(s3, bucket=SIF_BUCKET, prefix='indra_db_sif_dump',
                     with_dt=True)
    # Find all pickles
    keys = [key for key in tree.gets('key') if key[0].endswith('.pkl')]
    # Sort newest first
    keys.sort(key=lambda t: t[1], reverse=True)
    # Find newest set of files' directory
    latest_dir = keys[0][0].split('/')[-2]
    # Get keys of those pickles
    keys_in_latest_dir = [k[0] for k in keys if latest_dir in k[0] and
                          any(nf in k[0] for nf in necc_files)]
    # Map key to resource
    necc_keys = {}
    for n in necc_files:
        for k in keys_in_latest_dir:
            if n in k:
                necc_keys[n] = k
    df = _load_pickle_from_s3(s3, key=necc_keys['db_dump_df'],
                              bucket=SIF_BUCKET)
    sev = _load_pickle_from_s3(s3, key=necc_keys['strat_ev'],
                              bucket=SIF_BUCKET)
    bd = _load_pickle_from_s3(s3, key=necc_keys['belief_dict'],
                              bucket=SIF_BUCKET)
    return df, sev, bd


def _load_pickle_from_s3(s3, key, bucket):
    try:
        res = s3.get_object(Key=key, Bucket=bucket)
        pyobj = pickle.loads(res['Body'].read())
    except Exception as err:
        logger.error('Someting went wrong while loading, reading or '
                         'unpickling the object from s3')
        raise err
    return pyobj


def _dump_network_to_s3(name, indranet_graph_object):
    s3 = get_s3_client(unsigned=False)
    key = 'indra_db_files/' + name
    s3.put_object(Bucket=NET_BUCKET, Key=key,
                  Body=pickle.dumps(obj=indranet_graph_object))


def dump_new_nets(mdg=None, dg=None, sg=None, dump_to_s3=False, verbosity=0):
    """Main script function for dumping new networks from latest db dumps"""
    df, sev, bd = _get_latest_files()
    options = {'df': df,
               'belief_dict': bd,
               'strat_ev_dict': sev,
               'include_entity_hierarchies': True,
               'verbosity': verbosity}

    if mdg:
        network = nf.sif_dump_df_to_nx_digraph(graph_type='multi', **options)
        dump_it_to_pickle(INDRA_MDG_CACHE, network)
        if dump_to_s3:
            _dump_network_to_s3(INDRA_MDG, network)
    if dg:
        network = nf.sif_dump_df_to_nx_digraph(**options)
        dump_it_to_pickle(INDRA_DG_CACHE, network)
        if dump_to_s3:
            _dump_network_to_s3(INDRA_DG, network)
    if sg:
        network, isng = nf.sif_dump_df_to_nx_digraph(graph_type='signed',
                                                     **options)
        dump_it_to_pickle(INDRA_SEG_CACHE, network)
        dump_it_to_pickle(INDRA_SNG_CACHE, isng)
        if dump_to_s3:
            _dump_network_to_s3(INDRA_SEG, network)
            _dump_network_to_s3(INDRA_SNG, network)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Dump new networks')
    parser.add_argument('--mdg', help='Dump a new MultiDiGraph')
    parser.add_argument('--dg', help='Dump a new DiGraph')
    parser.add_argument('--sg', help='Dump new signed edge and node graphs')
    parser.add_argument('--s3', help='Also upload the new graphs to s3')
    args = parser.parse_args()
    dump_new_nets(mdg=args.mdg, dg=args.dg, sg=args.sg, dump_to_s3=args.s3)
