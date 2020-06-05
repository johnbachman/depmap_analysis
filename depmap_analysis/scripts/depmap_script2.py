"""The DepMap script

The script matches observed correlations with different graph represenations
if the IndraNetwork

# Inputs:
#   1. Pre-processed correlation matrix with only (gene, gene, z-score)
#   2. nx.DiGraph (or nx.MultiDiGraph?) of IndraNetwork containing at least
#      agA/B: (name, ns, id), hash, type, belief, sign

# Processing:
#   0a. If signed graph: match edge sign with correlation sign
#   0b. If pybel graph: get node representation and match edge sign with
#       correlation sign
#   1. Find direct links both ways
#   2. Find A-X-B links: both ways (2), common target, common regulator
#   3. Find famplex links (have common parent)

# Questions:
#   Q1. Where would the cutting down to specific SD ranges be done?
#   A: Probably outside match correlations, somewhere inside or after
#      preprocessing. Better to do it all at once for one dump of the data

# Output:
#   An object of a new class that wraps a dataframe that can generate
#   different explanations statistics
"""
import inspect
import logging
import argparse
import multiprocessing as mp
from time import time
from math import floor
from pathlib import Path
from itertools import product
from collections import defaultdict
from datetime import datetime

import numpy as np
import pandas as pd
import networkx as nx
from pybel.dsl.node_classes import CentralDogma

from depmap_analysis.util.io_functions import pickle_open, dump_it_to_pickle
from depmap_analysis.network_functions.net_functions import \
    INT_MINUS, INT_PLUS, ns_id_from_name, get_hgnc_node_mapping
from depmap_analysis.network_functions.famplex_functions import common_parent
from depmap_analysis.network_functions.depmap_network_functions import \
    corr_matrix_to_generator, iter_chunker
from depmap_analysis.util.statistics import DepMapExplainer
from depmap_analysis.scripts.depmap_preprocessing import run_corr_merge

logger = logging.getLogger('DepMap Script')
logger.setLevel(logging.DEBUG)

indranet = None
hgnc_node_mapping = None
output_list = None


def _match_correlation_body(corr_iter, expl_types, stats_columns,
                            expl_columns, bool_columns, min_columns,
                            explained_set, _type):
    # Separate out this part

    stats_dict = {k: [] for k in stats_columns}
    expl_dict = {k: [] for k in expl_columns}

    for gA, gB, zsc in corr_iter:
        # Initialize current iteration stats
        stats = {k: False for k in bool_columns}

        # Append to stats_dict
        stats_dict['agA'].append(gA)
        stats_dict['agB'].append(gB)
        stats_dict['z-score'].append(zsc)

        # Skip if A or B not in graph or if type is pybel and no node
        # mapping exists
        if _type == 'pybel' and \
                not (hgnc_node_mapping[gA] and hgnc_node_mapping[gB]) or \
                _type != 'pybel' and \
                (gA not in indranet.nodes or gB not in indranet.nodes):
            for k in set(stats_dict.keys()).difference(set(min_columns)):
                if k == 'not in graph':
                    # Flag not in graph
                    stats_dict[k].append(True)
                else:
                    # All columns are NaN's
                    stats_dict[k].append(np.nan)
            continue

        if _type == 'pybel':
            # Get ns, id
            a_pb_node = list(hgnc_node_mapping[gA])[0]
            a_ns, a_id = a_pb_node.namespace, a_pb_node.identifier
            b_pb_node = list(hgnc_node_mapping[gB])[0]
            b_ns, b_id = b_pb_node.namespace, b_pb_node.identifier
        else:
            a_ns, a_id, b_ns, b_id = get_ns_id(gA, gB, indranet)

        # Append to stats dict
        stats_dict['agA_ns'].append(a_ns)
        stats_dict['agB_ns'].append(b_ns)
        stats_dict['agA_id'].append(a_id)
        stats_dict['agB_id'].append(b_id)

        # If in expl set, skip other explanations
        if explained_set:
            if gA in explained_set and gB in explained_set:
                # Set explained set = True
                stats_dict['explained set'].append(True)

                # Set overall explained = True
                stats_dict['explained'].append(True)

                # All other columns to False
                for k in set(bool_columns).difference(
                        {'explained set', 'explained'}):
                    stats_dict[k].append(False)

                # Set explanation type and data
                # Append to expl_dict
                expl_dict['agA'].append(gA)
                expl_dict['agB'].append(gB)
                expl_dict['z-score'].append(zsc)
                expl_dict['expl type'].append('explained set')
                expl_dict['expl data'].append(np.nan)

                # And skip the rest of explanations
                continue

        # Create iterator for pairs
        expl_iter = product(hgnc_node_mapping[gA], hgnc_node_mapping[gB]) \
            if _type == 'pybel' else [(gA, gB)]

        expl_iterations = defaultdict(list)
        for A, B in expl_iter:
            # Loop expl functions
            for expl_type, expl_func in expl_types.items():
                # Function signature: s, o, corr, net, graph_type, **kwargs
                # Function should return what will be kept in the 'expl_data'
                # column of the expl_df

                # Skip if 'explained set', which is caught above
                if expl_type == 'explained set':
                    continue

                # Some functions reverses A, B hence the s, o assignment
                s, o, expl_data = expl_func(A, B, zsc, indranet, _type)
                if expl_data:
                    expl_dict['agA'].append(s)
                    expl_dict['agB'].append(o)
                    expl_dict['z-score'].append(zsc)
                    expl_dict['expl type'].append(expl_type)
                    expl_dict['expl data'].append(expl_data)

                    # Append to expl_iterations
                    expl_iterations[expl_type].append(expl_data)

        # Check which ones got explained
        for expl_type_, expl_data_ in expl_iterations.items():
            stats[expl_type_] = bool(expl_data_)

        # Set explained column
        stats['explained'] = any([b for b in stats.values()])

        # Add stats to stats_dict
        for expl_tp in stats:
            stats_dict[expl_tp].append(stats[expl_tp])

        # Assert that all columns are the same length
        if not all(len(ls) for ls in stats_dict.values()):
            raise IndexError('Unequal column lengths in stats_dict after '
                             'iteration')
    return stats_dict, expl_dict


def match_correlations(corr_z, sd_range, **kwargs):
    """The main loop for matching correlations with INDRA explanations

    Parameters
    ----------
    corr_z : pd.DataFrame
        The pre-processed correlation matrix. No more processing of the
        matrix should have to be done here, i.e. it should already have
        filtered the correlations to the proper SD ranges and removed the
        genes that are not applicable for this explanation,
        self correlations should also have been removed.
    indranet : nx.DiGraph
        The graph representation of the indra network. Each edge should
        have an attribute named 'statements' containing a list of sources
        supporting that edge. If signed search, indranet is expected to be an
        nx.MultiDiGraph with edges keyes by (gene, gene, sign) tuples.
    sd_range : tuple[float]
        The SD ranges that the corr_z is filtered to

    Returns
    -------
    depmap_explainer
        An instance of the DepMapExplainer class containing the explanations
        for the correlations.
    """
    min_columns = ('agA', 'agB', 'z-score')
    id_columns = min_columns + ('agA_ns', 'agA_id', 'agB_ns', 'agB_id')
    # Map each expl type to a function that handles that explanation
    expl_types = {'a-b': expl_ab,
                  'b-a': expl_ba,
                  'common parent': find_cp,
                  'explained set': explained,  # a priori explained
                  'a-x-b': expl_axb,
                  'b-x-a': expl_bxa,
                  'shared regulator': get_sr,
                  'shared target': get_st,
                  }
    bool_columns = ('not in graph', 'explained') + tuple(expl_types.keys())
    stats_columns = id_columns + bool_columns
    expl_columns = ('agA', 'agB', 'z-score', 'expl type', 'expl data')
    explained_set = kwargs.get('explained_set', {})

    _type = kwargs.get('graph_type', 'unsigned')
    logger.info(f'Doing correlation matching with {_type} graph')
    ymd_now = datetime.now().strftime('%Y%m%d')
    indra_date = kwargs['indra_date'] if kwargs.get('indra_date') \
        else ymd_now
    depmap_date = kwargs['depmap_date'] if kwargs.get('depmap_date') \
        else ymd_now

    bool_matrix = np.invert(np.isnan(corr_z.values))
    estim_pairs = floor((bool_matrix.sum() - bool_matrix.diagonal().sum())/2)
    print(f'Starting workers at {datetime.now().strftime("%H:%M:%S")} with '
          f'about {estim_pairs} pairs to check')
    tstart = time()

    with mp.Pool() as pool:
        MAX_SUB = 512
        n_sub = min(kwargs.get('n-chunks', 256), MAX_SUB)
        chunksize = max(estim_pairs // n_sub, 1)

        # Pick one more so we don't do more than MAX_SUB
        chunksize += 1 if n_sub == MAX_SUB else 0
        chunk_iter = iter_chunker(n=chunksize,
                                  iterable=corr_matrix_to_generator(corr_z))
        for chunk in chunk_iter:
            pool.apply_async(func=_match_correlation_body,
                             # corr_iter, expl_types, stats_columns,
                             # expl_columns, bool_columns, min_columns,
                             # explained_set, _type
                             args=(
                                 chunk,
                                 expl_types,
                                 stats_columns,
                                 expl_columns,
                                 bool_columns,
                                 min_columns,
                                 explained_set,
                                 _type
                             ),
                             callback=success_callback)

        logger.info('Done submitting work to pool workers')
        pool.close()
        pool.join()

    print(f'Execution time: {time() - tstart} seconds')
    print(f'Done at {datetime.now().strftime("%H:%M:%S")}')

    # Here initialize a DepMapExplainer and append the result fro the
    # different processes
    explainer = DepMapExplainer(stats_columns=stats_columns,
                                expl_columns=expl_columns,
                                info={'indra_network_date': indra_date,
                                      'depmap_date': depmap_date,
                                      'sd_range': sd_range,
                                      'graph_type': _type,
                                      **kwargs.get('info', {})},
                                )

    logger.info(f'Generating DepMapExplainer with output from '
                f'{len(output_list)} results')
    for stats_dict, expl_dict in output_list:
        explainer.stats_df = explainer.stats_df.append(other=pd.DataFrame(
            data=stats_dict))
        explainer.expl_df = explainer.expl_df.append(other=pd.DataFrame(
            data=expl_dict))

    explainer.has_data = True
    return explainer


def explained(s, o, corr, net, graph_type, **kwargs):
    # This function is used for a priori explained relationships
    return s, o, 'explained_set'


def find_cp(s, o, corr, net, _type, **kwargs):
    # This function is the same for all graph_types
    s_ns, s_id, o_ns, o_id = get_ns_id(s, o, net)

    if not s_id:
        s_ns, s_id = ns_id_from_name(s)
    if not o_id:
        o_ns, o_id = ns_id_from_name(o)

    if s_id and o_id:
        parents = list(common_parent(ns1=s_ns, id1=s_id, ns2=o_ns, id2=o_id))
        if parents:
            return s, o, parents

    return s, o, None


def expl_axb(s, o, corr, net, _type, **kwargs):
    x_set = set(net.succ[s]) & set(net.pred[o])
    if _type in {'signed', 'pybel'}:
        x_nodes = _get_signed_interm(s, o, corr, net, x_set)
    else:
        x_nodes = x_set

    if x_nodes:
        return s, o, list(x_nodes)
    else:
        return s, o, None


def expl_bxa(s, o, corr, net, _type, **kwargs):
    return expl_axb(o, s, corr, net, _type, **kwargs)


# Shared regulator: A<-X->B
def get_sr(s, o, corr, net, _type, **kwargs):
    x_set = set(net.pred[s]) & set(net.pred[o])

    if _type in {'signed', 'pybel'}:
        x_nodes = _get_signed_interm(s, o, corr, net, x_set)
    else:
        x_nodes = x_set

    if x_nodes:
        return s, o, list(x_nodes)
    else:
        return s, o, None


# Shared target: A->X<-B
def get_st(s, o, corr, net, _type, **kwargs):
    x_set = set(net.succ[s]) & set(net.succ[o])

    if _type in {'signed', 'pybel'}:
        x_nodes = _get_signed_interm(s, o, corr, net, x_set)
    else:
        x_nodes = x_set

    if x_nodes:
        return s, o, list(x_nodes)
    else:
        return s, o, None


def expl_ab(s, o, corr, net, _type, **kwargs):
    edge_dict = get_edge_statements(s, o, corr, net, _type, **kwargs)
    if edge_dict:
        return s, o, edge_dict.get('stmt_hash') if _type == 'pybel' else \
            edge_dict.get('statements')
    return s, o, None


def expl_ba(s, o, corr, net, _type, **kwargs):
    # Reverse order call to expl_ab
    return expl_ab(o, s, corr, net, _type, **kwargs)


def get_edge_statements(s, o, corr, net, _type, **kwargs):
    if _type in {'signed', 'pybel'}:
        int_sign = INT_PLUS if corr >= 0 else INT_MINUS
        return net.edges.get((s, o, int_sign), None)
    else:
        return net.edges.get((s, o))


def _get_signed_interm(s, o, corr, sign_edge_net, x_set):
    # Make sure we have the right sign type
    int_sign = INT_PLUS if corr >= 0 else INT_MINUS

    # ax and xb sign need to match correlation sign
    x_approved = set()
    for x in x_set:
        ax_plus = sign_edge_net.edges.get((s, x, INT_PLUS), {})
        ax_minus = sign_edge_net.edges.get((s, x, INT_MINUS), {})
        xb_plus = sign_edge_net.edges.get((x, o, INT_PLUS), {})
        xb_minus = sign_edge_net.edges.get((x, o, INT_MINUS), {})

        if int_sign == INT_PLUS:
            if ax_plus and xb_plus or ax_minus and xb_minus:
                x_approved.add(x)
        if int_sign == INT_MINUS:
            if ax_plus and xb_minus or ax_minus and xb_plus:
                x_approved.add(x)
    return x_approved


def get_ns_id(subj, obj, net):
    """Get ns:id for both subj and obj

    Parameters
    ----------

    subj : str
        The subject node
    obj : str
        The source node
    net : nx.Graph
        A networkx graph object that at least contains node entries.

    Returns
    -------
    tuple
        A tuple with four entries:
        (subj namespace, subj id, obj namespace, obj id)
    """
    if isinstance(subj, CentralDogma):
        s_pbn = list(hgnc_node_mapping[subj])[0]
        s_ns, s_id = s_pbn.namespace, s_pbn.identifier
        o_pbn = list(hgnc_node_mapping[obj])[0]
        o_ns, o_id = o_pbn.namespace, o_pbn.identifier

    else:
        s_ns = net.nodes[subj]['ns'] if net.nodes.get(subj) else None
        s_id = net.nodes[subj]['id'] if net.nodes.get(subj) else None
        o_ns = net.nodes[obj]['ns'] if net.nodes.get(obj) else None
        o_id = net.nodes[obj]['id'] if net.nodes.get(obj) else None

    return s_ns, s_id, o_ns, o_id


def success_callback(res):
    logger.info('Appending a result')
    output_list.append(res)


def graph_types(types):
    """Types is a set of strings with names of the allowed graph types"""
    def types_check(_type):
        """Check the input graph type

        Parameters
        ----------
        _type : str
            The input graph type

        Returns
        -------
        str
            Returns the lowercase of the input string representing the graph
            type
        """
        if _type.lower() not in types:
            raise argparse.ArgumentError(f'Provided graph type {_type} not '
                                         f'allowed. Have to be one of {types}')
        return _type.lower()
    return types_check


def file_path():
    """Checks if file at provided path exists"""
    def check_path(fpath):
        p = Path(fpath)
        if not p.is_file():
            raise argparse.ArgumentError(f'File {fpath} does not exist')
        return fpath
    return check_path


def main(indra_net, sd_range, outname, graph_type, z_score=None,
         raw_data=None, raw_corr=None, pb_model=None,
         pb_node_mapping=None, n_chunks=256, ignore_list=None, info=None,
         indra_date=None, depmap_date=None):
    """Set up correlation mtaching of depmap data with an indranet graph

    Parameters
    ----------
    indra_net : nx.DiGraph|nx.MultiDiGraph
    sd_range : tuple(float|None)
    outname : str
    graph_type : str
    z_score : pd.DataFrame
    raw_data : str
    raw_corr : str
    pb_model : PyBEL network
    pb_node_mapping : dict
    n_chunks : int
    ignore_list : list
    info : dict
    indra_date : str
    depmap_date : str

    Returns
    -------
    depmap_analysis.util.statistics.DepMapExplainer
    """
    global indranet, hgnc_node_mapping, output_list
    indranet = indra_net
    hgnc_node_mapping = pb_node_mapping

    run_options = {'n-chunks': n_chunks}

    # 1 Check options
    sd_l, sd_u = sd_range if len(sd_range) == 2 else \
        ((sd_range[0], None) if len(sd_range) == 1 else (None, None))

    if not sd_l and not sd_u:
        raise ValueError('Must specify at least a lower bound for the SD '
                         'range')

    if graph_type == 'pybel' and not pb_model:
        raise ValueError('Must provide PyBEL model with option pybel_model '
                         'if graph type is pybel')

    # Get the z-score corr matrix
    if not (bool(z_score is not None and len(z_score)) ^
            bool(raw_data or raw_corr)):
        raise ValueError('Must provide either z_score XOR either of raw_data '
                         'or raw_corr')

    outname = outname if outname.endswith('.pkl') else \
        outname + '.pkl'
    outpath = Path(outname)

    if z_score is not None:
        z_corr = z_score
    else:
        z_sc_options = {
            'crispr_raw': raw_data[0],
            'rnai_raw': raw_data[1],
            'crispr_corr': raw_corr[0],
            'rnai_corr': raw_corr[1]
        }
        z_corr = run_corr_merge(**z_sc_options)

    graph_type = graph_type
    run_options['graph_type'] = graph_type

    # Get mapping of correlation names to pybel nodes
    if graph_type == 'pybel':
        hgnc_node_mapping = get_hgnc_node_mapping(
            hgnc_names=z_corr.columns.values,
            pb_model=pb_model,
            pb_signed_edge_graph=indranet
        )

    # 2. Filter to SD range
    if sd_l and sd_u:
        logger.info(f'Filtering correlations to {sd_l} - {sd_u} SD')
        z_corr = z_corr[((z_corr > sd_l) & (z_corr < sd_u)) |
                        ((z_corr < -sd_l) & (z_corr > -sd_u))]
    elif sd_l and not sd_u:
        logger.info(f'Filtering correlations to {sd_l}+ SD')
        z_corr = z_corr[(z_corr > sd_l) | (z_corr < -sd_l)]
    run_options['corr_z'] = z_corr
    run_options['sd_range'] = (sd_l, sd_u) if sd_u else (sd_l, None)

    # 3. Ignore list as file

    if ignore_list:
        run_options['explained_set'] = set(ignore_list)

    # 4. Add meta data
    info_dict = {}
    if info:
        info_dict['info'] = info
    if depmap_date:
        info_dict['depmap_date'] = depmap_date
    if indra_date:
        info_dict['indra_date'] = indra_date
    run_options['info'] = info_dict

    # Create output list in global scope
    output_list = []
    explanations = match_correlations(**run_options)

    # mkdir in case it  doesn't exist
    outpath.parent.mkdir(parents=True, exist_ok=True)
    dump_it_to_pickle(fname=outpath.absolute().resolve().as_posix(),
                      pyobj=explanations)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('DepMap Explainer main script')
    #   1a Load depmap data from scratch | load crispr/rnai raw corr | z-score
    corr_group = parser.add_mutually_exclusive_group(required=True)
    corr_group.add_argument(
        '--raw-data', nargs=2, type=file_path(),
        help='File paths to CRISPR raw data and RNAi raw data from the '
             'DepMap Portal. The CRISPR file name should match '
             '*gene_effect.csv. The RNAi file name should match '
             '*gene_dep_scores.csv'
    )
    corr_group.add_argument(
        '--raw-corr', nargs=2, type=file_path(),
        help='File paths to raw correlation data (before z-score conversion) '
             'containing hdf compressed correlation data. These files '
             'contain the result of running `raw_df.corr()`'
    )
    corr_group.add_argument(
        '--z-score', type=file_path(),
        help='The file path to the fully merged correlation matrix '
             'containing z-scores.')

    #   1b Load indranet
    parser.add_argument(
        '--indranet', type=file_path(), required=True,
        help='The indra network to use for explanations. Should be either a '
             'DiGraph or signed DiGraph (a MultiDiGraph with max two edges '
             'per node pair, one edge per sign).'
    )

    #   1c Optionally provide PyBEL model
    parser.add_argument(
        '--pybel-model', type=file_path(),
        help='If graph type is "pybel", use this argument to provide the '
             'necessary pybel model used to precompute the pybel node mapping.'
    )

    #   1d Provide graph type
    allowed_types = {'unsigned', 'signed', 'pybel'}
    parser.add_argument(
        '--graph-type', type=graph_types(allowed_types),
        default='unsigned',
        help='Specify the graph type used. Allowed values are '
             f'{allowed_types}'
    )

    #   2. Filter to SD range
    parser.add_argument('--sd-range', nargs='+', type=float, required=True,
                        help='SD range to filter to')
    #   3. Ignore list as file
    parser.add_argument(
        '--ignore-list', type=str,
        help='Provide a text file with one gene name per line to skip in the'
             'explanations')

    # 4 output
    parser.add_argument(
        '--outname', required=True, type=str,
        help='The output name (could contain a path as well) of the pickle '
             'dump of the explainer object')

    # 5 Pick number of jobs
    parser.add_argument(
        '--n-chunks', type=int, default=256,
        help='Pick the number of slices to split the work into. Does not '
             'have to be equal to the amount of CPUs.'
    )

    # 6 Extra info
    parser.add_argument('--indra-date',
                        help='Provide a date for the dump from which the '
                             'indra network was created')
    parser.add_argument('--depmap-date',
                        help='Provide the release date of the depmap data '
                             'used.')

    args = parser.parse_args()
    arg_dict = vars(args)

    # Load z_corr, indranet and optionally pybel_model
    inet_graph = pickle_open(args.indranet)
    arg_dict['indra_net'] = inet_graph
    if arg_dict.get('z_score'):
        corr_matrix = pd.read_hdf(arg_dict['z_score'])
        arg_dict['z_score'] = corr_matrix
    else:
        arg_dict['raw_data'] = arg_dict.get('raw_data')
        arg_dict['corr_data'] = arg_dict.get('corr_data')

    if arg_dict.get('graph_type') == 'pybel' and \
            not arg_dict.get('pybel_model'):
        raise ValueError('Must provide PyBEL model with option pybel_model '
                         'if graph type is pybel')
    if arg_dict.get('pybel_model'):
        arg_dict['pb_model'] = pickle_open(arg_dict['pybel_model'])

    main_keys = inspect.signature(main).parameters.keys()
    kwargs = {k: v for k, v in arg_dict.items() if k in main_keys}

    main(**kwargs)
