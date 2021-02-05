"""The DepMap script

The script matches observed correlations with different graph
representations if the IndraNetwork

# Inputs:
1. Pre-processed correlation matrix with only (gene, gene, z-score)
2. nx.DiGraph (or nx.MultiDiGraph?) of IndraNetwork containing at least
   agA/B: (name, ns, id), hash, type, belief, sign

# Processing:
0a. If signed graph: match edge sign with correlation sign
0b. If pybel graph: get node representation and match edge sign with
    correlation sign
1. Find direct links both ways
2. Find A-X-B links: both ways (2), common target, common regulator
3. Find famplex links (have common parent)

# Questions:
Q1. Where would the cutting down to specific SD ranges be done?
A: Probably outside match correlations, somewhere inside or after
   preprocessing. Better to do it all at once for one dump of the data

# Output:
An instance of the DepMapExplainer class that wraps dataframes that can
generate different explanations statistics
"""
import pickle
import inspect
import logging
import argparse
import multiprocessing as mp
from time import time
from copy import deepcopy
from typing import Union, List, Dict, Iterable, Tuple, Optional, Generator, \
    Callable, Set, Hashable, Any
from pathlib import Path
from itertools import product
from collections import defaultdict
from datetime import datetime

import numpy as np
import pandas as pd
import networkx as nx

from indra.util import batch_iter
from indra.util.multiprocessing_traceback import WrapException
from indra_db.util.s3_path import S3Path
from depmap_analysis.util.aws import get_s3_client
from depmap_analysis.util.io_functions import file_opener, \
    dump_it_to_pickle, allowed_types, file_path, strip_out_date
from depmap_analysis.network_functions.net_functions import \
    pybel_node_name_mapping
from depmap_analysis.network_functions.depmap_network_functions import \
    corr_matrix_to_generator, down_sampl_size
from depmap_analysis.util.statistics import DepMapExplainer, min_columns, \
    id_columns
from depmap_analysis.preprocessing import *
from depmap_analysis.scripts.depmap_script_expl_funcs import *

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

mito_file_name = 'Human.MitoCarta3.0.xls'
mito_file = Path(__file__).absolute().parent.parent\
    .joinpath('resources', mito_file_name).as_posix()

indranet: nx.DiGraph = nx.DiGraph()
hgnc_node_mapping: Dict[str, Set] = dict()
output_list: List[Tuple[Dict[str, List], Dict[str, List]]] = []


def _match_correlation_body(corr_iter: Generator[Tuple[str, str, float],
                                                 None, None],
                            expl_types: Dict[str, Callable],
                            stats_columns: Tuple[str],
                            expl_columns: Tuple[str],
                            bool_columns: Tuple[str],
                            expl_mapping: Optional[Dict[str, str]],
                            _type: str,
                            allowed_ns: Optional[Set[str]] = None,
                            allowed_sources: Optional[Set[str]] = None,
                            is_a_part_of: Optional[List[str]] = None,
                            immediate_only: Optional[bool] = False,
                            strict_intermediates: Optional[bool] = False):
    try:
        global indranet

        stats_dict = {k: [] for k in stats_columns}
        expl_dict = {k: [] for k in expl_columns}
        options = {'immediate_only': immediate_only,
                   'strict_intermediates': strict_intermediates}
        if is_a_part_of:
            options['is_a_part_of'] = is_a_part_of
        if allowed_ns:
            options['ns_set'] = allowed_ns
        if allowed_sources:
            options['src_set'] = allowed_sources
        if expl_mapping:
            options['expl_mapping'] = expl_mapping

        for tup in corr_iter:
            # Break loop when batch_iter
            if tup is None:
                break
            gA, gB, zsc = tup
            # Initialize current iteration stats
            stats = {k: False for k in bool_columns}

            # Append to stats_dict; These assignments need to be before the
            # checks that will skip current iteration
            stats_dict['agA'].append(gA)
            stats_dict['agB'].append(gB)
            stats_dict['z-score'].append(zsc)

            # Skip if A or B not in graph or (if type is pybel) no node
            # mapping exists for either A or B
            if _type == 'pybel' and (gA not in hgnc_node_mapping or
                                     gB not in hgnc_node_mapping) or \
                    _type != 'pybel' and (gA not in indranet.nodes or
                                          gB not in indranet.nodes):
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
                a_ns, a_id = get_ns_id_pybel_node(gA, tuple(hgnc_node_mapping[gA]))
                b_ns, b_id = get_ns_id_pybel_node(gB, tuple(hgnc_node_mapping[gB]))
            else:
                a_ns, a_id, b_ns, b_id = get_ns_id(gA, gB, indranet)

            # Append to stats dict
            stats_dict['agA_ns'].append(a_ns)
            stats_dict['agB_ns'].append(b_ns)
            stats_dict['agA_id'].append(a_id)
            stats_dict['agB_id'].append(b_id)

            # Create iterator for pairs
            expl_iter = product(hgnc_node_mapping[gA], hgnc_node_mapping[gB]) \
                if _type == 'pybel' else [(gA, gB)]

            # Add hgnc symbol name to expl kwargs if pybel
            if _type == 'pybel':
                options['s_name'] = gA
                options['o_name'] = gB

            expl_iterations = defaultdict(list)
            for A, B in expl_iter:
                # Loop expl functions. Args:
                # s, o, corr, net, graph_type, **kwargs
                expl_args = (A, B, zsc, indranet, _type)
                for expl_type, expl_func in expl_types.items():
                    # Function should return what will be kept in the
                    # 'expl_data' column of the expl_df
                    # Some functions reverses A, B hence the s, o assignment
                    s, o, expl_data = expl_func(*expl_args, **options)
                    if expl_data:
                        # Use original name
                        s_name = s.name if _type == 'pybel' else s
                        o_name = o.name if _type == 'pybel' else o
                        expl_dict['agA'].append(s_name)
                        expl_dict['agB'].append(o_name)
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
    except Exception as exc:
        raise WrapException()


def match_correlations(corr_z: pd.DataFrame,
                       sd_range: Tuple[float, Union[float, None]],
                       script_settings: Dict[str, Union[str, int, float]],
                       **kwargs):
    """The main loop for matching correlations with INDRA explanations

    Note that indranet is assumed to be a global variable that needs to be
    set outside of this function and be set to global

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
        nx.MultiDiGraph with edges keys by (gene, gene, sign) tuples.
    sd_range : tuple[float]
        The SD ranges that the corr_z is filtered to
    script_settings : Dict[str, Union[str, int, float]]
        Dictionary with script settings for the purpose of book keeping

    Returns
    -------
    depmap_analysis.util.statistics.DepMapExplainer
        An instance of the DepMapExplainer class containing the explanations
        for the correlations.
    """
    # Map each expl type to a function that handles that explanation
    if not kwargs.get('expl_funcs'):
        # No function names provided, use all explanation functions
        logger.info('All explanation types used')
        expl_types = {funcname_to_colname[func_name]: func
                      for func_name, func in expl_functions.items()}
    else:
        # Map function names to functions, check if
        expl_types = {}
        for func_name in kwargs['expl_funcs']:
            if func_name not in expl_functions:
                logger.warning(f'{func_name} does not map to a registered '
                               f'explanation function. Allowed functions '
                               f'{", ".join(expl_functions.keys())}')
            else:
                expl_types[funcname_to_colname[func_name]] = \
                    expl_functions[func_name]

    if not len(expl_types):
        raise ValueError('No explanation functions provided')

    bool_columns = ('not in graph', 'explained') + tuple(expl_types.keys())
    stats_columns = id_columns + bool_columns
    expl_columns = ('agA', 'agB', 'z-score', 'expl type', 'expl data')
    expl_mapping = kwargs.get('expl_mapping', {})

    _type = kwargs.get('graph_type', 'unsigned')
    logger.info(f'Doing correlation matching with {_type} graph')

    # Get options
    if kwargs.get('allowed_ns'):
        allowed_ns = {n.lower() for n in kwargs['allowed_ns']}
        logger.info('Only allowing the following namespaces: %s' %
                    ', '.join(allowed_ns))
    else:
        allowed_ns = None
    if kwargs.get('allowed_sources'):
        allowed_sources = {s.lower() for s in kwargs['allowed_sources']}
        logger.info('Only allowing the following sources: %s' %
                    ', '.join(allowed_sources))
    else:
        allowed_sources = None
    is_a_part_of = kwargs.get('is_a_part_of')
    immediate_only = kwargs.get('immediate_only', False)
    strict_intermediates = kwargs.get('strict_intermediates', False)

    # Try to get dates of files from file names and file info
    ymd_now = datetime.now().strftime('%Y%m%d')
    inet_file = script_settings['indranet']
    indra_date = kwargs['indra_date'] if kwargs.get('indra_date') \
        else (strip_out_date(inet_file, r'\d{8}')
              if strip_out_date(inet_file, r'\d{8}') else ymd_now)
    dm_file = script_settings['z_score']
    depmap_date = kwargs['depmap_date'] if kwargs.get('depmap_date') \
        else (strip_out_date(dm_file, r'\d{8}')
              if strip_out_date(dm_file, r'\d{8}') else ymd_now)

    logger.info('Calculating number of pairs to check...')
    estim_pairs = _get_pairs(corr_z)
    logger.info(f'Starting workers at {datetime.now().strftime("%H:%M:%S")} '
                f'with about {estim_pairs} pairs to check')
    tstart = time()

    with mp.Pool() as pool:
        MAX_SUB = 512
        n_sub = min(kwargs.get('n-chunks', 256), MAX_SUB)
        chunksize = max(estim_pairs // n_sub, 1)

        # Pick one more so we don't do more than MAX_SUB
        chunksize += 1 if n_sub == MAX_SUB else 0
        chunk_iter = batch_iter(iterator=corr_matrix_to_generator(corr_z),
                                batch_size=chunksize, return_func=list)
        for chunk in chunk_iter:
            pool.apply_async(func=_match_correlation_body,
                             # args should match the args for func
                             args=(
                                 chunk,
                                 expl_types,
                                 stats_columns,
                                 expl_columns,
                                 bool_columns,
                                 expl_mapping,
                                 _type,
                                 allowed_ns,
                                 allowed_sources,
                                 is_a_part_of,
                                 immediate_only,
                                 strict_intermediates
                             ),
                             callback=success_callback,
                             error_callback=error_callback)

        logger.info('Done submitting work to pool workers')
        pool.close()
        pool.join()

    print(f'Execution time: {time() - tstart} seconds')
    print(f'Done at {datetime.now().strftime("%H:%M:%S")}')

    # Here initialize a DepMapExplainer and append the result for the
    # different processes
    explainer = DepMapExplainer(stats_columns=stats_columns,
                                expl_columns=expl_columns,
                                info={'indra_network_date': indra_date,
                                      'depmap_date': depmap_date,
                                      'sd_range': sd_range,
                                      'graph_type': _type,
                                      **kwargs.get('info', {})},
                                script_settings=script_settings,
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


def _get_pairs(corr_z: pd.DataFrame):
    # Kudos to https://stackoverflow.com/a/45631406/10478812
    return corr_z.mask(
        np.triu(np.ones(corr_z.shape)).astype(bool)
    ).notna().sum().sum()


def _down_sample_df(z_corr: pd.DataFrame, sample_size: int) -> pd.DataFrame:
    # Do small initial downsampling
    row_samples = len(z_corr) - 1
    n_pairs = _get_pairs(corr_z=z_corr)
    while n_pairs > int(1.1 * sample_size):
        logger.info(f'Down sampling from {n_pairs}')
        z_corr = z_corr.sample(row_samples, axis=0)
        z_corr = z_corr.filter(list(z_corr.index), axis=1)

        # Update n_pairs and row_samples
        n_pairs = _get_pairs(corr_z=z_corr)
        mm = max(row_samples - int(np.ceil(0.05 * row_samples))
                 if n_pairs - sample_size < np.ceil(0.1 * sample_size)
                 else down_sampl_size(n_pairs, len(z_corr), sample_size),
                 1)
        row_samples = row_samples - 1 if mm >= row_samples else mm
    return z_corr


def success_callback(res):
    logger.info('Appending a result')
    output_list.append(res)


def error_callback(err):
    logger.error(f'The following exception occurred in process '
                 f'{mp.current_process().pid} (print of traceback):')
    logger.exception(err)


def main(indra_net: Union[nx.DiGraph, nx.MultiDiGraph],
         outname: str,
         graph_type: str,
         sd_range: Tuple[float, Union[None, float]],
         random: bool = False,
         z_score: Optional[Union[pd.DataFrame, str]] = None,
         z_score_file: Optional[str] = None,
         raw_data: Optional[List[str]] = None,
         raw_corr: Optional[List[str]] = None,
         expl_funcs: Optional[List[str]] = None,
         pb_node_mapping: Optional[Dict[str, Set]] = None,
         n_chunks: Optional[int] = 256,
         expl_mapping: Optional[Dict[str, str]] = None,
         is_a_part_of: Optional[List[str]] = None,
         immediate_only: Optional[bool] = False,
         strict_intermediates: Optional[bool] = False,
         allowed_ns: Optional[List[str]] = None,
         allowed_sources: Optional[List[str]] = None,
         info: Optional[Dict[Hashable, Any]] = None,
         indra_date: Optional[str] = None,
         indra_net_file: Optional[str] = None,
         depmap_date: Optional[str] = None,
         sample_size: Optional[int] = None,
         shuffle: Optional[bool] = False,
         overwrite: Optional[bool] = False,
         normalize_names: Optional[bool] = False,
         argparse_dict: Optional[Dict[str, Union[str, float, int,
                                                 List[str]]]] = None):
    """Set up correlation matching of depmap data with an indranet graph

    Parameters
    ----------
    indra_net : Union[nx.DiGraph, nx.MultiDiGraph]
        The graph representation of the indra network. Each edge should
        have an attribute named 'statements' containing a list of sources
        supporting that edge. If signed search, indranet is expected to be an
        nx.MultiDiGraph with edges keyed by (gene, gene, sign) tuples.
    outname : str
        A file path (can be an S3 url) to where to store the final pickle
        file containing the DepmapExplainer
    graph_type : str
        The graph type of the graph used for the explanations. Can be one of
        'unsigned', 'signed', 'pybel'.
    sd_range : Tuple[float, Union[float, None]]
        A tuple of the lower and optionally the upper bound of the z-score
        range to use when getting correlations
    random : bool
        Whether to do a random sampling or not. If True do a random sample
        instead of cutting the correlations of to the given SD range.
    z_score : Union[pd.DataFrame, str]
        The correlation matrix as a pandas DataFrame or a path to the
        dataframe.
    z_score_file : str
        The path to the correlation DataFrame. This argument is only used
        for storing the info about the location. To provide the location for
        loading the DataFrame, use `z_score`.
    raw_data : Optional[List[str]]
        File paths to CRISPR raw data and RNAi raw data from the DepMap Portal
    raw_corr : Optional[List[str]]
        File paths to raw correlation data (before z-score conversion)
        containing hdf compressed correlation data. These files contain the
        result of running `raw_df.corr()`.
    expl_funcs : Optional[Iterable[str]]
        Provide a list of explanation functions to apply. Default: All
        functions are applied.
    pb_node_mapping : Optional[Union[Dict, Set[Any]]]
        If graph type is "pybel", use this argument to provide a mapping
        from HGNC symbols to pybel nodes in the pybel model
    n_chunks : Optional[int]
        How many chunks to split the data into in the multiprocessing part
        of the script
    expl_mapping : Optional[Dict[str, str]]
        A mapping from entity to a string with information about the entity
        that explains why it is considered 'explained'
    is_a_part_of : Optional[Iterable]
        A set of identifiers to look for when applying the common parent
        explanation between a pair of correlating nodes.
    immediate_only : Optional[bool]
        Only look for immediate parents. This option might limit the number
        of results that are returned. Default: False.
    strict_intermediates : Optional[bool]
        If True: return explanation only when there is a set intersection of
        nodes up- or downstream of A, B for shared regulators and shared
        targets, otherwise return if there are successors for either of A or
        B. Default: False.
    allowed_ns : Optional[List[str]]
        A list of allowed name spaces for explanations involving
        intermediary nodes. Default: Any namespace.
    allowed_sources : Optional[List[str]]
        The allowed sources for edges. This will not affect subsequent edges
        in explanations involving 2 or more edges. Default: all sources are
        allowed.
    info : Optional[Dict[Hashable, Any]]
        An optional dict in which to save meta data about this run
    indra_date : Optional[str]
        The date of the sif dump used to create the graph
    indra_net_file : Optional[str]
        The file path to the graph used
    depmap_date : Optional[str]
        The date (usually a quarter e.g. 19Q4) the depmap data was published
        on depmap.org
    sample_size : Optional[int]
        Number of correlation pairs to approximately get out of the
        correlation matrix after down sampling it
    shuffle : Optional[bool]
        If True, shuffle the correlation matrix. This is good to do in case
        the input data have some sort of structure that could lead to large
        discrepancies between compute times for the different processes.
        Default: False.
    overwrite : Optional[bool]
        If True, overwrite any output files. Default: False.
    normalize_names : Optional[bool]
        If True, try to normalize the names in the correlation matrix that
        are not found in the provided graph. Default: False.
    argparse_dict : Optional[Dict[str, Union[str, float, int, List[str]]]]
        Provide the argparse options from running this file as a script
    """
    global indranet, hgnc_node_mapping, output_list
    indranet = indra_net

    assert expl_funcs is None or isinstance(expl_funcs, (list, tuple, set))
    run_options = {'n-chunks': n_chunks, 'expl_funcs': expl_funcs or None}

    # 1 Check options
    sd_l, sd_u = sd_range if sd_range and len(sd_range) == 2 else \
        ((sd_range[0], None) if sd_range and len(sd_range) == 1 else
         (None, None))

    if not random and not sd_l and not sd_u:
        raise ValueError('Must specify at least a lower bound for the SD '
                         'range or flag run for random explanation')

    if graph_type == 'pybel' and not pb_node_mapping:
        raise ValueError('Must provide PyBEL node mapping with option '
                         'pb_node_mapping if graph type is "pybel"')

    # Get the z-score corr matrix
    if not (bool(z_score is not None and len(z_score)) ^
            bool(raw_data or raw_corr)):
        raise ValueError('Must provide either z_score XOR either of raw_data '
                         'or raw_corr')

    if expl_mapping:
        run_options['expl_mapping'] = expl_mapping
        logger.info(f'Using explained set with '
                    f'{len(expl_mapping)} genes')

    outname = outname if outname.endswith('.pkl') else \
        outname + '.pkl'
    if not overwrite:
        if outname.startswith('s3://'):
            s3 = get_s3_client(unsigned=False)
            if S3Path.from_string(outname).exists(s3):
                raise FileExistsError(f'File {str(outname)} already exists!')
        elif Path(outname).is_file():
            raise FileExistsError(f'File {str(outname)} already exists!')

    if z_score is not None:
        if isinstance(z_score, str) and Path(z_score).is_file():
            z_corr = pd.read_hdf(z_score)
        elif isinstance(z_score, pd.DataFrame):
            z_corr = z_score
        else:
            raise ValueError(f'Unknown type {z_score.__class__} of provided '
                             f'correlation matrix.')
    else:
        z_sc_options = {
            'crispr_raw': raw_data[0],
            'rnai_raw': raw_data[1],
            'crispr_corr': raw_corr[0],
            'rnai_corr': raw_corr[1]
        }
        z_corr = run_corr_merge(**z_sc_options)

    run_options['graph_type'] = graph_type
    # Add optional options
    run_options['immediate_only'] = immediate_only
    run_options['strict_intermediates'] = strict_intermediates
    if allowed_ns:
        run_options['allowed_ns'] = allowed_ns
    if allowed_sources:
        run_options['allowed_sources'] = allowed_sources
    if is_a_part_of:
        run_options['is_a_part_of'] = is_a_part_of

    # Get mapping of correlation names to pybel nodes
    if graph_type == 'pybel':
        if isinstance(pb_node_mapping, dict):
            hgnc_node_mapping = pb_node_mapping
        elif isinstance(pb_node_mapping, str) and \
                Path(pb_node_mapping).is_file():
            hgnc_node_mapping = file_opener(pb_node_mapping)
        else:
            raise ValueError('Could not load pybel node mapping')

    # 2. Filter to SD range OR run random sampling
    if random:
        logger.info('Doing random sampling through df.sample')
        z_corr = z_corr.sample(142, axis=0)
        z_corr = z_corr.filter(list(z_corr.index), axis=1)
        # Remove correlation values to not confuse with real data
        z_corr.loc[:, :] = 0
    else:
        if sd_l and sd_u:
            logger.info(f'Filtering correlations to {sd_l} - {sd_u} SD')
            z_corr = z_corr[((z_corr > sd_l) & (z_corr < sd_u)) |
                            ((z_corr < -sd_l) & (z_corr > -sd_u))]
        elif isinstance(sd_l, (int, float)) and sd_l and not sd_u:
            logger.info(f'Filtering correlations to {sd_l}+ SD')
            z_corr = z_corr[(z_corr > sd_l) | (z_corr < -sd_l)]

    run_options['sd_range'] = (sd_l, sd_u) if sd_u else (sd_l, None)

    # Pick a sample
    if sample_size is not None and not random:
        logger.info(f'Reducing correlation matrix to a random approximately '
                    f'{sample_size} correlation pairs.')
        z_corr = _down_sample_df(z_corr, sample_size)

    # Shuffle corr matrix without removing items
    elif shuffle and not random:
        logger.info('Shuffling correlation matrix...')
        z_corr = z_corr.sample(frac=1, axis=0)
        z_corr = z_corr.filter(list(z_corr.index), axis=1)

    if normalize_names:
        logger.info('Normalizing correlation matrix column names')
        z_corr = normalize_corr_names(z_corr, indra_net)
    else:
        logger.info('Leaving correlation matrix column names as is')

    run_options['corr_z'] = z_corr

    # 4. Add meta data
    info_dict = {}
    if info:
        info_dict['info'] = info
    if depmap_date:
        info_dict['depmap_date'] = depmap_date
    if indra_date:
        run_options['indra_date'] = indra_date
    run_options['info'] = info_dict

    # Set the script_settings
    run_options['script_settings'] = {
        'raw_data': raw_data,
        'raw_corr': raw_corr,
        'z_score': z_score if isinstance(z_score, str) else (z_score_file or
                                                             'no info'),
        'random': random,
        'indranet': indra_net_file or 'no info',
        'shuffle': shuffle,
        'sample_size': sample_size,
        'n_chunks': n_chunks,
        'outname': outname,
        'apriori_explained': apriori_explained if isinstance(
            apriori_explained, str) else 'no info',
        'graph_type': graph_type,
        'pybel_node_mapping': pb_node_mapping if isinstance(
            pb_node_mapping, str) else 'no info',
        'argparse_info': argparse_dict
    }

    # Create output list in global scope
    output_list = []
    explanations = match_correlations(**run_options)
    if outname.startswith('s3://'):
        try:
            logger.info(f'Uploading results to s3: {outname}')
            s3 = get_s3_client(unsigned=False)
            s3outpath = S3Path.from_string(outname)
            s3outpath.upload(s3=s3, body=pickle.dumps(explanations))
            logger.info('Finished uploading results to s3')
        except Exception:
            new_path = Path(outname.replace('s3://', ''))
            logger.warning(f'Something went wrong in s3 upload, trying to '
                           f'save locally instead to {new_path}')
            new_path.parent.mkdir(parents=True, exist_ok=True)
            dump_it_to_pickle(fname=new_path.absolute().resolve().as_posix(),
                              pyobj=explanations, overwrite=overwrite)

    else:
        # mkdir in case it doesn't exist
        outpath = Path(outname)
        logger.info(f'Dumping results to {outpath}')
        outpath.parent.mkdir(parents=True, exist_ok=True)
        dump_it_to_pickle(fname=outpath.absolute().resolve().as_posix(),
                          pyobj=explanations, overwrite=overwrite)
    logger.info('Script finished')


if __name__ == '__main__':
    parser = argparse.ArgumentParser('DepMap Explainer main script')
    #   1a Load depmap data from scratch | load crispr/rnai raw corr | z-score
    corr_group = parser.add_mutually_exclusive_group(required=True)
    range_group = parser.add_mutually_exclusive_group(required=True)
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
             'containing z-scores for gene-gene correlations.')
    corr_group.add_argument(
        '--raw-drugs', nargs=2, type=file_path(),
        help='File paths to the raw drug data from the DepMap portal\'s '
             'PRISM Repurposing data. The file names should match '
             '`primary-screen-replicate-collapsed-logfold-change.csv` and '
             '`primary-screen-replicate-collapsed-treatment-info.csv`'
    )

    #   1b Load indranet
    parser.add_argument(
        '--indranet', type=file_path(), required=True,
        help='The indra network to use for explanations. Should be either a '
             'DiGraph or signed DiGraph (a MultiDiGraph with max two edges '
             'per node pair, one edge per sign).'
    )

    #   1c-1 Optionally provide PyBEL model
    parser.add_argument(
        '--pybel-model', type=file_path(),
        help='If graph type is "pybel", use this argument to provide the '
             'necessary pybel model used to precompute the pybel node mapping.'
    )

    #   1c-2 Optionally provide node mapping for hgnc symbol - pybel node
    parser.add_argument(
        '--pybel-node-mapping', type=file_path(),
        help='If graph type is "pybel", use this argument to optionally '
             'provide a mapping from HGNC symbols to pybel nodes in the '
             'pybel model'
    )

    #   1d Provide graph type
    allowed_graph_types = {'unsigned', 'signed', 'pybel'}
    parser.add_argument(
        '--graph-type', type=allowed_types(allowed_graph_types),
        default='unsigned',
        help=f'Specify the graph type used. Allowed values are {allowed_types}'
    )

    #   1e Provide allowed_ns
    parser.add_argument(
        '--allowed-ns', nargs='+',
        help='Specify the allowed namespaces to be used in the graph for '
             'intermediate nodes. Default: all namespaces are allowed.'
    )

    #   1f Provide sources to filter to for edges
    parser.add_argument(
        '--allowed-sources', nargs='+',
        help='Specify the allowed sources for edges. This will not affect '
             'subsequent edges in explanations involving 2 or more edges. '
             'Default: all sources are allowed.'
    )

    #   1g Provide expl function names
    parser.add_argument(
        '--expl-funcs', nargs='+',
        type=allowed_types(set(expl_functions.keys())),
        help='Provide explainer function names to be used in the explanation '
             'loop'
    )

    #   2a. Filter to SD range
    range_group.add_argument('--sd-range', nargs='+', type=float,
                             help='SD range to filter to')
    # OR
    #   2b. do random sampling from correlation matrix genes
    range_group.add_argument('--random', action='store_true',
                             help='Check the explanation rate for randomly '
                                  'sampled pairs of genes from the full '
                                  'correlation matrix')

    #   3. A-priori explained entities
    parser.add_argument(
        '--apriori-explained', type=file_path(), nargs='?',
        const=mito_file,
        help='A-priori explained entities must be in a file that can be '
             'parsed as CSV/TSV with column names "name" for entity name and '
             '"description" for explanation why the entity is explained. '
             'This argument can be used as a flag as well: by only providing '
             '`--apriori-explained` (without any value) the default resource '
             'file containing MitoCarta 3.0 data is used.')

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

    # Sampling
    parser.add_argument(
        '--sample-size', type=int,
        help='If provided, down sample the correlation matrix so this many '
             'pairs (approximately) are picked at random.'
    )

    # 6 Extra info and options
    parser.add_argument('--indra-date',
                        help='Provide a date for the dump from which the '
                             'indra network was created')
    parser.add_argument('--depmap-date',
                        help='Provide the release date of the depmap data '
                             'used.')
    parser.add_argument('--shuffle', action='store_true',
                        help='Shuffle the correlation matrix before running '
                             'matching loop.')
    parser.add_argument('--is-a-part-of', nargs='+',
                        help='Identifiers that are considered to explain '
                             'pair connections in common parent search in '
                             'ontology.')
    parser.add_argument('--immediate-only', action='store_true',
                        help='Only look in immediate parents in common '
                             'parent search.')
    parser.add_argument('--strict-intermediates', action='store_true',
                        help='For shared target and shared regulators: only '
                             'return explanation if there are shared nodes '
                             'up- or downstream of A-B pair.')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite any output files that already exist.')
    parser.add_argument('--normalize-names', action='store_true',
                        help='Try to normalize the names in the correlation '
                             'matrix if they are not found in the provided '
                             'graph')

    args = parser.parse_args()
    arg_dict = vars(args)
    arg_dict['argparse_dict'] = deepcopy(arg_dict)

    # Load z_corr, indranet and optionally pybel_model
    inet_graph = file_opener(args.indranet)
    arg_dict['indra_net'] = inet_graph
    arg_dict['indra_net_file'] = args.indranet
    if arg_dict.get('z_score'):
        corr_matrix = pd.read_hdf(arg_dict['z_score'])
        arg_dict['z_score_file'] = arg_dict['z_score']
    elif arg_dict.get('raw_drugs'):
        raw_drug_data, raw_drug_info = arg_dict['raw_drugs']
        corr_matrix = drugs_to_corr_matrix(raw_drug_data, raw_drug_info)
    else:
        arg_dict['raw_data'] = arg_dict.get('raw_data')
        arg_dict['corr_data'] = arg_dict.get('corr_data')
        z_sc_options = {
            'crispr_raw': arg_dict['raw_data'][0],
            'rnai_raw': arg_dict['raw_data'][1],
            'crispr_corr': arg_dict['corr_data'][0],
            'rnai_corr': arg_dict['corr_data'][1]
        }
        corr_matrix = run_corr_merge(**z_sc_options)

    arg_dict['z_score'] = corr_matrix
    hgnc_names = corr_matrix.columns.values
    # Get hgnc node mapping
    if arg_dict.get('graph_type') == 'pybel' and \
            not arg_dict.get('pybel_node_mapping') and \
            not arg_dict.get('pybel_model'):
        raise ValueError('Must provide PyBEL model with option pybel_model '
                         'or provide node mapping with option '
                         'if graph type is pybel')
    # Only model provided: create mapping
    if arg_dict.get('pybel_model') and not arg_dict.get('pybel_node_mapping'):
        mapping = pybel_node_name_mapping(
            node_names=hgnc_names, node_ns='HGNC',
            pb_model=file_opener(arg_dict['pybel_model'])
        )
        arg_dict['pb_node_mapping'] = mapping
    # Mapping is provided: load the mapping
    elif arg_dict.get('pybel_node_mapping'):
        if arg_dict['pybel_node_mapping'].endswith('.pkl'):
            arg_dict['pb_node_mapping'] = \
                file_opener(arg_dict['pybel_node_mapping'])
        elif arg_dict['pybel_node_mapping'].endswith('.json'):
            arg_dict['pb_node_mapping'] = \
                file_opener(arg_dict['pybel_node_mapping'])
        else:
            raise ValueError('Unknown file type %s' %
                             arg_dict['pybel_node_mapping'].split('.')[-1])

    # Get ignore list
    apriori_explained = args.apriori_explained
    if apriori_explained:
        # Check if it's the default file
        if mito_file_name in apriori_explained:
            expl_map = get_mitocarta_info(apriori_explained)
        else:
            # Hope it's a csv
            try:
                expl_df = pd.read_csv(apriori_explained)
                expl_map = {e: d for e, d in zip(expl_df.name,
                                                 expl_df.description)}
            except Exception as err:
                raise ValueError('A-priori explained entities must be in a '
                                 'file that can be parsed as CSV/TSV with '
                                 'column names "name" for entity name and '
                                 '"description" for explanation why the '
                                 'entity is explained.') \
                    from err
        arg_dict['expl_mapping'] = expl_map

    main_keys = inspect.signature(main).parameters.keys()
    kwargs = {k: v for k, v in arg_dict.items() if k in main_keys}

    main(**kwargs)
