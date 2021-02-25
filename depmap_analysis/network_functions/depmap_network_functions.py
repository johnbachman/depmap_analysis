import os
import re
import sys
import math
import logging
import itertools as itt
from random import choices
from math import ceil, log10
from typing import Iterable
from collections import Mapping, OrderedDict, defaultdict

import numpy as np
import pandas as pd
import networkx as nx
from scipy.special import gamma, hyp2f1
from sqlalchemy.exc import StatementError
from scipy import interpolate as interpol
from scipy.optimize import curve_fit as opt_curve_fit
from pandas.core.series import Series as pd_Series_class

from indra_db import util as dbu
from indra_db import client as dbc
from indra.tools import assemble_corpus as ac
from indra.statements import Statement
from indra.sources.indra_db_rest import api as db_api
from indra.sources.indra_db_rest.exceptions import IndraDBRestAPIError
import depmap_analysis.util.io_functions as io
import depmap_analysis.network_functions.famplex_functions as ff

db_prim = dbu.get_primary_db()
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def entry_exist_dict(nest_dict, outer_key, inner_key):
    if nest_dict.get(outer_key) and nest_dict.get(outer_key).get(inner_key):
        return True
    else:
        return False


def _entry_exist_corr_matrix(corr, val1, val2):
    try:
        _ = corr[np.in1d(corr.index.get_level_values(0), val1)][val2]
        return True
    except KeyError:
        return False


def create_nested_dict():
    """Returns a nested dictionary of arbitrary depth

    Returns
    -------
    defaultdict(create_nested_dict)
    """
    return defaultdict(create_nested_dict)


def mean_z_score(mu1, sig1, c1, mu2, sig2, c2):
    return 0.5 * _z_sc(num=c1, mu=mu1, sigma=sig1) + \
        0.5 * _z_sc(num=c2, mu=mu2, sigma=sig2)


def comb_z_sc_gen(crispr_corr, rnai_corr, stats_dict):
    """Generate random samples of combined_z_sc from the input matrices

    Parameters
    crispr_corr : pd.DataFrame
        Correlation matrix with gene-gene-correlation score
    rnai_corr : pd.DataFrame
        Correlation matrix with gene-gene-correlation score
    stats_dict : dict
        Nested dictionary containing the mean and standard deviations of the
        input correlation matrices. Expected structure:
            {
                'crispr': {
                    'mu': <mean>,
                    'sigma': <standard deviation>
                },
                'rnai': {
                    'mu': <mean>,
                    'sigma': <standard deviation>
                }
            }

    Yields
    ------
    float
        A random sample of the combined z-score from the two correlations
        matrices
    """
    cmu = stats_dict['crispr']['mu']
    csig = stats_dict['crispr']['sigma']
    rmu = stats_dict['rnai']['mu']
    rsig = stats_dict['crispr']['sigma']
    comb_gene_tuples = tuple(
        set(crispr_corr.columns.values).intersection(rnai_corr.columns.values)
    )
    while True:
        g1, g2 = choices(comb_gene_tuples, k=2)
        while g1 == g2:
            g1, g2 = choices(comb_gene_tuples, k=2)
        cc = crispr_corr.loc[g1, g2]
        rc = rnai_corr.loc[g1, g2]
        yield mean_z_score(mu1=cmu, sig1=csig, c1=cc,
                           mu2=rmu, sig2=rsig, c2=rc)


def csv_file_to_generator(fname, column_list):
    # todo consider passing a file stream object that can be read line by line
    """Return a tuple generator given a csv file with
    gene-gene-correlation pairs and specified columns

    fname : str
        File path of csv file
    column_list : list[str]
        List of column names

    Returns
    -------
    tuple_generator : generator object
         A generator that returns a tuple of each row
    """
    pair_corr_file = pd.read_csv(fname, names=column_list)
    return (tuple(line[1]) for line in pair_corr_file.iterrows())


def corr_matrix_to_generator(corrrelation_df_matrix, max_pairs=None):
    """Return a tuple generator given a correlation matrix
    
    The function takes a correlation matrix and returns a consumable tuple 
    generator object. Once consumed, the object is exhausted and a new
    generator needs to be produced.

    corrrelation_df_matrix : pandas.DataFrame
        A correlation matrix as a pandas dataframe

    Returns
    -------
    tuple_generator : generator object
        A generator that returns a tuple of each row
    """
    # Sample at random: get a random sample of the correlation matrix that has
    # enough non-nan values to exhaustively generate at least max_pair
    all_pairs = corrrelation_df_matrix.notna().sum().sum()
    if all_pairs == 0:
        logger.warning('Correlation matrix is empty')
        raise ValueError('Script aborted due to empty correlation matrix')

    corr_df_sample = pd.DataFrame()
    if max_pairs:
        if max_pairs >= all_pairs:
            logger.info('The requested number of correlation pairs is '
                            'larger than the available number of pairs. '
                            'Resetting `max_pairs` to %i' % all_pairs)
            corr_df_sample = corrrelation_df_matrix

        elif max_pairs < all_pairs:
            n = int(np.floor(np.sqrt(max_pairs))/2 - 1)
            corr_df_sample = corrrelation_df_matrix.sample(
                n, axis=0).sample(n, axis=1)

            # Increase sample until number of extractable pairs exceed
            # max_pairs
            while corr_df_sample.notna().sum().sum() <= max_pairs:
                n += 1
                corr_df_sample = corrrelation_df_matrix.sample(
                    n, axis=0).sample(n, axis=1)

        logger.info('Created a random sample of the correlation matrix '
                        'with %i extractable correlation pairs.'
                    % corr_df_sample.notna().sum().sum())

    # max_pairs == None: no sampling, get all non-NaN correlations;
    else:
        corr_df_sample = corrrelation_df_matrix

    corr_value_matrix = corr_df_sample.values
    gene_name_array = corr_df_sample.index.values
    if isinstance(gene_name_array[0], tuple):
        gene_name_array = [n[0] for n in gene_name_array]
    tr_up_indices = np.triu_indices(n=len(corr_value_matrix), k=1)
    # Only get HGNC symbols (first in tuple) since we're gonna compare to
    # INDRA statements, which is currently done with HGNC symbols
    # todo change to output HGNC IDs instead when switching to HGNC id in
    #  stmts
    return ((gene_name_array[i], gene_name_array[j],
            float(corr_value_matrix[i, j]))
            for i, j in zip(*tr_up_indices)
            if not np.isnan(corr_value_matrix[i, j]))


def _dump_master_corr_dict_to_pairs_in_csv(fname, nest_dict):
    logger.info('Dumping master dict to pairs in %s' % fname)
    pairs = 0
    with open(fname, 'w') as fo:
        fo.write('gene1,gene2,correlation_list\n')
        for ok, od in nest_dict.items():
            for ik, id_ in od.items():
                vals = []
                for name in id_:
                    vals.append(id_[name])
                fo.write('%s,%s,"%s"\n' % (ok, ik, str(vals)))
                pairs += 1
    return pairs


def nx_undir_to_neighbor_lookup_json(expl_undir_graph, outbasename,
                                     neighbor_dir='neighbor_lookup/'):
    """Dumps one json dict per node in a undirected nx graph where entry is a
    json array containing all neighbors

    expl_undir_graph : nx.Graph
    neighbor_dir : directory to put jsons in. Must be specified or the default
    'neighbor_lookup/' is used.

    """
    path = '/'.join(
        outbasename.replace('//', '/').split('/')[:-1]
    ) + '/' + neighbor_dir
    if '/' == path[0]:
        path = path[1:]
    if not os.path.isdir(path):
        logger.info('Could not find path "%s", creating new directory.' %
                    path)
        os.makedirs(path)

    logger.info('Dumping node neighbor dicts to "%s'
                    '/neighbors_to_NODENAME.json"' % path)
    for node in expl_undir_graph.nodes:
        nnnl = []
        for other_node in expl_undir_graph[node]:
            inner_dict = expl_undir_graph[node][other_node]
            nnnl.append([other_node, inner_dict['attr_dict']['correlation']])
        io.dump_it_to_json(fname=path + '/neighbors_to_%s.json' % node,
                        pyobj=nnnl)
    logger.info('Finished dumping node neighbor dicts to %s' % path)


def nested_stmt_dict_to_nx_multidigraph(nest_d):
    """Returns a directed multigraph where each edge links a statement with
    u=subj, v=obj, edge_key=stmt,

    nest_d : defaultdict(dict)
        Nested dict of statements: nest_d[subj][obj]

    Returns
    -------
    nx_muldigraph : nx.MultiDiGraph
        An nx directed multigraph linking agents with statements
    """

    logger.info('Building directed multigraph from nested dict of '
                    'statements')
    nx_muldir = nx.MultiDiGraph()

    for subj in nest_d:
        if nest_d.get(subj):
            for obj in nest_d[subj]:
                # Check if subj-obj connection exists in dict
                if subj is not obj and nest_d.get(subj).get(obj):
                    # Get list of statements
                    dds_list = nest_d[subj][obj]
                    for stmt in dds_list:
                        # One edge per statement
                        # Could instead add stmt attributes like
                        # evidence.text, supported by, suppors, uuic, etc
                        nx_muldir.add_edge(u_for_edge=subj, v_for_edge=obj,
                                           attr_dict={'stmt': stmt})
    return nx_muldir


def _merge_belief(sif_df, belief_dict):
    if isinstance(belief_dict, str):
        belief_dict = io.file_opener(belief_dict)
    elif isinstance(belief_dict, dict):
        belief_dict = belief_dict

    hashes = []
    beliefs = []
    for k, v in belief_dict.items():
        hashes.append(k)
        beliefs.append(v)

    sif_df = sif_df.merge(
        right=pd.DataFrame(data={'stmt_hash': hashes, 'belief': beliefs}),
        how='left',
        on='stmt_hash'
    )
    # Check for missing hashes
    if sif_df['belief'].isna().sum() > 0:
        logger.warning('%d rows with missing belief score found' %
                       sif_df['belief'].isna().sum())
        logger.info('Setting missing belief scores to 1/evidence count')
    return sif_df


def sif_dump_df_to_nest_d(sif_df_in, belief_dict=None):
    """Convert a sif dump df to a nested dict

    Paramters
    ---------
    sif_df_in : pd.DataFrame|str
        A pd.DataFrame with at least the columns 'agA_name', 'agB_name',
        'stmt_type', 'stmt_hash'. Any other columns will be part of the
        list of dictionaries in the innermost entry.
    belief_dict : str|dict
        The file path to a pickled dict or a dict object keyed by statement
        hash containing the belief score for the corresponding statements.
        The hashes should correspond to the hashes in the loaded dataframe.

    Returns
    -------
    nest_d : dict(dict(())
    """

    if isinstance(sif_df_in, str):
        sif_df = io.file_opener(sif_df_in)
    elif isinstance(sif_df_in, pd.DataFrame):
        sif_df = sif_df_in
    else:
        raise ValueError('sif_df_in must be of type str or pd.DataFrame')

    mandatory_columns = ('agA_name', 'agB_name', 'stmt_type', 'stmt_hash')
    if not set(mandatory_columns).issubset(set(sif_df.columns.values)):
        raise ValueError('Data frame must at least have the columns %s' %
                         '"' + '", "'.join(mandatory_columns) + '"')

    if belief_dict:
        sif_df = _merge_belief(sif_df, belief_dict)

    logger.info('Producing nested dict of stmts from sif dump df')
    nest_d = {}
    for index, row in sif_df.iterrows():
        rd = row.to_dict()
        a_name = rd.pop('agA_name')
        b_name = rd.pop('agB_name')

        if nest_d.get(a_name):
            if nest_d[a_name].get(b_name):
                # a-b relation already exists and should be a list
                try:
                    nest_d[a_name][b_name].append(rd)
                except KeyError:
                    logger.error('KeyError when trying to append new '
                                     'statement for relation %s %s relation'
                                 % (a_name, b_name))
                except AttributeError:
                    logger.error('AttributeError for %s %s' %
                                 (a_name, b_name))
            else:
                # First a-b relation, add inner dict as list
                nest_d[a_name][b_name] = [rd]
        else:
            # First addition of this a_name agent
            nest_d[a_name] = {b_name: [rd]}
    return nest_d


def nested_stmt_dict_to_nx_digraph(nest_d, belief_dict=None):
    """Returns a directed graph from a two layered nested dictionary

    Nested dictionary

        d[subj][obj] = [stmts/stmt hashes]

    nest_d : defaultdict(dict)
        Nested dict of statements: nest_d[subj][obj]
    belief_dict : dict()
        dict with {hash: belief score} as key: value pairs

    Returns
    -------
    nx_digraph : nx.DiGraph
        An nx directed graph with statements
    """

    logger.info('Building directed simple graph from two layered nested '
                    'dict.')
    nx_dir_g = nx.DiGraph()

    if not belief_dict:
        logger.info('No Belief Score dictionary provided')
        logger.warning('API belief score checkup is not implemented yet')
        pass  # ToDo Download from sif dump on s3

    else:
        if isinstance(next(iter(belief_dict)), str):
            k_int = False
        elif isinstance(next(iter(belief_dict)), int):
            k_int = True
        else:
            logger.warning('Belief dict seems to have keys that are not '
                               'str or int.')
            k_int = None

    # Flag to check if the statement dict has belief score in it or not
    has_belief = False
    for k1, d1 in nest_d.items():
        for k2, v in d1.items():
            for tups in v:
                # Has to be (type, hash) or (type, hash, belief score)
                assert len(tups) == 2 or len(tups) == 3
                if len(tups) == 3:
                    has_belief = True
                break

    error_count = 0
    for subj in nest_d:
        if nest_d.get(subj):
            for obj in nest_d[subj]:
                # Check if subj-obj connection exists in dict
                if subj is not obj and nest_d.get(subj).get(obj):
                    # Add edge
                    inner_obj = nest_d[subj][obj]
                    inner_obj_b = []
                    if not has_belief:
                        for typ, hsh in inner_obj:
                            try:
                                bs = belief_dict[int(hsh)] if k_int else \
                                    belief_dict[str(hsh)]
                            except KeyError:
                                bs = 0.1  # Todo get min belief for any stmt
                                error_count += 1
                                if error_count < 10:
                                    logger.warning('No entry found in '
                                                       'belief dict for hash '
                                                       '%s' % str(hsh))
                                if error_count == 10:
                                    logger.warning('Muting belief dict '
                                                       'KeyError...')
                            t = (typ, hsh, bs)
                            inner_obj_b.append(t)
                    else:
                        inner_obj_b = inner_obj
                    nx_dir_g.add_edge(u_of_edge=subj,
                                      v_of_edge=obj,
                                      attr_dict={'stmts': inner_obj_b})
    if error_count > 100:
        logger.warning('%i hashes did not have a mathcing key in belief '
                           'dict' % error_count)
    return nx_dir_g


def nested_stmt_explained_dict_nx_digraph(nest_d):
    """Returns a directed graph from a three layered nested dictionary

    Form of nested dictionary

        d[subj][obj] = {correlation: float,
                        directed: [stmts/stmt hashes],
                        undirected: [stmts/stmt hashes],
                        x_is_intermediary: [X],
                        x_is_downstream: [X],
                        x_is_upstream: [X]}

    nest_d : defaultdict(dict)
        Nested dict of statements: nest_d[subj][obj] = {attr_dict}

    Returns
    -------
    nx_digraph : nx.DiGraph
        An nx directed graph with statements and/or connecting nodes as edges
    """

    logger.info('Building directed simple graph from three layered nested '
                    'dict.')
    nx_dir_g = nx.DiGraph()

    for subj in nest_d:
        if nest_d.get(subj):
            for obj in nest_d[subj]:
                # Check if subj-obj connection exists in dict
                if subj is not obj and nest_d.get(subj).get(obj):
                    # Add edge
                    inner_obj = nest_d[subj][obj]
                    nx_dir_g.add_edge(u_of_edge=subj,
                                      v_of_edge=obj,
                                      attr_dict=inner_obj)
    return nx_dir_g


def nested_stmt_dict_to_nx_graph(nest_d):
    """Returns an undirected graph built from a nested dict of statements

    Use this function to build a simple undirected graph. Suitable when the
    goal is to generate a node-edge graph for plain visualization.

    nest_d : defaultdict(dict)
        A nested dict with two or more layers

    Returns
    -------
    nx_undir : nx.Graph
        An undirected, unweighted networkx graph
    """

    nx_undir = nx.Graph()

    # Create queue from nested dict
    ndq = list(nest_d.items())

    logger.info('Building undirected graph from nested dict of statements')
    # Run until queue is empty
    while ndq:
        # get node u and dict d from top of queue
        u, d = ndq.pop()
        # Loop nodes nd and (possible) dicts nd of dict d
        for v, nd in d.items():
            # Add edge u-v if it's not a self-loop
            if u is not v:
                nx_undir.add_edge(u, v)
            # If nd has deeper layers, put that to the queue
            if isinstance(nd, Mapping):
                ndq.append((v, nd))

    return nx_undir


def pd_to_nx_graph(corr_sr, source='id1', target='id2',
                   edge_attr='correlation', use_abs_corr=True):
    """Return a graph from a pandas series containing correlaton between gene
    A and gene B, using the correlation as edge weight

    Parameters
    ----------
    corr_sr : pandas Series or DataFrame
        Pandas Series/DataFrame containing A, B corr
    source : str
        which column to identify as source node (output is still undirected)
    target : str
        which column to identify as target node (output is still undirected)
    edge_attr : int or str
        Column to use for edge attributes
    use_abs_corr : Bool
        Use absolute value as edge weight. Otherwise magnitude is used.

    Returns
    -------
    corr_weight_graph : nx.Graph
        An undirected, weighted, networkx graph
    """

    # check if corr_sr is series or dataframe
    if type(corr_sr) == pd_Series_class:
        logger.info('Converting Pandas Series to Pandas DataFrame')
        corr_df = pd.DataFrame(corr_sr).reset_index()
    else:
        corr_df = corr_sr

    corr_df = corr_df.rename(
        columns={'level_0': source, 'level_1': target, 0: edge_attr})

    if use_abs_corr:
        logger.info('Using absolute correlation values')
        corr_df = corr_df.apply(lambda c: c.abs() if np.issubdtype(
            c.dtype, np.number) else c)

    if type(edge_attr) is list or type(edge_attr) is bool:
        logger.warning('More than one attribute might be added to edges. '
                           'Resulting networkx graph might not be usable as '
                           'simple weighted graph.')
    logger.info('Creating weighted undirected graph from network data')
    corr_weight_graph = nx.from_pandas_dataframe(df=corr_df,
                                                 source=source,
                                                 target=target,
                                                 edge_attr=edge_attr)
    return corr_weight_graph


def nx_graph_from_corr_tuple_list(corr_list, use_abs_corr=False):
    """Return a graph from a list of edges, using the correlation as weight

    Parameters
    ----------
    corr_list : list or iterator
        Edge tuples
    use_abs_corr : Bool
        Use absolute value as edge weight. Otherwise magnitude is used.

    Returns
    -------
    corr_weight_graph : nx.Graph
        An undirected, weighted, networkx graph
    """
    corr_weight_graph = nx.Graph()

    if use_abs_corr:
        logger.info('Using absolute correlation values')
        corr_list = map(lambda t: (t[0], t[1], abs(t[2])), corr_list)

    logger.info('Converting tuples to an edge bunch')
    edge_bunch = map(lambda t: (t[0], t[1], {'weight': t[2]}), corr_list)

    logger.info('Creating weighted undirected graph from network data')
    corr_weight_graph.add_edges_from(ebunch_to_add=edge_bunch)

    return corr_weight_graph


def histogram_from_tuple_generator(tuple_gen, binsize, first,
                                   number_of_bins=None):
    """Returns a histogram for large data sets represented as tuple generators

    tuple_gen : generator object
        tuple_generator object that generates A, B, value tuples
    number_of_bins: int
        the number fo bins to use
    binsize : float
        the size of bins
    first : float
        The left most (min(x)) edge of the bin edges

    Returns
    -------
    home_brewed_histo: np.array
        A histrogram of the data in fpath according to number of bins,
        binsize and first.
    """
    if number_of_bins is None:
        number_of_bins = int(2*abs(first) / binsize)
    home_brewed_histo = np.zeros(number_of_bins, dtype=int)
    for g1, g1, flt in tuple_gen:
        home_brewed_histo[io.map2index(start=first, binsize=binsize,
                                       value=flt)] += 1
    return home_brewed_histo


def _my_gauss(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


def _pdf_bivariate_normal(r, rho, n):
    """PDF for the sample correlation coefficient r of a normal bivariate

    See:
    ( 'https://en.wikipedia.org/wiki/Pearson_correlation_coefficient'
      '#Using_the_exact_distribution' )
    and,
    ( 'https://stats.stackexchange.com/questions/191937/what-is-the'
      '-distribution-of-sample-correlation-coefficients-between-two-uncorrel' )

    """
    # RHO is mean of PDF?
    # n is number of cell lines?
    # gamma = scipy/reference/generated/scipy.special.gamma.html
    # hyp2f1 = scipy/reference/generated/scipy.special.hyp2f1.html
    gamma_n_1 = gamma(n-1)  # Gamma(n-1)
    gamma_n_1_2 = gamma(n-0.5)  # Gamma(n-1.2)
    gauss_hyperg = hyp2f1(0.5, 0.5, (2*n-1)/2, (rho*r+1)/1)

    denom = (n-2)*gamma_n_1*(1-rho**2)**((n-1)/2)*(1-r)**((n-4)/2)  # upper
    numer = np.sqrt(2*np.pi)*gamma_n_1_2*(1-rho*r)**(n-3/2)  # lower

    return gauss_hyperg*denom/numer


def get_stats(tuple_generator):
    """Get mean and standard deviation from large file with A,B-value pairs.

    tuple_generator: generator object
        tuple_generator object that generates A, B, value tuples

    Returns
    -------
        mean, standard_deviation
    """
    # See:
    # 'https://math.stackexchange.com/questions/102978/
    # incremental-computation-of-standard-deviation#103025'
    # OR
    # 'https://math.stackexchange.com/questions/374881/
    # recursive-formula-for-variance#375022'
    t1, t2, skip, m = 0, 0, 0, 0
    for m, (gene1, gene2, c) in enumerate(tuple_generator):
        corr = float(c)
        if corr == 1.0:
            skip += 1
            continue  # Ignore self correlations if they exist
        else:
            t1 += corr
            t2 += corr**2

    assert m != 0
    t0 = m + 1 - skip
    mean = t1 / t0
    std = np.sqrt(t0 * t2 - t1**2) / t0
    return mean, std


def get_gaussian_stats(bin_edges, hist):
    """Assuming a gaussian histogram, return scale (a), mean (mu) and sigma.

    We assume the guassian has the form:

        f(x,A,mu,sigma) = A*np.exp(-(x-mu)**2/(2*sigma**2))

    bin_edges: np.array
         edges for the bins as numpy array of length n+1
    hist: np.array
        histogram as numpy array of length n

    Returns
    -------
    a, mu, sigma: tuple(floats)
        a: scaling of the gaussian
        mu: mean of the distribution
        sigma: standard deviation
    """
    a0 = max(hist) / 2.0
    sigma0 = 0.125*(bin_edges[-1]-bin_edges[0])  # (bin distance)/8
    bin_positions = (bin_edges[:-1] + bin_edges[1:]) / 2.0
    mu0 = bin_positions[np.argmax(bin_positions)]  # position of guess max
    coeff, var_matrix = opt_curve_fit(_my_gauss, bin_positions, hist,
                                      [a0, mu0, sigma0])
    a = coeff[0]
    mu = coeff[1]
    sigma = coeff[2]
    return a, mu, sigma


def _get_partial_gaussian_stats(bin_edges, hist):
    bin_positions = (bin_edges[:-1] + bin_edges[1:]) / 2.0
    # Get log of nonzero x, y pairs
    log_hist = []
    saved_positions = []
    for x, y in zip(bin_positions, hist):
        if y > 0:
            log_hist.append(np.log(y))
            saved_positions.append(x)

    log_hist = np.array(log_hist)
    saved_positions = np.array(saved_positions)

    interp_log_gaussian = interpol.interp1d(x=saved_positions,
                                            y=log_hist,
                                            kind='quadratic')
    interp_gaussian = np.exp(interp_log_gaussian(bin_positions))

    return get_gaussian_stats(bin_edges, interp_gaussian)


def get_gene_gene_corr_dict(tuple_generator):
    """Returns a gene-gene correlation nested dict from a gene correlation
    generator

    tuple_generator : generator object
        A generator object returning (gene, gene, correlation) tuples

    Returns
    -------
    nested_dict : defaultdict(dict)
        Dict with gene-gene-correlation
    """
    corr_nest_dict = create_nested_dict()
    logger.info('Generating correlation lookup')
    skip = 0
    doublets = 0
    for count, (gene1, gene2, c) in enumerate(tuple_generator):
        # Self correlations should be filtered at this point but as a backup
        if gene1 == gene2:
            skip += 1
            continue
        else:
            corr = float(c)
            if entry_exist_dict(corr_nest_dict, gene1, gene2):
                doublets += 1
            corr_nest_dict[gene1][gene2] = corr
    count += 1
    logger.info('Created correlation dictionary of length %i, skipped %i, '
                    'found %i doublets' % (count, skip, doublets))
    corr_nest_dict.default_factory = None
    return corr_nest_dict


def merge_correlation_data(correlation_dicts_list, settings):
    """Merge multiple correlation data sets to one single iterable of
    (gene, gene, correlation_dict)

    correlation_dicts_list: list[(gene_set_name, dataset_dict, sigma_dict)]
        List of name-dataset_dict-distr_stats_dict tuples
    settings:

    Returns
    -------
    merged_corr_dict : defaultdict(dict)
        A master correlation dict containing the set of (A,B) correlations
        that exist in all sets. The structure is d[g][g] = {set_name: corr, ...}
    """

    # Since we're merging to an intersection, we only need to loop one of the
    # dictionaries and compare it to the other. Check which one is the
    # smallest and loop that one (assuming len(dict) gives the nested dict
    # with the shortest iteration).

    # We have three cases:
    # 1.    Only one data set
    # 2.    Two data sets
    # 3.    Three or more data sets (not yet implemented)
    #
    # For 1: Skip the merger and just return whatever data set came in.
    # For 2: Merge the two dictionaries.
    # For 3: Merge two of the data sets and then call a recursive

    # Create the return dict
    merged_corr_dict = create_nested_dict()

    # Case 1
    if len(correlation_dicts_list) == 1:
        only_name, only_dict, only_sigma_dict = correlation_dicts_list.pop()

        npairs = 0
        for o_gene, d in only_dict.items():
            for i_gene, corr in d.items():
                if o_gene is not i_gene:
                    merged_corr_dict[o_gene][i_gene][only_name] = corr
                    npairs += 1

    # Case 2
    elif len(correlation_dicts_list) == 2:
        # Get shortest dict from tuple list
        name_dict_sigma_tuple = next(d for d in correlation_dicts_list if len(
            d[1]) == min([len(i[1]) for i in correlation_dicts_list]))
        # Remove the tuple with the shortest dict  from the list
        correlation_dicts_list.remove(name_dict_sigma_tuple)
        set_name, shortest_dict, sigma_dict = name_dict_sigma_tuple
        # Pop the list to get the oter tuple to merge with
        other_name, other_dict, other_sigma_dict = correlation_dicts_list.pop()
        logger.info('Merging correlation dicts %s and %s' %
                    (set_name, other_name))

        # Loop shortest correlation lookup dict
        npairs = 0
        for o_gene, d in shortest_dict.items():
            for i_gene, corr in d.items():
                if o_gene is not i_gene and \
                        not entry_exist_dict(merged_corr_dict, o_gene, i_gene):
                    # Check both directions
                    other_corr = None
                    if entry_exist_dict(other_dict, o_gene, i_gene):
                        other_corr = other_dict[o_gene][i_gene]
                    elif entry_exist_dict(other_dict, i_gene, o_gene):
                        other_corr = other_dict[i_gene][o_gene]

                    if other_corr and pass_filter(
                            corr1=corr, mu1=sigma_dict['mean'],
                            sigma1=sigma_dict['sigma'],
                            corr2=other_corr, mu2=other_sigma_dict['mean'],
                            sigma2=other_sigma_dict['sigma'],
                            margin=settings['margin'],
                            filter_type=settings['filter_type']):
                        merged_corr_dict[o_gene][i_gene][set_name] = corr
                        merged_corr_dict[o_gene][i_gene][other_name] = \
                            other_corr
                        if 'z_score_mean' == settings['filter_type']:
                            merged_corr_dict[o_gene][i_gene][
                                'combined_z_score'] = \
                                0.5 * _z_sc(num=corr,
                                            mu=sigma_dict['mean'],
                                            sigma=sigma_dict['sigma']) + \
                                0.5 * _z_sc(num=other_corr,
                                            mu=other_sigma_dict['mean'],
                                            sigma=other_sigma_dict['sigma'])
                        assert merged_corr_dict[o_gene][i_gene].get(
                            set_name, None) is not None
                        assert merged_corr_dict[o_gene][i_gene].get(
                            other_name, None) is not None
                        assert merged_corr_dict[o_gene].get(o_gene, None)\
                            is None
                        assert merged_corr_dict[i_gene].get(i_gene, None)\
                            is None
                        npairs += 1

                    # Did not pass filter
                    else:
                        continue

                # Entry already exists
                else:
                    continue

    # todo: create recursive data set merger for 3 or more data sets
    # Case 3 (not yet implemented)
    # else:
    #     return merge_correlation_dicts_recursive(
    #         correlation_dicts_list.append(('master', merged_corr_dict)))

    merged_corr_dict.default_factory = None
    return merged_corr_dict, npairs


def merge_correlation_dicts_recursive(correlation_dicts_list):
    """This should be the recursive version of correlation_dicts_list(). Call
    this function from correlation_dicts_list()

    :param correlation_dicts_list:
    :return:
    """
    pass


def get_combined_correlations(dict_of_data_sets, filter_settings,
                              output_settings):
    """Return a combined dict of correlations given multiple gene data sets

    The input data [needs to be a collection of the gene expression data set.
    Need to ask user to list all sources]

    The input data set dict has the following format:

        dict_of_data_sets[gene_set_name] = dataset_dict

        dataset_dict = {data: str (depmap filepath),
                        corr: str (depmap corr file),
                        outbasename: str (base name for all output files),
                        ll: float (lower limit for correlation),
                        ul: float (upper limit for correlation),
                        max_pairs: int (max number of sampled pairs from corr)
                        mean: float (mean of correlation distr),
                        sigma: float (st-dev of correlation distr),
                        filter_margin: float (st-dev diff for filtering distr),
                        dump_unique_pairs: Bool (Output unique corr pairs),
                        }

    The filter settings should contain the following:

        filter_settings = {strict: Bool - If True, both A and B both have to
                            be in `gene_set_filter`.
                           gene_set_filter: list[genes to filter on]
                           cell_line_filter: list - cell line names in DepMap
                            ID format
                           margin: float - diff in terms of standard
                            deviations between correlations
                           filter_type: str - Type of filtering (Default: None)
                          }

    The output settings should contain the following:

        output_settings = {dump_unique_pairs: Bool - If True, dump a list of
                                                     all the unique pairs of
                                                     genes that have been
                                                     looped over.
                           outbasename: str - The output base name to be used
                                              in any file dump.
                           }

    The returned master correlation dict has the following format:

        d[gene][gene] = {gene_set1: correlation,
                         gene_set2: correlation,
                         ...}

    dict_of_data_sets: dict
        Dictionary containing the filepaths and settings for the data set
    filter_settings: dict
        Dictionary with filter settings

    Returns
    -------
    master_corr_dict: defaultdict(dict(...))
        A nested dict containing a lookup of the filtered set of
        gene-gene-correlations
    gene_set_intersection: set()
        The set of HGNC gene names in the master correlation lookup
    stats_dict: dict(dict)
    """
    name_dict_stats_list = []
    gene_set_intersection = set()
    stats_dict = dict()

    outbasename = output_settings['outbasename']

    for gene_set_name, dataset_dict in dict_of_data_sets.items():
        logger.info('-' * 37)
        logger.info(' > > > Processing set "%s" < < < ' % gene_set_name)
        logger.info('-' * 37)
        logger.info('Loading gene data...')
        gene_data = pd.read_csv(dataset_dict['data'], index_col=0, header=0)
        rows, cols = gene_data.shape
        if rows > cols:
            logger.info('Transposing data...')
            gene_data = gene_data.T

        # If filtering on cell lines, check if cell line IDs need to be
        # translated to DepMap ID (happens for RNAi)
        if filter_settings.get('cell_line_filter') and not \
                re.match('ACH-[0-9][0-9][0-9][0-9][0-9][0-9]',
                         gene_data.index.values[0]):
            assert filter_settings['cell_line_translation_dict'] is not None
            logger.info('Translating cell line names to DepMap ID')
            gene_data.rename(
                filter_settings['cell_line_translation_dict']['CCLE_Name'],
                inplace=True
            )

        if filter_settings.get('cell_line_filter'):
            logger.info('Filtering to provided cell lines')
            gene_data = gene_data[gene_data.index.isin(
                filter_settings['cell_line_filter'])]
            assert len(gene_data) > 0
            logger.info('Calculating Pearson correlation matrix from '
                            'cell line filtered data')
            full_corr_matrix = gene_data.corr()
            full_corr_matrix.to_hdf(outbasename +
                '_%s_cell_line_filtered_correlations.h5'
                % gene_set_name, 'correlations')
        else:
            if dataset_dict.get('corr'):
                logger.info('Reading pre-calculated correlation file.')
                full_corr_matrix = pd.read_hdf(
                    dataset_dict['corr'], 'correlations'
                )
            else:
                logger.info('No correlation file provided calculating '
                                'new Pearson correlation matrix...'
                            )
                full_corr_matrix = gene_data.corr()
                full_corr_matrix.to_hdf(outbasename +
                    '_%s_all_correlations.h5' % gene_set_name, 'correlations')

        logger.info('Removing self correlations for set %s' % gene_set_name)
        full_corr_matrix = full_corr_matrix[full_corr_matrix != 1.0]

        if full_corr_matrix.notna().sum().sum() == 0:
            logger.warning('Correlation matrix is empty')
            sys.exit('Script aborted due to empty correlation matrix')

        if dataset_dict.get('sigma'):
            logger.info('Using provided sigma of %f for set %s' %
                        (dataset_dict['sigma'], gene_set_name))
            sigma_dict = {'mean': dataset_dict['mean'],
                          'sigma': dataset_dict['sigma']}
        else:
            logger.info('Calculating mean and standard deviation for set %s'
                            ' from %s' % (gene_set_name, dataset_dict['data']))
            mu, si = get_stats(corr_matrix_to_generator(full_corr_matrix))
            logger.info('Set %s mean: %f, st dev: %f' %
                        (gene_set_name, mu, si))
            sigma_dict = {'mean': mu, 'sigma': si}

        # Get corr matrix and the accompanied set of genes
        filtered_corr_matrix, set_hgnc_syms, set_hgnc_ids,\
            sym2id_dict, id2sym_dict = get_correlations(
                depmap_data=gene_data,
                filter_gene_set=filter_settings['gene_set_filter'],
                pd_corr_matrix=full_corr_matrix,
                strict=filter_settings['strict'],
                dump_unique_pairs=output_settings['dump_unique_pairs'],
                outbasename=outbasename,
                sigma_dict=sigma_dict,
                lower_limit=dataset_dict['ll'],
                upper_limit=dataset_dict['ul']
            )
        logger.info('Created tuple generator with %i unique genes from '
                        'set "%s"' % (len(set_hgnc_syms), gene_set_name))
        if filtered_corr_matrix.notna().sum().sum() == 0:
            logger.warning('Correlation matrix is empty')
            sys.exit('Script aborted due to empty correlation matrix')

        logger.info('Dumping json HGNC symbol/id dictionaries...')
        io.dump_it_to_json(outbasename + '_%s_sym2id_dict.json' % gene_set_name,
                        sym2id_dict)
        io.dump_it_to_json(outbasename + '_%s_id2sym_dict.json' % gene_set_name,
                        id2sym_dict)

        # Generate correlation dict
        corr_dict = get_gene_gene_corr_dict(
            tuple_generator=corr_matrix_to_generator(
                corrrelation_df_matrix=filtered_corr_matrix,
                max_pairs=dataset_dict['max_pairs']
            )
        )

        # Append correlation dict and stats to list
        stats_dict[gene_set_name] = sigma_dict
        name_dict_stats_list.append((gene_set_name, corr_dict, sigma_dict))
        if len(gene_set_intersection) == 0:
            gene_set_intersection = set_hgnc_syms
        else:
            gene_set_intersection.intersection_update(set_hgnc_syms)

    if len(name_dict_stats_list) > 1:
        logger.info('---------------------')
        logger.info('Merging the data sets')
        logger.info('---------------------')

    # Merge the dictionaries
    master_corr_dict, npairs = merge_correlation_data(
        correlation_dicts_list=name_dict_stats_list,
        settings=filter_settings
    )
    logger.info('Created gene correlation master dictionary of length %i' %
                npairs)

    return master_corr_dict, gene_set_intersection, stats_dict


def get_correlations(depmap_data, filter_gene_set, pd_corr_matrix,
                     strict, dump_unique_pairs, outbasename, sigma_dict,
                     lower_limit=1.0, upper_limit=None):
    # todo make function take data dict as input or use args* + kwargs**
    """Return correlation data, filtered gene sets and gene id translation
    dictionaries given a depmap gene data file

    depmap_data: str
        Filepath to depmap data file to process
    filter_gene_set: str
        Filepath to a geneset to filter data to.
    pd_corr_matrix: str
        Filepath to pre-calculated correlations of depmap_data.
    strict: Bool
        If True, both genes in the gene pair have to exist in filter_gene_set
    outbasename: str
        Basename to use for output files
    unique_pair_corr_file: str
        Filepath to csvfile with unique tuples of gene,gene,corr.
    recalc: Bool
        If True, recalculate correlations (has to be True if pd_corr_matrix
        is None).
    lower_limit: float
        Smallest distance from correlation mean to consider
    upper_limit: float or None
        Largest distance from correlation mean to consider (Good for picking a
        sample in the middle of the correlation distribution). If None,
        there is no upper bound, i.e. all pairs with correlation == 1.0
        are considered.

    Returns
    -------
    filtered_correlation_matrix: pandas.DataFrame
        Correlation matrix as pd.DataFrame
    all_hgnc_ids: set()
        The set of all HGNC IDs in the correlation matrix
    all_hgnc_symb: set()
        The set of all HGNC symbols in the correlation matrix
    sym2id_dict: dict
        Dictionary translating HGNC symbols to HGNC IDs in the data set
    id2sym_dict: dict
        Dictionary translating HGNC IDs to HGNC symbols in the data set
    """

    filtered_correlation_matrix, sym2id_dict, id2sym_dict = _get_corr_df(
        depmap_data=depmap_data, corr_matrix=pd_corr_matrix,
        filter_gene_set=filter_gene_set, strict=strict,
        lower_limit=lower_limit, upper_limit=upper_limit,
        sigma_dict=sigma_dict
    )

    if filtered_correlation_matrix.notna().sum().sum() == 0:
        logger.warning('Correlation matrix is empty')
        raise ValueError('Script aborted due to empty correlation matrix')

    all_hgnc_symb = set(t[0] for t in filtered_correlation_matrix.index.values)
    all_hgnc_ids = set(t[1] for t in filtered_correlation_matrix.index.values)

    if dump_unique_pairs:
        if lower_limit == 0.0 and (upper_limit is None or upper_limit >= (1.0 -
        sigma_dict['mean']) / sigma_dict['sigma']):
            fname = outbasename + '_all_unique_correlation_pairs.csv'
        elif lower_limit > 0.0 and (upper_limit is None or upper_limit >= (
                1.0 - sigma_dict['mean']) / sigma_dict['sigma']):
            fname = outbasename + '_unique_correlation_pairs_ll%s.csv' % \
                    (str(lower_limit).replace('.', ''))

        else:
            fname = outbasename + '_unique_correlation_pairs_ll%s_ul%s.csv' % \
                    (str(lower_limit).replace('.', ''),
                     str(upper_limit).replace('.', ''))
        logger.info('Saving unique correlation pairs to %s. '
                        '(May take a while)' % fname)
        io.dump_it_to_csv(fname, corr_matrix_to_generator(
            filtered_correlation_matrix))

    return filtered_correlation_matrix,\
        all_hgnc_symb, all_hgnc_ids, \
        sym2id_dict, id2sym_dict


def _get_corr_df(depmap_data, corr_matrix, filter_gene_set,
                 strict, lower_limit, upper_limit, sigma_dict):
    # todo make function take data dict as input or use args* + kwargs**
    multi_index_data = pd.MultiIndex.from_tuples(
        tuples=[
            (t[0], t[1].strip('(').strip(')')) for t in [
                s.split() for s in depmap_data.columns.values
            ]
        ],
        names=['HGNCsymbol', 'HGNCid'])
    depmap_data.set_axis(axis=1, labels=multi_index_data, inplace=True)

    if len(corr_matrix.index.values[0].split()) == 2:
        # split 'HGNCsymb (HGNCid)' to 'HGNCsymb' '(HGNCid)' multiindexing:
        # pandas.pydata.org/pandas-docs/stable/advanced.html
        logger.info('Performing multi indexing of correlation matrix')

        # Get new indices
        hgnc_sym2id, hgnc_id2sym = {}, {}
        tuple_list = []
        for mystr in corr_matrix.index.values:
            hgnc_symb, hgnc_id = mystr.split()
            hgnc_id_ = hgnc_id.strip('(').strip(')')
            hgnc_sym2id[hgnc_symb] = hgnc_id_
            hgnc_id2sym[hgnc_id_] = hgnc_symb
            tuple_list.append((hgnc_symb, hgnc_id_))

        multi_index_corr = pd.MultiIndex.from_tuples(
            tuples=tuple_list,
            names=['HGNCsymbol', 'HGNCid']
        )

        # Add multi-index inplace as index and column name
        corr_matrix.set_axis(axis=0, labels=multi_index_corr, inplace=True)
        corr_matrix.set_axis(axis=1, labels=multi_index_corr, inplace=True)

    elif len(corr_matrix.index.values[0].split()) == 1:
        # leave intact? Check if there are IDs? Warning that you don't
        # have IDs/symbols but proceed?
        logger.warning('Only one identifier found in index column. '
                           'Assuming it is HGNC symbol.')
    else:
        logger.warning('Uknown index column. Output dictionaries will '
                           'likely be affected.')

    if filter_gene_set:
        # Read gene set to look at
        gene_filter_list = io.read_gene_set_file(
            gf=filter_gene_set, data=depmap_data
        )
        if len(gene_filter_list) == 0:
            logger.warning('Gene filter empty, continuing without filter')
            gene_filter_list = []
    else:
        gene_filter_list = []  # Evaluates to False

    row, col = corr_matrix.shape
    assert row > 0

    corr_matrix_df = None
    # 1. No gene set file, leave 'corr_matrix' intact
    if not filter_gene_set:
        logger.info('No gene filtering')
        corr_matrix_df = corr_matrix

    # 2. loaded gene list but not strict: filter correlation matrix to one of
    # the two genes in the pair being the
    elif filter_gene_set and not strict:
        try:
            # Try to split first item: raises AttributeError if tuple
            corr_matrix.index.values[0].split()
            corr_matrix_df = corr_matrix[gene_filter_list]
            logger.info('Non-strict gene filtering')
        except AttributeError:
            logger.info('Non-strict multi index gene filtering')
            corr_matrix_df = corr_matrix[np.in1d(
                corr_matrix.index.get_level_values(0),
                gene_filter_list)]

    # 3. Strict: both genes in interaction must be from loaded set;
    #    Filter data, then calculate correlations and then unstack
    elif filter_gene_set and strict:
        logger.info('Strict gene filtering')
        corr_matrix_df = corr_matrix_df[np.in1d(
            corr_matrix_df.index.get_level_values(0),
            gene_filter_list)]

    if corr_matrix_df is None or corr_matrix_df.notna().sum().sum() == 0:
        logger.warning('Correlation matrix is empty')
        sys.exit('Script aborted due to empty correlation matrix')

    # No filtering
    if lower_limit == 0.0 and (upper_limit is None or upper_limit >= (1.0 -
            sigma_dict['mean']) / sigma_dict['sigma']):
        logger.warning('No correlation filtering is performed. Be aware '
                           'of large RAM '
                          'usage.')
        return corr_matrix_df, hgnc_sym2id, hgnc_id2sym
    # Filter correlations
    else:
        return corr_limit_filtering(
            corr_matrix_df=corr_matrix_df,
            lower_limit=lower_limit,
            upper_limit=upper_limit,
            mu=sigma_dict['mean'],
            sigma=sigma_dict['sigma']
        ), hgnc_sym2id, hgnc_id2sym


def corr_limit_filtering(corr_matrix_df, lower_limit, upper_limit, mu, sigma):
    """Filters a correlation matrix to values in (lower_limit, upper_limit)

    corr_matrix_df: pandas.DataFrame
        A pandas correlation matrix as a pandas data frame
    lower_limit: float
        Smallest distance (measured in SD) from mean correlation to consider
    upper_limit: float
        Largest distance (measured in SD) from mean correlation to consider
        (good for picking a sample in the middle of the correlation
         distribution)

    Returns
    -------
    corr_matrix_df: pandas.DataFrame
        A filtered correlation dataframe matrix
    """
    # Filter by number of SD from mean
    if lower_limit > 0.0 and upper_limit and upper_limit < (1.0 - mu) / sigma:
        logger.info('Filtering correlations to range %.2f < abs(C-mu)/SD < '
                        '%.2f' % (lower_limit, upper_limit))
        corr_matrix_df = corr_matrix_df[
            abs(corr_matrix_df - mu) / sigma > lower_limit
        ]
        corr_matrix_df = corr_matrix_df[
            abs(corr_matrix_df - mu) / sigma < upper_limit
        ]
    elif lower_limit == 0.0 and upper_limit and upper_limit < (1.0 - mu) / \
            sigma:
        logger.info('Filtering correlations to range 0.0 <= abs(C-mu)/SD < '
                        '%.2f' % upper_limit)
        corr_matrix_df = corr_matrix_df[
            abs(corr_matrix_df - mu) / sigma < upper_limit
        ]
    elif lower_limit > 0.0 and not upper_limit:
        logger.info('Filtering correlations to range %.2f < abs(C-mu)/SD'
                    % lower_limit)
        corr_matrix_df = corr_matrix_df[
            abs(corr_matrix_df - mu) / sigma > lower_limit
        ]
    elif lower_limit == 0.0 and not upper_limit:
        logger.info('Not filtering correlations.')

    return corr_matrix_df


def pass_filter(corr1, mu1, sigma1, corr2, mu2, sigma2, margin=None,
                filter_type='z_scor_diff'):
    """Filter for passing correlation scores based on their difference in
    standard deviation

    corr1: float
        Correlation of first dataset
    mu1: float
        Mean of the correlations of the first dataset
    sigma1: float
        Standard deviation of the correlations of the first dataset
    n2: float
        Correlation from second data set
    mu2: float
        Mean of the correlations of the second dataset
    sigma2: float
        Standard deviation of the correlations of the second dataset
    margin: float
        Number to set as pass cutoff for some of the filters:
        - z-score-mean: mean of the z-scores is smaller than the margin
        - z-score-diff: difference of the z-scores is smaller than the margin
        - z-score-product: product of the z-scores is greater than the margin
    filter_type:
        The filter type to use

    Returns
    -------
    bool
        If True, the correlations are similar enough as measured by their
        difference in their distance from the mean standard deviation.
    """
    if filter_type in {'z_score_mean', 'z_score_diff', 'z_score_product'} \
            and margin is None:
        raise ValueError('Margin needs to be set for filter type %' %
                         filter_type)
    if filter_type == 'z_score_mean':
        return _z_score_mean_co(corr1, mu1, sigma1, corr2, mu2, sigma2, margin)
    elif filter_type == 'z_score_diff':
        return _z_score_diff(corr1, mu1, sigma1, corr2, mu2, sigma2, margin)
    elif filter_type == 'z_score_product':
        return _z_score_product(corr1, mu1, sigma1, corr2, mu2, sigma2, margin)
    elif filter_type == 'sign':
        return same_sign(corr1, corr2)
    # No filter/filter not recognized:
    else:
        return True


def _z_sc(num, mu, sigma):
    """Return z-score of given num drawn from given distribution"""
    return (num - mu) / sigma


def _z_score_mean_co(corr1, mu1, sigma1, corr2, mu2, sigma2, margin):
    """Pass if the mean of the z-scores is greater than margin"""
    return 0.5 * _z_sc(corr1, mu1, sigma1) + \
        0.5 * _z_sc(corr2, mu2, sigma2) > margin


def _z_score_diff(corr1, mu1, sigma1, corr2, mu2, sigma2, margin):
    """Pass if the difference in the z-score is smaller than margin
    """
    return abs(_z_sc(corr1, mu1, sigma1) - _z_sc(corr2, mu2, sigma2)) < margin


def _z_score_product(corr1, mu1, sigma1, corr2, mu2, sigma2, margin):
    """Pass if the product of the z-scores is greater than margin"""
    return _z_sc(corr1, mu1, sigma1) * _z_sc(corr2, mu2, sigma2) > margin


def get_sign(num):
    """Get the sign of num

    Parameters
    ----------
    num : number
        Any valid number as int, float or str

    Returns
    -------
    int
        -1 or +1
    """
    if isinstance(num, str):
        try:
            num = float(num)
        except ValueError as err:
            logger.warning('object is not recognized as a number')
            raise err
    return math.copysign(1, num)


def same_sign(n1, n2):
    """Return True if n1 and n2 have the same sign

    Special cases:
        zeros:
            1. if n1 == 0 AND n2 == 0, return True
            2. if n1 == 0 XOR n2 == 0, return False

        inf and NaN:
            if not math.isfinite(n1) OR not math.isfinite(n2), return False

        Can't be interpreted as number, return False
    """
    # Catch non-numeric correlations
    try:
        if isinstance(n1, str):
            n1 = float(n1)
        if isinstance(n2, str):
            n2 = float(n2)
    except ValueError:
        logger.warning('Correlation could not be interpreted as numeric. '
                           'Skipping...')
        return False

    # Catch nan and inf
    if not math.isfinite(n1) or not math.isfinite(n2):
        logger.warning('Correlation is undefined. Skipping...')
        return False

    # Both zero
    if n1 == 0 and n2 == 0:
        return True

    # XOR: if (n1==0 or n2==0) and not (n1==0 and n2==0)
    elif (n1 == 0) ^ (n2 == 0):
        return False

    return math.copysign(1, n1) == math.copysign(1, n2)


def get_directed(stmts, undirected_types=None):
    """Given a statement list, sort statements based on directionality.

    The statements can be either of regular INDRA statements or statement
    type-statement hash pairs or statement JSONs.

    stmts : list[stmts]
        A list of INDRA statements.
    undirected : [statement types]
        A list of name strings considered to be undirected.
        Default: ['Complex', 'SelfModification', 'parent']

    Returns
    -------
    dir_stmts, undir_stmts : ([stmts], [stmts])
        Two lists of statements, one containiing all undirected statements
        and one contaning all directed statements.
    """

    dir_stmts, undir_stmts = [], []
    if not undirected_types:
        undirected_types = ['Complex', 'SelfModification', 'parent']

    # check type of stmt list
    if stmts:
        # if (type, hash)
        if type(stmts[0]) == tuple:
            dir_stmts, undir_stmts = get_directed_type_hash(stmts,
                                                            undirected_types)
        # if normal statements
        elif isinstance(stmts[0], Statement):
            dir_stmts, undir_stmts = \
                get_directed_actual_statements(stmts, undirected_types)
        # if json statements
        elif isinstance(stmts[0], OrderedDict):
            dir_stmts, undir_stmts = get_directed_json(stmts, undirected_types)

    return dir_stmts, undir_stmts


def get_directed_type_hash(stmts, undirected_types):
    """Given a list of type-hash tuples, sort statements on directionality.

    stmts : [(type, hash)] or [(type, hash, belief)]
         A list of type-string, hash tuples.
    undirected : [statement types]
        A list of name strings considered to be undirected.

    Returns
    -------
    dir_stmts, undir_stmts : ([stmts], [stmts])
        Two lists of statements, one containiing all undirected statements
        and one contaning all directed statements.
    """

    dir_stmts, undir_stmts = [], []

    for stmt in stmts:
        if len(stmt) == 2:
            ctype, hsh = stmt
            hash_string = str(hsh)
            tup = (ctype, hash_string)
        elif len(stmt) == 3:
            ctype, hsh, bs = stmt
            hash_string = str(hsh)
            tup = (ctype, hash_string, bs)

        if ctype in undirected_types:
            undir_stmts.append(tup)
        else:
            dir_stmts.append(tup)

    return dir_stmts, undir_stmts


def get_directed_actual_statements(stmts, undirected_types):
    """Given a list of INDRA statements, sort statements on directionality.

    stmts: list[:py:class:`indra.statements.Statement`]
        List of INDRA statements
    undirected: list[statement types]
        A list of name strings considered to be undirected.

    Returns
    -------
    dir_stmts, undir_stmts : [stmts], [stmts]
        Two lists of statements, one containiing all undirected statements
        and one contaning all directed statements.
    """

    dir_stmts, undir_stmts = [], []

    for stmt in stmts:
        if stmt.to_json()['type'] in undirected_types:
            undir_stmts.append(stmt)
        else:
            dir_stmts.append(stmt)

    return dir_stmts, undir_stmts


def get_directed_json(stmts, undirected_types):
    """Given a list of json statements, sort statements on directionality.

    stmts: list[json statements]
        A list of INDRA statements in json format
    undirected: [statement types]
        A list of name strings considered to be undirected.

    Returns
    -------
    dir_stmts, undir_stmts : ([stmts], [stmts])
        Two lists of statements, one containiing all undirected statements
        and one contaning all directed statements.
    """

    dir_stmts, undir_stmts = [], []

    for stmt in stmts:
        if stmt['type'] in undirected_types:
            undir_stmts.append(stmt)
        else:
            dir_stmts.append(stmt)

    return dir_stmts, undir_stmts


def agent_name_set(stmt):
    """Return the list of agent names in a statement.

    stmt: :py:class:`indra.statements.Statement`
        INDRA statement

    Returns
    -------
    ags: list[agent names]
        List of agent names in statement
    """
    ags = []
    try:
        ags.update(list(map(lambda ag: ag.name, stmt.agent_list())))
    except AttributeError:
        for ag in stmt.agent_list():
            if ag is None:
                pass
            else:
                ags.append(ag.name)
    return ags


def nested_hash_dict_from_pd_dataframe(hash_pair_dataframe):
    """Returns a nested dict of

    hash_pair_dataframe: pandas.DataFrame
        A padnas dataframe containing indra statements/statement hashes

    Returns
    -------
    nest_hash_dict: defaultdict(dict)
        A nested dict of d[subj][obj] = [statement hashes]
    """
    nest_hash_dict = defaultdict(dict)

    # Row should be a mini dataframe with keys:
    # agent_1=subj, agent_2=obj, type, hash
    for index, row in hash_pair_dataframe.iterrows():
        if len(row) == 4:
            subj, obj, stmt_type, stmt_hash = row
        else:
            subj, obj, stmt_type, stmt_hash = row.agA_name, row.agB_name, \
                                              row.stmt_type, row.stmt_hash
        if entry_exist_dict(nest_hash_dict, subj, obj):
            # Entry subj-obj already exists and should be a list
            nest_hash_dict[subj][obj].append((stmt_type, stmt_hash))
        else:
            nest_hash_dict[subj][obj] = [(stmt_type, stmt_hash)]
        # # Add reverse direction if type is complex or selfmodification
        # if type.lower() in ['complex', 'selfmodification']:
        #     nest_hash_dict[obj][subj][stmt_type] = stmt_hash
    return nest_hash_dict


def nested_dict_of_stmts(stmts, belief_dict=None):
    """Generates a nested dict of the form dict[key1][key2] = [statement list]
    from INDRA statements.

    stmts :  list[:py:class:`indra.statements.Statement`]
        List or set of INDRA statements to find connections in

    Returns
    -------
    stmts_dict : collections.defaultdict
         dict of the form dict[subj][obj] = list[stmts]
    """

    if belief_dict is None:
        logger.warning('No belief score dict is provided! Please provide a '
                           'belief score dict through the `-b ('
                           '--belief-score-dict)` argument.')
        raise ValueError

    nested_stmt_dicts = defaultdict(dict)

    for st in stmts:
        # NOTE: If statement is complex, it migth have more than two agents
        # and the agents won't be distinguishable as subject,object
        try:
            bs = belief_dict[str(st)]
        except KeyError:
            bs = 1
        agent_list = agent_name_set(stmt=st)

        # It takes two agents to tango
        if len(agent_list) > 1:

            # Is complex or selfmodification
            if st.to_json()['type'].lower() in ['complex', 'selfmodification']:
                for agent, other_agent in itt.permutations(agent_list, r=2):
                    try:
                        nested_stmt_dicts[agent][other_agent].append((st, bs))
                    except KeyError:  # If pair does not exist yet
                        nested_stmt_dicts[agent][other_agent] = [(st, bs)]

            # Non-complex interaction
            else:
                subj = agent_list[0]
                obj = agent_list[1]

                try:
                    nested_stmt_dicts[subj][obj].append((st, bs))
                except KeyError:
                    nested_stmt_dicts[subj][obj] = [(st, bs)]

            # Check common parent (same family or complex)
            for agent, other_agent in itt.permutations(agent_list, r=2):
                if ff.has_common_parent(id1=agent, id2=other_agent):
                    bs = None
                    try:
                        if 'parent' not in \
                                nested_stmt_dicts[agent][other_agent]:
                            nested_stmt_dicts[agent][other_agent].append(
                                ('parent', bs))
                    except KeyError:
                        nested_stmt_dicts[agent][other_agent] = [('parent', bs)]

        # Ignore when we only have one agent
        else:
            continue

    logger.info('Created nested dict of length %i from %i statements.' %
                (len(nested_stmt_dicts), len(stmts)))
    return nested_stmt_dicts


def dedupl_nested_dict_gen(stmts, belief_dict):
    """De-duplicate a list of statments

    stmts: list[:py:class:`indra.statements.Statement`]
        List of INDRA statements
    belief_dict : dict()
        dict with {hash: belief score} as key: value pairs

    Returns
    -------
    nested_dict_stmts: defaultdict(dict)
        A nested dict of deduplicated statement
    """
    ddstmt_list = deduplicate_stmt_list(stmts=stmts, ignore_str='parent')
    nested_dict_stmts = nested_dict_of_stmts(ddstmt_list, belief_dict)

    return nested_dict_stmts


def deduplicate_stmt_list(stmts, ignore_str):
    """Takes a list of statements list[stmts] and runs
    indra.preassembler.Preassembler.combine_duplicate_stmts() while also
    taking care of non-statements

    stmts : list[:py:class:`indra.statements.Statement`]
        List of INDRA statements

    Returns
    -------
    stmts : list[:py:class:`indra.statements.Statement`]
         List of preassembled statments possibly including a non-statements
    """
    # subjects should be the outer keys and objects should be the inner

    if ignore_str in stmts:
        logger.info('Deduplicating statements and accounting for custom '
                        'string %s' % ignore_str)
        only_stmt_list = [s for s in stmts if type(s) is not str]
        stmts = pa_filter_unique_evidence(only_stmt_list)
        stmts += [ignore_str]
    else:
        stmts = pa_filter_unique_evidence(stmts)
    return stmts


def pa_filter_unique_evidence(stmts):
    """Wrapper function for chaining preassembly statements meant to reduce
    the number of statements.

    stmts : list[:py:class:`indra.statements.Statement`]
        List of INDRA statements

    Returns
    -------
    stmts : list[:py:class:`indra.statements.Statement`]
        List of preassembled indra statements
    """

    # Ground statemtens:
    grounded_stmts = ac.map_grounding(stmts)

    # Use curated site information to standardize modification sites in stmts
    ms_stmts = ac.map_sequence(grounded_stmts)

    # Compiles together raw statements to one statement per type
    opa_stmts = ac.run_preassembly(ms_stmts, return_toplevel=False)
    return opa_stmts


def latex_output(subj, obj, corr, stmts, ev_len_fltr, ignore_str='parent'):
    """Compiles information about statements and returns a LaTeX
    formatted string that can be written to a file.

    subj : str
        indra statement subject
    obj : str
        indra statment object
    corr : float
        Correlation between subject and object
    stmts : list[:py:class:`indra.statements.Statement`]
        List of indra statements
    ignore_str : str
        String to ignore if it appears in the list of indra statements

    Returns
    -------
    output : str
        string with information about the statements that connect subject and
        object formatted for printing or for writing to a text file.
    """

    output = ''

    # Build up a string that shows explanations for each connection
    # This string is put in the script instead
    # output = r'\section{{{}, {}: {}}}'.format(subj, obj, corr) + '\n' + \
    #          r'See correlation plot \href{{' \
    #          r'https://depmap.org/portal/interactive/?xDataset=Avana' \
    #          '&xFeature' \
    #          '={}&yDataset=Avana&yFeature={}&colorDataset=lineage' \
    #          '&colorFeature=all&filterDataset=context&filterFeature=' \
    #          '&regressionLine=false&statisticsTable=false&associationTable=' \
    #          'true&plotOnly=false}}{{here}}'.format(subj, obj) + '\n\n'

    pa_stmts = stmts.copy()

    # HERE: insert subsection A->B
    output += r'\subsection{{{A} $\rightarrow$ {B}}}'.format(A=subj, B=obj)+'\n'

    # Sort stmts by evidence length
    # stmts_dict = dict()
    # ev_lens = []
    # for stmt in pa_stmts:
    #     stmts_dict[str(stmt)] = stmt
    #     ev_text_list = list(set(['N/A' if ev.text is None else ev.text for ev in
    #                        stmt.evidence]))
    #     if 'N/A' in ev_text_list:
    #         ev_text_list.remove('N/A')
    #     # Save tuple (len, str(stmt))
    #     ev_lens.append((len(ev_text_list), str(stmt)))
    # ev_lens.sort(key=lambda tup: tup[0], reverse=True)

    # HERE: itemize per statement type in stmt for loop
    output += r'\begin{itemize}'+'\n'
    for stmt in pa_stmts:
    # for lene, st_key in ev_lens:
    #     stmt = stmts_dict[st_key]
        if type(stmt) is str and str(stmt) == ignore_str:
                output += r'\item {s} and {o} are in the same complex ' \
                          r'or family'.format(s=subj, o=obj) + '\n'
        else:
            # Remove duplicate evidence text
            ev_text_set = set(['N/A' if ev.text is None else ev.text for ev in
                               stmt.evidence])
            ev_text_list = list(ev_text_set)
            if 'N/A' in ev_text_list:
                ev_text_list.remove('N/A')

            # assert lene == len(ev_text_list)

            if len(ev_text_list) >= ev_len_fltr:

                output += r'\item {nstmts} instances found of statement '\
                          r'{stmt}; supports count: {supc}; Supported by ' \
                          r'count: {supbc}'.format(stmt=str(stmt),
                                                   nstmts=len(ev_text_list),
                                                   supc=len(stmt.supports),
                                                   supbc=len(
                                                       stmt.supported_by))+'\n'

                # There are statements with zero length evidence lists
                if len(ev_text_list) > 0:
                    # HERE enumerate evidence text
                    output += r'\begin{enumerate}'+'\n'

                    max_ev = 25  # Dont ouput more than 25 evidences
                    for count, ev_text in enumerate(ev_text_list):
                        output += r'\item Evidence: \texttt{' + ev_text + \
                                  r'}' + '\n'
                        if count+1 == max_ev:
                            break

                    output += r'\end{enumerate}'+'\n'
            else:
                output += 'Evidence count below threshold of {}.\n'\
                    .format(ev_len_fltr)

    output += r'\end{itemize}' + '\n'

    # Don't forget to escape latex characters in 'ev_text'
    output = output.replace('_', '\_')\
        .replace('%', '\%')\
        .replace('&', '\&')\
        .replace('^', '\^')\
        .replace('~', '\~')

    return output.encode('ascii', 'ignore').decode('ascii')


def dbc_load_statements(hgnc_syms):
    """Get statements where hgnc symbol is subject/object from indra.db.client

    Parameters
    ----------
    hgnc_syms : Iterable
        An iterable containing HGNC symbols

    Returns
    -------
    stmts : set{:py:class:`indra.statements.Statement`}
        A set of all retrieved INDRA statemetents containing HGNC symbols
    """
    stmts = set()
    counter = 0
    n_hgnc_ids = len(hgnc_syms)
    try:
        for hgnc_id in hgnc_syms:
            stmts.update(dbc.get_statements_by_gene_role_type(agent_id=hgnc_id,
                                                              db=db_prim,
                                                              preassembled=
                                                              False,
                                                              fix_refs=False))
            counter += 1
            if counter % max([10, 10 ** ceil(log10(n_hgnc_ids)) // 100]) == 0:
                logger.info(' : : : Finished %i queries out of %i '
                                ': : :' % (counter, n_hgnc_ids))

    except KeyboardInterrupt as e:
        db_prim.session.rollback()
        raise e
    except StatementError as e:
        db_prim.session.rollback()
        raise e
    return stmts


def direct_relation(id1, id2, long_stmts=set()):
    """Returns a list of INDRA statements

    Parameters
    ----------
    id1/id2 : str
        Strings of the two ids to check a direct relation between.
    long_stmts : set[:py:class:`indra.statements.Statement`]
        (Optional) List or set of INDRA statements to find connections in

    Returns
    -------
    stmts : list[:py:class:`indra.statements.Statement`]
        List of INDRA statements that directly relate id1 and id2
    """
    if not long_stmts:
        stmts = direct_relation_from_api(id1=id1, id2=id2)
    else:
        stmts = direct_relation_from_stmts(id1=id1, id2=id2,
                                           stmts_in=long_stmts)
    return stmts


def direct_relation_from_api(id1, id2, on_limit='sample'):
    """Queries INDRA DB for Statements linking two genes and returns a list
    containing the matching statements.

    Parameters
    ----------
    id1/id2 : str
        Strings of the two ids to check a direct relation between.
    on_limit : str
        There are four options for handling the a query that is to large:
        `sample` - (default) take a sample of statements from the result,
        `truncate` - simply return the first 10,000 statements of the result,
        `error` - raise an error if the query is too large, or
        `persist` - perform as many queries as needed to get all the statements.
        Note that this last option generally takes much much longer to execute

    Returns
    -------
    stmts : list[:py:class:`indra.statements.Statement`]
        A list of INDRA Statement instances.
    """
    try:
        stmts = db_api.get_statements(subject=id1,
                                      object=id2,
                                      on_limit=on_limit)
        stmts + db_api.get_statements(subject=id2,
                                      object=id1,
                                      on_limit=on_limit)
    except IndraDBRestAPIError:
        stmts = db_api.get_statements(subject=id1 + '@TEXT',
                                      object=id2 + '@TEXT',
                                      on_limit=on_limit)
        stmts + db_api.get_statements(subject=id2 + '@TEXT',
                                      object=id1 + '@TEXT',
                                      on_limit=on_limit)
    return stmts


def direct_relation_from_stmts(id1, id2, stmts_in):
    """Returns a list of INDRA statements that connect id1 and id2 queried
    from a provided list of statements,

    Parameters
    ----------
    id1/id2 : str
        Strings of the two ids to check a direct relation between.
    stmts_in : set[:py:class:`indra.statements.Statement`]
        List of INDRA statements to find connections in.

    Returns
    -------
    stmts_out : list[:py:class:`indra.statements.Statement`]
        List of INDRA statements that directly relate id1 and id2
    """
    target_ag = {id1, id2}
    stmts_out = []
    for stms in stmts_in:
        s_agents = agent_name_set(stms)
        if target_ag.issubset(s_agents):
            stmts_out.append(stms)
    return stmts_out


def relation_type(indra_stmt):
    """Return the statement type in an INDRA statement as a string.

    Parameters
    ----------
    indra_stmt : :py:class:`indra.statements.Statement`
        INDRA statment

    Returns
    -------
    relation type : str
        A string containing an INDRA relation type
    """
    return indra_stmt.to_json()['type']


def relation_types(stmts):
    """Returns the corresponding list of INDRA Statement types associated
    with a list of Statements.

    Parameters
    ----------
    stmts : list[:py:class:`indra.statements.Statement`]
        A list of INDRA Statement instances

    Returns
    -------
    types : list[INDRA statement types]
        A list of strings containing the INDRA statement types
    """
    types = []
    for stmt in stmts:
        types.append(relation_type(stmt))
    return types


def has_direct_relation(id1, id2, long_stmts=set()):
    """Indicates whether two genes are linked by Statements in the INDRA data
    bases.

    Parameters
    ----------
    id1/id2 : str
        HGNC names for the two genes.

    Returns
    -------
    bool
        True if the HGNC symbols has a direct relation found in the
        indra.sources.indra_db_rest.client_api databases.
    """
    return bool(direct_relation(id1, id2, long_stmts=long_stmts))


def are_connected(id1, id2, long_stmts=set()):
    """Indicates whether two genes have a connection either through a direct
    relation or a through a common parent.

    Parameters
    ----------
    id1/i2 : str
        HGNC symbol

    Returns
    -------
    bool
        True if the two HGNC symbols either have a common parent or if they
        have a direct relation found in the
        indra.sources.indra_db_rest.client_api databases.
    """
    return ff.has_common_parent(ns1='HGNC', id1=id1,
                                ns2='HGNC', id2=id2) or \
           has_direct_relation(id1=id1, id2=id2, long_stmts=long_stmts)


def connection_types(id1, id2, long_stmts=set()):
    """Returns a list of the connection types linking two genes.

    Parameters
    ----------
    id1/i2 : str
        HGNC symbol

    Returns
    -------
    ctypes : list[type]
        Returns a list of connection types.
        `[]` - empty list if no connections.
        Type is any of:
        `INDRA statement` - Any INDRA statement type
        `parent` - id1 and id2 are connected through common parent(s)
    """

    ctypes = relation_types(direct_relation(id1=id1, id2=id2,
                                            long_stmts=long_stmts))
    if ff.has_common_parent(id1=id1, id2=id2):
        ctypes += ['parent']
    return ctypes


def down_sampl_size(available_pairs, size_of_matrix, wanted_pairs,
                    buffer_factor=2):
    """Return a sample size that would make a new square dataframe contain
    close to but above the number of wanted pairs

    Assuming that the number of extractable pairs is N = k*L^2 (k is some
    constant), we want to pick the fraction p of the matrix size L that
    provides the number of rows to sample such that the resulting number of
    extractable pairs from the matrix comes close to s, the number of
    desired samples.

    Parameters
    ----------
    available_pairs : int
        Number of pairs in correlation matrix as counted by
        corr_z.mask(
            np.triu(np.ones(corr_z.shape)).astype(bool)
        ).notna().sum().sum()
    size_of_matrix : int
        The number of rows/columns of the correlation matrix as counted by
        len(corr)
    wanted_pairs : int
        The target number of pairs to the get out from the matrix
    buffer_factor : float
        The buffer factor will be multiplied with the wanted_pairs to
        safeguard that the output sample size does not yield a matrix
        containing fewer than the wanted number of extractable pairs.

    Returns
    -------
    int
    """
    # Set vars
    L = size_of_matrix
    N = available_pairs
    s = wanted_pairs

    # Calculate fraction
    p = np.sqrt(buffer_factor*s/N)

    # Get number of rows to sample
    return int(np.ceil(p*L))


def get_pairs(corr_z: pd.DataFrame) -> int:
    """Count the number of extractable pairs from a pandas correlation matrix

    Count the number of pairs that can be looped over from the DataFrame
    correlation matrix from the upper triangle of the matrix (since the
    matrix is assumed to be symmetric) with NaN's and potential values on
    the diagonal ignored.

    Parameters
    ----------
    corr_z : pd.DataFrame
        A DataFrame with correlations obtained from pandas.DataFrame.corr().
        The DataFrame is assumed to be pre-filtered such that values
        filtered out are NaN's.

    Returns
    -------
    int
        The count of pairs that can be looped over
    """
    # Expect corr_z to be filtered to the values of interest and that the
    # values that are filtered out are NaN's

    # Map to boolean with NaN => False, else True
    bm: pd.DataFrame = (~corr_z.isna())

    # Mask lower triangle and diagonal with zeroes
    ma = bm.mask(np.tril(np.ones(bm.shape).astype(bool)), other=0)

    # Return sum over full matrix
    return int(ma.sum().sum())


def get_chunk_size(n_chunks: int, total_items: int) -> int:
    """Find n such that `(n-1)*n_chunks < total_items <= n*n_chunks`

    Parameters
    ----------
    n_chunks : int
        The number of chunks of iterables wanted
    total_items : int
        The total number of items in the original iterable

    Returns
    -------
    int
        Return n in `(n-1)*n_chunks < total_items <= n*n_chunks`
    """
    # How many pairs does a chunk need to contain to get chunks_wanted chunks?
    return max(int(np.ceil(total_items / n_chunks)), 1)
