import os
import csv
import json
import logging
import numpy as np
import pandas as pd
import networkx as nx
import itertools as itt
from math import ceil, log10
from collections import Mapping
from collections import defaultdict
from collections import OrderedDict
from sqlalchemy.exc import StatementError
from scipy import interpolate as interpol
from scipy.optimize import curve_fit as opt_curve_fit
from pandas.core.series import Series as pd_Series_class
from indra_db import util as dbu
from indra_db import client as dbc
from indra.statements import Statement
from indra.tools import assemble_corpus as ac
from indra.preassembler import hierarchy_manager as hm
from indra.sources.indra_db_rest import client_api as capi
from indra.sources.indra_db_rest.client_api import IndraDBRestError

import pdb  # ToDo remove import before final commit

db_prim = dbu.get_primary_db()
dnf_logger = logging.getLogger('DepMapFunctionsLogger')


def _entry_exist(nest_dict, outer_key, inner_key):
    if nest_dict.get(outer_key) and nest_dict.get(outer_key).get(inner_key):
        return True
    else:
        return False


def create_nested_dict():
    """Returns a nested dictionary of arbitrary depth

    Returns
    -------
    defaultdict(create_nested_dict)
    """
    return defaultdict(create_nested_dict)


def csv_file_to_generator(fname, column_list):
    """Return a tuple generator given a csv file and specified columns

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


def corr_matrix_to_generator(corrrelation_df_matrix):
    """Return a tuple generator given a correlation matrix

    corrrelation_df_matrix : pandas.DataFrame
        A pandas correlation matrix as a dataframe

    Returns
    -------
    tuple_generator : generator object
        A generator that returns a tuple of each row
    """
    corr_value_matrix = corrrelation_df_matrix.values
    gene_name_array = corrrelation_df_matrix.index.values
    tr_up_indices = np.triu_indices(n=len(corr_value_matrix), k=1)
    return ((gene_name_array[i], gene_name_array[j],
            str(corr_value_matrix[i, j]))
            for i, j in zip(*tr_up_indices)
            if not np.isnan(corr_value_matrix[i, j]))


def _dump_it_to_json(fname, pyobj):
    with open(fname, 'w') as json_out:
        json.dump(pyobj, json_out)


def _dump_it_to_csv(fname, iterable, separator=','):
    with open(fname, 'w', newline='') as csvf:
        wrtr = csv.writer(csvf, delimiter=separator)
        wrtr.writerows(iterable)


def nx_undir_to_neighbor_lookup_json(expl_undir_graph,
                                     path_prefix='neighbor_lookup/'):
    """Dumps one json dict per node in a undirected nx graph where entry is a
    json array containing all neighbors

    expl_undir_graph : nx.Graph
    path_prefix : directory to put jsons in. Must be specified or the default
    'neighbor_lookup/' is used.

    """
    if not os.path.isdir(path_prefix):
        dnf_logger.info('Could not find path "%s", creating new directory.' %
                        path_prefix)
        os.mkdir(path_prefix)

    dnf_logger.info('Dumping node neighbor dicts to "%s'
                    'neighbors_to_NODENAME.json"' % path_prefix)
    for node in expl_undir_graph.nodes:
        nnnl = []
        for other_node in expl_undir_graph[node]:
            inner_dict = expl_undir_graph[node][other_node]
            nnnl.append([other_node, inner_dict['attr_dict']['correlation']])
        _dump_it_to_json(fname=path_prefix+'neighbors_to_%s.json' % node,
                         pyobj=nnnl)
    dnf_logger.info('Finished dumping node neighbor dicts to %s' % path_prefix)


def _filter_corr_data(corr, clusters, cl_limit):
    # To get gene names where clustering coefficient is above limit.
    # Clusters has to come from an already produced graph.
    # The filtered correlations can then be used to produce a new graph.
    # The usage of this would be to get a nicer plot that conecntrates on the
    # clusters instead of plotting _everything_, making it hard to see the
    # forest for all the trees.

    filtered_genes = [k for k in clusters if clusters[k] > cl_limit]
    filtered_correlations = corr[filtered_genes].unstack()
    return filtered_correlations


def rank_nodes(node_list, nested_dict_stmts, gene_a, gene_b, x_type):
    """Returns a list of tuples of nodes and their rank score

    The provided node list should contain the set of nodes that connects subj
    and obj through an intermediate node found in nested_dict_stmts.

    nested_dict_stmts

        d[subj][obj] = [stmts/stmt hashes]

    node_list : list[nodes]
    nested_dict_stmts : defaultdict(dict)
        Nested dict of statements: nest_d[subj][obj]
    gene_a : str
        HGNC name of gene A in an A-X-B connection
    gene_b : str
        HGNC name of gene B in an A-X-B connection
    x_type : str
        One of 'x_is_intermediary', 'x_is_downstream' or 'x_is_upstream'

    -------
    Returns
    dir_path_nodes_wb : list[(node, rank)]
        A list of node, rank tuples.
    """

    def _calc_rank(nest_dict_stmts, subj_ax, obj_ax, subj_xb, obj_xb):
        ax_stmts = nest_dict_stmts[subj_ax][obj_ax]
        xb_stmts = nest_dict_stmts[subj_xb][obj_xb]
        ax_score_list = []
        xb_score_list = []

        # The statment with the highest belief score should
        # represent the edge (potentially multiple stmts per edge)
        for typ, hsh_a, bs in ax_stmts:
            ax_score_list.append(bs)
        for typ, hsh_b, bs in xb_stmts:
            xb_score_list.append(bs)

        # Rank by multiplying the best two belief scores for each edge
        rank = max(ax_score_list) * max(xb_score_list)

        # No belief score should be zero, thus rank should never be zero
        try:
            assert rank != 0
        except AssertionError:
            dnf_logger.warning('Combined rank == 0 for hashes %s and %s' %
                               (hsh_a, hsh_b))
            pdb.set_trace()
        return rank

    dir_path_nodes_wb = []
    if x_type is 'x_is_intermediary':  # A->X->B or A<-X<-B
        for gene_x in node_list:
            x_rank = _calc_rank(nest_dict_stmts=nested_dict_stmts,
                                subj_ax=gene_a, obj_ax=gene_x,
                                subj_xb=gene_x, obj_xb=gene_b)
            dir_path_nodes_wb.append((gene_x, x_rank))

    elif x_type is 'x_is_downstream':  # A->X<-B
        for gene_x in node_list:
            x_rank = _calc_rank(nest_dict_stmts=nested_dict_stmts,
                                subj_ax=gene_a, obj_ax=gene_x,
                                subj_xb=gene_b, obj_xb=gene_x)
            dir_path_nodes_wb.append((gene_x, x_rank))
    elif x_type is 'x_is_upstream':  # A<-X->B

        for gene_x in node_list:
            x_rank = _calc_rank(nest_dict_stmts=nested_dict_stmts,
                                subj_ax=gene_x, obj_ax=gene_a,
                                subj_xb=gene_x, obj_xb=gene_b)
            dir_path_nodes_wb.append((gene_x, x_rank))

    return dir_path_nodes_wb


def nx_directed_multigraph_from_nested_dict(nest_d):
    """Returns a directed multigraph where each edge links a statement with
    u=subj, v=obj, edge_key=stmt,

    nest_d : defaultdict(dict)
        Nested dict of statements: nest_d[subj][obj]

    Returns
    -------
    nx_muldigraph : nx.MultiDiGraph
        An nx directed multigraph linking agents with statements
    """

    dnf_logger.info('Building directed multigraph from nested dict of '
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


def nx_directed_graph_from_nested_dict_2layer(nest_d, belief_dict=None):
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

    dnf_logger.info('Building directed simple graph from two layered nested '
                    'dict.')
    nx_dir_g = nx.DiGraph()

    if not belief_dict:
        dnf_logger.info('No Belief Score dictionary provided')
        dnf_logger.warning('API belief score checkup is not implemented yet')
        pass  # ToDo connect to API, calculate belief or use stmt belief

    # Flag to check if the statement dict has belief score in it or not
    has_belief = False
    for k1, d1 in nest_d.items():
        for k2, v in d1.items():
            for tups in v:
                # Has to be (type, hash) or (type, hash, belief score)
                try:
                    assert len(tups) == 2 or len(tups) == 3
                except AssertionError:
                    pdb.set_trace()  # Check what tups is
                if len(tups) == 3:
                    has_belief = True
                break

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
                                bs = belief_dict[str(hsh)]
                            except KeyError:
                                dnf_logger.warning('No entry found in belief '
                                                   'dict for hash %s' %
                                                   str(hsh))
                                bs = 1
                            t = (typ, hsh, bs)
                            inner_obj_b.append(t)
                    else:
                        inner_obj_b = inner_obj
                    nx_dir_g.add_edge(u_of_edge=subj,
                                      v_of_edge=obj,
                                      attr_dict={'stmts': inner_obj_b})
    return nx_dir_g


def nx_directed_graph_from_nested_dict_3layer(nest_d):
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

    dnf_logger.info('Building directed simple graph from three layered nested '
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


def nx_undirected_graph_from_nested_dict(nest_d):
    """Returns an undirected graph built from a nested dict of statements

    Use this function to build a simple undirected graph. Suitable when the
    goal is to generate a node-edge graph for plain visualization.

    nest_d : defaultdict(dict)
        A nested dict with two or more layers

    Returns
    -------
    nx_undir : networkx.classes.graph.Graph
        An undirected, unweighted networkx graph
    """

    nx_undir = nx.Graph()

    # Create queue from nested dict
    ndq = list(nest_d.items())

    dnf_logger.info('Building undirected graph from nested dict of statements')
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


def nx_graph_from_corr_pd_series(corr_sr, source='id1', target='id2',
                                 edge_attr='correlation', use_abs_corr=False):
    """Return a graph from a pandas sereis containing correlaton between gene A
    and gene B, using the correlation as edge weight

    corr_sr : pandas Series or DataFrame
        Pandas Series/DataFrame containing A, B corr

    source : str
        which column to identify as source node (output is still undirected)

    target : str
        which column to identify as target nodes (output is still undirected)

    edge_attr : int or str
        Column to use for edge attributes
    absolute : Bool
        Use absolute value as edge weight. Otherwise magnitude is used.

    Returns
    -------
    corr_weight_graph : nx.Graph
        An undirected, weighted, networkx graph
    """

    # check if corr_sr is series or dataframe
    if type(corr_sr) == pd_Series_class:
        dnf_logger.info('Converting Pandas Series to Pandas DataFrame')
        corr_df = pd.DataFrame(corr_sr).reset_index()
    else:
        corr_df = corr_sr

    corr_df = corr_df.rename(
        columns={'level_0': source, 'level_1': target, 0: edge_attr})

    if use_abs_corr:
        dnf_logger.info('Using absolute correlation values')
        corr_df = corr_df.apply(lambda c: c.abs() if np.issubdtype(
            c.dtype, np.number) else c)

    if type(edge_attr) is list or type(edge_attr) is bool:
        dnf_logger.warning('More than one attribute might be added to edges. '
                           'Resulting networkx graph might not be usable as '
                           'simple weighted graph.')
    dnf_logger.info('Creating weighted undirected graph from network data')
    corr_weight_graph = nx.from_pandas_dataframe(df=corr_df,
                                                 source=source,
                                                 target=target,
                                                 edge_attr=edge_attr)
    return corr_weight_graph


def nx_graph_from_corr_tuple_list(corr_list, use_abs_corr=False):
    """Return a graph from a list of edges, using the correlation as weight

    corr_list : list or iterator
        Edge tuples

    absolute : Bool
        Use absolute value as edge weight. Otherwise magnitude is used.

    Returns
    -------
    corr_weight_graph : nx.Graph
        An undirected, weighted, networkx graph
    """
    corr_weight_graph = nx.Graph()

    if use_abs_corr:
        dnf_logger.info('Using absolute correlation values')
        corr_list = map(lambda t: (t[0], t[1], abs(t[2])), corr_list)

    dnf_logger.info('Converting tuples to an edge bunch')
    edge_bunch = map(lambda t: (t[0], t[1], {'weight': t[2]}), corr_list)

    dnf_logger.info('Creating weighted undirected graph from network data')
    corr_weight_graph.add_edges_from(ebunch_to_add=edge_bunch)

    return corr_weight_graph


def _read_gene_set_file(gf, data):
    gset = []
    with open(gf, 'rt') as f:
        for g in f.readlines():
            gn = g.upper().strip()
            if gn in data:
                gset.append(gn)
    return gset


def _map2index(start, binsize, value):
    offset = int(abs(start//binsize))
    return offset + int(float(value) // binsize)


def histogram_for_large_files(fpath, number_of_bins, binsize, first):
    """Returns a histogram for very large files

    fpath: str(filename)
        filepath to file with data to be binned
    number_of_bins: int
        the number fo bins to use
    binsize: float
        the size of bins
    first: float
        The left most (min(x)) edge of the bin edges

    Returns
    -------
    home_brewed_histo: np.array
        A histrogram of the data in fpath according to number of bins,
        binsize and first.
    """
    home_brewed_histo = np.zeros(number_of_bins, dtype=int)
    with open(file=fpath) as fo:
        for line in fo:
            flt = line.strip()
            home_brewed_histo[_map2index(start=first, binsize=binsize,
                                         value=flt)] += 1
    return home_brewed_histo


def _manually_add_to_histo(hist, start, binsize, value):
    hist[_map2index(start, binsize, value)] += 1


def _my_gauss(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


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


def get_gene_gene_corr_dict_wstats(tuple_generator, nbins, binsize, hist_range):
    """Returns a gene-gene correlation nested dict given a gene-gene nested dict

    tuple_generator : generator object
        A generator object
    nbins: int
        number of bins in histogram
    binsize: float
        size of bins in histogram
    hist_range: tuple(float, float)
        Tuple of floats with start-stop

    Returns
    -------
    nested_dict : defaultdict(dict)
        Dict with gene-gene-correlation
    """
    corr_nest_dict = create_nested_dict()
    home_brewed_histo = np.zeros(nbins, dtype=int)
    edges = np.arange(hist_range[0], hist_range[1] + binsize, step=binsize)
    dnf_logger.info('Generating correlation lookup and statistics')
    mean = 0
    for gene1, gene2, c in tuple_generator:
        corr = float(c)

        mean = (mean + corr)/2
        corr_nest_dict[gene1][gene2] = corr
        home_brewed_histo[_map2index(hist_range[0], binsize, corr)] += 1

    try:
        assert sum(home_brewed_histo) > 0
    except AssertionError:
        pdb.set_trace()  # Why is sum(home_brewed_histo) == 0
    a, mu, sigma = get_gaussian_stats(bin_edges=edges, hist=home_brewed_histo)
    sigma_dict = {'a': a,
                  'gauss_mean': mu,
                  'gauss_sigma': sigma,
                  'real_mean': mean}
    return corr_nest_dict, sigma_dict


def merge_correlation_dicts(correlation_dicts_list, settings):
    """Merge multiple correlation data sets to one single iterable of
    (gene, gene, correlation_dict)

    correlation_dicts_list: list[(gene_set_name, corr_nest_dict)]
        List of name-corr_dict-distr_stats_dict tuples
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

    # Get shortest dict from tuple list
    name_dict_sigma_tuple = next(d for d in correlation_dicts_list if len(
        d[1]) == min([len(i[1]) for i in correlation_dicts_list]))
    # Remove the tuple with the shortest dict  from the list
    correlation_dicts_list.remove(name_dict_sigma_tuple)
    set_name, shortest_dict, sigma_dict = name_dict_sigma_tuple
    # Pop the list to get the oter tuple to merge with
    other_name, other_dict, other_sigma_dict = correlation_dicts_list.pop()

    merged_corr_dict = create_nested_dict()
    dnf_logger.info('Merging correlation dicts %s and %s' %
                    (set_name, other_name))

    # Loop shortest correlation lookup dict
    for o_gene, d in shortest_dict.items():
        for i_gene, corr in d.items():
            if not _entry_exist(merged_corr_dict, o_gene, i_gene):
                # Check both directions
                other_corr = None
                if _entry_exist(other_dict, o_gene, i_gene):
                    other_corr = other_dict[o_gene][i_gene]
                elif _entry_exist(other_dict, i_gene, o_gene):
                    other_corr = other_dict[i_gene][o_gene]

                if other_corr and pass_filter(
                        corr1=corr,
                        sigma1=sigma_dict['gauss_sigma'],
                        corr2=other_corr,
                        sigma2=other_sigma_dict['gauss_sigma'],
                        margin=settings['margin'],
                        filter_type=settings['filter_type']):
                    merged_corr_dict[o_gene][i_gene][set_name] = corr
                    merged_corr_dict[o_gene][i_gene][other_name] = other_corr
                    assert merged_corr_dict[o_gene][i_gene].get(set_name) is \
                        not None
                    assert merged_corr_dict[o_gene][i_gene].get(other_name) is \
                        not None
                else:
                    continue

    # recursive case: more than 2 dicts. Merge the shortest and another,
    # call same function with and the merged dict plus the rest of the list.
    # if correlation_dicts_list:
    #     return merge_correlation_dicts_recursive(
    #         correlation_dicts_list.append(('master', merged_corr_dict)))

    # base case: 2 dicts in list. Return merged dict.
    # else:
    return merged_corr_dict


# def merge_correlation_dicts_recursive(correlation_dicts_list):
#     """ This should be the recursive version of correlation_dicts_list(). Call
#     this function from correlation_dicts_list()
#
#     :param correlation_dicts_list:
#     :return: pass
#     """
#     pass


def get_combined_correlations(dict_of_data_sets, filter_settings):
    """Return a combined dict of correlations given multiple gene data sets

    The input data [needs to be a collection of the gene expression data set.
    Need to ask user to list all sources]

    The input data set dict has the following format:

        dict_of_data_sets[gene_set1] = dataset1_dict

        dataset1_dict = {data: (depmap filepath),
                         corr: (depmap corr file),
                         filter_gene_set: list[genes to filter on]
                         unique_pair_corr_file: filepath,
                         ll: lower_limit_for_correlation,
                         ul: upper_limit_for_correlation,
                         sigma: st-dev of correlation distr,
                         filter_margin: float,
                         merge_filter_type: str,
                         outbasename: str,
                         strict: Bool,
                         recalc: Bool}

    The filter settings should contain the following:

        filter_settings = {'margin':      diff in terms of standard deviations
                                          between correlations,
                           'filter_type': Type of filtering
                                          (Default: 'sigma-diff'),
                           'nbins':       Number of bins to use when binning
                                          the correlation data,
                           'binsize':     Bin size when binning the
                                          correlation data,
                           'hist_range':  Upper and lower edges when
                                          binning the correlation data
                                          (Default: -1.0, 1.0)}

    The returned master correlation dict has the following format:

        d[gene][gene] = {gene_set1: correlation,
                         gene_set2: correlation,
                         ...}

    dict_of_data_sets: dict()
        Dictionary containing the filepaths and settings for the data set
    filter_settings: dict
        Dictionary with filter settings

    Returns
    -------
    master_corr_dict: defaultdict(dict)
        A nested dict containing a lookup of the filtered set of
        gene-gene-correlations
     gene_set_intersection: set()
        The set of HGNC gene names in the master correlation lookup
    """
    corr_dicts_list = []
    gene_set_intersection = set()

    for gene_set_name, dataset_dict in dict_of_data_sets.items():
        dnf_logger.info('Processing set "%s"' % gene_set_name)

        # Get tuple generator and the accompanied set of genes
        tuple_generator, set_of_genes = get_correlations(
            depmap_data_file=dataset_dict['data'],
            geneset_file=dataset_dict['filter_gene_set'],  # [] for no set
            corr_file=dataset_dict['corr'],
            strict=dataset_dict['strict'],
            outbasename=dataset_dict['outbasename'],
            unique_pair_corr_file=dataset_dict['unique_pair_corr_file'],
            recalc=dataset_dict['recalc'],
            lower_limit=dataset_dict['ll'],
            upper_limit=dataset_dict['ul'])
        dnf_logger.info('Created tuple generator with %i unique genes from '
                        'set "%s"' % (len(set_of_genes), gene_set_name))

        # Generate correlation dict and get the statistics of the distribution
        corr_dict, sigma_dict = get_gene_gene_corr_dict_wstats(
            tuple_generator=tuple_generator,
            nbins=filter_settings['nbins'],
            binsize=filter_settings['binsize'],
            hist_range=filter_settings['hist_range'])
        dnf_logger.info('Created correlation dictionary of length %i for set '
                        '"%s"' % (len(corr_dict), gene_set_name))

        # Append correlation dict and stats to list
        corr_dicts_list.append((gene_set_name, corr_dict, sigma_dict))
        gene_set_intersection.intersection_update(set_of_genes)

    # Merge the dictionaries and the set of genes
    master_corr_dict = merge_correlation_dicts(corr_dicts_list,
                                               settings=filter_settings)
    dnf_logger.info('Merged gene sets to master dictionary of length %i' %
                    len(master_corr_dict))

    return master_corr_dict, gene_set_intersection


def get_correlations(depmap_data_file, geneset_file, corr_file, strict,
                     outbasename, unique_pair_corr_file, recalc=False,
                     lower_limit=0.2, upper_limit=1.0):
    # todo make function take gene set data dict as input?
    """ given a gene-feature data matrix in csv format.

    depmap_data_file: str
        Filepath to depmap data file to process
    geneset_file: str
        Filepath to a geneset to filter data to.
    corr_file: str
        Filepath to pre-calculated correlations of depmap_data_file.
    strict: Bool
        If True, all genes in correlations have to exist in geneset_file
    outbasename: str
        Basename to use for output files
    unique_pair_corr_file: str
        Filepath to csvfile with unique tuples of gene,gene,corr.
    recalc: Bool
        If True, recalculate correlations (has to be True if corr_file is None).
    lower_limit: float
        Lowest correlation magnitude to consider
    upper_limit: float
        Highest correlation magnitude to consider (good for picking a sample
        in the middle of the correlation distribution).

    Returns
    -------
    uniq_pair_gen: generator
        Generator of gene,gene,correlation tuples from file or correlation
        calculation
    all_hgnc_ids: set()
        The set of all genes in the correlation
    """

    filtered_correlation_matrix = _get_corr_df(
        depmap_data_file=depmap_data_file,
        corr_file=corr_file,
        geneset_file=geneset_file,
        strict=strict,
        recalc=recalc,
        outbasename=outbasename,
        lower_limit=lower_limit,
        upper_limit=upper_limit)

    all_hgnc_ids = set(filtered_correlation_matrix.index.values)

    # Return a generator object from either a loaded file or a pandas dataframe
    if unique_pair_corr_file:
        dnf_logger.info('Loading unique correlation pairs from %s' %
                        unique_pair_corr_file)
        uniq_pair_gen = csv_file_to_generator(unique_pair_corr_file,
                                              ['gene1', 'gene2', 'corr'])
    else:
        if lower_limit == 0:
            fname = outbasename + '_all_unique_correlation_pairs.csv'
        else:
            fname = outbasename + '_unique_correlation_pairs_ll%s.csv' % \
                    str(lower_limit).replace('.', '')

        uniq_pair_gen = corr_matrix_to_generator(filtered_correlation_matrix)
        uniq_pair_gen_forcsv = corr_matrix_to_generator(
            filtered_correlation_matrix)
        dnf_logger.info('Saving unique correlation pairs to %s. '
                        '(May take a while)' % fname)
        _dump_it_to_csv(fname, uniq_pair_gen_forcsv)
    return uniq_pair_gen, all_hgnc_ids


def _get_corr_df(depmap_data_file, corr_file, geneset_file, strict, recalc,
                 outbasename, lower_limit, upper_limit):

    # Open data file
    dnf_logger.info('Reading DepMap data from %s' % depmap_data_file)
    data = pd.read_csv(depmap_data_file, index_col=0, header=0)
    data = data.T

    if geneset_file:
        # Read gene set to look at
        gene_filter_list = _read_gene_set_file(gf=geneset_file, data=data)
    else:
        gene_filter_list = []  # Evaluates to False

    # 1. no loaded gene list OR 2. loaded gene list but not strict -> data.corr
    corr_matrix_df = None
    if not geneset_file or (geneset_file and not strict):
        # Calculate the full correlations, or load from cached
        if recalc:
            dnf_logger.info('Calculating correlations (may take a long time)')
            corr = data.corr()
            corr.to_hdf(outbasename+'correlations.h5', 'correlations')
        else:
            if corr_file:
                dnf_logger.info('Loading correlations from %s' % corr_file)
                corr = pd.read_hdf(corr_file, 'correlations')
            else:
                dnf_logger.error('No correlation file provdided or calculated!')
                raise FileNotFoundError
        # No gene set file, leave 'corr' intact
        if not geneset_file:
            corr_matrix_df = corr

        # Gene set file present: filter and unstack
        elif geneset_file and not strict:
            corr_matrix_df = corr[gene_filter_list]

    # 3. Strict: both genes in interaction must be from loaded set;
    #    Filter data, then calculate correlations and then unstack
    elif geneset_file and strict:
        corr_matrix_df = data[gene_filter_list].corr()

    assert corr_matrix_df is not None

    # Remove self correlation, correlations below ll, sort on magnitude,
    # leave correlation intact
    dnf_logger.info('Removing self correlations')
    corr_matrix_df = corr_matrix_df[corr_matrix_df != 1.0]  # Self correlations

    # Filter low correlations
    if lower_limit > 0:
        return corr_limit_filtering(corr_matrix_df, lower_limit, upper_limit)
    # No filtering
    else:
        dnf_logger.warning('No filtering requested. Be aware of large RAM '
                          'usage.')
        return corr_matrix_df


def corr_limit_filtering(corr_matrix_df, lower_limit, upper_limit):
    """Filters a correlation matrix to values in (lower_limit, upper_limit)

    corr_matrix_df: pandas.DataFrame
        A pandas correlation matrix as a pandas data frame
    lower_limit: float
        Lowest correlation magnitude to consider
    upper_limit: float
        Highest correlation magnitude to consider (good for picking a sample
        in the middle of the correlation distribution)

    Returns
    -------
    corr_matrix_df: pandas.DataFrame
        A filtered correlation dataframe matrix
    """
    dnf_logger.info('Filtering correlations to %.1f < C < %.1f' %
                    (lower_limit, upper_limit))
    corr_matrix_df = corr_matrix_df[corr_matrix_df.abs() > lower_limit]
    if upper_limit < 1.0:
        corr_matrix_df = corr_matrix_df[corr_matrix_df.abs() < upper_limit]

    return corr_matrix_df


def pass_filter(corr1, sigma1, corr2, sigma2, margin, filter_type='sigma-diff'):
    """Filter for passing correlation scores based on their difference in
    standard deviation

    corr1: float
        Correlation from first data set
    sigma1: float
        Standard deviation of first dataset
    corr2: float
        Correlation from second data set
    sigma2: float
        Standard deviation of second dataset
    margin: float
        How far off the correlations can be to pass as "similar"
    filter_type:
        The filter type to use ("sigma-diff" is currently the only one)

    Returns
    -------
    bool
        If True, the correlations are similar enough as measured by their
        similarity in the standard deviation.
    """
    if filter_type == 'sigma-diff':
        return abs(corr1/sigma1 - corr2/sigma2) > margin
    else:
        return True


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
        (subj, obj, stmt_type, stmt_hash) = row
        if _entry_exist(nest_hash_dict, subj, obj):
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
        dnf_logger.warning('No belief score dict is provided! Please provide a '
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

            # Check common parent (same familiy or complex)
            for agent, other_agent in itt.permutations(agent_list, r=2):
                if has_common_parent(id1=agent, id2=other_agent):
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

    dnf_logger.info('Created nested dict of length %i from %i statements.' %
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
        dnf_logger.info('Deduplicating statements and accounting for custom '
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


def _old_str_output(subj, obj, corr, stmts, ignore_str='parent'):

    # Build up a string that shows explanations for each connection
    output = 'subj: %s; obj: %s; corr: %f \n' % (subj, obj, corr)

    if ignore_str in stmts:
        pure_list = [s for s in stmts if type(s) is not str]
        types = relation_types(pure_list)
        types += [ignore_str] * stmts.count(ignore_str)
    else:
        types = relation_types(stmts)

    dedupl_stmts = stmts.copy()

    types_set = set(types)
    types_sstmt = []
    for tp in types_set:
        for st in dedupl_stmts:
            if type(st) is not str:
                if st.to_json()['type'] == tp:
                    types_sstmt.append((tp, str(st)))
            elif type(st) is str and tp is ignore_str:
                types_sstmt.append((tp, str(st)))

    for tp, str_stmt in types_sstmt:
        if tp is not ignore_str:
            output += '- - - - - - - - - - - - - - - - - - - - - - - - - - - -'
            output += '\nInstances found of statement %s: %i\n' % \
                      (str_stmt, types.count(tp))
        for stmt in stmts:
            if type(stmt) is str and str(stmt) == ignore_str:
                output += '%s and %s are in the same complex or family\n' % \
                          (subj, obj)
            elif type(stmt) is not str and stmt.to_json()['type'] == tp:
                output += 'Evidence for uuid %s: ' % stmt.uuid
                ev = stmt.evidence[0].text
                output += ('N/A' if ev is None else ev)+'\n'
            else:
                continue

    # Add separator between each connection
    output += '\n\n#### #### #### #### #### ####\n'
    return output


def str_output(subj, obj, corr, stmts, ignore_str='parent'):
    """Formats information about statements and returns a string.

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
    output += 'subj: %s; obj: %s; corr: %f \n' % (subj, obj, corr) + \
              'https://depmap.org/portal/interactive/?xDataset=Avana&xFeature' \
              '={}&yDataset=Avana&yFeature={}&colorDataset=lineage' \
              '&colorFeature=all&filterDataset=context&filterFeature=' \
              '&regressionLine=false&statisticsTable=false&associationTable=' \
              'true&plotOnly=false\n'.format(subj, obj)

    pa_stmts = stmts.copy()

    for stmt in pa_stmts:
        output += '- - - - - - - - - - - - - - - - - - - - - - - - - - - -\n'
        if type(stmt) is str and str(stmt) == ignore_str:
            output += '%s and %s are in the same complex or family\n' % \
                      (subj, obj)
        else:
            # Remove duplicate evidence text
            ev_text_set = set(['N/A' if ev.text is None else ev.text for ev in
                               stmt.evidence])
            ev_text_list = list(ev_text_set)
            if 'N/A' in ev_text_list:
                ev_text_list.remove('N/A')

            output += '\nInstances found of statement %s: %i; supports ' \
                      'count: %i; Supported by count: %i\n' % \
                      (str(stmt), len(ev_text_list), len(stmt.supports),
                       len(stmt.supported_by))
            for count, ev_text in enumerate(ev_text_list):
                output += 'Evidence %i: ' % (count+1) + ev_text + '\n'

    # Add separator between each connection
    output += '\n\n#### #### #### #### #### ####\n'
    return output


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
    #     # pdb.set_trace()
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


def dbc_load_statements(hgnc_ids):
    """Load statements where hgnc id is subject or object from indra.db.client

    Parameters
    ----------
    hgnc_ids : iterable
        An iterable containing HGNC ids

    Returns
    -------
    stmts : set{:py:class:`indra.statements.Statement`}
        A set of all retrieved INDRA statemetents containing HGNC id
    """
    stmts = set()
    counter = 0
    n_hgnc_ids = len(hgnc_ids)
    try:
        for hgnc_id in hgnc_ids:
            stmts.update(dbc.get_statements_by_gene_role_type(agent_id=hgnc_id,
                                                              db=db_prim,
                                                              preassembled=
                                                              False,
                                                              fix_refs=False))
            counter += 1
            if counter % max([10, 10 ** ceil(log10(n_hgnc_ids)) // 100]) == 0:
                dnf_logger.info(' : : : Finished %i queries out of %i '
                                ': : :' % (counter, n_hgnc_ids))

    except KeyboardInterrupt as e:
        db_prim.session.rollback()
        raise e
    except StatementError as e:
        db_prim.session.rollback()
        raise e
    return stmts


def find_parent(ho=hm.hierarchies['entity'], ns='HGNC',
                id_=None, type_='all'):
    """A wrapper function for he.get_parents to make the functionality more
    clear.

    Parameters
    ----------
    ho : HierarchyManager object
        A HierarchyManager object. Default: entity hierarchy object
    ns : str
        namespace id. Default: HGNC
    id_ : str
        id to check parents for. Default: None
    type_ : str
        'all': (Default) return all parents irrespective of level;
        'immediate': return only the immediate parents;
        'top': return only the highest level parents

    Returns
    -------
    set
        set of parents of database id in namespace ns
    """
    return ho.get_parents(ho.get_uri(ns, id_), type_)


def common_parent(ho=hm.hierarchies['entity'], ns1='HGNC',
                  id1=None, ns2='HGNC', id2=None, type_='all'):
    """Returns the set of common parents.

    Parameters
    ----------
    ho : HierarchyManager object
        A HierarchyManager object. Default: entity hierarchy object
    ns1 : str
        namespace id. Default: HGNC
    id1 : str
        First id to check parents for. Default: None
    ns2 : str
        namespace id. Default: HGNC
    id2 : str
        Second id to check parents for. Default: None
    type_ : str
        'all': (Default) return all parents irrespective of level;
        'immediate': return only the immediate parents;
        'top': return only the highest level parents

    Returns
    -------
    set
        set of common parents in uri format
    """
    return find_parent(ho, ns1, id1, type_) & find_parent(ho, ns2, id2, type_)


def has_common_parent(ho=hm.hierarchies['entity'], ns1='HGNC', id1=None,
                      ns2='HGNC', id2=None, type='all'):

    """Returns True if id1 and id2 has at least one common parent.

    Parameters
    ----------
    ho : HierarchyManager object
        A HierarchyManager object. Default: entity hierarchy object
    ns1 : str
        namespace id. Default: HGNC
    id1 : str
        First id to check parents for. Default: None
    ns2 : str
        namespace id. Default: HGNC
    id2 : str
        Second id to check parents for. Default: None
    type : str
        'all': return all parents irrespective of level;
        'immediate': return only the immediate parents;
        'top': return only the highest level parents

    Returns
    -------
    bool
        True if hgnc1 and hgnc2 has one or more common parents.
    """
    return bool(common_parent(ho, ns1, id1, ns2, id2, type))


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
        stmts = capi.get_statements(subject=id1, object=id2, on_limit=on_limit)
        stmts + capi.get_statements(subject=id2, object=id1, on_limit=on_limit)
    except IndraDBRestError:
        stmts = capi.get_statements(subject=id1 + '@TEXT', object=id2 + '@TEXT',
                                    on_limit=on_limit)
        stmts + capi.get_statements(subject=id2 + '@TEXT', object=id1 + '@TEXT',
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
        True if the HGNC ids has a direct relation found in the
        indra.sources.indra_db_rest.client_api databases.
    """
    return bool(direct_relation(id1, id2, long_stmts=long_stmts))


def are_connected(id1, id2, long_stmts=set()):
    """Indicates whether two genes have a connection either through a direct
    relation or a through a common parent.

    Parameters
    ----------
    id1/i2 : str
        HGNC id

    Returns
    -------
    bool
        True if the two HGNC ids either have a common parent or if they have a
        direct relation found in the indra.sources.indra_db_rest.client_api
        databases.
    """
    return has_common_parent(ns1='HGCN', id1=id1, ns2='HGCN', id2=id2) or \
        has_direct_relation(id1=id1, id2=id2, long_stmts=long_stmts)


def connection_types(id1, id2, long_stmts=set()):
    """Returns a list of the connection types linking two genes.

    Parameters
    ----------
    id1/i2 : str
        HGNC id

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
    if has_common_parent(id1=id1, id2=id2):
        ctypes += ['parent']
    return ctypes
