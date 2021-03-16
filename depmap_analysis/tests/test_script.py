from typing import Dict, Callable

import numpy as np
import pandas as pd
import networkx as nx

from indra.util import batch_iter
from indra.databases.hgnc_client import uniprot_ids, hgnc_names
from depmap_analysis.util.io_functions import file_opener
from depmap_analysis.post_processing import *
from depmap_analysis.network_functions.depmap_network_functions import \
    corr_matrix_to_generator, get_pairs, get_chunk_size, down_sample_df
from depmap_analysis.scripts.depmap_script2 import _match_correlation_body, \
    expl_columns, id_columns
from depmap_analysis.scripts.depmap_script_expl_funcs import *
from . import *

reverse_uniprot = {v: k for k, v in uniprot_ids.items()}


def _gen_sym_df(size):
    # Get square, symmetric matrix in dataframe
    m = np.random.rand(size, size)
    m = (m + m.T) / 2
    np.fill_diagonal(m, 1.)
    return pd.DataFrame(m)


def _get_off_diag_pair(max_index: int):
    if max_index == 0:
        raise ValueError('Cannot have max_index == 0')
    r = np.random.randint(0, max_index)
    c = np.random.randint(0, max_index)
    while r == c:
        c = np.random.randint(0, max_index)
    return r, c


def test_df_pair_calc():
    size = 4
    a = _gen_sym_df(size)

    # Put NaN's in two symmetrical, but random positions
    row, col = _get_off_diag_pair(size)
    a.iloc[row, col] = np.nan
    a.iloc[col, row] = np.nan

    n_pairs = get_pairs(a)
    assert n_pairs == int((size**2 - size - 2) / 2)

    # Test having some of the diagonal elements gone as well, this should
    # not affect the count since diagonal elements should be ignored
    d = np.random.randint(0, size)
    a.iloc[d, d] = np.nan
    n_pairs = get_pairs(a)
    assert n_pairs == int((size**2 - size - 2) / 2)


def test_down_sampling():
    size = 40
    a = _gen_sym_df(size)

    pairs = set()
    for n in range(5):
        row, col = _get_off_diag_pair(size)
        while (row, col) in pairs:
            row, col = _get_off_diag_pair(size)
        a.iloc[row, col] = np.nan
        a.iloc[col, row] = np.nan
        pairs.add((row, col))
        pairs.add((col, row))

    goal_pairs = 10
    a = down_sample_df(z_corr=a, sample_size=goal_pairs)
    assert goal_pairs <= get_pairs(a) <= 1.1 * goal_pairs, get_pairs(a)


def test_iterator_slicing():
    size = 50
    a = _gen_sym_df(size)

    pairs = set()
    n = 0
    for n in range(size):
        k = 0
        row, col = _get_off_diag_pair(size)
        while (row, col) in pairs:
            row, col = _get_off_diag_pair(size)
            k += 1
            if k > 1000:
                print('Too many while iterations, breaking')
                break
        if k > 1000:
            break
        a.iloc[row, col] = np.nan
        a.iloc[col, row] = np.nan
        pairs.add((row, col))
        pairs.add((col, row))

    pairs_removed = n + 1

    # Assert that we're correct so far

    # Get total pairs available:
    total_pairs = get_pairs(a)

    # all items - diagonal - all removed items off diagonal
    assert (size**2 - size - 2*pairs_removed) / 2 == total_pairs

    # Check that the iterator slicing for multiprocessing runs through all
    # the pairs

    # Chunks wanted
    chunks_wanted = 10

    chunksize = get_chunk_size(chunks_wanted, total_pairs)

    chunk_iter = batch_iter(iterator=corr_matrix_to_generator(a),
                            batch_size=chunksize, return_func=list)

    pair_count = 0
    chunk_ix = 0
    for chunk_ix, list_of_pairs in enumerate(chunk_iter):
        pair_count += len([(t[0][0], t[0][1], t[1]) for t in
                           list_of_pairs if t is not None])

    # Were all pairs looped?
    assert pair_count == total_pairs, \
        f'pair_count={pair_count} total_pairs={total_pairs}'
    # Does the number of loop iterations correspond to the number of chunks
    # wanted?
    assert chunk_ix + 1 == chunks_wanted, \
        f'chunk_ix+1={chunk_ix + 1}, chunks_wanted={chunks_wanted}'

    # Redo the same with subset of names
    name_subset = list(np.random.choice(a.columns.values,
                                        size=size // 3,
                                        replace=False))
    # Add a name that does not exist in the original df
    name_subset.append(size+2)

    # Get total pairs available
    total_pairs_permute = get_pairs(a, subset_list=name_subset)

    # Chunks wanted
    chunks_wanted = 10

    chunksize = get_chunk_size(chunks_wanted, total_pairs_permute)

    chunk_iter = batch_iter(
        iterator=corr_matrix_to_generator(a, subset_list=name_subset),
        batch_size=chunksize,
        return_func=list
    )

    pair_count = 0
    chunk_ix = 0
    for chunk_ix, list_of_pairs in enumerate(chunk_iter):
        pair_count += len([(t[0][0], t[0][1], t[1]) for t in
                           list_of_pairs if t is not None])

    # Were all pairs looped?
    assert pair_count == total_pairs_permute, \
        f'pair_count={pair_count} total_pairs={total_pairs_permute}'

    # Does the number of loop iterations correspond to the number of chunks
    # wanted?
    assert chunk_ix + 1 == chunks_wanted, \
        f'chunk_ix+1={chunk_ix + 1}, chunks_wanted={chunks_wanted}'


def test_sampling():
    size = 50
    a = _gen_sym_df(size)

    pairs = set()
    n = 0
    for n in range(size):
        k = 0
        row, col = _get_off_diag_pair(size)
        while (row, col) in pairs:
            row, col = _get_off_diag_pair(size)
            k += 1
            if k > 1000:
                print('Too many while iterations, breaking')
                break
        if k > 1000:
            break
        a.iloc[row, col] = np.nan
        a.iloc[col, row] = np.nan
        pairs.add((row, col))
        pairs.add((col, row))

    pairs_removed = n + 1
    # Assert that we're correct so far

    # Get total pairs available:
    total_pairs = get_pairs(a)

    # all items - diagonal - all removed items off diagonal
    assert (size**2 - size - 2*pairs_removed) / 2 == total_pairs

    # Check that the iterator slicing for multiprocessing runs through all
    # the pairs

    # Test max pairs
    max_pairs = total_pairs//2
    total_pairs = max_pairs
    chunks_wanted = 10

    chunksize = get_chunk_size(chunks_wanted, total_pairs)

    chunk_iter = batch_iter(
        iterator=corr_matrix_to_generator(a, max_pairs=total_pairs),
        batch_size=chunksize, return_func=list
    )

    pair_count = 0
    chunk_ix = 0
    for chunk_ix, list_of_pairs in enumerate(chunk_iter):
        pair_count += len([(t[0][0], t[0][1], t[1]) for t in
                           list_of_pairs if t is not None])

    # Were all pairs looped?
    assert pair_count == total_pairs, \
        f'pair_count={pair_count} total_pairs={total_pairs}'
    # Does the number of loop iterations correspond to the number of chunks
    # wanted?
    assert chunk_ix + 1 == chunks_wanted, \
        f'chunk_ix+1={chunk_ix + 1}, chunks_wanted={chunks_wanted}'

    # Redo the same with subset of names
    name_subset = list(np.random.choice(a.columns.values,
                                        size=size // 3,
                                        replace=False))
    # Add a name that does not exist in the original df
    name_subset.append(size+2)

    # Get total pairs available
    total_pairs_permute = get_pairs(a, subset_list=name_subset)
    max_pairs = total_pairs_permute // 4
    total_pairs_permute = max_pairs

    # Chunks wanted
    chunks_wanted = 10

    chunksize = get_chunk_size(chunks_wanted, total_pairs_permute)

    chunk_iter = batch_iter(
        iterator=corr_matrix_to_generator(a,
                                          max_pairs=total_pairs_permute,
                                          subset_list=name_subset),
        batch_size=chunksize,
        return_func=list
    )

    pair_count = 0
    chunk_ix = 0
    for chunk_ix, list_of_pairs in enumerate(chunk_iter):
        pair_count += len([(t[0][0], t[0][1], t[1]) for t in
                           list_of_pairs if t is not None])

    # Were all pairs looped?
    assert pair_count == total_pairs_permute, \
        f'pair_count={pair_count} total_pairs={total_pairs_permute}'

    # Does the number of loop iterations correspond to the number of chunks
    # wanted?
    assert chunk_ix + 1 == chunks_wanted, \
        f'chunk_ix+1={chunk_ix + 1}, chunks_wanted={chunks_wanted}'


def test_depmap_script():
    up2path, _, pathid2pathname = file_opener(
        's3://depmap-analysis/misc_files/reactome_pathways.pkl')
    reactome_dict = {'uniprot_mapping': up2path,
                     'pathid_name_mapping': pathid2pathname}
    df = get_df()
    idg = get_dg()

    up1 = 'A0A075B6P5'
    up2 = 'A5LHX3'
    hgnc_id1 = reverse_uniprot[up1]
    hgnc_name1 = hgnc_names[hgnc_id1]
    hgnc_id2 = reverse_uniprot[up2]
    hgnc_name2 = hgnc_names[hgnc_id2]

    idg.add_node(hgnc_name1, ns='HGNC', id=hgnc_id1)
    idg.add_node(hgnc_name2, ns='HGNC', id=hgnc_id2)
    not_in_graph = 'not_in_graph'

    # Make correlation matrix with all combinations from the df pairs
    all_names = list(set(df.agA_name.values) | set(df.agB_name.values)) + \
        [not_in_graph, hgnc_name1, hgnc_name2]
    all_names.sort()
    corr_m = _gen_sym_df(len(all_names))
    corr_m.columns = all_names
    corr_m.index = all_names

    func_names = ['expl_ab', 'expl_ba', 'expl_axb', 'expl_bxa', 'get_sr',
                  'get_st', react_funcname]

    func_map = {funcname_to_colname[fname]: expl_functions[fname]
                for fname in func_names}
    bool_columns = ('not_in_graph', 'explained') + tuple(func_map.keys())
    stats_columns = id_columns + bool_columns

    _type = 'unsigned'

    corr_pairs = corr_matrix_to_generator(corr_m)
    stats_dict, expl_dict = _match_correlation_body(
        corr_iter=corr_pairs, expl_types=func_map,
        stats_columns=stats_columns, expl_cols=expl_columns,
        bool_columns=bool_columns, _type=_type,
        return_unexplained=False, reactome_dict=reactome_dict,
        local_indranet=idg, apriori_explained=None
    )

    assert set(stats_columns) == set(stats_dict.keys())
    assert set(expl_columns) == set(expl_dict.keys())

    expl_df = pd.DataFrame(expl_dict)
    stats_df = pd.DataFrame(stats_dict)

    # Test content
    # Any connection with not_in_graph should be
    assert all(b for b in
               stats_df[(stats_df.agA == not_in_graph) |
                        (stats_df.agB == not_in_graph)].not_in_graph), \
        str([b for b in stats_df[(stats_df.agA == not_in_graph) |
                                 (stats_df.agB == not_in_graph)].not_in_graph])

    assert all(np.isnan(b) for b in
               stats_df[(stats_df.agA == not_in_graph) |
                        (stats_df.agB == not_in_graph)].explained)

    expected = {'not_in_graph': False,
                'explained': True,
                ab_colname: False,
                ba_colname: False,
                axb_colname: False,  # Not True, as pairs go alphabetically
                bxa_colname: True,  # True from testing Y2,Z2
                sr_colname: False,
                st_colname: False,
                react_colname: False}
    p = 'Y2_Z2'
    res = stats_df[list(bool_columns)][stats_df.pair == p].to_dict(
        orient='records')[0]
    for k, b in res.items():
        assert b == expected[k]

    assert expl_df[
               (expl_df.pair == p) & (expl_df.expl_type == bxa_colname)
    ].expl_data.values[0] == ['X2']
    assert len(expl_df[(expl_df.pair == p) &
                       (expl_df.expl_type == sr_colname)]) == 0
    assert len(expl_df[(expl_df.pair == p) &
                       (expl_df.expl_type == st_colname)]) == 0

    expected = {'not_in_graph': False,
                'explained': True,
                ab_colname: False,
                ba_colname: False,
                axb_colname: False,
                bxa_colname: False,
                sr_colname: True,
                st_colname: True,
                react_colname: False}
    p = 'X1_X2'
    res: Dict = stats_df[list(bool_columns)][stats_df.pair == p].to_dict(
        orient='records')[0]
    for k, b in res.items():
        assert b == expected[k]

    assert expl_df[
               (expl_df.pair == p) & (expl_df.expl_type == sr_colname)
    ].expl_data.values[0][2] == ['Z2']
    assert expl_df[
               (expl_df.pair == p) & (expl_df.expl_type == st_colname)
    ].expl_data.values[0][2] == ['Z1']

    # Check that reactome is explained, and not counted as among the explained
    len_react = len(stats_df[stats_df[react_colname] == True])
    assert len_react == 1, len_react
    len_react = len(stats_df[(stats_df[react_colname] == True) &
                             (stats_df.explained == False)])
    assert len_react == 1, len_react

    # Test getting interesting df
    interesting_df = get_non_reactome_axb_expl_df(graph=idg, stats_df=stats_df,
                                                  expl_df=expl_df,
                                                  z_corr=corr_m)
    assert len(interesting_df) == 5
    assert set(interesting_df.pair) == {'X1_X2', 'Y1_Z2', 'Y2_Z2', 'Z1_Z2'}


def test_reactome_expl():
    up2path, _, pathid2pathname = file_opener(
        's3://depmap-analysis/misc_files/reactome_pathways.pkl')
    reactome_dict = {'uniprot_mapping': up2path,
                     'pathid_name_mapping': pathid2pathname}

    react_func: Callable = expl_functions[react_funcname]
    up1 = 'A0A075B6P5'
    up2 = 'A5LHX3'
    res = {'R-HSA-2871837'}
    descr = ['FCERI mediated NF-kB activation']
    assert res == set(up2path[up1]) & set(up2path[up2])

    hgnc_id1 = reverse_uniprot[up1]
    hgnc_name1 = hgnc_names[hgnc_id1]
    hgnc_id2 = reverse_uniprot[up2]
    hgnc_name2 = hgnc_names[hgnc_id2]
    func_args = (hgnc_name1, hgnc_name2, 0.0, nx.DiGraph(), 'unsigned',
                 reactome_dict)
    s, o, explained, data = react_func(*func_args)
    assert explained
    assert data == descr, str(data or 'None returned')
