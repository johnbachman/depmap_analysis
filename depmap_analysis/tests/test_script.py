import numpy as np
import pandas as pd
import networkx as nx
from typing import Optional

from indra.util import batch_iter
from depmap_analysis.network_functions.depmap_network_functions import \
    corr_matrix_to_generator, get_pairs, get_chunk_size
from depmap_analysis.scripts.depmap_script2 import _down_sample_df, \
    _match_correlation_body, expl_columns, id_columns
from depmap_analysis.scripts.depmap_script_expl_funcs import *
from . import *


indranet: Optional[nx.DiGraph] = None


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
    a = _down_sample_df(z_corr=a, sample_size=goal_pairs)
    assert goal_pairs <= get_pairs(a) <= 1.1 * goal_pairs, get_pairs(a)


def test_iterator_slicing():
    size = 50
    a = _gen_sym_df(size)

    pairs = set()
    n = 0
    for n in range(50):
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
    assert (size**2 - size - 2*pairs_removed) / 2 == get_pairs(a)

    # Check that the iterator slicing for multiprocessing runs through all
    # the pairs

    # Get total pairs available
    total_pairs = get_pairs(a)

    # Chunks wanted
    chunks_wanted = 10

    chunksize = get_chunk_size(chunks_wanted, total_pairs)

    chunk_iter = batch_iter(iterator=corr_matrix_to_generator(a),
                            batch_size=chunksize, return_func=list)

    pair_count = 0
    chunk_ix = 0
    for chunk_ix, list_of_pairs in enumerate(chunk_iter):
        pair_count += len([t for t in list_of_pairs if t is not None])

    # Were all pairs looped?
    assert pair_count == total_pairs, \
        f'pair_count={pair_count} total_pairs={total_pairs}'
    # Does the number of loop iterations correspond to the number of chunks
    # wanted?
    assert chunk_ix + 1 == chunks_wanted, \
        f'chunk_ix+1={chunk_ix + 1}, chunks_wanted={chunks_wanted}'


def test_depmap_script():
    df = get_df()
    idg = get_dg()

    not_in_graph = 'not_in_graph'
    # Make correlation matrix with all combinations from the df pairs
    all_names = list(set(df.agA_name.values) | set(df.agB_name.values)) + \
        [not_in_graph]
    all_names.sort()
    corr_m = _gen_sym_df(len(all_names))
    corr_m.columns = all_names
    corr_m.index = all_names

    func_names = ['expl_ab', 'expl_ba', 'expl_axb', 'expl_bxa', 'get_sr',
                  'get_st']
    bool_columns = ('not_in_graph', 'explained') + tuple(func_names)
    stats_columns = id_columns + bool_columns
    _type = 'unsigned'

    func_map = {fname: expl_functions[fname] for fname in func_names}

    corr_pairs = corr_matrix_to_generator(corr_m)
    stats_dict, expl_dict = _match_correlation_body(
        corr_iter=corr_pairs, expl_types=func_map,
        stats_columns=stats_columns, expl_cols=expl_columns,
        bool_columns=bool_columns, expl_mapping={}, _type=_type,
        strict_intermediates=True, local_indranet=idg
    )

    assert set(stats_columns) == set(stats_dict.keys())
    assert set(expl_columns) == set(expl_dict.keys())

    expl_df = pd.DataFrame(expl_dict)
    stats_df = pd.DataFrame(stats_dict)

    assert all(b for b in
               stats_df[(stats_df.agA == not_in_graph) |
                        (stats_df.agB == not_in_graph)].not_in_graph), \
        str([b for b in stats_df[(stats_df.agA == not_in_graph) |
                                 (stats_df.agB == not_in_graph)].not_in_graph])
    assert all(np.isnan(b) for b in
               stats_df[(stats_df.agA == not_in_graph) |
                        (stats_df.agB == not_in_graph)].explained)
