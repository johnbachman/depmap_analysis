import numpy as np
import pandas as pd

from depmap_analysis.scripts.depmap_script2 import _get_pairs, _down_sample_df


def _gen_sym_matrix(size):
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
    a = _gen_sym_matrix(size)

    # Put NaN's in two symmetrical, but random positions
    row, col = _get_off_diag_pair(size)
    a.iloc[row, col] = np.nan
    a.iloc[col, row] = np.nan

    n_pairs = _get_pairs(a)
    assert n_pairs == int((size**2 - size - 2) / 2)

    # Test having some of the diagonal elements gone as well, this should
    # not affect the count since diagonal elements should be ignored
    d = np.random.randint(0, size)
    a.iloc[d, d] = np.nan
    n_pairs = _get_pairs(a)
    assert n_pairs == int((size**2 - size - 2) / 2)


def test_down_sampling():
    size = 40
    a = _gen_sym_matrix(size)

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
    assert goal_pairs <= _get_pairs(a) <= 1.1 * goal_pairs, _get_pairs(a)
