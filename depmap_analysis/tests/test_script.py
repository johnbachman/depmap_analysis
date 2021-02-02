import numpy as np
import pandas as pd

from depmap_analysis.scripts.depmap_script2 import _get_pairs


def test_df_pair_calc():
    def _get_off_diag_pair(max_index: int):
        if max_index == 0:
            raise ValueError('Cannot have max_index == 0')
        r = np.random.randint(0, max_index)
        c = np.random.randint(0, max_index)
        while r == c:
            c = np.random.randn(0, max_index)
        return r, c

    # Get square, symmetric matrix
    size = 4
    m = np.random.rand(size, size)
    m = (m + m.T) / 2
    np.fill_diagonal(m, 1.)
    a = pd.DataFrame(m)

    # Put NaN's in two symmetrical, but random positions
    row, col = _get_off_diag_pair(size)
    a.iloc[row, col] = np.nan
    a.iloc[col, row] = np.nan

    n_pairs = _get_pairs(a)
    assert n_pairs == int((size**2 - size - 2) / 2)
