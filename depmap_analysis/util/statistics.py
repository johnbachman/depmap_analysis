import logging
from time import time

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import erfcinv, ndtri_exp

logger = logging.getLogger(__name__)


__all__ = ['get_z', 'get_logp', 'get_n']


def logerfcinv(logp):
    """Calculate inverse of complementary error function given log of argument

    """
    return np.where(logp > -10, erfcinv(np.exp(logp)),
                    ndtri_exp(logp))


def norminv_logcdf(logp):
    """Inverse CDF corresponding to log of p value"""
    return -np.sqrt(2) * logerfcinv(np.log(2) + logp)


def get_logp(recalculate, data_n, data_corr, filepath, method='beta'):
    if method not in ('t', 'beta'):
        raise ValueError('Method must be t or beta')
    start = time()
    if recalculate:
        # T-statistic method
        # See https://stackoverflow.com/a/24469099
        # See https://support.minitab.com/en-us/minitab-express/1/help-and-how-to/basic-statistics/inference/supporting-topics/basics/manually-calculate-a-p-value/
        if method == 't':
            t = data_corr * np.sqrt((data_n - 2)/(1 - data_corr * data_corr))
            logp = np.log(2) + stats.t.logsf(t.abs(), data_n-2)
        # Beta-distribution method
        # https://github.com/scipy/scipy/blob/v1.6.2/scipy/stats/stats.py#L3781-L3962
        elif method == 'beta':
            ab = data_n/2 - 1
            logp = np.log(2) + stats.beta.logcdf(-abs(data_corr), ab, ab, loc=-1, scale=2)
        # Make dataframe
        data_logp = pd.DataFrame(logp, columns=data_corr.columns, index=data_corr.index)
        data_logp.to_hdf('%s.h5' % filepath, filepath.split('/')[-1])
    else:
        data_logp = pd.read_hdf('%s.h5' % filepath)
    elapsed = time() - start
    print(elapsed, "sec")
    return data_logp


def get_z(recalculate, data_logp, data_corr, filepath):
    start = time()
    if recalculate:
        #z_mat = stats.norm.ppf(1 - np.exp(data_logp) / 2)
        #z_mat = -norminv_logcdf(data_logp - np.log(2))
        z_mat = abs(norminv_logcdf(data_logp - np.log(2)))
        data_sign = data_corr.copy()
        data_sign[data_sign < 0] = -1
        data_sign[data_sign > 0] = 1
        data_z = data_sign * pd.DataFrame(z_mat, index=data_logp.columns,
                                          columns=data_logp.columns)
        data_z.to_hdf('%s.h5' % filepath, filepath.split('/')[-1])
    else:
        data_z = pd.read_hdf('%s.h5' % filepath)
    elapsed = time() - start
    print(elapsed, "sec")
    return data_z


def get_n(recalculate, data_df, data_corr, filepath):
    start = time()
    if recalculate:
        num_cols = data_df.shape[1]
        data_mat = data_df._get_numeric_data().to_numpy(dtype=float,
                                                   na_value=np.nan, copy=False)
        n_mat = np.zeros((num_cols, num_cols))
        group_start = time()
        for a_ix in range(num_cols):
            if a_ix % 100 == 0:
                print(a_ix, '%.2f' % (time() - group_start), 'sec per round,',
                      int(time() - start), 'sec total')
                group_start = time()
            n_mat[a_ix, a_ix] = (~np.isnan(data_mat[:, a_ix])).sum()
            for b_ix in range(a_ix + 1, num_cols):
                n = (~np.isnan(data_mat[:, a_ix]) &
                     ~np.isnan(data_mat[:, b_ix])).sum()
                n_mat[a_ix, b_ix] = n
                n_mat[b_ix, a_ix] = n
        data_n = pd.DataFrame(n_mat, index=data_df.columns, columns=data_df.columns)
        data_n.to_hdf('%s.h5' % filepath, filepath.split('/')[-1])
    else:
        data_n = pd.read_hdf('%s.h5' % filepath)
    elapsed = time() - start
    print(int(elapsed), 'sec')
    return data_n