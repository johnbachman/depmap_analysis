import logging
from time import time
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import ndtri_exp

logger = logging.getLogger(__name__)


__all__ = ['get_z', 'get_logp', 'get_n']


def get_logp(recalculate: bool, data_n: pd.DataFrame,
             data_corr: pd.DataFrame, filepath: Optional[str] = None,
             method: Optional[str] = 'beta') -> pd.DataFrame:
    """Get the log of p values

    Parameters
    ----------
    recalculate :
        If True, recalculate the log of the p-values
    data_n :
        A dataframe with sampling size values
    data_corr :
        A dataframe with correlation values
    filepath :
        If `recalculate==True`: read the logp values from this file.
        If `recalculate==False`: write the logp values to this file.
        If not provided, run the calculation and return the logp values
        without writing them to a file.
    method :
        Provided the method by which to calculate the log of the p-values.
        Default: 'beta'.

    Returns
    -------
    :
        The logp values calculated using the provided method or read from
        the filepath provided.
    """
    if method not in ('t', 'beta'):
        raise ValueError('Method must be t or beta')
    start = time()
    if recalculate or filepath is None:
        # T-statistic method
        # See https://stackoverflow.com/a/24469099
        # See https://support.minitab.com/en-us/minitab-express/1/help-and-how-to/basic-statistics/inference/supporting-topics/basics/manually-calculate-a-p-value/
        if method == 't':
            logger.info('Getting p values using t statistic method')
            t = data_corr * np.sqrt((data_n - 2)/(1 - data_corr * data_corr))
            logp = np.log(2) + stats.t.logsf(t.abs(), data_n-2)
        # Beta-distribution method
        # https://github.com/scipy/scipy/blob/v1.6.2/scipy/stats/stats.py#L3781-L3962
        else:
            logger.info('Getting p values using beta distribution method')
            ab = data_n/2 - 1
            logp = np.log(2) + \
                stats.beta.logcdf(-abs(data_corr), ab, ab, loc=-1, scale=2)
        # Make dataframe
        data_logp = pd.DataFrame(logp, columns=data_corr.columns,
                                 index=data_corr.index)
        if filepath is not None:
            logger.info(f'Saving logp dataframe to {filepath}')
            data_logp.to_hdf('%s.h5' % filepath, filepath.split('/')[-1])
    else:
        logger.info(f'Reading logp dataframe from file: {filepath}')
        data_logp = pd.read_hdf('%s.h5' % filepath)
    elapsed = time() - start
    print(elapsed, "sec")
    return data_logp


def get_z(recalculate: bool, data_logp: pd.DataFrame, data_corr: pd.DataFrame,
          filepath: Optional[str] = None) -> pd.DataFrame:
    """Get the z-score based on p-values of the correlation matrix

    Parameters
    ----------
    recalculate :
        If True, recalculate the z-scores
    data_logp :
        The logp values
    data_corr :
        The correlation matrix of entity-entity correlations.
    filepath :
        If `recalculate==True`: read the z-score values from this file.
        If `recalculate==False`: write the z-score values to this file.
        If not provided, run the calculation and return the z-score dataframe
        without writing it to a file.

    Returns
    -------
    :
        A dataframe with the z-scores
    """
    start = time()
    if recalculate or filepath is None:
        # z_mat = stats.norm.ppf(1 - np.exp(data_logp) / 2)
        # z_mat = -norminv_logcdf(data_logp - np.log(2))
        z_mat = abs(ndtri_exp(data_logp - np.log(2)))
        data_sign = data_corr.copy()
        data_sign[data_sign < 0] = -1
        data_sign[data_sign > 0] = 1
        data_z = data_sign * pd.DataFrame(z_mat, index=data_logp.columns,
                                          columns=data_logp.columns)
        if filepath is not None:
            logger.info(f'Saving z score dataframe to {filepath}')
            data_z.to_hdf('%s.h5' % filepath, filepath.split('/')[-1])
    else:
        logger.info(f'Reading z-score dataframe from {filepath}')
        data_z = pd.read_hdf('%s.h5' % filepath)
    elapsed = time() - start
    print(elapsed, "sec")
    return data_z


def get_n(recalculate: bool, data_df: pd.DataFrame,
          filepath: Optional[str] = None) -> pd.DataFrame:
    """Get sample sizes

    Parameters
    ----------
    recalculate :
        If True, recalculate the sample sizes
    data_df :
        Correlation data as a dataframe
    filepath :
        If `recalculate==True`: read the correlation values from this file.
        If `recalculate==False`: write the correlation values to this file.
        If not provided, run the calculation and return the correlation data
        without writing it to a file.

    Returns
    -------
    :
        A dataframe holding the sample sizes
    """
    start = time()
    if recalculate or filepath is None:
        logger.info('Calculating sampling values')
        num_cols = data_df.shape[1]
        data_mat = data_df._get_numeric_data().to_numpy(dtype=float,
                                                        na_value=np.nan,
                                                        copy=False)
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
        data_n = pd.DataFrame(n_mat, index=data_df.columns,
                              columns=data_df.columns)
        if filepath is not None:
            data_n.to_hdf('%s.h5' % filepath, filepath.split('/')[-1])
    else:
        logger.info(f'Reading sampling values from file {filepath}')
        data_n = pd.read_hdf('%s.h5' % filepath)
    elapsed = time() - start
    print(int(elapsed), 'sec')
    return data_n
