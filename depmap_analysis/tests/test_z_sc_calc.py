import numpy as np
import pandas as pd
from depmap_analysis.util.statistics import *
from depmap_analysis.preprocessing.depmap_preprocessing import run_corr_merge
from . import *


def test_merge_corr():
    size = 50
    a, nan_count_a = _get_df_w_nan(size=size)
    b, nan_count_b = _get_df_w_nan(size=size)

    # Get samples
    an = get_n(recalculate=True, data_df=a)
    bn = get_n(recalculate=True, data_df=b)

    # Get logp
    alog = get_logp(recalculate=True, data_corr=a, data_n=an)
    blog = get_logp(recalculate=True, data_corr=b, data_n=bn)

    # Get z
    az = get_z(recalculate=True, data_logp=alog, data_corr=a)
    bz = get_z(recalculate=True, data_logp=blog, data_corr=b)

    stouffer_merged: pd.DataFrame = (az + bz) / np.sqrt(2)
    merged: pd.DataFrame = run_corr_merge(crispr_corr=a, rnai_corr=b,
                                          output_dir='temp')

    # Check inf count
    assert (stouffer_merged.abs() == np.inf).sum().sum() == \
        (merged.abs() == np.inf).sum().sum()

    # Check NaN count
    assert pd.isna(stouffer_merged).sum().sum() == \
        pd.isna(merged).sum().sum()

    # Are they the same?
    assert stouffer_merged.equals(merged)
