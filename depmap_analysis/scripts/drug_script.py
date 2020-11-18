import sys
import numpy as np
import logging
from typing import Dict
from collections import Counter
from pandas import DataFrame
from depmap_analysis.util.io_functions import file_opener
from depmap_analysis.util.statistics import DepMapExplainer


logger = logging.getLogger(__name__)


def get_rankings(expl_df: DataFrame, sampl_size: int = None)\
        -> Dict[str, Counter]:
    """Get the count of next nearest neighborhood

    Parameters
    ----------
    expl_df : DataFrame
        Explanation dataframe
    sampl_size : int
        Sampling size.

    Returns
    -------
    dict
        A dict of Counters for each examined agent
    """
    # Get the agents that have any shared downstream explanations
    all_agents_sampled = \
        set(expl_df[expl_df['expl type'] == 'shared downstream'].agA.values) |\
        set(expl_df[expl_df['expl type'] == 'shared downstream'].agB.values)
    if sampl_size and sampl_size < len(all_agents_sampled):
        # Sample rows
        logger.info(f'Downsizing to {sampl_size}')
        all_agents_sampled = np.random.choice(list(all_agents_sampled),
                                              sampl_size, replace=False)

    nnn_counters = {}
    expl_df = expl_df[expl_df['expl type'] == 'shared downstream']
    for ag_name in all_agents_sampled:
        ll = []
        for lk in expl_df['expl data'][(expl_df.agA == ag_name) |
                                       (expl_df.agB == ag_name)].values:
            ll += lk
        nnn_counters[ag_name] = Counter(ll)
    return nnn_counters


if __name__ == '__main__':
    drug_file = sys.argv[1]
    try:
        sample_size = sys.argv[2]
    except IndexError:
        sample_size = None
    drug_expl = file_opener(drug_file)
    assert isinstance(drug_expl, DepMapExplainer)
    rankings = get_rankings(drug_expl.expl_df)

    # Sum up the counters to get a full counter
    overall_ranking = Counter()
    for c in rankings.values():
        overall_ranking += c
