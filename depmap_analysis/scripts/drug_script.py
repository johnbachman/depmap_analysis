import sys
import numpy as np
from collections import Counter
from depmap_analysis.util.io_functions import file_opener


def get_rankings(expl_df, sampl_size=None):
    """Get the count of next nearest neighborhood

    Parameters
    ----------
    expl_df : pd.DataFrame
        Explanation dataframe
    sampl_size : int
        Sampling size.

    Returns
    -------
    dict
        A dict of Counters for each examined agent
    """
    # Get the agents that have any shared downstream explanations
    all_agents_sampled = set(
        set(expl_df[expl_df['expl type'] == 'shared downstream'].agA.values) &
        set(expl_df[expl_df['expl type'] == 'shared downstream'].agB.values)
    )
    if sampl_size and sampl_size < len(all_agents_sampled):
        # Sample rows
        all_agents_sampled = np.random.choice(list(all_agents_sampled),
                                              sampl_size, replace=False)

    nnn_expl = {}
    expl_df = expl_df[expl_df['expl type'] == 'shared downstream']
    for ag_name in all_agents_sampled:
        ll = []
        for lk in expl_df['expl data'][(expl_df.agA == ag_name) |
                                       (expl_df.agB == ag_name)].values:
            ll += lk
        nnn_expl[ag_name] = Counter(ll)
    return nnn_expl


if __name__ == '__main__':
    drug_file = sys.argv[1]
    try:
        sample_size = sys.argv[2]
    except IndexError:
        sample_size = None
    drug_expl = file_opener(drug_file)
    rankings = get_rankings(drug_expl)
