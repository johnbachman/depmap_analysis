import sys
import numpy as np
import logging
from typing import Dict, Tuple, List
from collections import Counter, defaultdict
from pandas import DataFrame
from depmap_analysis.util.io_functions import file_opener
from depmap_analysis.util.statistics import DepMapExplainer


logger = logging.getLogger(__name__)


def get_rankings(expl_df: DataFrame, stats_df: DataFrame,
                 sampl_size: int = None) -> \
        Tuple[Dict[str, Counter],
              Dict[str, List[Tuple[int, int, float]]],
              List[Tuple[str, str, float, float, float, int,
                         int, int, int, int, int, int, int]
              ]]:
    """Get the count of next nearest neighborhood

    Parameters
    ----------
    expl_df : DataFrame
        Explanation dataframe
    stats_df : DataFrame
        Statistics dataframe
    sampl_size : int
        Sampling size

    Returns
    -------
    dict
        A dict of Counters for each examined agent
    """
    jaccard_ranks = []

    # agA, agB, z-score, agA_ns, agA_id, agB_ns, agB_id, not in graph,
    # explained, shared downstream, shared target
    for (agA, agB, corr, _, _, _, _, _, _, _, _) in \
            stats_df[stats_df.explained == True].itertuples(index=False):
        if agA != agB:
            st_a_succ, st_b_succ, st_int, st_uni = ([],)*4
            sd_a_succ, sd_b_succ, sd_int, sd_uni = ([],)*4

            for n, (expl_type, expl_data) in enumerate(
                    expl_df[['expl type', 'expl data']][
                        ((expl_df.agA == agA) & (expl_df.agB == agB) |
                         (expl_df.agA == agB) & (expl_df.agB == agA))
                    ].itertuples(index=False)):
                if expl_data is not None:
                    if expl_type == 'shared target':
                        st_a_succ, st_b_succ, st_int, st_uni = expl_data
                    elif expl_type == 'shared downstream':
                        sd_a_succ, sd_b_succ, sd_int, sd_uni = expl_data
                    else:
                        continue
                if n > 1:
                    raise IndexError('Should not have more than one '
                                     'explanation per (A,B) pair per '
                                     'category ("shared downstream" and '
                                     '"shared target" should be only '
                                     'explanations)')
            # Save:
            # A, B, corr, st JI, sd JI,
            # n_a_st, n_b_st, st_int, st_uni,
            # n_a_sd, n_b_sd, sd_int, sd_uni
            st_ji = len(st_int)/len(st_uni) if len(st_uni) else 0
            sd_ji = len(sd_int)/len(sd_uni) if len(sd_uni) else 0
            jaccard_ranks.append(
                (agA, agB, corr, st_ji, sd_ji,
                 len(st_a_succ), len(st_b_succ), len(st_int), len(st_uni),
                 len(sd_a_succ), len(sd_b_succ), len(sd_int), len(sd_uni))
            )

    # Get the agents that have any shared downstream explanations
    all_agents_sampled = \
        set(expl_df[expl_df['expl type'] == 'shared downstream'].agA.values) |\
        set(expl_df[expl_df['expl type'] == 'shared downstream'].agB.values)
    if sampl_size and sampl_size < len(all_agents_sampled):
        # Sample rows
        logger.info(f'Downsizing to {sampl_size}')
        all_agents_sampled = np.random.choice(list(all_agents_sampled),
                                              sampl_size, replace=False)

    nnn_counters = {}  # Counters of the intersection of downstream targets
    jaccard_index = defaultdict(list)
    expl_df = expl_df[expl_df['expl type'] == 'shared downstream']
    for ag_name in all_agents_sampled:
        ll = []
        for ins, uni in expl_df['expl data'][(expl_df.agA == ag_name) |
                                             (expl_df.agB == ag_name)].values:
            ll += ins
            jaccard_index[ag_name].append((len(ins), len(uni),
                                           len(ins)/len(uni)))
        nnn_counters[ag_name] = Counter(ll)
    return nnn_counters, jaccard_index, jaccard_ranks


if __name__ == '__main__':
    drug_file = sys.argv[1]
    try:
        sample_size = sys.argv[2]
    except IndexError:
        sample_size = None
    drug_expl = file_opener(drug_file)
    assert isinstance(drug_expl, DepMapExplainer)
    rankings, ji, jr = get_rankings(drug_expl.expl_df, drug_expl.stats_df)

    # Sum up the counters to get a full counter
    overall_ranking = Counter()
    for c in rankings.values():
        overall_ranking += c

    # Get average Jaccard index per drug
    jaccard_ranking = []
    for name, jvs in ji.items():
        li, lu, ljr = list(zip(*jvs))
        jaccard_ranking.append((name,
                                sum(ljr)/len(ljr),
                                sum(li)/len(li),
                                sum(lu)/len(lu)))
    jaccard_ranking.sort(key=lambda t: t[1], reverse=True)
    jaccard_df_per_drug = DataFrame(
        data=jaccard_ranking,
        columns=['drug', 'jaccard_index', 'n_intersection', 'n_union']
    )
    jaccard_df_per_pair = DataFrame(
        data=jr,
        # A, B, corr, st JI, sd JI,
        # n_a_st, n_b_st, st_int, st_uni,
        # n_a_sd, n_b_sd, sd_int, sd_uni
        columns=[
            'drugA', 'drugB', 'corr', 'st_jaccard_index', 'sd_jaccard_index',
            'succ_a_st', 'succ_b_st', 'n_st_intersection', 'n_st_union',
            'succ_a_sd', 'succ_b_sd', 'n_sd_intersection', 'n_sd_union'
        ]
    ).sort_values(by=['jaccard_index', 'corr'], ascending=[False, True])
    logger.info('Done with script, results are in variables '
                '`overall_ranking`, `jaccard_df_per_drug` and '
                '`jaccard_df_per_pair`')
