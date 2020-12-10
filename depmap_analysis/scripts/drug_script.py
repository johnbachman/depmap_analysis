import sys
import numpy as np
import pandas as pd
import logging
from typing import Tuple, List, Optional
from networkx import DiGraph
from collections import Counter, defaultdict
from depmap_analysis.util.io_functions import file_opener
from depmap_analysis.util.statistics import DepMapExplainer

logger = logging.getLogger(__name__)


def _src_count(succ_list: List[Tuple[str, List[str]]], src: List[str]):
    return sum([any([s in srcl for s in src]) for _, srcl in succ_list])


def _percent_in_tas_or_db_A(row):
    if len(row.succ_a_st):
        count = _src_count(row.succ_a_st, ['tas', 'drugbank'])
        return count / len(row.succ_a_st)
    else:
        return 0


def _percent_in_tas_or_db_B(row):
    if len(row.succ_b_st):
        count = _src_count(row.succ_b_st, ['tas', 'drugbank'])
        return count / len(row.succ_b_st)
    else:
        return 0


def get_jaccard_rankings_per_pair(expl_df: pd.DataFrame,
                                  stats_df: pd.DataFrame,
                                  graph: Optional[DiGraph] = None) \
        -> pd.DataFrame:
    """

    Parameters
    ----------
    expl_df : DataFrame
        Explanation dataframe
    stats_df : DataFrame
        Statistics dataframe
    graph: Optional[DiGraph]
        pass

    Returns
    -------
    DataFrame
    """

    def _get_sources(g: DiGraph, s: str, x_list: List[str]):
        """Return a list of tuples with (x, bool) where bool refers to if
        the edges a-x or b-x has any source from the provided sources"""
        w_srcs = []
        for x in x_list:
            stmt_dict_list = g.edges.get((s, x), {}).get('statements', [])
            x_srcs = set()
            for stmt_dict in stmt_dict_list:
                for src in stmt_dict['source_counts']:
                    x_srcs.add(src)
            w_srcs.append((x, list(x_srcs)))
        return w_srcs

    if graph is None:
        logger.info('No graph provided, will skip getting sources of targets')
    jaccard_ranks = []

    # agA, agB, z-score, agA_ns, agA_id, agB_ns, agB_id, not in graph,
    # explained, shared downstream, shared target
    for (agA, agB, corr, _, _, _, _, _, _, _, _) in \
            stats_df[stats_df.explained == True].itertuples(index=False):
        if agA != agB:
            st_a_succ, st_b_succ, st_int, st_uni = ([],) * 4
            sd_a_succ, sd_b_succ, sd_int, sd_uni = ([],) * 4

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
                    # raise IndexError('Should not have more than one '
                    #                  'explanation per (A,B) pair per '
                    #                  'category ("shared downstream" and '
                    #                  '"shared target" should be only '
                    #                  'explanations)')
                    # Todo: handle misnamed entities
                    logger.warning(f'Skipping unexpected data from {agA}, '
                                   f'{agB}')
                    continue
            # Save:
            # A, B, corr, st JI, sd JI,
            # a_st, b_st, st_int, st_uni,
            # a_sd, b_sd, sd_int, sd_uni
            st_ji = len(st_int) / len(st_uni) if len(st_uni) else 0
            sd_ji = len(sd_int) / len(sd_uni) if len(sd_uni) else 0
            st_a_succ_srcs = _get_sources(g=graph, s=agA, x_list=st_a_succ) \
                if graph is not None else st_a_succ
            st_b_succ_srcs = _get_sources(g=graph, s=agB, x_list=st_b_succ) \
                if graph is not None else st_b_succ
            jaccard_ranks.append(
                (agA, agB, corr, st_ji, sd_ji,
                 st_a_succ_srcs, st_b_succ_srcs, st_int, st_uni,
                 sd_a_succ, sd_b_succ, sd_int, sd_uni)
            )

    # A, B, corr, st JI, sd JI,
    # n_a_st, n_b_st, st_int, st_uni,
    # n_a_sd, n_b_sd, sd_int, sd_uni
    output_cols = (
        'drugA', 'drugB', 'corr', 'st_jaccard_index', 'sd_jaccard_index',
        'succ_a_st', 'succ_b_st', 'st_intersection', 'st_union',
        'succ_a_sd', 'succ_b_sd', 'sd_intersection', 'sd_union'
    )
    df = pd.DataFrame(data=jaccard_ranks, columns=output_cols)

    # Add some columns
    if graph:
        df['agA_percent_in_tas_db'] = df.apply(_percent_in_tas_or_db_A, axis=1)
        df['agB_percent_in_tas_db'] = df.apply(_percent_in_tas_or_db_B, axis=1)
        df['percent_in_tas_db'] = 0.5*(df.agB_percent_in_tas_db +
                                       df.agA_percent_in_tas_db)
    return df


def get_rankings_per_drug(expl_df: pd.DataFrame, sampl_size: int = None) -> \
        Tuple[Counter, pd.DataFrame]:
    """Get the count of next nearest neighborhood

    Parameters
    ----------
    expl_df : DataFrame
        Explanation dataframe
    sampl_size : int
        Sampling size

    Returns
    -------
    dict
        A dict of Counters for each examined agent
    """
    # Get the agents that have any shared downstream explanations
    all_agents_sampled = \
        set(expl_df[expl_df['expl type'] == 'shared downstream'].agA.values)\
        | \
        set(expl_df[expl_df['expl type'] == 'shared downstream'].agB.values)
    if sampl_size and sampl_size < len(all_agents_sampled):
        # Sample rows
        logger.info(f'Downsizing to {sampl_size}')
        all_agents_sampled = np.random.choice(list(all_agents_sampled),
                                              sampl_size, replace=False)

    nnn_counters = {}  # Counters of the intersection of downstream targets
    jaccard_index = defaultdict(list)
    for ag_name in all_agents_sampled:
        ll = []
        for sd_a_succ, sd_b_succ, ins, uni in \
                expl_df['expl data'][
                    ((expl_df.agA == ag_name) | (expl_df.agB == ag_name)) &
                    (expl_df['expl type'] == 'shared downstream')
                ].values:
            ll += ins
            jaccard_index[ag_name].append((len(ins), len(uni),
                                           len(ins) / len(uni)))
        nnn_counters[ag_name] = Counter(ll)

    # Sum up the counters to get a full counter
    global_ranking = Counter()
    for c in nnn_counters.values():
        global_ranking += c

    # Get average Jaccard index per drug
    jaccard_ranking = []
    for name, jvs in jaccard_index.items():
        li, lu, ljr = list(zip(*jvs))
        jaccard_ranking.append((name,
                                sum(ljr) / len(ljr),
                                sum(li) / len(li),
                                sum(lu) / len(lu)))
    jaccard_ranking.sort(key=lambda t: t[1], reverse=True)
    df = pd.DataFrame(
        data=jaccard_ranking,
        columns=['drug', 'jaccard_index', 'n_intersection', 'n_union']
    )

    return global_ranking, df


if __name__ == '__main__':
    drug_file = sys.argv[1]
    try:
        sample_size = sys.argv[2]
    except IndexError:
        sample_size = None
    drug_expl = file_opener(drug_file)
    assert isinstance(drug_expl, DepMapExplainer)
    overall_ranking, jaccard_df_per_drug = \
        get_rankings_per_drug(drug_expl.expl_df)
    jaccard_df_per_pair = get_jaccard_rankings_per_pair(drug_expl.expl_df,
                                                        drug_expl.stats_df)

    logger.info('Done with script, results are in variables '
                '`overall_ranking`, `jaccard_df_per_drug` and '
                '`jaccard_df_per_pair`')
