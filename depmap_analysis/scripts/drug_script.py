"""
Post processing script after running the depmap script on the drug data
"""
import sys
import numpy as np
import pandas as pd
import logging
from typing import Tuple, List, Optional, Set, Dict
from networkx import DiGraph
from collections import Counter, defaultdict
from depmap_analysis.util.io_functions import file_opener
from depmap_analysis.util.statistics import DepMapExplainer

logger = logging.getLogger(__name__)


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

    def _get_sources(g: DiGraph, s: str, x_list: List[str],
                     allowed_srcs: Set[str]) -> List[str]:
        """Return a list of tuples (x, L) where L=[source] is a list of
        sources for the edge s-x. If allowed_srcs is provided, the list will
        be filtered to only allow x if the s-x edge has support from the
        allowed sources"""
        w_srcs = []
        for x in x_list:
            found = False
            stmt_dict_list: List[Dict[str, Dict[str, int]]] = \
                g.edges.get((s, x), {}).get('statements', [])
            for stmt_dict in stmt_dict_list:
                for src in stmt_dict['source_counts']:
                    if src.lower() in allowed_srcs:
                        w_srcs.append(x)
                        found = True
                        break  # Break source dict loop
                if found:
                    break  # Break stmt dict loop

        return w_srcs

    if graph is None:
        logger.info('No graph provided, will skip getting sources of targets')
    jaccard_ranks = []

    # agA, agB, z-score, agA_ns, agA_id, agB_ns, agB_id, not in graph,
    # explained, shared downstream, shared target
    for (agA, agB, corr, _, _, _, _, _, _, _, _) in \
            stats_df[stats_df.explained == True].itertuples(index=False):
        if agA != agB:
            l2_a_succ, l2_b_succ, l2_int, l2_uni = ([],) * 4
            l3_a_succ, l3_b_succ, l3_int, l3_uni = ([],) * 4

            for n, (expl_type, expl_data) in enumerate(
                    expl_df[['expl type', 'expl data']][
                        ((expl_df.agA == agA) & (expl_df.agB == agB) |
                         (expl_df.agA == agB) & (expl_df.agB == agA))
                    ].itertuples(index=False)):
                if n >= 2:
                    # raise IndexError('Should not have more than one '
                    #                  'explanation per (A,B) pair per '
                    #                  'category ("shared downstream" and '
                    #                  '"shared target" should be only '
                    #                  'explanations)')
                    # Todo: handle misnamed entities
                    if n == 2:
                        logger.warning(f'Unexpected data from {agA}, {agB}, '
                                       f'skipping further data')
                    continue

                if expl_data is not None:
                    if expl_type == 'shared target':
                        l2_a_succ, l2_b_succ, l2_int, l2_uni = expl_data
                    elif expl_type == 'shared downstream':
                        l3_a_succ, l3_b_succ, l3_int, l3_uni = expl_data
                    else:
                        continue
            # Save:
            # A, B, corr, st JI, sd JI,
            # l2_a, l2_b, l2_int, l2_uni,
            # l3_a, l3_b, l3_int, l3_uni
            l2_ji = len(l2_int) / len(l2_uni) if len(l2_uni) else 0
            l3_ji = len(l3_int) / len(l3_uni) if len(l3_uni) else 0

            # Get the l1 and l2 columns
            if graph:
                # Get filtered successors
                # L0
                l0_a_succ = _get_sources(graph, agA, l2_a_succ, {'drugbank'})
                l0_b_succ = _get_sources(graph, agB, l2_b_succ, {'drugbank'})
                l0_int = set() & set()
                l0_uni = set() | set()
                l0_ji = len(l0_int)/len(l0_uni) if len(l0_uni) else 0
                # L1
                l1_a_succ = _get_sources(
                    graph, agA, l2_a_succ, {'tas', 'drugbank'})
                l1_b_succ = _get_sources(
                    graph, agB, l2_b_succ, {'tas', 'drugbank'})
                l1_int = set(l1_a_succ) & set(l1_b_succ)
                l1_uni = set(l1_a_succ) | set(l1_b_succ)
                l1_ji = len(l1_int)/len(l1_uni) if len(l1_uni) else 0
            # If no graph, just set None
            else:
                l0_a_succ, l0_b_succ, l0_int, l0_uni, l0_ji = (None,) * 5
                l1_a_succ, l1_b_succ, l1_int, l1_uni, l1_ji = (None,) * 5

            jaccard_ranks.append(
                (agA, agB, corr,
                 l0_a_succ, l0_b_succ, l0_int, l0_uni, l0_ji,
                 l1_a_succ, l1_b_succ, l1_int, l1_uni, l1_ji,
                 l2_a_succ, l2_b_succ, l2_int, l2_uni, l2_ji,
                 l3_a_succ, l3_b_succ, l3_int, l3_uni, l3_ji)
            )
    output_cols = (
        'drugA', 'drugB', 'correlation',

        'L0_succ_a', 'L0_succ_b', 'L0_intersection',
        'L0_union', 'L0_jaccard_index',

        'L1_succ_a', 'L1_succ_b', 'L1_intersection',
        'L1_union', 'L1_jaccard_index',

        'L2_succ_a', 'L2_succ_b', 'L2_intersection',
        'L2_union', 'L2_jaccard_index',

        'L3_succ_a', 'L3_succ_b', 'L3_intersection',
        'L3_union', 'L3_jaccard_index'
    )
    return pd.DataFrame(data=jaccard_ranks, columns=output_cols)


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
