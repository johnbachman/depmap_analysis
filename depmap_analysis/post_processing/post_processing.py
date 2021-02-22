import logging
import pandas as pd
from typing import Dict, List, Union, Tuple, Optional
from networkx import DiGraph, MultiDiGraph
from collections import defaultdict, Counter

from depmap_analysis.network_functions.net_functions import ag_belief_score
from depmap_analysis.scripts.depmap_script_expl_funcs import *

logger = logging.getLogger(__name__)

__all__ = ['get_non_reactome_axb_expl_df']


class NotInGraph(Exception):
    """Raise for missing nodes or edges in graph"""


def _get_edge_data(G: Union[DiGraph, MultiDiGraph], a: str,
                   b: str, corr_sign: int = None) \
        -> Dict[str, Union[float, List[Tuple[str, str, float]]]]:
    # Get:
    # 'agg_belief': float
    # 'types': Dict[stmt_type, Tuple[(hash, belief)]]
    edge = (a, b) if corr_sign is None else (a, b, corr_sign)

    # Check if edge not in graph
    if edge not in G.edges:
        raise NotInGraph

    # Get belief
    agg_belief: float = ag_belief_score(float(sd['belief']) for sd in
                                        G.edges[edge]['statements'])

    # Get List[Tuple[stmt type, hash, belief]]
    types_list = []
    for sd in G.edges[edge].get('statements', {}):
        types_list.append((sd['stmt_type'], sd['stmt_hash'], sd['belief']))

    # If no edge data was collected, raise NotInGraph
    if not types_list:
        raise NotInGraph

    return {'agg_belief': agg_belief,
            'types': types_list}


def _check_hashes(a: str, x: str, b: str, ab_corr: float,
                  G: Union[DiGraph, MultiDiGraph], expl_type: str) \
        -> Dict[str, Union[float, str, List[Tuple[str, str, float]]]]:
    # Get edge data:
    # 'agX_ns', 'agX_id', 'ax_belief', 'bx_belief', 'ax_data', bx_data
    if expl_type in {axb_colname, bxa_colname}:
        # a->x; x->b  a/b already flipped in the expl_df where the pairs
        # come from
        s, z1, z2, o = a, x, x, b
    elif expl_type == st_colname:  # a->x; b->x
        s, z1, z2, o = a, x, b, x
    else:
        raise ValueError(f'Unrecognized explanation type {expl_type}')

    x_ns, x_id = '', ''
    # If signed
    if G.is_multigraph():
        p, m = 0, 1

        if ab_corr > 0:
            # Grab edges of equal sign for positive correlation
            sign_tups = [(p, p), (m, m)]
        else:
            # Grab edges of opposite sign for negative correlation
            sign_tups = [(p, m), (m, p)]

        # If any edges passed, get their associated stmt types and hashes
        ax_edge_data = []
        bx_edge_data = []
        ax_beliefs = []
        bx_beliefs = []
        for s1, s2 in sign_tups:
            try:
                # 'agg_belief': float
                # 'types': Dict[stmt_type, Tuple[(hash, belief)]]
                ax_dict = _get_edge_data(G=G, a=s, b=z1, corr_sign=s1)
                ax_beliefs.append(ax_dict['agg_belief'])
                ax_edge_data.append(ax_dict['types'])

                bx_dict = _get_edge_data(G=G, a=z2, b=o, corr_sign=s2)
                bx_beliefs.append(ax_dict['agg_belief'])
                bx_edge_data.append(bx_dict['types'])
                x_ns: str = G.nodes[x]['ns']
                x_id: str = G.nodes[x]['id']

            except NotInGraph:
                continue

        # Take the better belief if two
        ax_belief = max(ax_beliefs) if ax_beliefs else None
        bx_belief = max(bx_beliefs) if bx_beliefs else None

    # If unsigned, just grab the edges
    else:
        try:
            ax_dict = _get_edge_data(G=G, a=s, b=z1)
            bx_dict = _get_edge_data(G=G, a=z2, b=o)
            ax_belief = ax_dict['agg_belief']
            bx_belief = bx_dict['agg_belief']
            ax_edge_data = ax_dict['types']
            bx_edge_data = bx_dict['types']
            x_ns: str = G.nodes[x]['ns']
            x_id: str = G.nodes[x]['id']
        except NotInGraph:
            ax_edge_data = {}
            bx_edge_data = {}
            ax_belief = None
            bx_belief = None

    if ax_edge_data and bx_edge_data:
        return {'agX_ns': x_ns, 'agX_id': x_id,
                'ax_belief': ax_belief, 'bx_belief': bx_belief,
                'ax_data': ax_edge_data, 'bx_data': bx_edge_data}
    return {}


def _get_df_per_key(key: str, stats_df: pd.DataFrame, expl_df: pd.DataFrame,
                    corr_zsc_df: pd.DataFrame, graph: DiGraph,
                    columns: Tuple[str, ...]) \
        -> Tuple[Dict[str, List], Counter]:
    # Ignored expl types
    ign_types = set(funcname_to_colname.values())\
        .difference({st_colname, axb_colname, bxa_colname})

    # Initialize rows
    rows_dict = defaultdict(list)

    # Initialize skip counter
    skips = Counter({'ValueError': 0, 'KeyError': 0})

    # Get stats row for current key
    stats_row: pd.Series = stats_df[stats_df.pair == key]

    # Loop expl rows for current key
    for ix, expl_row in expl_df[expl_df.pair == key].iterrows():
        if expl_row.expl_type in ign_types:
            continue

        # 'expl_type', 'agX', 'agX_ns', 'agX_id', 'ax_corr',
        # 'bx_corr', 'ax_belief', 'bx_belief', 'ax_data', 'bx_data'
        x_iter = expl_row.expl_data[2] if expl_row.expl_type == st_colname \
            else expl_row.expl_data

        for x in x_iter:
            # Get edge data:
            # 'agX_ns', 'agX_id', 'ax_belief',
            # 'bx_belief', 'ax_data', 'bx_data'
            try:
                edge_dict = _check_hashes(a=expl_row.agA, x=x, b=expl_row.agB,
                                          ab_corr=expl_row.z_score, G=graph,
                                          expl_type=expl_row.expl_type)
            except ValueError:
                skips['ValueError'] += 1
                continue

            # Get 'ax_corr', 'bx_corr'
            try:
                ax_corr = corr_zsc_df.loc[expl_row.agA, x]
                bx_corr = corr_zsc_df.loc[expl_row.agB, x]
            except KeyError:
                skips['KeyError'] += 1
                continue

            # Append new row
            rows_dict['pair'].append(stats_row.pair.values[0])
            rows_dict['agA'].append(stats_row.agA.values[0])
            rows_dict['agB'].append(stats_row.agB.values[0])
            rows_dict['z_score'].append(stats_row.z_score.values[0])
            rows_dict['agA_ns'].append(stats_row.agA_ns.values[0])
            rows_dict['agA_id'].append(stats_row.agA_id.values[0])
            rows_dict['agB_ns'].append(stats_row.agB_ns.values[0])
            rows_dict['agB_id'].append(stats_row.agB_id.values[0])
            rows_dict['expl_type'].append(expl_row.expl_type)
            rows_dict['agX'].append(x)
            rows_dict['agX_ns'].append(edge_dict['agX_ns'])
            rows_dict['agX_id'].append(edge_dict['agX_id'])
            rows_dict['ax_corr'].append(ax_corr)
            rows_dict['bx_corr'].append(bx_corr)
            rows_dict['ax_belief'].append(edge_dict['ax_belief'])
            rows_dict['bx_belief'].append(edge_dict['bx_belief'])
            rows_dict['ax_data'].append(edge_dict['ax_data'])
            rows_dict['bx_data'].append(edge_dict['bx_data'])

    if rows_dict:
        # Check that all columns are in the dict
        assert set(columns) == set(rows_dict.keys())
        return rows_dict, skips
    return {}, skips


def get_non_reactome_axb_expl_df(graph: Union[DiGraph, MultiDiGraph],
                                 stats_df: pd.DataFrame,
                                 expl_df: pd.DataFrame,
                                 z_corr: pd.DataFrame) -> pd.DataFrame:
    """Generate an exploded dataframe on the edge data level for filtered expl

    The resulting data frame will have the following columns:
    'pair', 'agA', 'agB', 'z_score', 'agA_ns', 'agA_id', 'agB_ns', 'agB_id',
    'expl_type', 'agX', 'agX_ns', 'agX_id', 'ax_corr', 'bx_corr',
    'ax_belief', 'bx_belief', 'ax_data', 'bx_data'

    'pair' is the unique key identifying a group of explanations per A, B, corr
    'ax_data'/'bx_data' are a collection of tuples, each one containing
        (statement type, statement hash, belief)

    Parameters
    ----------
    graph: Union[DiGraph, MultiDiGraph]
        The graph used to generate the explanations in stats_df and expl_df
    stats_df: pd.DataFrame
        The stats data frame
    expl_df: pd.DataFrame
        The expl data frame
    z_corr: pd.DataFrame
        The square correlation data frame used to generate the expl_df and
        stats_df data frames

    Returns
    -------
    pd.DataFrame
        The exploded data frame
    """
    # Get interesting keys
    df = stats_df[stats_df.not_in_graph == False]
    df = df[((df[st_colname]) |
             (df[axb_colname]) |
             (df[bxa_colname])) &
            (df[apriori_colname] == False) &
            (df[ab_colname] == False) &
            (df[ba_colname] == False) &
            (df[react_colname] == False)]
    ab_keys = df.pair.values

    # Loop AB given from outside, then collect the columns:
    columns = ('pair', 'agA', 'agB', 'z_score', 'agA_ns', 'agA_id',
               'agB_ns', 'agB_id', 'expl_type', 'agX', 'agX_ns', 'agX_id',
               'ax_corr', 'bx_corr', 'ax_belief', 'bx_belief', 'ax_data',
               'bx_data')
    results: Dict[str, List] = {c: [] for c in columns}
    counters = []
    for key in ab_keys:
        # Get df data per x:
        rows_data, skips_counter = \
            _get_df_per_key(key=key, stats_df=stats_df, expl_df=expl_df,
                            corr_zsc_df=z_corr, graph=graph, columns=columns)

        # Append the data to its list
        if rows_data:
            for k, dl in results.items():
                dl.extend(rows_data[k])

    all_skip = sum(counters, Counter())
    if any(all_skip.values()):
        for skip, count in all_skip.items():
            if count > 0:
                logger.warning(f'Skipped {skip} {count} times')

    return pd.DataFrame(data=results)
