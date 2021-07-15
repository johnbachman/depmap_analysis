from typing import Tuple

import numpy as np
import pandas as pd
from networkx import DiGraph, MultiDiGraph
from datetime import datetime
from depmap_analysis.network_functions.net_functions import \
    sif_dump_df_to_digraph


__all__ = ['get_df', 'get_dg', '_gen_sym_df', '_get_off_diag_pair',
           '_get_df_w_nan']


def get_df() -> pd.DataFrame:
    # Add input
    agA_names = ['X1', 'X2', 'X1', 'X2', 'Z2', 'Z2']
    agA_ns_list = ['nsX1', 'nsX2', 'nsX1', 'nsX2', 'nsZ2', 'nsZ2']
    agA_ids = ['idX1', 'idX2', 'idX1', 'idX2', 'idZ2', 'idZ2']
    agB_names = ['Y1', 'Y2', 'Z1', 'Z1', 'X1', 'X2']
    agB_ns_list = ['nsY1', 'nsY2', 'nsZ1', 'nsZ1', 'nsX1', 'nsX2']
    agB_ids = ['idY1', 'idY2', 'idZ1', 'idZ1', 'idX1', 'idX2']
    h1, h2, h3, h4, h5, h6 = 1234657890, 9876543210, 1212121212, \
                             5454545454, 7878787878, 9191919191
    hashes = [h1, h2, h3, h4, h5, h6]
    bd = [0.685, 0.95, 0.64, 0.897, 0.486, 0.684]
    stmt_types = ['Activation', 'Complex', 'Activation', 'Phosphorylation',
                  'IncreaseAmount', 'DecreaseAmount']
    ev_counts = [7, 13, 13, 13, 13, 13]
    src = [{'srcA': 2, 'srcB': 5}, {'srcA': 5, 'srcB': 8},
           {'srcA': 5, 'srcB': 8},
           {'srcA': 2, 'srcB': 5}, {'srcA': 5, 'srcB': 8},
           {'srcA': 5, 'srcB': 8}]

    sif_dict = {'agA_name': agA_names, 'agA_ns': agA_ns_list,
                'agA_id': agA_ids, 'agB_name': agB_names,
                'agB_ns': agB_ns_list,
                'agB_id': agB_ids, 'stmt_type': stmt_types,
                'evidence_count': ev_counts, 'stmt_hash': hashes,
                'source_counts': src, 'belief': bd}

    sif_df = pd.DataFrame(sif_dict)
    return sif_df


def get_dg() -> DiGraph:
    sif_df = get_df()
    date = datetime.utcnow().strftime('%Y-%m-%d')
    idg: DiGraph = sif_dump_df_to_digraph(df=sif_df, date=date,
                                          graph_type='digraph',
                                          include_entity_hierarchies=False)
    return idg


def _gen_sym_df(size):
    # Get square, symmetric matrix in dataframe
    m = np.random.rand(size, size)
    m = (m + m.T) / 2
    np.fill_diagonal(m, 1.)
    return pd.DataFrame(m)


def _get_off_diag_pair(max_index: int):
    if max_index == 0:
        raise ValueError('Cannot have max_index == 0')
    r = np.random.randint(0, max_index)
    c = np.random.randint(0, max_index)
    while r == c:
        c = np.random.randint(0, max_index)
    return r, c


def _get_df_w_nan(size: int = 50,
                  nan_count: int = 50) -> Tuple[pd.DataFrame, int]:
    """Return symmetric dataframe with some NaNs"""
    # Cannot remove more than a full off-diagonal triangle of the matrix
    remove = min(nan_count, (size**2 - size) / 2)
    a = _gen_sym_df(size)

    pairs = set()
    n = 0
    for n in range(remove):
        k = 0
        row, col = _get_off_diag_pair(size)
        while (row, col) in pairs:
            row, col = _get_off_diag_pair(size)
            k += 1
            if k > 1000:
                print('Too many while iterations, breaking')
                break
        if k > 1000:
            break
        a.iloc[row, col] = np.nan
        a.iloc[col, row] = np.nan
        pairs.add((row, col))
        pairs.add((col, row))

    return a, n+1
