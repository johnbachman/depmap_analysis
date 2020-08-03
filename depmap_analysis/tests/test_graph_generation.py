import pandas as pd
from depmap_analysis.network_functions.net_functions import \
    sif_dump_df_to_digraph

# Add input
hashes = [1234657890, 9876543210]
bd = {1234657890: 0.685, 9876543210: 0.95}
src = {1234657890: {'srcA': 2, 'srcB': 5}, 9876543210: {'srcA': 5, 'srcB': 8}}
sif = {'agA_name': ['nameX1', 'nameX2'], 'agA_ns': ['nsX1', 'nsX2'],
       'agA_id': ['idX1', 'idX2'], 'agB_name': ['nameY1', 'nameY2'],
       'agB_ns': ['nsY1', 'nsY2'], 'agB_id': ['idY1', 'idY2'],
       'stmt_type': ['Activation', 'Complex'], 'evidence_count': [7, 13],
       'stmt_hash': hashes}


def test_sif_df():
    sif_df = pd.DataFrame(sif)
    idg = sif_dump_df_to_digraph(df=sif_df, strat_ev_dict=src,
                                 belief_dict=bd, graph_type='digraph')
    assert idg.graph.get('edge_by_hash')
    assert idg.graph['edge_by_hash'][1234657890] == ('nameX1', 'nameY1')
    assert idg.graph['edge_by_hash'][9876543210] == ('nameX2', 'nameY2')
    assert idg.edges.get(('nameX1', 'nameY1'))
    assert isinstance(idg.edges[('nameX1', 'nameY1')]['statements'], list)
    assert idg.edges[('nameX1', 'nameY1')]['statements'][0]['stmt_hash'] == \
           1234657890
