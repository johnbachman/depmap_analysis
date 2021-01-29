from datetime import datetime
import pandas as pd
from depmap_analysis.network_functions.net_functions import \
    sif_dump_df_to_digraph

# Add input
h1, h2 = 1234657890, 9876543210
hashes = [h1, h2]
bd = [0.685, 0.95]
src = [{'srcA': 2, 'srcB': 5}, {'srcA': 5, 'srcB': 8}]
sif = {'agA_name': ['nameX1', 'nameX2'], 'agA_ns': ['nsX1', 'nsX2'],
       'agA_id': ['idX1', 'idX2'], 'agB_name': ['nameY1', 'nameY2'],
       'agB_ns': ['nsY1', 'nsY2'], 'agB_id': ['idY1', 'idY2'],
       'stmt_type': ['Activation', 'Complex'], 'evidence_count': [7, 13],
       'stmt_hash': hashes, 'source_counts': src, 'belief': bd}


def test_digraph_dump():
    sif_df = pd.DataFrame(sif)
    idg = sif_dump_df_to_digraph(df=sif_df,
                                 date=datetime.utcnow().strftime('%Y-%m-%d'),
                                 graph_type='digraph',
                                 include_entity_hierarchies=False)
    assert idg.graph.get('edge_by_hash')
    assert idg.graph['edge_by_hash'][h1] == ('nameX1', 'nameY1')
    assert idg.graph['edge_by_hash'][h2] == ('nameX2', 'nameY2')
    assert idg.edges.get(('nameX1', 'nameY1'))
    assert isinstance(idg.edges[('nameX1', 'nameY1')]['statements'], list)
    assert idg.edges[('nameX1', 'nameY1')]['statements'][0]['stmt_hash'] == \
           h1


def test_digraph_signed_types_dump():
    pass


def test_signed_graph_dump():
    pass


def test_expanded_signed_graph_dump():
    pass
