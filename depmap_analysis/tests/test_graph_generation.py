import pandas as pd
from depmap_analysis.network_functions.net_functions import \
    sif_dump_df_to_digraph

# Add input
bd = {}
src = {}
sif = {}


def test_sif_df():
    sif_df = pd.DataFrame(sif)
    idg = sif_dump_df_to_digraph(df=sif_df, strat_ev_dict=src,
                                 belief_dict=bd, graph_type='digraph')
    assert idg.graph.get('edge_by_hash')
