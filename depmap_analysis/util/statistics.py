import pandas as pd


class DepMapExplainer:
    stats_columns = ('agA', 'agA_ns', 'agA_id', 'agB', 'agB_ns', 'agB_id',
                     'z-score', 'explained', 'direct', 'common parent',
                     'explained set', 'a-x-b', 'b-x-a', 'shared regulator',
                     'shared target', )
    expl_columns = ('agA', 'agB', 'z-score', 'expl_type', 'expl_data')

    def __init__(self, indra_network_date, tag=None, network_type='digraph'):
        self.tag = tag if tag else None
        self.indra_network_date = indra_network_date
        self.network_type = network_type
        self.stats_df = pd.DataFrame(columns=self.stats_columns)
        self.expl_df = pd.DataFrame(columns=self.expl_columns)
