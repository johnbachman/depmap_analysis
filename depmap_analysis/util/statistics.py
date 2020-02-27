import pandas as pd


class DepMapExplainer:

    def __init__(self, stats_columns, expl_columns, indra_network_date,
                 tag=None, network_type='digraph'):
        self.tag = tag if tag else None
        self.indra_network_date = indra_network_date
        self.network_type = network_type
        self.stats_df = pd.DataFrame(columns=stats_columns)
        self.expl_df = pd.DataFrame(columns=expl_columns)
        self.has_data = False
        self.is_signed = True if network_type in {'signed', 'pybel'} else False

    # todo add methods for statistics output
