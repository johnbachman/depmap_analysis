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
        self.summary = {}
        self.summary_str = ''

    def __str__(self):
        return self.summary_str

    def summarize(self):
        if not self.summary_str:
            self.summary_str = self.get_summary_str()
        print(self.summary_str)

    def get_summary(self):
        if not self.summary:
            # Total pairs checked
            self.summary['total checked'] = len(self.stats_df)
            # Not in graph
            self.summary['not in graph'] = sum(self.stats_df['not in graph'])
            # unexplained
            self.summary['unexplained'] = \
                sum(self.stats_df['explained'] == False)
            # explained
            self.summary['explained'] = self.stats_df['explained'].sum()
            # count common parent
            self.summary['common parent'] = \
                self.stats_df['common parent'].sum()
            # count "explained set"
            self.summary['explained set'] = \
                self.stats_df['explained set'].sum()
            # count "complex or direct"
            self.summary['complex or direct'] = \
                sum(self.stats_df['a-b'] | self.stats_df['b-a'])
            # count directed a-x-b: a->x->b or b->x->a
            self.summary['x intermediate'] = \
                sum(self.stats_df['a-x-b'] | self.stats_df['b-x-a'])
            # count shared target: a->x<-b
            self.summary['shared regulator'] = \
                self.stats_df['shared regulator'].sum()
            # count shared regulator: a<-x->b
            self.summary['shared target'] = \
                self.stats_df['shared target'].sum()
            # count shared regulator as only expl
            self.summary['sr only'] = self._get_sr_only()
            # explained - (shared regulator as only expl)
            self.summary['explained (excl sr)'] = \
                self.summary['explained'] - self.summary['sr only']

        return self.summary

    def get_summary_str(self):
        if not self.summary_str:
            self.summary_str = ''
        return self.summary_str

    def save_summary(self, fname):
        """Save summary to a file"""
        pass

    def _get_sr_only(self):
        # Select rows that match the following conditions
        indices = self.stats_df[
            (self.stats_df['shared regulator'] == True) &
            (self.stats_df['a-b'] == False) &
            (self.stats_df['b-a'] == False) &
            (self.stats_df['common parent'] == False) &
            (self.stats_df['explained set'] == False) &
            (self.stats_df['a-x-b'] == False) &
            (self.stats_df['b-x-a'] == False) &
            (self.stats_df['shared target'] == False) &
            (self.stats_df['not in graph'] == False)
        ].index
        return len(indices)
