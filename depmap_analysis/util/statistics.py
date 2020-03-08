import logging
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from depmap_analysis.scripts.corr_stats_axb import main as axb_stats

logger = logging.getLogger(__name__)


class DepMapExplainer:

    def __init__(self, stats_columns, expl_columns, info, tag=None,
                 network_type='digraph'):
        self.tag = tag
        self.indra_network_date = info['indra_network_date']
        self.depmap_date = info['depmap_date']
        self.sd_range = info['sd_range']
        self.network_type = network_type
        self.stats_df = pd.DataFrame(columns=stats_columns)
        self.expl_df = pd.DataFrame(columns=expl_columns)
        self._has_data = False
        self.is_signed = True if network_type in {'signed', 'pybel'} else False
        self.summary = {}
        self.summary_str = ''
        self.corr_stats_axb = {}

    def __str__(self):
        return self.summary_str if self.has_data() else \
            'DepMapExplainer is empty'

    def has_data(self):
        if len(self.stats_df) > 0 or len(self.expl_df) > 0:
            self._has_data = True
        else:
            self._has_data = False
        return self._has_data

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
            for expl in ['total checked', 'not in graph', 'explained',
                         'explained (excl sr)', 'unexplained',
                         'explained set', 'common parent',
                         'complex or direct', 'x intermediate',
                         'shared regulator', 'shared target', 'sr only']:
                summary = self.get_summary()
                self.summary_str +=\
                    (expl +": ").ljust(22) + str(summary[expl]) + '\n'
        return self.summary_str

    def save_summary(self, fname):
        """Save summary to a file"""
        summary = self.get_summary()
        with open(fname, 'w') as f:
            f.write('explanation,count\n')
            for e, c in summary.items():
                f.write(f'{e},{c}\n')

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

    def get_corr_stats_axb(self, z_corr=None):
        """Get statistics of the correlations associated with different
        explanation types

        Parameters
        ----------
        z_corr : pd.DataFrame
            A pd.DataFrame containing the correlation z scores used to
            create the statistics in this object

        Returns
        -------
        dict
            A Dict containing correlation data for different explanations
        """
        if not self.corr_stats_axb:
            if not z_corr:
                raise ValueError('The z score correlation matrix must be '
                                 'provided when running get_corr_stats_axb '
                                 'for the first time.')
            if isinstance(z_corr, str):
                z_corr = pd.read_hdf(z_corr)
            self.corr_stats_axb = axb_stats(self.expl_df, z_corr)
        return self.corr_stats_axb

    def plot_corr_stats(self, outdir, z_corr=None):
        od = Path(outdir)
        if not od.is_dir():
            od.mkdir(parents=True, exist_ok=True)
        corr_stats = self.get_corr_stats_axb(z_corr=z_corr)
        sd = f'{self.sd_range[0]} - {self.sd_range[1]} SD' \
            if self.sd_range[1] else f'{self.sd_range[0]}+ SD'
        for k, v in corr_stats.items():
            for plot_type in ['all_azb_corrs', 'azb_avg_corrs', 'all_x_corrs',
                              'avg_x_corrs', 'top_x_corrs']:
                if len(v[plot_type]) > 0:
                    if isinstance(v[plot_type][0], tuple):
                        data = [t[-1] for t in v[plot_type]]
                    else:
                        data = v[plot_type]
                    plt.hist(x=data, bins='auto')
                    plt.title('%s %s; %s' %
                              (plot_type.replace('_', ' ').capitalize(),
                               k.replace('_', ' '),
                               sd))
                    plt.xlabel('combined z-score')
                    plt.ylabel('count')
                    plt.savefig(
                        od.joinpath('%s_%s.png' % (plot_type, k)).as_posix(),
                        format='png'
                    )
                    plt.show()
                else:
                    logger.warning('Empty result for %s (%s) in range %s'
                                   % (k, plot_type, sd))

