import logging
from math import floor
from io import BytesIO
from pathlib import Path
from datetime import datetime

import boto3
import pandas as pd
import matplotlib.pyplot as plt

from indra.util.aws import get_s3_client
from depmap_analysis.scripts.corr_stats_axb import main as axb_stats

logger = logging.getLogger(__name__)


class DepMapExplainer:
    """Contains the result of the matching of correlations and an indranet
    graph

    Attributes
    ----------
    tag : str
    indra_network_date : str
    depmap_date : str
    sd_range : tuple(float|None)
    info : dict
    network_type : str
    stats_df : pd.DataFrame
    expl_df : pd.DataFrame
    is_signed : Bool
    summary : dict
    summary_str : str
    corr_stats_axb : dict
    """

    def __init__(self, stats_columns, expl_columns, info, script_settings,
                 tag=None, network_type='digraph'):
        """
        Parameters
        ----------
        stats_columns : list[str]|tuple[str]
        expl_columns : list[str]|tuple[str]
        info : dict
        script_settings : dict
        tag : str
        network_type : str
        """
        self.tag = tag
        self.indra_network_date = info.pop('indra_network_date')
        self.depmap_date = info.pop('depmap_date')
        self.sd_range = info.pop('sd_range')
        self.info = info
        self.script_settings = script_settings
        self.network_type = network_type
        self.stats_df = pd.DataFrame(columns=stats_columns)
        self.expl_df = pd.DataFrame(columns=expl_columns)
        self._has_data = False
        self.is_signed = True if network_type in {'signed', 'pybel'} else False
        self.summary = {}
        self.summary_str = ''
        self.corr_stats_axb = {}

    def __str__(self):
        return self.get_summary_str() if self.has_data else \
            'DepMapExplainer is empty'

    def __len__(self):
        # Will return the number of pairs checked
        return len(self.stats_df)

    def has_data(self):
        """Check if any of the data frames have data in them

        Returns
        -------
        bool
        """
        if len(self.stats_df) > 0 or len(self.expl_df) > 0:
            self._has_data = True
        else:
            self._has_data = False
        return self._has_data

    def summarize(self):
        """Count explanations and print a summary count of them"""
        if not self.summary_str:
            self.summary_str = self.get_summary_str()
        print(self.summary_str)

    def get_summary(self):
        """Calculate and return a dict with the summary counts per explanation

        Returns
        -------
        dict
        """
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
        """Get the summary string or fill it out from the summary dictionary

        Returns
        -------
        str
        """
        if not self.summary_str:
            for expl in ['total checked', 'not in graph', 'explained',
                         'explained (excl sr)', 'unexplained',
                         'explained set', 'common parent',
                         'complex or direct', 'x intermediate',
                         'shared regulator', 'shared target', 'sr only']:
                summary = self.get_summary()
                self.summary_str +=\
                    (expl + ": ").ljust(22) + str(summary[expl]) + '\n'
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

    def get_corr_stats_axb(self, z_corr=None, max_proc=None, reactome=None,
                           max_so_pairs_size=10000, mp_pairs=True):
        """Get statistics of the correlations associated with different
        explanation types

        Parameters
        ----------
        z_corr : pd.DataFrame
            A pd.DataFrame containing the correlation z scores used to
            create the statistics in this object
        max_proc : int > 0
            The maximum number of processes to run in the multiprocessing
            in get_corr_stats_mp. Default: multiprocessing.cpu_count()
        reactome : tuple[dict]|list[dict]
            A tuple or list of dicts. The first dict is expected to contain
            mappings from UP IDs of genes to Reactome pathway IDs. The second
            dict is expected to contain the reverse mapping (i.e Reactome IDs
            to UP IDs). The third dict is expected to contain mappings from
            the Reactome IDs to their descriptions.
        max_so_pairs_size : int
            The maximum number of correlation pairs to process. If the
            number of eligble pairs is larger than this number, a random
            sample of max_so_pairs_size is used. Default: 10 000. If the
            number of pairs to check is smaller than 10 000, no sampling is
            done.
        mp_pairs : bool
            If True, get the pairs to process using multiprocessing if larger
            than 10 000. Default: True.

        Returns
        -------
        dict
            A Dict containing correlation data for different explanations
        """
        if not self.corr_stats_axb:
            if z_corr is None:
                raise ValueError('The z score correlation matrix must be '
                                 'provided when running get_corr_stats_axb '
                                 'for the first time.')
            if isinstance(z_corr, str):
                z_corr = pd.read_hdf(z_corr)
            self.corr_stats_axb = axb_stats(
                self.expl_df, z_corr=z_corr, reactome=reactome,
                eval_str=False, max_proc=max_proc,
                max_corr_pairs=max_so_pairs_size, do_mp_pairs=mp_pairs
            )
        return self.corr_stats_axb

    def plot_corr_stats(self, outdir, z_corr=None, reactome=None,
                        show_plot=False, max_proc=None, index_counter=None,
                        max_so_pairs_size=10000, mp_pairs=True):
        """Plot the results of running explainer.get_corr_stats_axb()

        Parameters
        ----------
        outdir : str
            The output directory to save the plots in. If string starts with
            's3:' upload to s3. outdir must then have the form
            's3:<bucket>/<sub_dir>' where <bucket> must be specified and
            <sub_dir> is optional and may contain subdirectories.
        z_corr : pd.DataFrame
            A pd.DataFrame containing the correlation z scores used to
            create the statistics in this object
        reactome : tuple[dict]|list[dict]
            A tuple or list of dicts. The first dict is expected to contain
            mappings from UP IDs of genes to Reactome pathway IDs. The second
            dict is expected to contain the reverse mapping (i.e Reactome IDs
            to UP IDs). The third dict is expected to contain mappings from
            the Reactome IDs to their descriptions.
        show_plot : bool
            If True also show plots
        max_proc : int > 0
            The maximum number of processes to run in the multiprocessing in
            get_corr_stats_mp. Default: multiprocessing.cpu_count()
        index_counter : generator
            An object which produces a new int by using 'next()' on it. The
            integers are used to separate the figures so as to not append
            new plots in the same figure.
        max_so_pairs_size : int
            The maximum number of correlation pairs to process. If the
            number of eligible pairs is larger than this number, a random
            sample of max_so_pairs_size is used. Default: 10000.
        mp_pairs : bool
            If True, get the pairs to process using multiprocessing if larger
            than 10 000. Default: True.
        """
        # Local file or s3
        if outdir.startswith('s3:'):
            full_path = outdir.replace('s3:', '').split('/')
            bucket = full_path[0]
            if not _bucket_exists(bucket):
                raise FileNotFoundError(f'The bucket {bucket} seems to not '
                                        f'exist on s3.')
            key_base = '/'.join(full_path[1:]) if len(full_path) > 1 else \
                'output_data_' + datetime.utcnow().strftime('%Y%m%d%H%M%S')
            od = None
        else:
            bucket = None
            key_base = None
            od = Path(outdir)
            if not od.is_dir():
                od.mkdir(parents=True, exist_ok=True)

        # Get corr stats
        corr_stats = self.get_corr_stats_axb(
            z_corr=z_corr, max_proc=max_proc, reactome=reactome,
            max_so_pairs_size=max_so_pairs_size, mp_pairs=mp_pairs
        )
        sd = f'{self.sd_range[0]} - {self.sd_range[1]} SD' \
            if self.sd_range[1] else f'{self.sd_range[0]}+ SD'
        for n, (k, v) in enumerate(corr_stats.items()):
            for m, plot_type in enumerate(['all_azb_corrs', 'azb_avg_corrs',
                                           'all_x_corrs', 'avg_x_corrs',
                                           'top_x_corrs',
                                           'all_reactome_corrs',
                                           'reactome_avg_corrs']):
                if len(v[plot_type]) > 0:
                    name = '%s_%s_%s.pdf' % \
                           (plot_type, k, self.script_settings['graph_type'])
                    if od is None:
                        fname = BytesIO()
                    else:
                        fname = od.joinpath(name).as_posix()
                    if isinstance(v[plot_type][0], tuple):
                        data = [t[-1] for t in v[plot_type]]
                    else:
                        data = v[plot_type]
                    fig_index = next(index_counter) if index_counter \
                        else int(f'{n}{m}')
                    plt.figure(fig_index)
                    plt.hist(x=data, bins='auto')
                    title = '%s %s; %s (%s)' % \
                            (plot_type.replace('_', ' ').capitalize(),
                             k.replace('_', ' '), sd,
                             self.script_settings['graph_type'])
                    plt.title(title)
                    plt.xlabel('combined z-score')
                    plt.ylabel('count')

                    # Save to file or ByteIO and S3
                    plt.savefig(fname, format='pdf')
                    if od is None:
                        # Reset pointer
                        fname.seek(0)
                        # Upload to s3
                        _upload_bytes_io_to_s3(bytes_io_obj=fname,
                                               bucket=bucket,
                                               key=key_base + '/' + name)

                    # Show plot
                    if show_plot:
                        plt.show()

                    # Close figure
                    plt.close(fig_index)
                else:
                    logger.warning(f'Empty result for {k} ({plot_type}) in '
                                   f'range {sd} for graph type '
                                   f'{self.script_settings["graph_type"]}')

    def plot_dists(self, outdir, z_corr=None, reactome=None,
                   show_plot=False, max_proc=None, index_counter=None,
                   max_so_pairs_size=10000, mp_pairs=True):
        """Compare the distributions of differently sampled A-X-B correlations

        Parameters
        ----------
        outdir : str
            The output directory to save the plots in. If string starts with
            's3:' upload to s3. outdir must then have the form
            's3:<bucket>/<sub_dir>' where <bucket> must be specified and
            <sub_dir> is optional and may contain subdirectories.
        z_corr : pd.DataFrame
            A pd.DataFrame containing the correlation z scores used to
            create the statistics in this object
        reactome : tuple[dict]|list[dict]
            A tuple or list of dicts. The first dict is expected to contain
            mappings from UP IDs of genes to Reactome pathway IDs. The second
            dict is expected to contain the reverse mapping (i.e Reactome IDs
            to UP IDs). The third dict is expected to contain mappings from
            the Reactome IDs to their descriptions.
        show_plot : bool
            If True also show plots
        max_proc : int > 0
            The maximum number of processes to run in the multiprocessing in
            get_corr_stats_mp. Default: multiprocessing.cpu_count()
        index_counter : generator
            An object which produces a new int by using 'next()' on it. The
            integers are used to separate the figures so as to not append
            new plots in the same figure.
        max_so_pairs_size : int
            The maximum number of correlation pairs to process. If the
            number of eligible pairs is larger than this number, a random
            sample of max_so_pairs_size is used. Default: 10000.
        mp_pairs : bool
            If True, get the pairs to process using multiprocessing if larger
            than 10 000. Default: True.
        """
        # Local file or s3
        if outdir.startswith('s3:'):
            full_path = outdir.replace('s3:', '').split('/')
            bucket = full_path[0]
            if not _bucket_exists(bucket):
                raise FileNotFoundError(f'The bucket {bucket} seems to not '
                                        f'exist on s3.')
            key_base = '/'.join(full_path[1:]) if len(full_path) > 1 else \
                'output_data_' + datetime.utcnow().strftime('%Y%m%d%H%M%S')
            od = None
        else:
            bucket = None
            key_base = None
            od = Path(outdir)
            if not od.is_dir():
                od.mkdir(parents=True, exist_ok=True)

        # Get corr stats
        corr_stats = self.get_corr_stats_axb(
            z_corr=z_corr, max_proc=max_proc, reactome=reactome,
            max_so_pairs_size=max_so_pairs_size, mp_pairs=mp_pairs
        )
        fig_index = next(index_counter) if index_counter \
            else floor(datetime.timestamp(datetime.utcnow()))
        plt.figure(fig_index)
        all_ind = corr_stats['axb_not_dir']
        plt.hist(all_ind['azb_avg_corrs'], bins='auto', density=True,
                 color='b', alpha=0.3)
        plt.hist(all_ind['avg_x_corrs'], bins='auto', density=True,
                 color='r', alpha=0.3)
        legend = ['A-X-B for all X', 'A-X-B for X in network']
        if len(all_ind['reactome_avg_corrs']):
            plt.hist(all_ind['reactome_avg_corrs'], bins='auto',
                     density=True, color='g', alpha=0.3)
            legend.append('A-X-B for X in reactome path')

        sd_str = f'{self.sd_range[0]} - {self.sd_range[1]} SD' \
            if self.sd_range[1] else f'{self.sd_range[0]}+ SD'
        title = 'avg X corrs %s (%s)' % (sd_str,
                                         self.script_settings['graph_type'])
        plt.title(title)
        plt.ylabel('Norm. Density')
        plt.xlabel('mean(abs(corr(a,x)), abs(corr(x,b))) (SD)')
        plt.legend()
        name = '%s_%s_axb_hist_comparison.pdf' % \
               (sd_str, self.script_settings['graph_type'])

        # Save to file or ByteIO and S3
        if od is None:
            fname = BytesIO()
        else:
            fname = od.joinpath(name).as_posix()
        plt.savefig(fname, format='pdf')
        if od is None:
            # Reset pointer
            fname.seek(0)
            # Upload to s3
            _upload_bytes_io_to_s3(bytes_io_obj=fname, bucket=bucket,
                                   key=key_base + '/' + name)

        # Show plot
        if show_plot:
            plt.show()

        # Close figure
        plt.close(fig_index)


def _upload_bytes_io_to_s3(bytes_io_obj, bucket, key):
    """Upload a BytesIO object to s3

    Parameters
    ----------
    bytes_io_obj : BytesIO
        Object to upload
    bucket : str
        S3 bucket to save object in
    key : str
        S3 key to use for object
    """
    bytes_io_obj.seek(0)  # Just in case
    s3 = get_s3_client(unsigned=False)
    s3.put_object(Body=bytes_io_obj, Bucket=bucket, Key=key)


def _bucket_exists(buck):
    s3 = boto3.resource('s3')
    return s3.Bucket(buck).creation_date is not None
