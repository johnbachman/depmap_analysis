import logging
from math import floor
from io import BytesIO
from typing import Tuple, Dict, Union, Hashable, List, Optional, Any
from pathlib import Path
from datetime import datetime

import boto3
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

from indra.util.aws import get_s3_client
from indra_db.util.s3_path import S3Path
from depmap_analysis.util.io_functions import file_opener
from depmap_analysis.scripts.depmap_script_expl_funcs import *
from depmap_analysis.scripts.corr_stats_axb import main as axb_stats
from depmap_analysis.scripts.corr_stats_data_functions import HistData, \
    Histogram, HistogramDirectory, Results, get_hist_dir
from depmap_analysis.post_processing import *

logger = logging.getLogger(__name__)


__all__ = ['min_columns', 'id_columns', 'expl_columns', 'DepMapExplainer']


min_columns = ('pair', 'agA', 'agB', 'z_score')
id_columns = min_columns + ('agA_ns', 'agA_id', 'agB_ns', 'agB_id')
expl_columns = min_columns + ('expl_type', 'expl_data')


class DepMapExplainer:
    """Contains the result of the matching of correlations and an indranet
    graph

    Attributes
    ----------
    tag : str
        A string that describes the data
    indra_network_date : str
        The date of the sif dump used to create the graph
    depmap_date : str
        The date (usually a quarter e.g. 19Q4) the depmap data was published
        on depmap.org
    sd_range : Tuple[float, Union[float, None]]
        A tuple of the lower and optionally the upper bound of the z-score
        range to use when getting correlations
    info : Dict[str, Any]
        Dictionary with meta data
    network_type : str
        The graph type used, e.g. unsigned, signed, pybel
    stats_df : pd.DataFrame
        The dataframe that per row contains which entity pairs where checked
        and which explanations were applicable to them.
    expl_df : pd.DataFrame
        The dataframe that per row contains one explanation type for an entity
        pair that has been explained and the data supporting the explanation.
    is_signed : bool
        True if the graph used was signed or not
    summary : Dict[str, int]
        A dict mapping different explanations to counts of the explanation
    summary_str : str
        A printable string that summarizes the data in terms of explanation
        count
    corr_stats_axb : Dict[Dict[str, int]]
        A Dict containing correlation data for different explanations

    Methods
    -------
    has_data : bool
        Checks and returns whether there is any data in any of the
        associated data frames
    summarize
        Prints a human readable summary of the results for a quick overview
    get_summary : Dict
        Generates counts of explanations and returns a dict
    get_summary_str : str
        Generates a string of the results from the dict returned from
        get_summary
    save_summary
        Save the summary dict from get_summary to a csv file
    get_corr_stats_axb : Dict
        Get correlation statistics
    plot_corr_stats
        Plot the results of get_corr_stats_axb
    plot_dists
        Compare the distributions of differently sampled A-X-B correlations
    """

    def __init__(self, stats_columns: Tuple[str],
                 expl_columns: Tuple[str],
                 info: Dict[Hashable, Any],
                 script_settings: Dict[str, Union[str, float, int, List[str]]],
                 tag: Optional[str] = None,
                 network_type: str = 'digraph'):
        """
        Parameters
        ----------
        stats_columns : Union[List[str], Tuple[str]]
            Columns for stats_df
        expl_columns : Union[List[str], Tuple[str]]
            Columns for expl_df
        info : Dict[Hashable, Any]
            Dictionary with meta data
        script_settings : Dict[str, Union[str, float, int, List[str]]]
            Dictionary containing the settings and input files used to run
            the script
        tag : str
            A string that describes the data
        network_type : str
            The graph type used, e.g. unsigned, signed, pybel
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
        self.expl_cols = list(set(stats_columns).difference(id_columns))
        self._has_data = False
        self.is_signed = True if network_type in {'signed', 'pybel'} else False
        self.summary = {}
        self.summary_str = ''
        self.corr_stats_axb: Optional[Results] = None
        self.loaded_corr: Optional[pd.DataFrame] = None
        self.loaded_graph: Optional[Union[nx.DiGraph, nx.MultiDiGraph]] = None
        self.loaded_reactome: Optional[Tuple[Dict[str, Any], ...]] = None
        self.histograms: Optional[HistogramDirectory] = None

    def __str__(self):
        return self.get_summary_str() if self.__len__() else \
            'DepMapExplainer is empty'

    def __len__(self):
        # Will return the number of pairs checked
        return len(self.stats_df)

    def load_graph(self) -> Union[nx.DiGraph, nx.MultiDiGraph]:
        """Load and return the graph used in script"""
        if self.loaded_graph is None:
            if self.script_settings.get('indranet'):
                indranet_file = self.script_settings['indranet']
            elif self.script_settings.get('argparse_info', {}).get('indranet'):
                indranet_file = \
                    self.script_settings['argparse_info']['indranet']
            else:
                raise FileNotFoundError('No graph file location seems to be '
                                        'present in script settings.')

            graph = file_opener(indranet_file)
            assert isinstance(graph, (nx.DiGraph, nx.MultiDiGraph))
            self.loaded_graph = graph
        return self.loaded_graph

    def load_z_corr(self, local_file_path: Optional[str] = None) \
            -> pd.DataFrame:
        """Load and return the correlation data frame used in script

        Note: Deprecate arg when pd.read_hdf can take S3 urls
        https://github.com/pandas-dev/pandas/issues/31902

        Parameters
        ----------
        local_file_path : str
            File path to the correlation matrix data frame. Provide it if the
            file path in the script settings attribute does not exist or is
            inaccessible.

        Returns
        -------
        pd.DataFrame
        """
        if self.loaded_corr is None:
            if local_file_path:
                z_corr_file = local_file_path
            else:
                if self.script_settings.get('z_score'):
                    z_corr_file = self.script_settings['z_score']
                elif self.script_settings.get('argparse_info', {}).get('z_score'):
                    z_corr_file = self.script_settings['argparse_info']['z_score']
                else:
                    raise FileNotFoundError('No file location seems to be '
                                            'present in script settings. Please '
                                            'provide a file using '
                                            '`local_file_path`')
            logger.info(f'Loading {z_corr_file}')
            z_corr = pd.read_hdf(z_corr_file)
            logger.info('Finished loading hdf file')
            assert isinstance(z_corr, pd.DataFrame)
            self.loaded_corr = z_corr
        return self.loaded_corr

    def load_reactome(self):
        """Load and return the reactome data used in script"""
        if self.loaded_reactome is None:
            if self.script_settings.get('argparse_info', {}).get(
                    'reactome_dict'):
                reactome_file = \
                    self.script_settings['argparse_info']['reactome_dict']
                reactome = file_opener(reactome_file)
                assert isinstance(reactome, (tuple, list)), \
                    f'{reactome_file} does not seem to contain tuple of ' \
                    f'(upid - pathway mapping, pathway - upid mapping, ' \
                    f'pathway id - pathway description).'
            else:
                raise FileNotFoundError('No reactome file location seems to '
                                        'be present in script settings.')
            self.loaded_reactome = reactome

        return self.loaded_reactome

    def has_data(self):
        """Check if any of the data frames have data in them

        Returns
        -------
        bool
        """
        return len(self.stats_df) > 0 or len(self.expl_df) > 0

    def _filter_stats_to_interesting(self) -> pd.DataFrame:
        """Filter to axb/bxa/shared target, excl direct, reactome, apriori"""
        return filter_to_interesting(self.stats_df)

    def summarize(self):
        """Count explanations and print a summary count of them"""
        if not self.summary_str:
            self.summary_str = self.get_summary_str()
        print(self.summary_str)

    def extend_stats(self):
        """Extend stats_df with the calculated booleans from self.summary"""
        pass

    def get_summary(self):
        """Return a dict with the summary counts

        Returns
        -------
        Dict
        """
        if not self.summary:
            # Get explanation column counts
            for expl_type in self.stats_df.columns:
                if expl_type in id_columns:
                    continue
                self.summary[expl_type] = self.stats_df[expl_type].sum()

            # Special Counts #
            # Total pairs checked
            self.summary['total checked'] = self.__len__()
            # unexplained: can be NaN, True, False
            self.summary['unexplained'] = \
                sum(self.stats_df['explained'] == False)
            # count "complex or direct"
            if ab_colname in self.stats_df.columns and \
                    ba_colname in self.stats_df.columns:
                self.summary['complex or direct'] = \
                    sum(self.stats_df[ab_colname] | self.stats_df[ba_colname])
            # count directed a-x-b: a->x->b or b->x->a
            if axb_colname in self.stats_df.columns and \
                    bxa_colname in self.stats_df.columns:
                self.summary['x intermediate'] = \
                    sum(self.stats_df[axb_colname] |
                        self.stats_df[bxa_colname])
            # count shared regulator as only expl
            if sr_colname in self.stats_df.columns:
                self.summary['sr only'] = self._get_sr_only()
                # explained - (shared regulator as only expl)
                self.summary['explained (excl sr)'] = \
                    self.summary['explained'] - self.summary['sr only']
            # Count axb type explanations that does not have reactome,
            # direct/complex or apriori explanations
            if all([cn in self.stats_df.colummns for cn in
                    [st_colname, axb_colname, bxa_colname, apriori_colname,
                     ab_colname, ba_colname, react_colname]]):
                self.summary['explained no reactome, direct, apriori'] = \
                    self._get_axb_type_no_react()

        return self.summary

    def get_summary_str(self):
        """Get the summary string or fill it out from the summary dictionary

        Returns
        -------
        str
        """
        # ToDo: Fix order of output
        if not self.summary_str:
            summary = self.get_summary()
            self.summary_str = \
                'Explanation'.ljust(22) + 'count\n' + '-'*len('Explanation') +\
                ' '*(22-len('Explanation')) + '-'*len('count') + '\n'
            for expl, count in summary.items():
                self.summary_str +=\
                    (expl + ": ").ljust(22) + str(count) + '\n'
        return self.summary_str

    def save_summary(self, fname):
        """Save summary to a file

        Parameters
        ----------
        fname : str
        """
        summary = self.get_summary()
        with open(fname, 'w') as f:
            f.write('explanation,count\n')
            for e, c in summary.items():
                f.write(f'{e},{c}\n')

    def _get_sr_only(self):
        # Get indices where 'shared regulator' is True
        sr_true = self.stats_df[
                self.stats_df[sr_colname] == True
            ].index.values
        # Exclude overall explained and shared regulator
        other_cols = [col for col in self.expl_cols if col not in
                      {sr_colname, 'explained'}]
        others_false = self.stats_df[
                self.stats_df[other_cols] == False
            ].index.values

        return len(set(sr_true).intersection(others_false))

    def _get_axb_type_no_react(self):
        df = self._filter_stats_to_interesting()
        return len(df)

    def get_filtered_triples_df(self,
                                z_corr_file_path: Optional[str] = None) \
            -> pd.DataFrame:
        """Generate a data frame containing a-x-b with their metadata

        The columns are:
        'pair', 'agA', 'agB', 'z_score', 'agA_ns', 'agA_id', 'agB_ns',
        'agB_id', 'expl_type', 'agX', 'agX_ns', 'agX_id', 'ax_corr',
        'xb_corr', 'ax_belief', 'xb_belief', 'hashes'

        - 'pair' is the unique key identifying a group of explanations per
           A, B, corr
        - 'ax_data'/'bx_data' are a collection of tuples, each one
           containing (statement type, statement hash, belief)

        Parameters
        ----------
        z_corr_file_path : str
            File path to the correlation matrix data frame. Provide this
            argument if the file path either does not exist in the script
            settings attribute or it is inaccessible.

        Returns
        -------
        pd.DataFrame
            A data frame of filtered pairs with merged data from stats_df,
            expl_df and meta data from the explaining edges in the graph.
        """
        # Load indra graph used
        graph: Union[nx.DiGraph, nx.MultiDiGraph] = self.load_graph()

        # Load the correlation matrix used
        z_corr: pd.DataFrame = \
            self.load_z_corr(local_file_path=z_corr_file_path)

        return get_non_reactome_axb_expl_df(graph=graph,
                                            stats_df=self.stats_df,
                                            expl_df=self.expl_df,
                                            z_corr=z_corr)

    def get_sd_str(self):
        """Construct a string """
        if self.sd_range[0] is None:
            return 'RND'

        if self.sd_range[1]:
            return f'{self.sd_range[0]}-{self.sd_range[1]}SD'
        else:
            return f'{self.sd_range[0]}+SD'

    def get_corr_stats_axb(self, z_corr=None, max_proc=None, reactome=None,
                           max_so_pairs_size=10000, mp_pairs=True,
                           run_linear=False) -> Results:
        """Get statistics of the correlations from different explanation types

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
            number of eligible pairs is larger than this number, a random
            sample of max_so_pairs_size is used. Default: 10 000. If the
            number of pairs to check is smaller than 10 000, no sampling is
            done.
        mp_pairs : bool
            If True, get the pairs to process using multiprocessing if larger
            than 10 000. Default: True.
        run_linear : bool
            If True, gather the data without multiprocessing. This option is
            good when debugging or if the environment for some reason does
            not support multiprocessing. Default: False.

        Returns
        -------
        Results
            A BaseModel containing correlation data for different explanations
        """
        # Todo
        #  1. Save histogram counts: Basemodel of s3url + np histogram/counts
        #  instead of running the loop everytime to get full data. Code to
        #  produce the same plot as plt.hist():
        #  hist, bin_edges = np.histogram(a=data, bins='auto', density=True)
        #  bin_width = bin_edges[1] - bin_edges[0]
        #  bar_placements = bin_edges[:-1] + bin_width/2
        #  plt.bar(bar_placements, hist, color='b', alpha=0.3, width=bin_width)
        #  plt.show()
        #  2. log along y axis
        if not self.corr_stats_axb:
            if z_corr is None:
                raise ValueError('The z score correlation matrix must be '
                                 'provided when running get_corr_stats_axb '
                                 'for the first time.')
            if isinstance(z_corr, str):
                z_corr = pd.read_hdf(z_corr)
            self.corr_stats_axb: Results = axb_stats(
                self.expl_df, self.stats_df, z_corr=z_corr, reactome=reactome,
                eval_str=False, max_proc=max_proc,
                max_corr_pairs=max_so_pairs_size, do_mp_pairs=mp_pairs,
                run_linear=run_linear
            )
        return self.corr_stats_axb

    def get_data(self,
                 z_corr: pd.DataFrame,
                 reactome: Tuple[Dict[str, Union[List[str], str]]],
                 max_corr_pairs: int = 10000,
                 eval_str: bool = False,
                 max_proc: Optional[int] = None,
                 do_mp_pairs: Optional[bool] = False,
                 run_linear: bool = False):
        """

        Parameters
        ----------
        z_corr
        reactome
        max_corr_pairs
        eval_str
        max_proc
        do_mp_pairs
        run_linear

        Returns
        -------

        """
        # Todo: Up- and download to/from S3
        if not self.histograms:
            results: Results = self.get_corr_stats_axb(
                z_corr=z_corr, max_proc=max_proc, reactome=reactome,
                max_so_pairs_size=max_corr_pairs, mp_pairs=do_mp_pairs,
                run_linear=run_linear
            )
            hist_dir = get_hist_dir(results)
            self.histograms = hist_dir
        return self.histograms

    def plot_corr_stats(self, outdir, z_corr=None, reactome=None,
                        show_plot=False, max_proc=None, index_counter=None,
                        max_so_pairs_size=10000, mp_pairs=True,
                        run_linear=False):
        """Plot the results of running explainer.get_corr_stats_axb()

        Parameters
        ----------
        outdir : str
            The output directory to save the plots in. If string starts with
            's3://' upload to s3. outdir must then have the form
            's3://<bucket>/<sub_dir>' where <bucket> must be specified and
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
        run_linear : bool
            If True, gather the data without multiprocessing. This option is
            good when debugging or if the environment for some reason does
            not support multiprocessing. Default: False.
        """
        # Local file or s3
        if outdir.startswith('s3://'):
            s3_path = S3Path.from_string(outdir)
            logger.info(f'Outdir path is on S3: {str(s3_path)}')
            od = None
        else:
            s3_path = None
            od = Path(outdir)
            if not od.is_dir():
                logger.info(f'Creating directory/ies for {od}')
                od.mkdir(parents=True, exist_ok=True)

        # Get corr stats
        corr_stats: Results = self.get_corr_stats_axb(
            z_corr=z_corr, max_proc=max_proc, reactome=reactome,
            max_so_pairs_size=max_so_pairs_size, mp_pairs=mp_pairs,
            run_linear=run_linear
        )
        sd_str = self.get_sd_str()
        for m, (plot_type, data) in enumerate(corr_stats.dict().items()):
            if len(data) > 0:
                name = f'{plot_type}_{self.script_settings["graph_type"]}.pdf'
                logger.info(f'Using file name {name}')
                if od is None:
                    fname = BytesIO()
                else:
                    fname = od.joinpath(name).as_posix()
                if isinstance(data[0], tuple):
                    data = [t[-1] for t in data]

                fig_index = next(index_counter) if index_counter else m
                plt.figure(fig_index)
                plt.hist(x=data, bins='auto')
                title = f'{plot_type.replace("_", " ").capitalize()}; '\
                        f'{sd_str} {self.script_settings["graph_type"]}'

                plt.title(title)
                plt.xlabel('combined z-score')
                plt.ylabel('count')

                # Save to file or ByteIO and S3
                plt.savefig(fname, format='pdf')
                if od is None:
                    # Reset pointer
                    fname.seek(0)
                    # Upload to s3
                    full_s3_path = _joinpath(s3_path, name)
                    _upload_bytes_io_to_s3(bytes_io_obj=fname,
                                           s3p=full_s3_path)

                # Show plot
                if show_plot:
                    plt.show()

                # Close figure
                plt.close(fig_index)
            else:
                logger.warning(f'Empty result for {plot_type} in '
                               f'range {sd_str} for graph type '
                               f'{self.script_settings["graph_type"]}')

    def plot_dists(self, outdir, z_corr: pd.DataFrame = None, reactome=None,
                   show_plot=False, max_proc=None, index_counter=None,
                   max_so_pairs_size=10000, mp_pairs=True, run_linear=False):
        """Compare the distributions of differently sampled A-X-B correlations

        Parameters
        ----------
        outdir : str
            The output directory to save the plots in. If string starts with
            's3://' upload to s3. outdir must then have the form
            's3://<bucket>/<sub_dir>' where <bucket> must be specified and
            <sub_dir> is optional and may contain subdirectories.
        z_corr : Optional[pd.DataFrame]
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
        run_linear : bool
            If True, gather the data without multiprocessing. This option is
            good when debugging or if the environment for some reason does
            not support multiprocessing. Default: False.
        """
        # Local file or s3
        if outdir.startswith('s3://'):
            s3_path = S3Path.from_string(outdir)
            od = None
        else:
            s3_path = None
            od = Path(outdir)
            if not od.is_dir():
                od.mkdir(parents=True, exist_ok=True)

        # Get corr stats
        corr_stats: Results = self.get_corr_stats_axb(
            z_corr=z_corr, max_proc=max_proc, reactome=reactome,
            max_so_pairs_size=max_so_pairs_size, mp_pairs=mp_pairs,
            run_linear=run_linear
        )
        fig_index = next(index_counter) if index_counter \
            else floor(datetime.timestamp(datetime.utcnow()))
        plt.figure(fig_index)
        legend = ['A-X-B for all X', 'A-X-B for X in network']
        # Plot A-Z-B
        plt.hist(corr_stats.azb_avg_corrs, bins='auto', density=True,
                 color='b', alpha=0.3)
        # Plot A-X-B
        plt.hist(corr_stats.avg_x_corrs, bins='auto', density=True,
                 color='r', alpha=0.3)
        # Plot reactome expl in
        if len(corr_stats.reactome_avg_corrs):
            plt.hist(corr_stats.reactome_avg_corrs, bins='auto',
                     density=True, color='g', alpha=0.3)
            legend.append('A-X-B for X in reactome path')

        sd_str = self.get_sd_str()
        title = 'avg X corrs %s (%s)' % (sd_str,
                                         self.script_settings['graph_type'])
        plt.title(title)
        plt.ylabel('Norm. Density')
        plt.xlabel('mean(abs(corr(a,x)), abs(corr(x,b))) (SD)')
        plt.legend(legend)
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
            full_s3_path = _joinpath(s3_path, name)
            _upload_bytes_io_to_s3(bytes_io_obj=fname,
                                   s3p=full_s3_path)

        # Show plot
        if show_plot:
            plt.show()

        # Close figure
        plt.close(fig_index)

    def plot_interesting(self, outdir, z_corr: pd.DataFrame = None,
                         reactome=None, show_plot=False, max_proc=None,
                         index_counter=None, max_so_pairs_size=10000,
                         mp_pairs=True, run_linear=False):
        """Plots the same type of plot as plot_dists, but filters A, B

        A, B are filtered to those that fulfill the following:
            - No a-b or b-a explanations
            - Not explained by apriori explanations
            - Without common reactome pathways
            - With a-x-b, b-x-a or shared target explanation

        Parameters
        ----------
        outdir : str
            The output directory to save the plots in. If string starts with
            's3://' upload to s3. outdir must then have the form
            's3://<bucket>/<sub_dir>' where <bucket> must be specified and
            <sub_dir> is optional and may contain subdirectories.
        z_corr : Optional[pd.DataFrame]
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
        run_linear : bool
            If True, gather the data without multiprocessing. This option is
            good when debugging or if the environment for some reason does
            not support multiprocessing. Default: False.
        """
        # Local file or s3
        if outdir.startswith('s3://'):
            s3_path = S3Path.from_string(outdir)
            od = None
        else:
            s3_path = None
            od = Path(outdir)
            if not od.is_dir():
                od.mkdir(parents=True, exist_ok=True)

        # Get corr stats
        corr_stats: Results = self.get_corr_stats_axb(
            z_corr=z_corr, max_proc=max_proc, reactome=reactome,
            max_so_pairs_size=max_so_pairs_size, mp_pairs=mp_pairs,
            run_linear=run_linear
        )
        fig_index = next(index_counter) if index_counter \
            else floor(datetime.timestamp(datetime.utcnow()))
        plt.figure(fig_index)
        plt.hist(corr_stats.azfb_avg_corrs, bins='auto', density=True,
                 color='b', alpha=0.3)
        plt.hist(corr_stats.avg_x_filtered_corrs, bins='auto', density=True,
                 color='r', alpha=0.3)
        legend = ['Filtered A-X-B for any X',
                  'Filtered A-X-B for X in network']

        sd_str = self.get_sd_str()
        title = f'avg X corrs, filtered {sd_str} ' \
                f'({self.script_settings["graph_type"]})'
        plt.title(title)
        plt.ylabel('Norm. Density')
        plt.xlabel('mean(abs(corr(a,x)), abs(corr(x,b))) (SD)')
        plt.legend(legend)
        name = '%s_%s_axb_filtered_hist_comparison.pdf' % \
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
            full_s3_path = _joinpath(s3_path, name)
            _upload_bytes_io_to_s3(bytes_io_obj=fname,
                                   s3p=full_s3_path)

        # Show plot
        if show_plot:
            plt.show()

        # Close figure
        plt.close(fig_index)


def _upload_bytes_io_to_s3(bytes_io_obj: BytesIO, s3p: S3Path):
    """Upload a BytesIO object to s3

    Parameters
    ----------
    bytes_io_obj : BytesIO
        Object to upload
    s3p : S3Path
        An S3Path instance of the full upload url
    """
    logger.info(f'Uploading BytesIO object to s3: {str(s3p)}')
    bytes_io_obj.seek(0)  # Just in case
    s3 = get_s3_client(unsigned=False)
    s3p.put(body=bytes_io_obj, s3=s3)


def _bucket_exists(buck):
    s3 = boto3.resource('s3')
    return s3.Bucket(buck).creation_date is not None


def _exists(fpath: Union[Path, S3Path]) -> bool:
    if isinstance(fpath, S3Path):
        s3 = boto3.client('s3')
        return fpath.exists(s3)
    else:
        return fpath.is_file()


def _joinpath(fpath: Union[S3Path, Path], other: str) -> Union[S3Path, Path]:
    if isinstance(fpath, Path):
        return fpath.joinpath(other).absolute()
    else:
        if fpath.to_string().endswith('/') and not other.startswith('/') or \
                not fpath.to_string().endswith('/') and other.startswith('/'):
            return S3Path.from_string(fpath.to_string() + other)
        elif fpath.to_string().endswith('/') and other.startswith('/'):
            return S3Path.from_string(fpath.to_string() + other[1:])
        elif not fpath.to_string().endswith('/') and not other.startswith('/'):
            return S3Path.from_string(fpath.to_string() + '/' + other)
        else:
            raise ValueError(f'Unable to join {fpath.to_string()} and '
                             f'{other} with "/"')
