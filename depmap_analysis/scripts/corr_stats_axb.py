"""Get the distribution of correlation z scores for pairs a-x or x-b where x
is an intermediate between gene pair a-b.
"""
import sys
import ast
import pickle
import logging
from os import path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from depmap_analysis.scripts.depmap_script_expl_funcs import axb_colname, \
    bxa_colname, ab_colname, ba_colname, st_colname
from .corr_stats_async import get_corr_stats_mp, GlobalVars, get_pairs_mp

logger = logging.getLogger('DepMap Corr Stats')
logger.setLevel(logging.INFO)


def main(expl_df: pd.DataFrame, stats_df: pd.DataFrame, z_corr: pd.DataFrame,
         reactome: Optional[Tuple[Dict[str, List[str]],
                                  Dict[str, List[str]],
                                  Dict[str, str]]] = None,
         eval_str: Optional[bool] = False,
         max_proc: Optional[int] = None,
         max_corr_pairs: Optional[int] = 10000,
         do_mp_pairs: Optional[bool] = True,
         run_linear: bool = False) \
        -> Dict[str, List[float]]:
    """Get statistics of the correlations associated with different
    explanation types

    Parameters
    ----------
    expl_df: pd.DataFrame
        A pd.DataFrame containing all available explanations for the pairs
        of genes in z_corr. Available in the DepmapExplainer as
        DepmapExplainer.expl_df.
    stats_df: pd.DataFrame
        A pd.DataFrame containing all checked A-B pairs and if they are
        explained or not. Available in the DepmapExplainer as
        DepmapExplainer.stats_df.
    z_corr : pd.DataFrame
        A pd.DataFrame of correlation z scores
    reactome : tuple[dict]|list[dict]
        A tuple or list of dicts. The first dict is expected to contain
        mappings from UP IDs of genes to Reactome pathway IDs. The second
        dict is expected to contain the reverse mapping (i.e Reactome IDs
        to UP IDs). The third dict is expected to contain mappings from the
        Reactome IDs to their descriptions.
    eval_str : bool
        If True, run ast.literal_eval() on the 'expl data' column of expl_df
    max_proc : int > 0
        The maximum number of processes to run in the multiprocessing in
        get_corr_stats_mp. Default: multiprocessing.cpu_count()
    max_corr_pairs : int
        The maximum number of correlation pairs to process. If the number of
        eligible pairs is larger than this number, a random sample of
        max_so_pairs_size is used. Default: 10 000. If the number of pairs
        to check is smaller than 1000, no sampling is done.
    do_mp_pairs : bool
        If True, get the pairs to process using multiprocessing if larger
        than 10 000. Default: True.
    run_linear : bool
        If True, run the script without multiprocessing. This option is good
        when debugging or if the environment for some reason does not
        support multiprocessing. Default: False.

    Returns
    -------
    dict
        A Dict containing correlation data for different explanations
    """
    # Limit to any a-x-b OR a-b expl (this COULD include explanations where
    # 'direct' and NOT 'pathway' is the explanation, but this should be a
    # very small set)
    logger.info("Filter expl_df to pathway, direct, shared_target")
    expl_df = expl_df[
        (expl_df['expl_type'] == axb_colname) |
        (expl_df['expl_type'] == bxa_colname) |
        (expl_df['expl_type'] == ab_colname) |
        (expl_df['expl_type'] == ba_colname) |
        (expl_df['expl_type'] == st_colname)
    ]

    # Re-map the columns containing string representations of objects
    if eval_str:
        expl_df['expl data'] = \
            expl_df['expl data'].apply(lambda x: ast.literal_eval(x))

    # Get all correlation pairs that were explained
    all_ab_corr_pairs = set(map(lambda p: tuple(p),
                                expl_df[['agA', 'agB']].values))

    gbv = GlobalVars(expl_df=expl_df, stats_df=stats_df, sampl=16)
    if not run_linear and do_mp_pairs and len(all_ab_corr_pairs) > 10000:
        # Do multiprocessing
        logger.info('Getting axb subj-obj pairs through multiprocessing')
        gbv.assert_global_vars({'expl_df', 'stats_df'})
        pairs_axb_only = get_pairs_mp(all_ab_corr_pairs, max_proc=max_proc,
                                      max_pairs=max_corr_pairs)
    else:
        logger.info('Assembling axb subj-obj pairs linearly')
        # Pairs where a-x-b AND a-b explanation exists
        pairs_axb_direct = set()

        # Pairs where a-x-b AND NOT a-b explanation exists
        pairs_axb_only = set()

        logger.info("Stratifying correlations by interaction type")
        for s, o in all_ab_corr_pairs:
            # Make sure we don't try to explain self-correlations
            if s == o:
                continue

            # Get all interaction types associated with s and o
            int_types = \
                set(expl_df['expl_type'][(expl_df['agA'] == s) &
                                         (expl_df['agB'] == o)].values)

            # Filter to a-x-b, b-x-a, st
            axb_types = \
                {axb_colname, bxa_colname, st_colname}.intersection(int_types)

            # Only allow pairs where we do NOT have ab or ba explanation
            if axb_types and \
                    ab_colname not in int_types and \
                    ba_colname not in int_types:
                pairs_axb_only.add((s, o))

    # Check for and remove self correlations
    if not np.isnan(z_corr.loc[z_corr.columns[0], z_corr.columns[0]]):
        logger.info('Removing self correlations')
        diag_val = z_corr.loc[z_corr.columns[0], z_corr.columns[0]]
        z_corr = z_corr[z_corr != diag_val]

    # a-x-b AND NOT direct
    logger.info("Getting correlations for a-x-b AND NOT direct")
    options = {'so_pairs': pairs_axb_only, 'run_linear': run_linear}
    if max_proc:
        options['max_proc'] = max_proc

    # Set and assert existence of global variables
    assert_vars = {'z_cm', 'expl_df', 'stats_df'}
    if reactome is not None:
        gbv.update_global_vars(z_cm=z_corr, reactome=reactome)
        assert_vars.add('reactome')
    else:
        logger.info('No reactome file provided')
        gbv.update_global_vars(z_cm=z_corr)
    if gbv.assert_global_vars(assert_vars):
        all_x_corrs_no_direct, avg_x_corrs_no_direct, top_x_corrs_no_direct, \
            all_azb_corrs_no_direct, azb_avg_corrs_no_direct, \
            all_reactome_corrs, reactome_avg_corrs, \
            axb_filtered_avg_corrs, all_axb_filtered_corrs = \
            get_corr_stats_mp(**options)
    else:
        raise ValueError('Global variables could not be set')

    return {'all_x_corrs': all_x_corrs_no_direct,
            'avg_x_corrs': avg_x_corrs_no_direct,
            'top_x_corrs': top_x_corrs_no_direct,
            'all_azb_corrs': all_azb_corrs_no_direct,
            'azb_avg_corrs': azb_avg_corrs_no_direct,
            'all_reactome_corrs': all_reactome_corrs,
            'reactome_avg_corrs': reactome_avg_corrs,
            'all_x_filtered_corrs': all_axb_filtered_corrs,
            'avg_x_filtered_corrs': axb_filtered_avg_corrs}


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <basepath_to_expls> "
              f"<path_to_combined_z_sc_corr_h5>")
        sys.exit(0)
    expl_pairs_csv = sys.argv[1]  # Path to output data folder
    z_corr_file = sys.argv[2]  # Path to merged z scored correlations file
    logger.info('Loading correlation matrix...')
    z_cm = pd.read_hdf(z_corr_file)
    names = [n.split()[0] for n in z_cm.columns.values]
    z_cm.columns = names
    z_cm.index = names

    #sds = ['1_2sd', '2_3sd', '3_4sd', '4_5sd', '5_sd', 'rnd']
    sds = ['3_4sd'] #, '4_5sd']
    results_by_sd = {}
    for sd in sds:
        expl_fname = path.join(expl_pairs_csv, sd, '_explanations_of_pairs.csv')
        expl_dir = path.dirname(expl_fname)
        results_file = path.join(expl_dir, '%s_results_dict.pkl' % sd)
        # Check if we already have the results
        if path.isfile(results_file):
            with open(results_file, 'rb') as f:
                results_by_sd[sd] = pickle.load(f)
        # If we don't already have the results, compute them
        else:
            # Make sure we have the explanations CSV
            if not path.isfile(expl_fname):
                logger.info('Skipping %s, file does not exist' % expl_fname)
            else:
                logger.info('Getting pairs from %s' % expl_fname)
                df = pd.read_csv(expl_fname, delimiter=',')
                results = main(expl_df=df, z_corr=z_cm)
                results_by_sd[sd] = results
                logger.info("Pickling results file %s" % results_file)
                with open(results_file, 'wb') as f:
                    pickle.dump(results, f)

    # Load _explanations_of_pairs.csv for each range
    #for sd in ['1_2sd', '2_3sd', '3_4sd', '4_5sd', '5_sd', 'rnd']:
    for sd in sds:
        try:
            results = results_by_sd[sd]
        except KeyError:
            logger.info("Results for %s not found in dict, skipping" % sd)
            continue
        # Loop the different sets:
        #   - axb_and_dir - subset where direct AND pathway explains
        #   - axb_not_dir - subset where pathway, NOT direct explans
        #   - all_axb - the distribution of the any explanation
        for k, v in results.items():
            # Plot:
            #   1: all_x_corrs - the distribution of all gathered a-x,
            #      x-b combined z-scores
            #   2: top_x_corrs - the strongest (over the a-x, x-b average)
            #      z-score per A-B. List contains (A, B, topx).
            #   3:
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
                    plt.savefig(path.join(expl_dir,
                                          '%s_%s.png' % (plot_type, k)),
                                format='png')
                    plt.show()
                else:
                    logger.warning('Empty result for %s (%s) in range %s'
                                   % (k, plot_type, sd))
