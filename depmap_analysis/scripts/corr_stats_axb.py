"""Get the distribution of correlation z scores for pairs a-x or x-b where x
is an intermediate between gene pair a-b.
"""
import sys
import ast
import pickle
import logging
from os import path, environ
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from depmap_analysis.util.io_functions import pickle_open


logger = logging.getLogger('DepMap Corr Stats')

# DM_INPUT_DIR = environ.get(
#     'DEPMAP_DIR',
#     '/home/klas/repos/depmap_analysis/input_data'
# )
# OUTPUT_PARDIR = environ.get(
#     'OUTDIR',
#     '/home/klas/repos/depmap_analysis/output_data/19Q4_hgnc_fplx'
# )
#
# if not path.isdir(DM_INPUT_DIR):
#     raise ValueError('Must set DEPMAP_DIR as environment variable pointing '
#                      'to the directory where depmap input data is stored. '
#                      'For CRISPR in a subdir of format YYQq (Y=last two '
#                      'digits of year, q=1,2,3,4) for data release. For RNAi '
#                      'in a subdir called demeter')
# if not path.isdir(OUTPUT_PARDIR):
#     raise ValueError('Must set OUTPUT_PARDIR as environment variable '
#                      'pointing to the parent directory where the output '
#                      'data for all the different SD ranges are in '
#                      'subdirectories.')


def get_corr_stats(df, z_sc_cm, so_pairs):
    # Check for and remove self correlations
    if not np.isnan(z_sc_cm.loc[z_sc_cm.columns[0], z_sc_cm.columns[0]]):
        logger.info('Removing self correlations')
        diag_val = z_sc_cm.loc[z_sc_cm.columns[0], z_sc_cm.columns[0]]
        z_sc_cm = z_sc_cm[z_sc_cm != diag_val]

    all_axb_corrs = []
    top_axb_corrs = []
    #ab_to_axb_avg = []
    azb_avg_corrs = []
    all_azb_corrs = []
    axb_avg_corrs = []
    for subj, obj in so_pairs:
        ab_avg_corrs = []
        path_rows = df[(df['agA'] == subj) &
                       (df['agB'] == obj) &
                       ((df['expl type'] == 'a-x-b') |
                        (df['expl type'] == 'b-x-a') |
                       (df['expl type'] == 'shared target'))]

        # Make sure we don't double-count Xs such that X is shared target
        # and also a pathway; also, don't include X if X = subj or obj
        x_set = set()
        for ix, path_row in path_rows.iterrows():
            x_set.update([x for x in path_row['expl data']
                          if x not in (subj, obj)])

        for x in x_set:
            if x in z_sc_cm.columns:
                ax_corr = z_sc_cm.loc[subj, x]
                xb_corr = z_sc_cm.loc[x, obj]
                all_axb_corrs += [ax_corr, xb_corr]
                avg_corr = 0.5 * abs(ax_corr) + 0.5 * abs(xb_corr)
                ab_avg_corrs.append(avg_corr)
                axb_avg_corrs.append(avg_corr)
            else:
                # if warn < 3:
                #     warn += 1
                #     logger.warning('%s does not exist in the crispr and/or '
                #                    'rnai correlation matrices.' % x)
                # else:
                #     logger.warning('Muting warnings...')
                continue

        # Get the set of possible correlation z scores
        for z in set(z_sc_cm.columns):
            if z == subj or z == obj:
                continue
            az_corr = z_sc_cm.loc[z, subj]
            bz_corr = z_sc_cm.loc[z, obj]
            if np.isnan(az_corr) or np.isnan(bz_corr):
                continue  # Is there a more efficient way of doing this?
            all_azb_corrs.extend([az_corr, bz_corr])
            azb_avg_corrs.append(0.5 * abs(az_corr) + 0.5 * abs(bz_corr))

        # if warn:
        #     logger.warning('%d missing X genes out of %d in correlation '
        #                    'matrices' % (warn, len(x_list)))
        if len(ab_avg_corrs) > 0:
            max_magn_avg = max(ab_avg_corrs)
            top_axb_corrs.append((subj, obj, max_magn_avg))
            #ab_to_axb_avg.append((subj, obj, comb_zsc, max_magn_avg))

    return (all_axb_corrs, axb_avg_corrs, top_axb_corrs, all_azb_corrs,
            azb_avg_corrs)


def main(expl_df, z_corr, eval_str=False):
    """Get statistics of the correlations associated with different
    explanation types

    Parameters
    ----------
    expl_df: pd.DataFrame
        A pd.DataFrame containing all available explanations for the pairs
        of genes in z_corr
    z_corr : pd.DataFrame
        A pd.DataFrame of correlation z scores
    eval_str : bool
        If True, run ast.literal_eval() on the 'expl data' column of expl_df

    Returns
    -------
    dict
        A Dict containing correlation data for different explanations
    """
    # Limit to any a-x-b OR a-b expl (this COULD include explanations where
    # 'direct' and NOT 'pathway' is the explanation, but this should be a
    # very small set)
    logger.info("Filter expl_df to pathway, direct, shared_target")
    expl_df = expl_df[(expl_df['expl type'] == 'a-x-b') |
                      (expl_df['expl type'] == 'b-x-a') |
                      (expl_df['expl type'] == 'a-b') |
                      (expl_df['expl type'] == 'b-a') |
                      (expl_df['expl type'] == 'shared target')]

    # Re-map the columns containing string representations of objects
    if eval_str:
        expl_df['expl data'] = \
            expl_df['expl data'].apply(lambda x: ast.literal_eval(x))

    # Get all correlation pairs
    all_ab_corr_pairs = set(map(lambda p: tuple(p),
                                expl_df[['agA', 'agB']].values))
    # Pairs where a-x-b AND a-b explanation exists
    pairs_axb_direct = set()

    # Pairs where a-x-b AND NOT a-b explanation exists
    pairs_axb_only = set()

    # all a-x-b "pathway" explanations, should be union of the above two
    pairs_any_axb = set()

    logger.info("Stratifying correlations by interaction type")
    for s, o in all_ab_corr_pairs:
        # Make sure we don't try to explain self-correlations
        if s == o:
            continue
        # Get all interaction types associated with given subject s and
        # object o
        int_types = set(expl_df['expl type'][(expl_df['agA'] == s) &
                                             (expl_df['agB'] == o)].values)
        # Check intersection of types
        axb_types = {'a-x-b', 'b-x-a', 'shared target'}.intersection(int_types)
        if axb_types:
            pairs_any_axb.add((s, o))
            if 'direct' in int_types:
                # Direct and pathway
                pairs_axb_direct.add((s, o))
            else:
                # Pathway and NOT direct
                pairs_axb_only.add((s, o))
        else:
            # Skip direct only
            continue

    # The union should be all pairs where a-x-b explanations exist
    ab_axb_union = pairs_axb_direct.union(pairs_axb_only)
    assert ab_axb_union == pairs_any_axb

    # a-x-b AND direct
    """
    logger.info("Getting correlations for a-x-b AND direct")
    all_x_corrs_direct, avg_x_corrs_direct, top_x_corrs_direct, \
            all_azb_corrs_direct, azb_avg_corrs_direct = \
        get_corr_stats(df=expl_df, crispr_cm=crispr_corr_matrix,
                       rnai_cm=rnai_corr_matrix, so_pairs=pairs_axb_direct)
    """

    # a-x-b AND NOT direct
    logger.info("Getting correlations for a-x-b AND NOT direct")
    all_x_corrs_no_direct, avg_x_corrs_no_direct, top_x_corrs_no_direct, \
        all_azb_corrs_no_direct, azb_avg_corrs_no_direct = \
        get_corr_stats(df=expl_df, z_sc_cm=z_corr, so_pairs=pairs_axb_only)

    """
    # a-x-b (with and without direct)
    logger.info("Getting correlations for all a-x-b (direct and indirect)")
    all_x_corrs_union, avg_x_corrs_union, top_x_corrs_union, \
            all_azb_corrs_union, azb_avg_corrs_union = \
        get_corr_stats(df=expl_df, crispr_cm=crispr_corr_matrix,
                       rnai_cm=rnai_corr_matrix, so_pairs=ab_axb_union)
    """

    # All corrs for range (all pairs regardless of explanation type)
    # all_x_corrs, top_x_corrs, corr_vs_maxavg = \
    #     get_corr_stats(df=expl_df, crispr_cm=crispr_corr_matrix,
    #                    rnai_cm=rnai_corr_matrix, so_pairs=all_ab_corr_pairs)
    """
    return {'axb_and_dir': {'all_x_corrs': all_x_corrs_direct,
                            'avg_x_corrs': avg_x_corrs_direct,
                            'top_x_corrs': top_x_corrs_direct,
                            'all_azb_corrs': all_azb_corrs_direct,
                            'azb_avg_corrs': azb_avg_corrs_direct},
            'axb_not_dir': {'all_x_corrs': all_x_corrs_no_direct,
                            'avg_x_corrs': avg_x_corrs_no_direct,
                            'top_x_corrs': top_x_corrs_no_direct,
                            'all_azb_corrs': all_azb_corrs_no_direct,
                            'azb_avg_corrs': azb_avg_corrs_no_direct},
            'all_axb': {'all_x_corrs': all_x_corrs_union,
                        'avg_x_corrs': avg_x_corrs_union,
                        'top_x_corrs': top_x_corrs_union,
                        'all_azb_corrs': all_azb_corrs_union,
                        'azb_avg_corrs': azb_avg_corrs_union}}
    """
    return {'axb_not_dir': {'all_x_corrs': all_x_corrs_no_direct,
                            'avg_x_corrs': avg_x_corrs_no_direct,
                            'top_x_corrs': top_x_corrs_no_direct,
                            'all_azb_corrs': all_azb_corrs_no_direct,
                            'azb_avg_corrs': azb_avg_corrs_no_direct}}


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
