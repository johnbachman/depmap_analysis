import ast
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from os import path, environ
from depmap_analysis.util.io_functions import pickle_open
from depmap_analysis.network_functions.depmap_network_functions import \
    _z_sc, pass_filter


logger = logging.getLogger('DepMap Corr Stats')


# Stats for 19Q4 release for crispr set and for version 5 of the RNAi set
cmu = 0.003076
csig = 0.056813
rmu = 0.006854
rsig = 0.077614

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


def mean_z_score(mu1, sig1, c1, mu2, sig2, c2):
    return 0.5 * _z_sc(num=c1, mu=mu1, sigma=sig1) + \
        0.5 * _z_sc(num=c2, mu=mu2, sigma=sig2)


def func1(d):
    return mean_z_score(mu1=cmu, sig1=csig, c1=d['crispr'],
                        mu2=rmu, sig2=rsig, c2=d['rnai'])


def get_corr_stats(df, crispr_cm, rnai_cm, so_pairs):
    all_ = []
    top = []
    ab_to_axb_avg = []
    for subj, obj in so_pairs:
        ab_avg_corrs = []
        path_row = df[(df['subj'] == subj) &
                      (df['obj'] == obj) &
                      (df['type'] == 'pathway')]
        x_list = [x[0] for x in path_row['X'].values[0]]
        comb_zsc = path_row['comb_zsc'].values[0]
        warn = 0
        for x in x_list:
            if x in crispr_cm.columns and x in rnai_cm.columns and \
                    pass_filter(corr1=crispr_cm.loc[subj, x],
                                mu1=cmu, sigma1=csig,
                                corr2=rnai_cm.loc[subj, x],
                                mu2=rmu, sigma2=rsig,
                                filter_type='sign'):
                ax_corr = mean_z_score(
                    mu1=cmu, sig1=csig, c1=crispr_cm.loc[subj, x],
                    mu2=rmu, sig2=rsig, c2=rnai_cm.loc[subj, x])
                xb_corr = mean_z_score(
                    mu1=cmu, sig1=csig, c1=crispr_cm.loc[x, obj],
                    mu2=rmu, sig2=rsig, c2=rnai_cm.loc[x, obj])
                all_ += [ax_corr, xb_corr]
                ab_avg_corrs.append(0.5 * ax_corr + 0.5 * xb_corr)
            else:
                # if warn < 3:
                #     warn += 1
                #     logger.warning('%s does not exist in the crispr and/or '
                #                    'rnai correlation matrices.' % x)
                # else:
                #     logger.warning('Muting warnings...')
                continue
        # if warn:
        #     logger.warning('%d missing X genes out of %d in correlation '
        #                    'matrices' % (warn, len(x_list)))
        if len(ab_avg_corrs) > 0:
            max_magn_avg = max(ab_avg_corrs) if max(ab_avg_corrs) > abs(min(
                ab_avg_corrs)) else min(ab_avg_corrs)
            top.append((subj, obj, max_magn_avg))
            ab_to_axb_avg.append((subj, obj, comb_zsc, max_magn_avg))

    return all_, top, ab_to_axb_avg


def main(expl_df, crispr_corr_matrix, rnai_corr_matrix):
    # Limit to any a-x-b OR a-b expl (this COULD include explanations where
    # 'direct' and NOT 'pathway' is the explanation, but this should be a
    # very small set)
    expl_df = expl_df[(expl_df['type'] == 'pathway') |
                      (expl_df['type'] == 'direct')]

    # Re-map the columns contatining string representations of objects
    expl_df.meta_data = expl_df['meta_data'].apply(lambda x:
                                                   ast.literal_eval(x))
    expl_df.X = expl_df['X'].apply(lambda x: ast.literal_eval(x))

    # Get combined z score
    expl_df['comb_zsc'] = None
    for index, row in expl_df.iterrows():
        cd = row.meta_data
        if not cd:
            cd = {'crispr': crispr_corr_matrix.loc[row.subj, row.obj],
                  'rnai': rnai_corr_matrix.loc[row.subj, row.obj]}
        row.comb_zsc = func1(cd)

    # Get all correlation pairs
    all_ab_corr_pairs = set(map(lambda p: tuple(p),
                                expl_df[['subj', 'obj']].values))
    # Pairs where a-x-b AND a-b explanation exists
    pairs_axb_direct = set()

    # Pairs where a-x-b AND NOT a-b explanation exists
    pairs_axb_only = set()

    # all a-x-b "pathway" explanations, should be union of the above two
    pairs_any_axb = set()

    for s, o in all_ab_corr_pairs:
        # Get all interaction types associated with given subject s and
        # object o
        int_types = set(expl_df['type'][(expl_df['subj'] == s) &
                                        (expl_df['obj'] == o)].values)
        if 'pathway' in int_types:
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
    all_x_corrs_direct, top_x_corrs_direct, corr_vs_maxavg_direct = \
        get_corr_stats(df=expl_df, crispr_cm=crispr_corr_matrix,
                       rnai_cm=rnai_corr_matrix, so_pairs=pairs_axb_direct)

    # a-x-b AND NOT direct
    all_x_corrs_no_direct, top_x_corrs_no_direct, corr_vs_maxavg_no_direct = \
        get_corr_stats(df=expl_df, crispr_cm=crispr_corr_matrix,
                       rnai_cm=rnai_corr_matrix, so_pairs=pairs_axb_only)

    # a-x-b (with and without direct)
    all_x_corrs_union, top_x_corrs_union, corr_vs_maxavg_union = \
        get_corr_stats(df=expl_df, crispr_cm=crispr_corr_matrix,
                       rnai_cm=rnai_corr_matrix, so_pairs=ab_axb_union)

    # All corrs for range (all pairs regardless of explanation type)
    # all_x_corrs, top_x_corrs, corr_vs_maxavg = \
    #     get_corr_stats(df=expl_df, crispr_cm=crispr_corr_matrix,
    #                    rnai_cm=rnai_corr_matrix, so_pairs=all_ab_corr_pairs)
    return {'axb_and_dir': {'all_x_corrs': all_x_corrs_direct,
                            'top_x_corrs': top_x_corrs_direct,
                            'corr_vs_maxavg': corr_vs_maxavg_direct},
            'axb_not_dir': {'all_x_corrs': all_x_corrs_no_direct,
                            'top_x_corrs': top_x_corrs_no_direct,
                            'corr_vs_maxavg': corr_vs_maxavg_no_direct},
            'all_axb': {'all_x_corrs': all_x_corrs_union,
                        'top_x_corrs': top_x_corrs_union,
                        'corr_vs_maxavg': corr_vs_maxavg_union}}


if __name__ == '__main__':
    expl_pairs_csv = '/home/klas/repos/depmap_analysis/output_data' \
                     '/19Q4_hgnc_fplx/{range}/_explanations_of_pairs.csv'
    crispr_corr = '/home/klas/repos/depmap_analysis/input_data/depmap' \
                  '/19Q4/_crispr_all_correlations.h5'
    rnai_corr = '/home/klas/repos/depmap_analysis/input_data/depmap' \
                '/demeter/_rnai_all_correlations.h5'
    # Load rnai and crispr correlation sets for correlation lookup
    # - could possibly merge them to cut time in half and save RAM
    logger.info('Loading correlation matrices...')
    ccorr_matrix = pd.read_hdf(crispr_corr)
    names = [n.split()[0] for n in ccorr_matrix.columns.values]
    ccorr_matrix.columns = names
    ccorr_matrix.index = names
    rcorr_matrix = pd.read_hdf(rnai_corr)
    names = [n.split()[0] for n in rcorr_matrix.columns.values]
    rcorr_matrix.columns = names
    rcorr_matrix.index = names

    # Load _explanations_of_pairs.csv for each range
    for sd in ['1_2sd', '2_3sd', '3_4sd', '4_5sd', '5_sd', 'rnd']:
        fname = expl_pairs_csv.format(range=sd)
        if not path.isfile(fname):
            logger.info('Skipping %s, file does not exist' % fname)
        else:
            logger.info('Getting pairs from %s' % fname)
            df = pd.read_csv(fname, delimiter=',')
            results = main(expl_df=df,
                           crispr_corr_matrix=ccorr_matrix,
                           rnai_corr_matrix=rcorr_matrix)
            expl_dir = path.dirname(expl_pairs_csv.format(range=sd))

            for k, v in results.items():
                # Plot the
                if len(v['all_x_corrs']) > 0:
                    plt.hist([v['all_x_corrs']], bins='auto')
                    plt.savefig(path.join(expl_dir, 'all_corrs_%s.png' % k),
                                format='png')
                    plt.show()
                else:
                    logger.warning('Empty result for %s in range %s'
                                   % (k, sd))
