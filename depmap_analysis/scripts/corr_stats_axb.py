import ast
import pandas as pd
from os import path, environ
from depmap_analysis.util.io_functions import pickle_open
from depmap_analysis.network_functions.depmap_network_functions import _z_sc


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

expl_csv = '/home/klas/repos/depmap_analysis/output_data/19Q4_hgnc_fplx/' \
           '{range}/_explanations_of_pairs.csv'
crispr_corr = '/home/klas/repos/depmap_analysis/input_data/depmap' \
              '/19Q4/_crispr_all_correlations.h5'
rnai_corr = '/home/klas/repos/depmap_analysis/input_data/depmap' \
            '/demeter/_rnai_all_correlations.h5'


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
        x_list = [x[0] for x in df['X'][
            (df['subj'] == s) & (df['obj'] == o) & (df['type'] == 'pathway')
            ].values[0]]
        comb_zsc = df['comb_zsc']
        for x in x_list:
            if x in crispr_cm.columns and x in rnai_cm.columns:
                ax_corr = mean_z_score(
                    mu1=cmu, sig1=csig, c1=crispr_cm.loc[subj, x],
                    mu2=rmu, sig2=rsig, c2=rnai_cm.loc[subj, x])
                xb_corr = mean_z_score(
                    mu1=cmu, sig1=csig, c1=crispr_cm.loc[x, obj],
                    mu2=rmu, sig2=rsig, c2=rnai_cm.loc[x, obj])
                all_ += [ax_corr, xb_corr]
                ab_avg_corrs.append(0.5 * ax_corr + 0.5 * xb_corr)
        top.append(max(ab_avg_corrs))
        ab_to_axb_avg.append((comb_zsc, max(ab_avg_corrs)))

    return all_, top, ab_to_axb_avg


# Load rnai and crispr correlation sets for correlation lookup
# - could possibly merge them to cut time in half and save RAM
crispr_corr_matrix = pd.read_hdf(crispr_corr)
names = [n[0] for n in crispr_corr_matrix.columns.values]
crispr_corr_matrix.columns = names
crispr_corr_matrix.index = names
rnai_corr_matrix = pd.read_hdf(rnai_corr)
names = [n[0] for n in rnai_corr_matrix.columns.values]
rnai_corr_matrix.columns = names
rnai_corr_matrix.index = names


# Load _explanations_of_pairs.csv for your range
sd_range = '5_sd'
expl_df = pd.read_csv(expl_csv.format(range=sd_range), delimiter=',')
# Limit to all a-x-b expl
# TodO We want to keep all rows (or at least 'direct'+'pathway) for the pairs
#  where 'pathway' is part of the explanations, this removes all but the
#  'pathway' rows
expl_df = expl_df[expl_df['type'] == 'pathway']

# Re-map the columns contatining string representations of objects
expl_df.meta_data = expl_df['meta_data'].apply(lambda x: ast.literal_eval(x))
expl_df.X = expl_df['X'].apply(lambda x: ast.literal_eval(x))

# Get combined z score
expl_df['comb_zsc'] = None
for index, row in expl_df.iterrows():
    sd = row.meta_data
    if not sd:
        sd = {'crispr': crispr_corr_matrix.loc[row.subj, row.obj],
              'rnai': rnai_corr_matrix.loc[row.subj, row.obj]}
    row.comb_zsc = func1(sd)


# Get all correlation pairs
all_ab_corr_pairs = set(map(lambda p: tuple(p),
                            expl_df[['subj', 'obj']].values))
# Get pairs where a-x-b AND a-b explanation exists
ab_axb_direct = set()
# Get pairs where a-x-b AND NOT a-b explanation exists
ab_axb_only = set()
for s, o in all_ab_corr_pairs:
    # Get all interaction types associated with given subject s and object o
    int_types = set(expl_df['type'][(expl_df['subj'] == s) &
                                    (expl_df['obj'] == o)].values)
    if 'pathway' in int_types:
        if 'direct' in int_types:
            ab_axb_direct.add((s, o))
        else:
            ab_axb_only.add((s, o))
# Get the combined set of pairs
ab_axb_union = ab_axb_direct.union(ab_axb_only)

# Sort by comb_zsc
# print(expl_df.sort_values(by='comb_zsc', ascending=False).head(10))

# a-x-b AND direct
all_x_corrs_direct, top_x_corrs_direct, corr_vs_maxavg_direct = \
    get_corr_stats(df=expl_df, crispr_cm=crispr_corr_matrix,
                   rnai_cm=rnai_corr_matrix, so_pairs=ab_axb_direct)

# a-x-b AND NOT direct
all_x_corrs_no_direct, top_x_corrs_no_direct, corr_vs_maxavg_no_direct = \
    get_corr_stats(df=expl_df, crispr_cm=crispr_corr_matrix,
                   rnai_cm=rnai_corr_matrix, so_pairs=ab_axb_only)

# a-x-b (with and without direct)
all_x_corrs_union, top_x_corrs_union, corr_vs_maxavg_union = \
    get_corr_stats(df=expl_df, crispr_cm=crispr_corr_matrix,
                   rnai_cm=rnai_corr_matrix, so_pairs=ab_axb_union)

# All corrs for range (all pairs regardless of explanation type)
all_x_corrs, top_x_corrs, corr_vs_maxavg = \
    get_corr_stats(df=expl_df, crispr_cm=crispr_corr_matrix,
                   rnai_cm=rnai_corr_matrix, so_pairs=all_ab_corr_pairs)


# ToDo
#  - Functionalize the setup that has to be repeated for each SD range (so
#    we dont reload the correlation lookup
#  - Get the underlying correlation distribution for each range (all
#    correlations with a-x-b explantion, regardless of direct)
#  - also get the
