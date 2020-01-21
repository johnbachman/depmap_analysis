import logging
import argparse
from os import path
import pandas as pd

from depmap_analysis.network_functions.depmap_network_functions import \
    merge_corr_df, raw_depmap_to_corr


logger = logging.getLogger('DepMap PreProcessing')


def run_corr_merge(crispr_raw=None, rnai_raw=None,
                   crispr_corr=None, rnai_corr=None, output_dir=None):
    if not crispr_raw and not crispr_corr:
        raise ValueError('Need to provide one of crispr_raw or '
                         'cripsr_corr')
    if not rnai_raw and not rnai_corr:
        raise ValueError('Need to provide one of rnai_raw or rnai_corr')

    # First check for correlation matrix, then get it if it doesn't exist
    if crispr_corr:
        crispr_corr_df = pd.read_hdf(crispr_corr)
    else:
        # Create new one, write to input file's directory
        crispr_corr_df = raw_depmap_to_corr(pd.read_csv(crispr_raw,
                                                        index_col=0))
        in_dir = output_dir if output_dir else path.dirname(
            crispr_raw)
        logger.info(f'Saving correlation matrix to {in_dir}')
        name = '_crispr_all_correlations.h5'
        crispr_corr_df.to_hdf(path.join(in_dir, name), name)

    if rnai_corr:
        rnai_corr_df = pd.read_hdf(rnai_corr)
    else:
        # Create new one, write to input file's directory
        rnai_corr_df = raw_depmap_to_corr(pd.read_csv(rnai_raw,
                                                      index_col=0))
        in_dir = output_dir if output_dir else path.dirname(
            rnai_raw)
        logger.info(f'Saving correlation matrix to {in_dir}')
        name = '_rnai_all_correlations.h5'
        rnai_corr_df.to_hdf(path.join(in_dir, name), name)

    # Merge the correlation matrices
    return merge_corr_df(crispr_corr_df, rnai_corr_df)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('DepMap Data pre-processing')
    # Input dirs
    parser.add_argument('--crispr-raw',
                        help='The raw CRISPR gene effect data. File name '
                             'should match *gene_effect.csv')
    parser.add_argument('--rnai-raw',
                        help='The raw RNAi gene dependency data. File name '
                             'should match *gene_dep_scores.csv')
    # Option to start from raw correlations matrix instead of running
    # correlation calculation directly
    parser.add_argument('--crispr-corr',
                        help='The file containing an hdf compressed '
                             'correlation data frame of the crispr data. If '
                             'this file is provided, the raw crispr data is '
                             'ignored.')
    parser.add_argument('--rnai-corr',
                        help='The file containing an hdf compressed '
                             'correlation data frame of the rnai data. If '
                             'this file is provided, the raw rnai data is '
                             'ignored.')
    # Output dirs
    parser.add_argument('--output-dir', '-o',
                        help='Optional. A directory where to put the '
                             'output of the script. If not provided the '
                             'ouput will go to the directory/ies where '
                             'the corresponding input came from, i.e. the '
                             'output from the RNAi input will be placed '
                             'in the RNAi input directory ')

    args = parser.parse_args()
    options = {
        'crispr_raw': args.crispr_raw,
        'rnai_raw': args.rnai_raw,
        'crispr_corr': args.crispr_corr,
        'rnai_corr': args.rnai_corr,
        'output_dir': args.output_dir
    }

    # Run main script
    z_corr = run_corr_merge(**options)

    # Write merged correlations combined z score
    outdir = args.output_dir if args.output_dir else (path.dirname(
        args.crispr_corr) if args.crispr_corr else path.dirname(
        args.crispr_raw))
    logger.info(f'Writing combined correlations to {outdir}')
    fname = 'combined_z_score.h5'
    z_corr.to_hdf(path.join(outdir, fname), fname)
