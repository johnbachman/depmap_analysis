import logging
import argparse
from pathlib import Path

import pandas as pd

from depmap_analysis.network_functions.depmap_network_functions import \
    merge_corr_df, raw_depmap_to_corr
from depmap_analysis.network_functions.indra_network import get_top_ranked_name

logger = logging.getLogger('DepMap PreProcessing')


def run_corr_merge(crispr_raw=None, rnai_raw=None,
                   crispr_corr=None, rnai_corr=None,
                   output_dir='correlation_output',
                   remove_self_corr=True, dropna=False, random_sampl=0):
    """Return a merged correlation matrix from DepMap data

    Start with with either the raw DepMap files or pre-calculated
    correlation matrices

    Parameters
    ----------
    crispr_raw : str|pd.DataFrame
        Path to the raw crispr data. This file is typically named
        'Achilles_gene_effect.csv' at the DepMap portal.
    rnai_raw : str|pd.DataFrame
        Path to the raw RNAi data. This file is typically named
        'D2_combined_gene_dep_scores.csv'
    crispr_corr : str|pd.DataFrame
        Path to the pre-calculated crispr data matrix. This data structure
        is the result from running `crispr_raw_df.corr()`.
    rnai_corr : str|pd.DataFrame
        Path to the pre-calculated rnai data matrix. This data structure
        is the result from running `rnai_raw_df.corr()`.
    output_dir : str
        If used, write the correlation matrices to this directory.
        Otherwise they will be written to the same directory as the raw
        input data.
    remove_self_corr : bool
        If True, remove self correlations from the resulting DataFrame.
        Default: True
    random_sampl : int
        If specified, provides the size of the final correlation matrix
        where the genes are picked at random from the intersection of genes
        from both the RNAI and CRISPR data sets.
    dropna : bool
        If True, return the result of
        corr_df.dropna(axis=0, how='all').dropna(axis=1, how='all')
        Default: False.

    Returns
    -------
    pd.DataFrame
        A data frame containing the combined z-score matrix with NaN's
        removed.
    """
    if crispr_raw is None and crispr_corr is None:
        raise ValueError('Need to provide one of crispr_raw or cripsr_corr')
    if rnai_raw is None and rnai_corr is None:
        raise ValueError('Need to provide one of rnai_raw or rnai_corr')

    # First check for correlation matrix, then get it if it doesn't exist
    if crispr_corr:
        if isinstance(crispr_corr, str):
            logger.info(f'Reading crispr correlations from file {crispr_corr}')
            crispr_corr_df = pd.read_hdf(crispr_corr)
        else:
            crispr_corr_df = crispr_corr
    else:
        # Create new one, write to input file's directory
        if isinstance(crispr_raw, str):
            logger.info(f'Reading raw DepMap data from {crispr_raw}')
            crispr_raw_df = pd.read_csv(crispr_raw, index_col=0)
        else:
            crispr_raw_df = crispr_raw
        crispr_corr_df = raw_depmap_to_corr(crispr_raw_df, dropna=True)

        crispr_fpath = Path(output_dir).joinpath('_crispr_all_correlations.h5')
        logger.info(f'Saving crispr correlation matrix to {crispr_fpath}')
        if not crispr_fpath.parent.is_dir():
            crispr_fpath.parent.mkdir(parents=True, exist_ok=True)
        crispr_corr_df.to_hdf(crispr_fpath.absolute(), 'corr')

    if rnai_corr:
        if isinstance(rnai_corr, str):
            logger.info(f'Reading rnai correlations from file {crispr_corr}')
            rnai_corr_df = pd.read_hdf(rnai_corr)
        else:
            rnai_corr_df = rnai_corr
    else:
        # Create new one, write to input file's directory
        if isinstance(rnai_raw, str):
            logger.info(f'Reading raw DepMap data from {rnai_raw}')
            rnai_raw_df = pd.read_csv(rnai_raw, index_col=0)
        else:
            rnai_raw_df = rnai_raw

        # Check if we need to transpose the df
        if len(set(crispr_corr_df.columns.values) &
               set([n.split()[0] for n in rnai_raw_df.columns])) == 0:
            logger.info('Transposing RNAi raw data dataframe...')
            rnai_raw_df = rnai_raw_df.T

        rnai_corr_df = raw_depmap_to_corr(rnai_raw_df, dropna=True)

        rnai_fpath = Path(output_dir).joinpath('_rnai_all_correlations.h5')
        if not rnai_fpath.parent.is_dir():
            rnai_fpath.mkdir(parents=True, exist_ok=True)
        logger.info(f'Saving rnai correlation matrix to {rnai_fpath}')
        rnai_corr_df.to_hdf(rnai_fpath.absolute().as_posix(), 'corr')

    # Merge the correlation matrices
    z_cm = merge_corr_df(crispr_corr_df, rnai_corr_df,
                         remove_self_corr, dropna)

    if random_sampl and random_sampl < len(z_cm.columns):
        # Get n random rows
        z_cm = z_cm.sample(n=random_sampl)

        # Make square
        z_cm = z_cm[list(z_cm.index.values)]

    assert z_cm.notna().sum().sum() > 0, \
        print(f'Correlation matrix is empty')

    return z_cm


def get_drug_corr_matrix(drug_resp_file, drug_info_file):
    def _get_drug_name(drug_id):
        drug_rec = drug_info_df.loc[drug_id]
        return drug_rec['name']
    if isinstance(drug_resp_file, (str, Path)):
        drug_resp_df = pd.read_csv(drug_resp_file, index_col=0)
    elif isinstance(drug_resp_file, pd.DataFrame):
        drug_resp_df = drug_resp_file
    if isinstance(drug_info_file, (str, Path)):
        drug_info_df = pd.read_csv(drug_info_file, index_col=0)
    elif isinstance(drug_info_file, pd.DataFrame):
        drug_info_df = drug_info_file

    # Translate ids to names
    col_names = [_get_drug_name(did) for did in drug_resp_df.columns]
    # Gildaify?
    col_names = [get_top_ranked_name(name)[-1] or name for name in col_names]
    drug_resp_df.columns = col_names

    # Drop duplicate columns
    drug_resp_df = \
        drug_resp_df.loc[:, ~drug_resp_df.columns.duplicated()]

    # Make correlation matrix
    corr_df = drug_resp_df.corr()

    return corr_df


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
                             'in the RNAi input directory. The combined '
                             'z-score matrix will be written to the crispr '
                             'input directory if this option is not '
                             'provided.')
    parser.add_argument('--random', '-r', type=int,
                        help='Optional. If specified, provide the size of '
                             'the final correlation matrix where the genes '
                             'are picked at random from the intersection of '
                             'genes from both the RNAI and CRISPR data sets.')
    parser.add_argument('--fname', help='A file name for the output '
                                        'correlations DataFrame.')

    args = parser.parse_args()
    options = {
        'crispr_raw': args.crispr_raw,
        'rnai_raw': args.rnai_raw,
        'crispr_corr': args.crispr_corr,
        'rnai_corr': args.rnai_corr,
        'output_dir': args.output_dir,
        'random_sampl': args.random
    }

    # Run main script
    z_corr = run_corr_merge(**options)

    # Write merged correlations combined z score
    outdir = args.output_dir if args.output_dir else (Path(
        args.crispr_corr).parent.as_posix() if args.crispr_corr else Path(
        args.crispr_raw).parent.as_posix())
    logger.info(f'Writing combined correlations to {outdir}')
    fname = args.fname if args.fname else 'combined_z_score.h5'
    z_corr.to_hdf(Path(outdir, fname), 'zsc')
