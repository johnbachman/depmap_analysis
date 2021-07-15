import logging
import argparse
from typing import Dict, Union, Optional
from pathlib import Path

import numpy as np
import pandas as pd

from depmap_analysis.util import io_functions as io
from depmap_analysis.network_functions.famplex_functions import ns_id_xref, \
    ns_id_to_name
from depmap_analysis.util.statistics import *

logger = logging.getLogger(__name__)
__all__ = ['run_corr_merge', 'drugs_to_corr_matrix', 'get_mitocarta_info']

MERGE_METHODS = ('average', 'stouffer')
Z_SC_METHODS = ('standard', 't', 'beta')

path_obj = Union[str, pd.DataFrame]


def run_corr_merge(crispr_raw: Optional[path_obj] = None,
                   rnai_raw: Optional[path_obj] = None,
                   crispr_corr: Optional[path_obj] = None,
                   rnai_corr: Optional[path_obj] = None,
                   output_dir: str = 'correlation_output',
                   remove_self_corr: bool = False,
                   dropna: bool = False,
                   random_sampl: int = 0,
                   save_corr_files: bool = False,
                   z_corr_path: Optional[str] = None):
    """Return a merged correlation matrix from DepMap data

    Start with with either the raw DepMap files or pre-calculated
    correlation matrices

    Parameters
    ----------
    crispr_raw :
        Path to the raw crispr data. This file is typically named
        'Achilles_gene_effect.csv' at the DepMap portal.
    rnai_raw :
        Path to the raw RNAi data. This file is typically named
        'D2_combined_gene_dep_scores.csv'
    crispr_corr :
        Path to the pre-calculated crispr data matrix. This data structure
        is the result from running `crispr_raw_df.corr()`.
    rnai_corr :
        Path to the pre-calculated rnai data matrix. This data structure
        is the result from running `rnai_raw_df.corr()`.
    output_dir :
        If used, write the correlation matrices to this directory.
        Otherwise they will be written to the same directory as the raw
        input data.
    remove_self_corr :
        If True, remove self correlations from the resulting DataFrame.
        Default: False
    dropna :
        If True, drop NaN values after merging the dataframes using
        `pd.DataFrame.dropna(axis=0, how='all').dropna(axis=1, how='all')`
    random_sampl :
        If specified, provides the size of the final correlation matrix
        where the genes are picked at random from the intersection of genes
        from both the RNAI and CRISPR data sets.
    save_corr_files :
        If True, save the intermediate correlation data frames for both
        crispr and rnai. Default: True.
    z_corr_path :
        If provided, save the final correlation dataframe here

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
            crispr_raw_df = io.file_opener(crispr_raw, index_col=0)
        else:
            crispr_raw_df = crispr_raw
        crispr_corr_df = raw_depmap_to_corr(crispr_raw_df, split_names=True,
                                            dropna=False)

        if save_corr_files:
            crispr_fpath = Path(output_dir).joinpath(
                '_crispr_all_correlations.h5')
            logger.info(f'Saving crispr correlation matrix to {crispr_fpath}')
            if not crispr_fpath.parent.is_dir():
                crispr_fpath.parent.mkdir(parents=True, exist_ok=True)
            crispr_corr_df.to_hdf(crispr_fpath.absolute().as_posix(), 'corr')

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
            rnai_raw_df = io.file_opener(rnai_raw, index_col=0)
        else:
            rnai_raw_df = rnai_raw

        # Check if we need to transpose the df
        if len(set(crispr_corr_df.columns.values) &
               set([n.split()[0] for n in rnai_raw_df.columns])) == 0:
            logger.info('Transposing RNAi raw data dataframe...')
            rnai_raw_df = rnai_raw_df.T

        rnai_corr_df = raw_depmap_to_corr(rnai_raw_df, split_names=True,
                                          dropna=False)

        if save_corr_files:
            rnai_fpath = Path(output_dir).joinpath('_rnai_all_correlations.h5')
            if not rnai_fpath.parent.is_dir():
                rnai_fpath.parent.mkdir(parents=True, exist_ok=True)
            logger.info(f'Saving rnai correlation matrix to {rnai_fpath}')
            rnai_corr_df.to_hdf(rnai_fpath.absolute().as_posix(), 'corr')

    # Merge the correlation matrices
    z_cm = merge_corr_df(crispr_corr_df, rnai_corr_df, remove_self_corr,
                         merge_method='stouffer', z_sc_method='beta',
                         dropna=dropna)

    if random_sampl and random_sampl < len(z_cm.columns):
        # Get n random rows
        z_cm = z_cm.sample(n=random_sampl)

        # Make square
        z_cm = z_cm[list(z_cm.index.values)]

    assert z_cm.notna().sum().sum() > 0, 'Correlation matrix is empty'

    if z_corr_path:
        zc_path = Path(z_corr_path)
        zc_path.parent.mkdir(parents=True, exist_ok=True)
        z_cm.to_hdf(zc_path.absolute().as_posix(), 'corr')

    return z_cm


def drugs_to_corr_matrix(raw_file: str, info_file: str):
    """Preprocess and create a correlation matrix from raw drug data

    Parameters
    ----------
    raw_file : str
        Path to DepMap PRISM drug repurposing data file. Should match
        primary-screen-replicate-collapsed-logfold-change.csv
    info_file : str
        Path to DepMap PRISM drug repurposing info file. Should match
        primary-screen-replicate-collapsed-treatment-info.csv
    """
    def _get_drug_name(drug_id):
        drug_rec = info_df.loc[drug_id]
        return drug_rec['name']

    raw_df: pd.DataFrame = io.file_opener(raw_file, index_col=0)
    info_df: pd.DataFrame = io.file_opener(info_file, index_col=0)
    col_names = [_get_drug_name(did) for did in raw_df.columns]
    raw_df.columns = col_names

    return raw_depmap_to_corr(raw_df)


def raw_depmap_to_corr(depmap_raw_df: pd.DataFrame,
                       split_names: bool = False,
                       dropna: bool = False):
    """Pre-process and create a correlation matrix

    Any multi indexing is removed. Duplicated columns are also removed.

    Parameters
    ----------
    depmap_raw_df : pd.DataFrame
        The raw data from the DepMap portal as a pd.DataFrame
    split_names : bool
        If True, check if column names contain whitespace and if the do,
        split the name and keep the first part.
    dropna : bool
        If True, drop nan columns (should be genes) before calculating the
        correlations

    Returns
    -------
    corr : pd.DataFrame
        A pd.DataFrame containing the pearson correlations of the raw data.
    """
    # Rename
    if split_names and len(depmap_raw_df.columns[0].split()) > 1:
        logger.info('Renaming columns to contain only first part of name')
        gene_names = [n.split()[0] for n in depmap_raw_df.columns]
        depmap_raw_df.columns = gene_names

    # Drop duplicates
    if sum(depmap_raw_df.columns.duplicated()) > 0:
        logger.info('Dropping duplicated columns')
        depmap_raw_df = \
            depmap_raw_df.loc[:, ~depmap_raw_df.columns.duplicated()]

    # Drop nan's
    if dropna:
        logger.info('Dropping nan columns (axis=1)')
        depmap_raw_df = depmap_raw_df.dropna(axis=1)

    # Calculate correlation
    logger.info('Calculating data correlation matrix. This can take up to '
                '10 min depending on the size of the dataframe.')
    corr = depmap_raw_df.corr()
    logger.info('Done calculating data correlation matrix.')
    return corr


def merge_corr_df(corr_df, other_corr_df, remove_self_corr=True,
                  merge_method: str = 'stouffer', z_sc_method: str = 'beta',
                  dropna: bool = False) -> pd.DataFrame:
    """Merge two correlation matrices containing their combined z-scores

    Parameters
    ----------
    corr_df :
        A square pandas DataFrame containing gene-gene correlation values
        to be merged with other_corr_df.
    other_corr_df :
        A square pandas DataFrame containing gene-gene correlation values
        to be merged with corr_df.
    remove_self_corr :
        If True, remove self correlations from the resulting DataFrame.
        Default: True
    dropna :
        If True, drop NaN values after merging the dataframes using
        `pd.DataFrame.dropna(axis=0, how='all').dropna(axis=1, how='all')`
    merge_method :
        The method used when merging the dataframes. Available methods:
        - 'average': z_merged = (z1 + z2) / 2
        - 'stouffer': z_merged = (z1 + z2) / sqrt(2)
    z_sc_method :
        The method to use when getting the z-scores. Available methods:
        - 'standard' ("standard score"): z = (x - mu) / sigma.
        - 'beta': based on getting the log of the p-values, here using the
           beta distribution method.
        - 't': based on getting the log of the p-values, here using the
           t statistic method.

    Returns
    -------
    :
        A merged correlation matrix containing the merged values a z-scores
    """

    def _rename(corr: pd.DataFrame) -> pd.DataFrame:
        gene_names = [n.split()[0] for n in corr.columns]
        corr.columns = gene_names
        corr.index = gene_names
        return corr

    # Rename columns/indices to gene name only
    if len(corr_df.columns[0].split()) > 1:
        corr_df = _rename(corr_df)
    if len(other_corr_df.columns[0].split()) > 1:
        other_corr_df = _rename(other_corr_df)

    # Get corresponding z-score matrices
    logger.info('Getting z-score matrix of first data frame.')
    corr_z = _z_scored(corr_df, method=z_sc_method)
    logger.info('Getting z-score matrix of second data frame.')
    other_z = _z_scored(other_corr_df, method=z_sc_method)

    # Merge
    dep_z = _merge_z_corr(zdf=corr_z, other_z_df=other_z, method=merge_method,
                          remove_self_corr=remove_self_corr, dropna=dropna)
    assert dep_z.notna().sum().sum() > 0, 'Correlation matrix is empty!'
    return dep_z


def _merge_z_corr(zdf: pd.DataFrame, other_z_df: pd.DataFrame,
                  remove_self_corr: bool, method: str = 'average',
                  dropna: bool = False) -> pd.DataFrame:
    if method not in MERGE_METHODS:
        raise ValueError(f'Unrecognized merge method {method}. Valid methods '
                         f'are {", ".join(MERGE_METHODS)}')

    def _average(z1: pd.DataFrame, z2: pd.DataFrame) -> pd.DataFrame:
        return (z1 + z2) / 2

    def _stouffer_z(z1: pd.DataFrame, z2: pd.DataFrame) -> pd.DataFrame:
        return (z1 + z2) / np.sqrt(2)

    if method == 'average':
        dep_z = _average(zdf, other_z_df)
    else:
        dep_z = _stouffer_z(zdf, other_z_df)

    if dropna:
        dep_z = dep_z.dropna(axis=0, how='all').dropna(axis=1, how='all')

    if remove_self_corr:
        # Assumes the max correlation ONLY occurs on the diagonal
        self_corr_value = dep_z.loc[dep_z.columns[0], dep_z.columns[0]]
        dep_z = dep_z[dep_z != self_corr_value]
    assert dep_z.notna().sum().sum() > 0, 'Correlation matrix is empty!'
    return dep_z


def _z_scored_pvals(corr_df: pd.DataFrame, method: str = 'beta') \
        -> pd.DataFrame:
    # Get z-scores based on getting correlation dataframe
    if method not in ('beta', 't'):
        raise ValueError(f'Unrecognized p-value z-score method {method}. '
                         f'Available methods are beta or t.')

    # Get sample sizes from get_n
    n_df = get_n(recalculate=True, data_df=corr_df)

    # Get logp from get_logp
    logp_df = get_logp(recalculate=True, data_n=n_df, data_corr=corr_df,
                       method=method)

    # Get z-scores from get_z
    z_df = get_z(recalculate=True, data_logp=logp_df, data_corr=corr_df)
    return z_df


def _z_scored(corr: pd.DataFrame, method: str = 'beta') -> pd.DataFrame:
    if method not in Z_SC_METHODS:
        raise ValueError(f'Unrecognized z-score method {method}. Valid '
                         f'methods are {", ".join(Z_SC_METHODS)}')
    logger.info(f'Getting z-scores using {method} method')
    if method == 'standard':
        mean = _get_mean(corr)
        sd = _get_sd(corr)
        logger.info('Mean value: %f; St dev: %f' % (mean, sd))
        z_df: pd.DataFrame = (corr - mean) / sd
    else:
        z_df = _z_scored_pvals(corr_df=corr, method=method)
    return z_df


def _get_sd(df: pd.DataFrame) -> float:
    ma = _mask_array(df)
    return ma.std()


def _get_mean(df: pd.DataFrame) -> float:
    ma = _mask_array(df)
    return ma.mean()


def _mask_array(df: pd.DataFrame) -> np.ma.MaskedArray:
    """Mask any NaN's in dataframe values"""
    logger.info('Masking DataFrame values')
    return np.ma.array(df.values, mask=np.isnan(df.values))


def get_mitocarta_info(mitocarta_file: str) -> Dict[str, str]:
    """Load mitocarta file and get mapping of gene to gene function

    Parameters
    ----------
    mitocarta_file : str
        Path as string to mitocarta file

    Returns
    -------
    Dict[str, str]

    """
    xls = pd.ExcelFile(mitocarta_file)
    # Sheet A is second sheet of MitoCarta 3.0 info
    sheet_a = xls.parse(xls.sheet_names[1])
    hgnc_expl = {}
    for eid, expl in zip(sheet_a.HumanGeneID.values,
                         sheet_a.Description.values):
        hgnc_tup = ns_id_xref(from_ns='EGID', from_id=eid, to_ns='HGNC')
        if hgnc_tup:
            hgnc_symb = ns_id_to_name(*hgnc_tup)
            if hgnc_symb:
                hgnc_expl[hgnc_symb] = expl
    return hgnc_expl


if __name__ == '__main__':
    parser = argparse.ArgumentParser('DepMap Data pre-processing')
    # Input dirs
    parser.add_argument('--crispr-raw', type=io.file_path('.csv'),
                        help='The raw CRISPR gene effect data. File name '
                             'is usually Achilles_gene_effect.csv')
    parser.add_argument('--rnai-raw', type=io.file_path('.csv'),
                        help='The raw RNAi gene dependency data. File name '
                             'is usually D2_combined_gene_dep_scores.csv')
    # Option to start from raw correlations matrix instead of running
    # correlation calculation directly
    parser.add_argument('--crispr-corr', type=io.file_path('.h5'),
                        help='The file containing an hdf compressed '
                             'correlation data frame of the crispr data. If '
                             'this file is provided, the raw crispr data is '
                             'ignored.')
    parser.add_argument('--rnai-corr', type=io.file_path('.h5'),
                        help='The file containing an hdf compressed '
                             'correlation data frame of the rnai data. If '
                             'this file is provided, the raw rnai data is '
                             'ignored.')
    # Output dirs
    parser.add_argument('--output-dir', '-o',
                        help='Optional. A directory where to put the '
                             'output of the script. If not provided the '
                             'output will go to the directory/ies where '
                             'the corresponding input data came from, '
                             'e.g. the output from the RNAi input will be '
                             'placed in the RNAi input directory. The '
                             'combined z-score matrix will be written to the '
                             'crispr input directory if this option is not '
                             'provided. Note that s3 URL are not allowed as '
                             'pd.DataFrame.to_hdf does not support urls or '
                             'buffers at the moment.')
    parser.add_argument('--random', '-r', type=int,
                        help='Optional. If specified, provide the size of '
                             'the final correlation matrix where the genes '
                             'are picked at random from the intersection of '
                             'genes from both the RNAI and CRISPR data sets.')
    parser.add_argument('--fname',
                        help='A file name for the output correlation '
                             'DataFrame.')
    parser.add_argument('--save-corr', action='store_true',
                        help='Also save the intermediate z-scored '
                             'correlations matrices from each of the input '
                             'data files.')

    args = parser.parse_args()

    # Run main script
    z_corr = run_corr_merge(crispr_raw=args.crispr_raw,
                            rnai_raw=args.rnai_raw,
                            crispr_corr=args.crispr_corr,
                            rnai_corr=args.rnai_corr,
                            output_dir=args.output_dir,
                            random_sampl=args.random,
                            remove_self_corr=False,
                            save_corr_files=args.save_corr)

    # Write merged correlations combined z score
    outdir: Path = Path(args.output_dir) if args.output_dir else (Path(
        args.crispr_corr).parent.as_posix() if args.crispr_corr else Path(
        args.crispr_raw).parent.as_posix())
    logger.info(f'Writing combined correlations to {outdir}')
    fname = args.fname if args.fname else 'combined_z_score.h5'
    outdir.mkdir(parents=True, exist_ok=True)
    z_corr.to_hdf(Path(outdir, fname).absolute().as_posix(), 'zsc')
