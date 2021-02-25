import pickle
import logging
import argparse
from math import floor
from typing import Union
from pathlib import Path
from itertools import count
from datetime import datetime

import boto3
import pandas as pd

from indra_db.util.s3_path import S3Path
from depmap_analysis.util.statistics import DepMapExplainer
from depmap_analysis.util.io_functions import file_opener, file_path, \
    is_dir_path, dump_it_to_pickle

logger = logging.getLogger(__name__)


def _exists(fpath: str) -> bool:
    if fpath.startswith('s3://'):
        return S3Path.from_string(fpath).exists(s3)
    else:
        return Path(fpath).is_file()


def _joinpath(fpath: Union[S3Path, Path], other: str) -> str:
    if isinstance(fpath, Path):
        return fpath.joinpath(other).absolute().as_posix()
    else:
        return fpath.to_string() + other


def _save(fpath: str, expl_inst: DepMapExplainer):
    if fpath.startswith('s3://'):
        try:
            _ = expl_inst.s3_location
        except AttributeError:
            # For backwards compatibility
            expl_inst.s3_location = fpath
        logger.info(f'Uploading to {expl_inst.s3_location}')
        if not dry:
            s3p = expl_inst.get_s3_path()
            s3p.upload(s3=s3, body=pickle.dumps(expl_inst))
    else:
        # Just dump to local pickle
        dump_it_to_pickle(fname=fpath, pyobj=expl_inst)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Corr script looper')
    parser.add_argument(
        '--z-corr', required=True, type=file_path('h5'),
        help='The correlation matrix that was used to create the explanations'
    )
    parser.add_argument(
        '--base-path', required=True, type=is_dir_path(),
        help='The path to the pickled explainer classes for each SD range.'
    )
    parser.add_argument(
        '--outdir',
        help='Path to an output directory. Will be created if it non '
             'existing. If not provided, output files will be put in a '
             'directory called "output" in the "--base-path". If path '
             'starts with "s3:" upload to s3 instead and must then have the '
             'form "s3:<bucket>/<sub_dir>" where <bucket> must be specified '
             'and <sub_dir> is optional and may contain subdirectories.'
    )
    parser.add_argument(
        '--dry', action='store_true',
        help='If flag is set, only print log messages of what will be done '
             'instead of running.'
    )
    parser.add_argument(
        '--max-proc', type=int,
        help='The maximum number of processes to run in the multiprocessing '
             'in get_corr_stats_mp. The number will automatically be floored '
             'if decimal. Default: multiprocessing.cpu_count()'
    )

    parser.add_argument(
        '--max-so-pairs', type=int, default=10000,
        help='The maximum number of correlation pairs to process. If the '
             'number of eligble pairs is larger than this number, a random '
             'sample of max_so_pairs_size is used. Default: 10 000. If the '
             'number of pairs to check is smaller than 1000, no sampling is '
             'done.'
    )

    parser.add_argument(
        '--mp-pairs', action='store_true',
        help='Perform multi processing when gathering the subj-obj axb '
             'pairs to process.'
    )

    parser.add_argument(
        '--single-proc', action='store_true',
        help='Run all scripts on a single process. This option is good when '
             'debugging or if the environment for some reason does not '
             'support multiprocessing. Default: False.'
    )

    args = parser.parse_args()
    base_path: str = args.base_path
    outdir: str = args.outdir
    single_proc: bool = args.single_proc

    s3 = boto3.client('s3')

    # Set input dir
    if base_path.startswith('s3://'):
        s3_base_path = S3Path.from_string(base_path)
        input_iter = \
            [s3p.to_string() for s3p in s3_base_path.list_objects(s3)
             if s3p.to_string().endswith('.pkl')]
    else:
        local_base_path = Path(base_path)
        input_iter = [f.absolute().as_posix()
                      for f in local_base_path.glob('*.pkl')]

    # Set output dir
    if outdir.startswith('s3://'):
        output_dir = S3Path.from_string(outdir)
    else:
        output_dir = Path(outdir)

    dry = args.dry
    if dry:
        logger.info('DRY RUN: no calculations are run but file errors are '
                    'raised for missing files')

    if args.max_proc:
        max_proc = floor(args.max_proc)
        if max_proc < 1:
            raise argparse.ArgumentTypeError(
                f'{max_proc} is not a valid positive integer'
            )
    else:
        max_proc = None

    logger.info(f'Loading correlation file {args.z_corr}')
    if not dry:
        z_corr = pd.read_hdf(args.z_corr)
    else:
        if not Path(args.z_corr).is_file():
            raise FileNotFoundError(f'{args.z_corr} was not found')
        z_corr = pd.DataFrame()

    # Create a global indexer to separate each figure
    indexer = count(0)
    processed_explainers = []
    for explainer_file in input_iter:
        logger.info(
            f'> > > > '
            f'Processing {explainer_file} '
            f'{datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")} (UTC)'
            f' < < < <'
        )
        explainer_out = _joinpath(
            output_dir, str(explainer_file).split('/')[-1].split('.')[0]
        )
        logger.info(f'Saving output to {explainer_out}')
        if not dry:
            # Load pickle
            explainer = file_opener(explainer_file)
            try:
                assert isinstance(explainer, DepMapExplainer)
            except AssertionError:
                logger.warning(f'File {explainer_file} is not '
                               f'DepMapExplainer, skipping...')
                continue
            # Backwards compatibility: check if s3_location attribute
            # exists, otherwise set it and then re-upload
            try:
                _ = explainer.s3_location
            except AttributeError:
                _save(fpath=explainer_file, expl_inst=explainer)

            # Run stuff
            explainer.plot_corr_stats(outdir=explainer_out,
                                      z_corr=z_corr,
                                      show_plot=False,
                                      max_proc=max_proc,
                                      index_counter=indexer,
                                      max_so_pairs_size=args.max_so_pairs,
                                      mp_pairs=args.mp_pairs,
                                      run_linear=single_proc)

            explainer.plot_dists(outdir=explainer_out,
                                 z_corr=None,
                                 show_plot=False,
                                 index_counter=indexer,
                                 max_so_pairs_size=args.max_so_pairs,
                                 mp_pairs=args.mp_pairs,
                                 run_linear=single_proc)
            explainer.plot_interesting(outdir=explainer_out,
                                       z_corr=None,
                                       show_plot=False,
                                       index_counter=indexer,
                                       max_so_pairs_size=args.max_so_pairs,
                                       mp_pairs=args.mp_pairs,
                                       run_linear=single_proc)
        else:
            if not _exists(explainer_file):
                raise FileNotFoundError(f'{explainer_file} does not exist')
        logger.info(f'Writing output to {explainer_out}/*.pdf')
