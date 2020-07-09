import logging
import argparse
from math import floor
from pathlib import Path
from itertools import count
from datetime import datetime

import pandas as pd

from depmap_analysis.util.io_functions import pickle_open, file_path, \
    is_dir_path

logger = logging.getLogger(__name__)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Corr script looper')
    parser.add_argument(
        '--z-corr', required=True,
        help='The correlation matrix that was used to create the explanations'
    )
    parser.add_argument(
        '--reactome', type=file_path('pkl'),
        help='A tuple or list of dicts. The first dict is expected to'
             'contain mappings from UP IDs of genes to Reactome pathway IDs.'
             'The second dict is expected to contain the reverse mapping '
             '(i.e Reactome IDs to UP IDs). The third dict is expected to '
             'contain mappings from Reactome IDs to their descriptions.'
    )
    parser.add_argument(
        '--base-path', required=True, type=is_dir_path(),
        help='The path to the pickled explainer classes for each SD range.'
    )
    parser.add_argument(
        '--outdir',
        help='Path to an output directory. Will be created if it non '
             'existing. If not provided, outut files will be put in a '
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
        '--max-so-pairs', type=int,
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

    args = parser.parse_args()
    base_path = Path(args.base_path)
    output_dir = Path(args.outdir) if args.outdir else \
        base_path.joinpath('output')
    dry = args.dry

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
        if args.reactome:
            logger.info(f'Loading reactome file {args.reactome}')
            reactome = pickle_open(args.reactome)
        else:
            reactome = None
    else:
        if not Path(args.z_corr).is_file():
            raise FileNotFoundError(f'{args.z_corr} was not found')
        z_corr = pd.DataFrame()
        reactome = None
    # Create a global indexer to separate each figure
    indexer = count(0)
    for explainer_file in base_path.glob('*.pkl'):
        print(f'> > > > Processing {explainer_file} '
              f'{datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")} (UTC) < < '
              f'< <')
        explainer_out = output_dir.joinpath(explainer_file.stem).as_posix()
        if not dry:
            # Load pickle
            explainer = pickle_open(explainer_file)
            # Run stuff
            explainer.plot_corr_stats(outdir=explainer_out,
                                      z_corr=z_corr,
                                      reactome=reactome,
                                      show_plot=False,
                                      max_proc=max_proc,
                                      index_counter=indexer,
                                      max_so_pairs_size=args.max_so_pairs,
                                      mp_pairs=args.mp_pairs)

            explainer.plot_dists(outdir=explainer_out,
                                 z_corr=None, show_plot=False,
                                 index_counter=indexer,
                                 max_so_pairs_size=args.max_so_pairs,
                                 mp_pairs=args.mp_pairs)
        else:
            if not Path(explainer_file).is_file():
                raise FileNotFoundError(f'{explainer_file} does not exist')
        logger.info(f'Writing output to {explainer_out}/*.pdf')
