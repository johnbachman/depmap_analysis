import logging
import argparse
import pandas as pd
from pathlib import Path
from depmap_analysis.util.io_functions import pickle_open

logger = logging.getLogger(__name__)

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Corr script looper')
    parser.add_argument(
        '--z-corr', required=True,
        help='The correlation matrix that was used to create the explanations'
    )
    parser.add_argument(
        '--base-path', required=True,
        help='The path to the pickled explainer classes for each SD range.'
    )
    parser.add_argument(
        '--outdir',
        help='Path to an output directory. Will be created if it non '
             'existing. If not provided, outut files will be put in a '
             'directory called "output" in the "--base-path".'
    )
    parser.add_argument(
        '--dry', action='store_true',
        help='If flag is set, only print log messages of what will be done '
             'instead of running.'
    )

    args = parser.parse_args()
    base_path = Path(args.base_path)
    output_dir = Path(args.outdir) if args.outdir else \
        base_path.joinpath('output')
    dry = args.dry

    logger.info(f'Loading correlation file {args.z_corr}')
    if not dry:
        z_corr = pd.read_hdf(args.z_corr)
    else:
        z_corr = pd.DataFrame()

    for explainer_file in base_path.glob('*.pkl'):
        logger.info(f'Processing {explainer_file}')
        # Set outdir
        explainer_out = output_dir.joinpath(explainer_file.stem)
        if not dry:
            # Load pickle
            explainer = pickle_open(explainer_file)
            # Run stuff
            explainer.plot_corr_stats(outdir=explainer_out,
                                      z_corr=z_corr, show_plot=False)
            explainer.plot_dists(outdir=explainer_out,
                                 z_corr=None, show_plot=False)
        logger.info(f'Writing output to {explainer_out.as_posix()}/*.pdf')
