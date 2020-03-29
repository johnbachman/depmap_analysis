"""This script should be run post DepMap script

It aggregates the explanation statisticcs for all the SD ranges after
running 'run_script.sh'
"""

import re
import logging
import argparse
from os import listdir, path
from datetime import datetime

import pandas as pd


patt1 = re.compile('(.*): ([0-9]+)\n?')
patt2 = re.compile('(.*) ([0-9]+)\n?')
dirpatt = re.compile('([0-9]+_([0-9]*)(sd|SD|sD|Sd)|rnd|RND)$')
sumfilepatt = re.compile('(.*)_script_summary.txt')

logger = logging.getLogger(__file__)
logging.basicConfig(level=logging.INFO)
logger.info('INFO working')


def process_summary_file(logfile):
    """Read log file and extract at least
    ['total_expl_excl_sr', 'common_parent', 'mitochondrial', 'direct',
    'axb_excl_sr']

    Parameters
    ----------
    logfile : str
        The file path to a log file

    Returns
    -------
    dict
        A dict of {metric: count}
    """
    metric_dict = {}
    with open(logfile, 'r') as f:
        line = f.readline()
        while line:
            if line.startswith('> '):
                metric_dict.update(_match_line(line))
            line = f.readline()
    return metric_dict


def _match_line(line):
    """Match the counts of the following lines:

    > Total number of correlation pairs checked: 3420
    > Total correlations unexplained: 1379
    > Total correlations explained: 2041
    > Total correlations explained, ignoring shared regulator: 1782
    > Total correlations explained, excluding shared regulator \
        (total - shared only): 1782
    >    0 correlations have an explanation involving a common parent
    >    33 gene pairs were considered explained as part of the\
            "explained set"
    >    277 explanations involving direct connection or complex
    >    277 correlations have a directed explanation involving an \
            intermediate node (A->X->B/A<-X<-B)
    >    1780 correlations have an explanation involving an intermediate \
            node excluding shared regulators
    >    1947 correlations have an explanation involving a shared regulator \
            (A<-X->B)
    >    259 correlations have shared regulator as only explanation

    Parameters
    ----------
    line : str
        A string from the logfile containing the count statistics of one of
        the ineresting metrics

    Returns
    -------
    {metric: count}
        A dict with matching metric name and its count as key: value
    """

    if 'Total number of correlation pairs checked' in line:
        metric = 'total_corr'
        count = int(patt1.search(line).groups()[1])

    elif 'Total correlations unexplained' in line:
        metric = 'unexplained'
        count = int(patt1.search(line).groups()[1])

    elif 'Total correlations explained:' in line:
        metric = 'total_expl'
        count = int(patt1.search(line).groups()[1])

    elif 'Total correlations explained, ignoring shared regulator:' in line:
        metric = 'total_expl_excl_sr'
        count = int(patt1.search(line).groups()[1])

    elif 'Total correlations explained, excluding shared regulator' in line:
        metric = 'total_expl_excl_sr_assert'
        count = int(patt2.search(line).groups()[1])

    elif 'correlations have an explanation involving a common' in line:
        metric = 'common_parent'
        count = int(patt2.search(line).groups()[1])

    elif 'as part of the "explained set"' in line:
        metric = 'mitochondrial'
        count = int(patt2.search(line).groups()[1])

    elif 'involving direct connection or complex' in line:
        metric = 'direct'
        count = int(patt2.search(line).groups()[1])

    elif 'directed explanation involving an intermediate node ' \
         '(A->X->B/A<-X<-B)' in line:
        metric = 'axb_dir'
        count = int(patt2.search(line).groups()[1])

    elif 'intermediate node excluding shared regulators' in line:
        metric = 'axb_excl_sr'
        count = int(patt2.search(line).groups()[1])

    elif 'involving a shared regulator (A<-X->B)' in line:
        metric = 'sr_all'
        count = int(patt2.search(line).groups()[1])

    elif 'shared regulator as only explanation' in line:
        metric = 'sr_only'
        count = int(patt2.search(line).groups()[1])

    else:
        metric = None
        count = None
    return {metric: [count]} if metric else {}


if __name__ == '__main__':
    parser = argparse.ArgumentParser('log script')
    parser.add_argument('--logdir', required=True,
                        help='The parent directory of all the logs. In this '
                             'directory, there should be sub directories '
                             'with names matching "ll_(ul)sd" or '
                             '"ll_(ul)SD" in which there should be files '
                             'matching the name "*_script_summary.txt"')
    parser.add_argument('--tag', default=datetime.now().strftime('%Y%m%d'),
                        help='Append a tag to the output filename of the '
                             'form "filename_tag.csv".')

    args = parser.parse_args()
    logdir = args.logdir
    if not logdir:
        raise ValueError('Must specify log directory: --logdir <logdir>')
    elif not path.isdir(logdir):
        raise ValueError('--logdir must specify a valid path: --logdir '
                         '<logdir>')

    # Initialize df
    df = pd.DataFrame()

    # Get subdirectories
    log_dirs = [d for d in listdir(logdir) if dirpatt.search(d)]
    for d in log_dirs:
        # There should only be one summary file
        files = [f for f in listdir(path.join(logdir, d)) if
                 sumfilepatt.search(f)]
        if len(files) == 0:
            logger.warning(f'Directory {d} does not contain any summary '
                           f'file. Skipping...')
            continue
        else:
            if len(files) > 1:
                logger.warning('Found more than one summary file, picking '
                               'first one.')
            summary_file = path.join(logdir, d, files[0])

        if path.isfile(summary_file):
            logger.info(f'Processing {summary_file}')
            metrics = process_summary_file(summary_file)
            metrics['range'] = [d.upper()]
            df = df.append(other=pd.DataFrame(data=metrics))
        else:
            logger.warning(f'{summary_file} is not a file!')
            continue
    logger.info(f'Writing dataframe to csv file at {logdir}')
    tag = '_' + args.tag
    df.to_csv(path.join(logdir, f'expl_{logdir.split("/")[-1]}{tag}.csv'),
              sep=',', index=False)
