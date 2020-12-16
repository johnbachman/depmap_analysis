import re
import csv
import json
import pandas as pd
import pickle
import logging
import platform
from io import StringIO
from os import path, stat
from typing import Iterable
from pathlib import Path
from datetime import datetime
from argparse import ArgumentError
from functools import wraps
from itertools import repeat, takewhile

import numpy as np

logger = logging.getLogger(__name__)

DT_YmdHMS_ = '%Y-%m-%d-%H-%M-%S'
DT_YmdHMS = '%Y%m%d%H%M%S'
DT_Ymd = '%Y%m%d'
RE_YmdHMS_ = r'\d{4}\-\d{2}\-\d{2}\-\d{2}\-\d{2}\-\d{2}'
RE_YYYYMMDD = r'\d{8}'


def file_opener(fname: str, **kwargs) -> object:
    """Open file based on file extension

    kwargs can be provided and are used for s3_file_opener (s3 calls) and
    pd.read_csv() (local files)

    Parameters
    ----------
    fname : str
        The filename. If an s3 url, load object directly from s3.

    Returns
    -------
    object
        Object stored in file fname
    """
    if fname.startswith('s3://'):
        return s3_file_opener(fname, **kwargs)
    if fname.endswith('pkl'):
        return pickle_open(fname)
    elif fname.endswith('json'):
        return json_open(fname)
    elif fname.endswith(('csv', 'tsv')):
        return pd.read_csv(fname, **kwargs)  # Can provide e.g. index_col=0
    else:
        raise ValueError(f'Unknown file extension for file {fname}')


def s3_file_opener(s3_url: str, unsigned: bool = False, **kwargs) -> object:
    """Open a file from s3 given a standard s3-path

    kwargs are only relevant for csv/tsv files and are used for pd.read_csv()

    Parameters
    ----------
    s3_url : str
        S3 url of the format 's3://<bucket>/<key>'. The key is assumed to
        also contain a file ending
    unsigned : bool
        If True, perform S3 calls unsigned. Default: False

    Returns
    -------
    object
    """
    from indra_db.util.s3_path import S3Path
    from .aws import load_pickle_from_s3, read_json_from_s3, get_s3_client
    logger.info(f'Loading {s3_url} from s3')
    s3_path = S3Path.from_string(s3_url)
    s3 = get_s3_client(unsigned=unsigned)
    bucket, key = s3_path.bucket, s3_path.key
    if key.endswith('.json'):
        return read_json_from_s3(s3=s3, key=key, bucket=bucket)
    elif key.endswith('.pkl'):
        return load_pickle_from_s3(s3=s3, key=key, bucket=bucket)
    elif key.endswith(('.csv', '.tsv')):
        fileio = S3Path.from_string(s3_url).get(s3=s3)
        csv_str = fileio['Body'].read().decode('utf-8')
        raw_file = StringIO(csv_str)
        return pd.read_csv(raw_file, **kwargs)
    else:
        return S3Path.from_string(s3_url).get(s3=s3)


def file_dump_wrapper(f):
    """Wrapper for any function that dumps a python object

    The wrapped functions must contain at least the args "fname" and
    "pyobj", the kwarg "overwrite" to flag for forcing overwriting of the
    file is strongly encouraged

    Parameters
    ----------
    f : function

    Returns
    -------
    decorator
    """

    @wraps(f)
    def decorator(*args, **kwargs):
        """Decorates file dumping functions

        Parameters
        ----------
        *args
            Must contain at least two arguments: a string for a filename and
            the python object to write to a file
        **kwargs
            Can contain the flag "overwrite" that is used to raise an error
            if the provided filepath already exists

        Returns
        -------
        function
        """
        overwrite = kwargs.get('overwrite', False)
        if len(args) == 0:
            try:
                fname = kwargs['fname']
            except KeyError:
                raise TypeError("f() missing a required positional argument: "
                                "'fname'")
            file = Path(fname)
        else:
            fname = args[0]
            file = Path(fname)

        # 1. Check if file exists if overwrite is False
        if not overwrite and file.is_file():
            raise FileExistsError(f'File {str(file)} already exists! Use '
                                  f'overwrite=True to overwrite current file.')

        # 2. Create parent directories if not present
        if not file.parent.is_dir():
            file.parent.mkdir(parents=True)

        return f(*args, **kwargs)
    return decorator


@file_dump_wrapper
def dump_it_to_pickle(fname, pyobj, overwrite=False):
    """Save pyobj to fname as pickle"""
    logger.info('Dumping to pickle file %s' % fname)
    with Path(fname).open('wb') as po:
        pickle.dump(obj=pyobj, file=po)
    logger.info('Finished dumping to pickle')


@file_dump_wrapper
def dump_it_to_json(fname, pyobj, overwrite=False):
    """Save pyobj to fname as json"""
    logger.info('Dumping to json file %s' % fname)
    with Path(fname).open('w') as json_out:
        json.dump(pyobj, json_out)
    logger.info('Finished dumping to json')


@file_dump_wrapper
def dump_it_to_csv(fname, pyobj, separator=',', header=None, overwrite=False):
    """Save pyobj to fname as csv file"""
    logger.info('Dumping to csv file %s' % fname)
    file = Path(fname)
    if header:
        logger.info('Writing csv header')
        with file.open('w') as fo:
            fo.write(','.join(header)+'\n')
    with file.open('a', newline='') as csvf:
        wrtr = csv.writer(csvf, delimiter=separator)
        wrtr.writerows(pyobj)
    logger.info('Finished dumping to csv')


def pickle_open(fname):
    """Open pickle fname and return the contained object"""
    logger.info('Loading pickle file %s' % fname)
    file = Path(fname)
    with file.open('rb') as pi:
        pkl = pickle.load(file=pi)
    logger.info('Finished loading pickle file')
    return pkl


def json_open(fname):
    """Open json fname and return the object"""
    logger.info('Loading json file %s' % fname)
    file = Path(fname)
    with file.open('r') as jo:
        js = json.load(fp=jo)
    logger.info('Finished loading json file')
    return js


def read_gene_set_file(gf, data):
    """Read HGNC symbols from df and match with data"""
    gset = []
    try:
        # Works if string is returned: we assume this is when we only have
        # HGNC symbols
        data.columns[0].split()
        dset = set(data.columns)
    except AttributeError:
        # multi index
        dset = set([t[0] for t in data.columns])
    file = Path(gf)
    with file.open('rt') as f:
        for g in f.readlines():
            gn = g.upper().strip()
            if gn in dset:
                gset.append(gn)
    return gset


def rawincount(filename):
    """Count lines in filename

    filename: str
        Path to file to count lines in

    Returns
    -------
    line_count: int
        The number of lines in the file 'filename'
    """
    f = open(filename, 'rb')
    bufgen = takewhile(lambda x: x, (f.read(1024*1024) for _ in repeat(None)))
    return sum(buf.count(b'\n') for buf in bufgen)


def map2index(start, binsize, value):
    """Get bin index in symmetric histogram for value given start, binsize"""
    offset = int(abs(start//binsize))
    return offset + int(float(value) // binsize)


def histogram_for_large_files(fpath, number_of_bins, binsize, first):
    """Returns a histogram for very large files

    fpath: str(filename)
        filepath to file with data to be binned
    number_of_bins: int
        the number fo bins to use
    binsize: float
        the size of bins
    first: float
        The left most (min(x)) edge of the bin edges

    Returns
    -------
    home_brewed_histo: np.array
        A histrogram of the data in fpath according to number of bins,
        binsize and first.
    """
    home_brewed_histo = np.zeros(number_of_bins, dtype=int)
    file = Path(fpath)
    with file.open('r') as fo:
        for line in fo:
            flt = line.strip()
            home_brewed_histo[map2index(start=first, binsize=binsize,
                                        value=flt)] += 1
    return home_brewed_histo


def _manually_add_to_histo(hist, start, binsize, value):
    hist[map2index(start, binsize, value)] += 1


def allowed_types(types: Iterable):
    """Types is a set of strings with names of the allowed types"""
    def types_check(_type: str) -> str:
        """Check the input type

        Parameters
        ----------
        _type : str
            The input type

        Returns
        -------
        str
            Returns the lowercase of the input string representing the type
        """
        if _type.lower() not in types:
            raise ArgumentError(f'Provided graph type {_type} not allowed. '
                                f'Have to be one of {types}')
        return _type.lower()
    return types_check


def file_path(file_ending: str = None):
    """Checks if file at provided path exists"""
    def check_path(fpath: str):
        if fpath.startswith('s3://'):
            if file_ending and not fpath.endswith(file_ending):
                raise ValueError(f'Unrecognized file type '
                                 f'{fpath.split("/")[-1]}')
            from indra_db.util.s3_path import S3Path
            from .aws import get_s3_client
            if not S3Path.from_string(fpath).exists(s3=get_s3_client(False)):
                raise ValueError(f'File {fpath} does not exist')
            return fpath
        p = Path(fpath)
        if not p.is_file():
            raise ArgumentError(f'File {fpath} does not exist')
        if file_ending and not p.name.endswith(file_ending):
            raise ArgumentError(f'Unrecognized file type '
                                f'{p.name.split(".")[-1]}')
        return fpath
    return check_path


def is_dir_path():
    """Checks if provided path exists"""
    def is_dir(path):
        dp = Path(path)
        if not dp.is_dir():
            raise ArgumentError(f'Path {path} does not exist')
        return path
    return is_dir


def todays_date(dt_fmt=DT_Ymd):
    return datetime.now().strftime(dt_fmt)


def get_earliest_date(file):
    """Returns creation or modification timestamp of file

    # todo: Add s3 option

    Parameters
    ----------
    file : str
        File path

    Returns
    -------
    float
        Timestamp in seconds with microseconds as a float
    """
    # https://stackoverflow.com/questions/237079/
    # how-to-get-file-creation-modification-date-times-in-python
    if platform.system().lower() == 'windows':
        return path.getctime(file)
    else:
        st = stat(file)
        try:
            return st.st_birthtime
        except AttributeError:
            return st.st_mtime


def get_date_from_str(date_str, dt_format):
    """Returns a datetime object from a datestring of format FORMAT"""
    return datetime.strptime(date_str, dt_format)


def strip_out_date(keystring, re_format):
    """Strips out datestring of format re_format from a keystring"""
    try:
        return re.search(re_format, keystring).group()
    except AttributeError:
        logger.warning('Can\'t parse string "%s" for date using regex pattern '
                       'r"%s"' % (keystring, re_format))
        return None
