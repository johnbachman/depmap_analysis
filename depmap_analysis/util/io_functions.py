import csv
import json
import pickle
import logging
from argparse import ArgumentError
from pathlib import Path
from functools import wraps
from itertools import repeat, takewhile

import numpy as np

logger = logging.getLogger(__name__)


def file_opener(fname):
    """Open file based on file extension

    Parameters
    ----------
    fname : str
        The filename

    Returns
    -------
    object
        Object stored in file fname
    """
    if fname.endswith('pkl'):
        return pickle_open(fname)
    elif fname.endswith('json'):
        return json_open(fname)
    else:
        raise ValueError(f'Unknown file extension for file {fname}')


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
    logger.info('Finished dumping to pickle')


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


def graph_types(types):
    """Types is a set of strings with names of the allowed graph types"""
    def types_check(_type):
        """Check the input graph type

        Parameters
        ----------
        _type : str
            The input graph type

        Returns
        -------
        str
            Returns the lowercase of the input string representing the graph
            type
        """
        if _type.lower() not in types:
            raise ArgumentError(f'Provided graph type {_type} not allowed. '
                                f'Have to be one of {types}')
        return _type.lower()
    return types_check


def file_path(file_ending=None):
    """Checks if file at provided path exists"""
    def check_path(fpath):
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
