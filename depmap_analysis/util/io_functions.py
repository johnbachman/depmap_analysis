import csv
import json
import pickle
import logging
import numpy as np
from itertools import repeat, takewhile

logger = logging.getLogger('dnf utils')


def dump_it_to_pickle(fname, pyobj):
    """Save pyobj to fname as pickle"""
    logger.info('Dumping to pickle file %s' % fname)
    with open(fname, 'wb') as po:
        pickle.dump(obj=pyobj, file=po)
    logger.info('Finished dumping to pickle')


def dump_it_to_json(fname, pyobj):
    """Save pyobj to fname as json"""
    logger.info('Dumping to json file %s' % fname)
    with open(fname, 'w') as json_out:
        json.dump(pyobj, json_out)
    logger.info('Finished dumping to pickle')


def dump_it_to_csv(fname, pyobj, separator=',', header=None):
    """Save pyobj to fname as csv file"""
    logger.info('Dumping to csv file %s' % fname)
    if header:
        logger.info('Writing csv header')
        with open(fname, 'w') as fo:
            fo.write(','.join(header)+'\n')
    with open(fname, 'a', newline='') as csvf:
        wrtr = csv.writer(csvf, delimiter=separator)
        wrtr.writerows(pyobj)
    logger.info('Finished dumping to csv')


def pickle_open(fname):
    """Open pickle fname and return the contianed object"""
    logger.info('Loading pickle file %s' % fname)
    with open(fname, 'rb') as pi:
        pkl = pickle.load(file=pi)
    logger.info('Finished loading pickle file')
    return pkl


def json_open(fname):
    """Open json fname and return the object"""
    logger.info('Loading json file %s' % fname)
    with open(fname, 'r') as jo:
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

    with open(gf, 'rt') as f:
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
    with open(file=fpath) as fo:
        for line in fo:
            flt = line.strip()
            home_brewed_histo[map2index(start=first, binsize=binsize,
                                        value=flt)] += 1
    return home_brewed_histo


def _manually_add_to_histo(hist, start, binsize, value):
    hist[map2index(start, binsize, value)] += 1
