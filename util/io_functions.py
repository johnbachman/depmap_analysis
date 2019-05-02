import csv
import json
import pickle
import logging

logger = logging.getLogger('dnf utils')


def _dump_it_to_pickle(fname, pyobj):
    logger.info('Dumping to pickle file %s' % fname)
    with open(fname, 'wb') as po:
        pickle.dump(obj=pyobj, file=po)
    logger.info('Finished dumping to pickle')


def _dump_it_to_json(fname, pyobj):
    logger.info('Dumping to json file %s' % fname)
    with open(fname, 'w') as json_out:
        json.dump(pyobj, json_out)
    logger.info('Finished dumping to pickle')


def _dump_it_to_csv(fname, pyobj, separator=',', header=None):
    logger.info('Dumping to csv file %s' % fname)
    if header:
        logger.info('Writing csv header')
        with open(fname, 'w') as fo:
            fo.write(','.join(header)+'\n')
    with open(fname, 'a', newline='') as csvf:
        wrtr = csv.writer(csvf, delimiter=separator)
        wrtr.writerows(pyobj)
    logger.info('Finished dumping to csv')


def _pickle_open(file_path_to_pickle):
    logger.info('Loading pickle file %s' % file_path_to_pickle)
    with open(file_path_to_pickle, 'rb') as pi:
        pkl = pickle.load(file=pi)
    logger.info('Finished loading pickle file')
    return pkl


def _json_open(file_path_to_json):
    logger.info('Loading json file %s' % file_path_to_json)
    with open(file_path_to_json, 'r') as jo:
        js = json.load(fp=jo)
    logger.info('Finished loading json file')
    return js
