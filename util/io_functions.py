import csv
import json
import pickle


def _dump_it_to_pickle(fname, pyobj):
    with open(fname, 'wb') as po:
        pickle.dump(obj=pyobj, file=po)


def _dump_it_to_json(fname, pyobj):
    with open(fname, 'w') as json_out:
        json.dump(pyobj, json_out)


def _dump_it_to_csv(fname, pyobj, separator=',', header=None):
    if header:
        with open(fname, 'w') as fo:
            fo.write(','.join(header)+'\n')
    with open(fname, 'a', newline='') as csvf:
        wrtr = csv.writer(csvf, delimiter=separator)
        wrtr.writerows(pyobj)


def _pickle_open(file_path_to_pickle):
    with open(file_path_to_pickle, 'rb') as pi:
        return pickle.load(file=pi)


def _json_open(file_path_to_json):
    with open(file_path_to_json, 'r') as jo:
        return json.load(fp=jo)
