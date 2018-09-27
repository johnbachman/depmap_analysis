import os
import json
import logging
import pickle as pkl
import argparse as ap
from time import time
import depmap_network_functions as dnf

logger = logging.getLogger('jsonDump')


def _dump_it_to_json(fname, pyobj):
    with open(fname, 'w') as json_out:
        json.dump(pyobj, json_out)


parser = ap.ArgumentParser()
parser.add_argument('-p', '--pickle-file', required=True,
                    help='Pickle file containing a nested dict of explained '
                         'network: d[subj][obj] = {"directed": [(type, '
                         'hash, belief)] ...}')
parser.add_argument('-o', '--output-name',
                    help='Output base name of json files. With no input, the '
                         'default is "./expl_output/<hgnc gene>_is_subj.json".',
                    default='./expl_output/<hgnc gene>_is_subj.json')
args = parser.parse_args()
stamp = int(time())

if args.output_name.endswith('.json'):
    outbasename = args.output_name[:-5]  # Removes .json from basename
else:
    outbasename = args.output_name

if '/' not in outbasename:
    os.makedirs('./expl_output', exist_ok=True)
    outbasename = './expl_output/' + outbasename
logger.info('Using basename %s' % outbasename)

with open(args.pickle_file, 'rb') as pr:
    nest_dict = pkl.load(file=pr)

for subj, d in nest_dict.items():
    _dump_it_to_json(fname=outbasename + '%s_is_subj.json' % subj, pyobj=d)
