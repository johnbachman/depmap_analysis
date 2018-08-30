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
                    help='Pickle file containing a nested dict '
                         'd[subj][obj][[type, hash], ...]')
parser.add_argument('-o', '--output-name',
                    help='Output base name of json files. With no input, the '
                         'default is "indra_db_<time stamp>".',
                    default='./output/indra_db_{}_'.format(int(time())))

args = parser.parse_args()
stamp = int(time())

if not args.output_name.endswith('.json'):
    outbasename = args.output_name[-5:]  # Removes .json from basename
else:
    outbasename = args.output_name

os.makedirs('./output', exist_ok=True)
outbasename = './output/'+outbasename
logger.info('Using baseame %s' % outbasename)

with open(args.pickle_file, 'rb') as pr:
    nest_dict = pkl.load(file=pr)

# Create nested dict
nest_dict_out = nest_dict()

# Convert hash to strings
for s, inner_d in nest_dict.items():
    for o in inner_d:
        type_hash_list = inner_d[o]
        t_h_list_out = []

        for tp, hsh in type_hash_list:
            hash_string = str(hsh)
            t_h_list_out.append((tp, hash_string))

        nest_dict_out[s][o] = t_h_list_out

# Output:
# 1. subj list
# 2. For each subj: dict[obj] -> [[type, hash], ...]
# 3. obj list
# 4. For each obj: (reverse lookup) dict[subj] -> [[type, hash], ...]

# 1. subj list
subj_list = list(set(nest_dict_out.keys()))
_dump_it_to_json(outbasename+'_subjects.json', subj_list)

obj_set = set()

rev_dict = {}

for subj, d in nest_dict_out.items():
    obj_set.update(set(d.keys()))
    # 2. dump each subj dict as json
    _dump_it_to_json(outbasename+'_%s_is_subj.json' % subj, d)

    # Build reverse dicts:
    for obj, entry in d.items():
        rev_dict.setdefault(obj, {}).update({subj: entry})

# 3. obj list
_dump_it_to_json(outbasename+'_objects.json', list(obj_set))

for obj, d in rev_dict.items():
    # 4. dump the reverse/obj dicts
    _dump_it_to_json(outbasename+'_%s_is_obj.json' % obj, d)
