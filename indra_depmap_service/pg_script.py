import pickle
import logging
from collections import Counter
from paths_graph import PathsGraph, get_reachable_sets, \
                        CombinedPathsGraph


print("Loading network")
with open('_cache/nx_dir_graph_db_dump_20190417.pkl', 'rb') as f:
    g = pickle.load(f)

print("Done loading network")

source = 'NCKAP1'
target = 'TEAD1'
max_depth = 5
cur_length = 3
num_samples = 20000

print("Getting reachable sets")
fwd_reach, back_reach = get_reachable_sets(g, source, target, max_depth,
                               signed=False)

print("Building PG")
pg_list = []
for cur_length in range(1, max_depth+1):
    print("Building paths graph for length %d" % cur_length)
    pg = PathsGraph.from_graph(g, source, target, cur_length,
                                fwd_reach,  back_reach, signed=False,
                                target_polarity=0)
    pg_list.append(pg)

print("Building combined paths graph")
cpg = CombinedPathsGraph(pg_list)

print("Sampling %d paths" % num_samples)
paths = cpg.sample_cf_paths(num_samples)
path_ctr = Counter(paths)
path_ctr = sorted([(k, v) for k, v in path_ctr.items()],
                  key=lambda x: x[1], reverse=True)

