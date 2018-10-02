import pylab
import numpy as np
import scipy as sp
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
from collections import defaultdict

# run in ~/repos/depmap_analysis


def nest_dict():
    """Returns a nested dictionary of arbitrary depth

    Returns
    -------
    defaultdict(nest_dict)
    """
    return defaultdict(nest_dict)


def _map2index(start, binsize, value):
    offset = int(abs(start//binsize))
    return offset + int(float(value.strip()) // binsize)


def _histogram_for_large_files(fpath, number_of_bins, binsize, first):
    home_brewed_histo = np.zeros(number_of_bins, dtype=int)
    with open(file=fpath) as fo:
        for line in fo:
            flt = line.strip()
            home_brewed_histo[_map2index(first, binsize, flt)] += 1
    return home_brewed_histo


def _my_gauss(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


def _pickle_open(file_path_to_pickle):
    with open(file_path_to_pickle, 'rb') as pi:
        return pkl.load(file=pi)


def _entry_exist(nested_dict, outer_key, inner_key):
    if nested_dict.get(outer_key) and nested_dict.get(outer_key).get(inner_key):
        return True
    else:
        return False


read_correlation_file = input('Do you want to read the full correlation file? ')

nested_explained_dict = _pickle_open(
    'Q3_depmap_INDRA_db_20180730_explained_nest_dict_belief.pkl')

corr_lookup = _pickle_open('input_data/depmap/correlation_lookup_ll02.pkl')

direct = 0
explained = 0
only_x_upstream = 0
has_intermeditate = 0
intermediates_only = 0
AB_correlation_list = []
only_upstream_pairs = []
x_is_downstream_count = []
x_is_intermediate_count = []
AB_correlation_list_direct = []
AB_correlation_list_intermediate = []
AXB_correlation_list_intermediate_ll02 = []
x_count_df = pd.DataFrame(columns=['subj', 'obj', 'x_type', 'count'])
x_count_df_list = []
for subj in nested_explained_dict:
    for obj in nested_explained_dict[subj]:
        directed = nested_explained_dict[subj][obj]['directed']
        undirected = nested_explained_dict[subj][obj]['undirected']
        x_is_downstream = nested_explained_dict[subj][obj]['x_is_downstream']
        x_is_upstream = nested_explained_dict[subj][obj]['x_is_upstream']
        x_is_intermediary = nested_explained_dict[subj][obj][
            'x_is_intermediary']
        # Check if all but x_is_upstream are empty
        if not directed and not undirected and not x_is_downstream and not \
                x_is_intermediary:
            if x_is_upstream:
                only_x_upstream += 1
                only_upstream_pairs.append((subj, obj, x_is_upstream))
            else:
                continue
        else:
            correlation = nested_explained_dict[subj][obj]['correlation']
            AB_correlation_list.append(correlation)
            explained += 1
            if directed or undirected:
                AB_correlation_list_direct.append(correlation)
                direct += 1
            if x_is_downstream or x_is_intermediary:
                if x_is_downstream:
                    for x, s in x_is_downstream:
                        # get A-X correlation
                        if _entry_exist(corr_lookup, subj, x):
                            AXB_correlation_list_intermediate_ll02.append(
                                corr_lookup[subj][x])
                        elif _entry_exist(corr_lookup, x, subj):
                            AXB_correlation_list_intermediate_ll02.append(
                                corr_lookup[x][subj])

                        # get X-B correlation
                        if _entry_exist(corr_lookup, x, obj):
                            AXB_correlation_list_intermediate_ll02.append(
                                corr_lookup[x][obj])
                        elif _entry_exist(corr_lookup, obj, x):
                            AXB_correlation_list_intermediate_ll02.append(
                                corr_lookup[obj][x])

                    x_ds_count = len(x_is_downstream)
                    x_is_downstream_count.append(x_ds_count)
                    x_count_df_list.append({'subj': subj,
                                            'obj': obj,
                                            'x_type': 'x_is_downstream',
                                            'count': x_ds_count})
                    if x_ds_count > 500:
                        print('%s and %s share more than 500 targets' % (
                            subj, obj))
                if x_is_intermediary:
                    for x, s in x_is_intermediary:
                        # get A-X correlation
                        if _entry_exist(corr_lookup, subj, x):
                            AXB_correlation_list_intermediate_ll02.append(
                                corr_lookup[subj][x])
                        elif _entry_exist(corr_lookup, x, subj):
                            AXB_correlation_list_intermediate_ll02.append(
                                corr_lookup[x][subj])

                        # get X-B correlation
                        if _entry_exist(corr_lookup, x, obj):
                            AXB_correlation_list_intermediate_ll02.append(
                                corr_lookup[x][obj])
                        elif _entry_exist(corr_lookup, obj, x):
                            AXB_correlation_list_intermediate_ll02.append(
                                corr_lookup[obj][x])
                    x_im_count = len(x_is_intermediary)
                    x_is_intermediate_count.append(x_im_count)
                    x_count_df_list.append({'subj': subj,
                                            'obj': obj,
                                            'x_type': 'x_is_intermediary',
                                            'count': x_im_count})
                    if x_im_count > 500:
                        print('%s and %s have more than 500 intermediaries' % (
                            subj, obj))
                AB_correlation_list_intermediate.append(correlation)
                has_intermeditate += 1
                if not directed and not undirected:
                    intermediates_only += 1
print('Total number of explanations excluding those with only shared '
      'regulator: %i' % explained)
print('Total number of cases where shared regulator is only explanation: %i'
      % only_x_upstream)
print('Total number of explanations with direct link, excluding those with '
      'only shared regulator: %i' % direct)
print('Total number of explanations with intermediates, excluding those with '
      'only shared regulator: %i' % has_intermeditate)
print('Total number of explanations with only intermediates as explanation, '
      'excluding those with only shared regulator: %i' % intermediates_only)

correlation_array = np.array(AB_correlation_list)
correlation_array_intermediate = np.array(AB_correlation_list_intermediate)
correlation_array_direct = np.array(AB_correlation_list_direct)
correlation_array_AXB = np.array(AXB_correlation_list_intermediate_ll02)

max_of_x = max(max(x_is_intermediate_count), max(x_is_downstream_count))
x_is_intermediate_count_hist, x_count_binedges = \
    np.histogram(a=x_is_intermediate_count,
        bins=max_of_x,
        range=(0.5, max_of_x+0.5))
x_is_downstream_count_hist, x_count_binedges = \
    np.histogram(a=x_is_downstream_count,
        bins=max_of_x,
        range=(0.5, max_of_x+0.5))

expl_correlation_hist, binedges = np.histogram(
    a=correlation_array, bins=200, range=(-1.0, 1.0))
expl_correlation_direct_hist, binedges = np.histogram(
    a=correlation_array_direct, bins=200, range=(-1.0, 1.0))
expl_correlation_intermediate_hist, binedges = np.histogram(
    a=correlation_array_intermediate, bins=200, range=(-1.0, 1.0))
expl_correlation_AXB_hist, binedges = np.histogram(
    a=correlation_array_AXB, bins=200, range=(-1, 1))


if read_correlation_file.strip().lower() == 'y':
    start = -1.0
    stop = 1.0
    step = 0.01
    bins = int((stop-start)/step)
    all_corr_edges = np.arange(start, stop+step, step)
    all_corr_hist = _histogram_for_large_files('input_data/corr_only.dat',
                                               bins, step, start)
    all_corr_hist_norm = all_corr_hist / sum(all_corr_hist)
    pylab.plot(all_corr_edges[:-1], all_corr_hist/2, label='all corr')

pylab.figure(1)
# all_corr_ll03 = pd.read_csv(
#    'input_data/depmap/_saved_Q3_depmap_INDRA_db_20180730_ll03_correlations'
#    '.csv', names=['gene1', 'gene2', 'corr'])
# all_corr_ll03_array = all_corr_ll03['corr'].values
# correlations_ll03_hist, binedges = np.histogram(a=all_corr_ll03_array,
#                                                bins=200, range=(-1.0, 1.0))
# correlations_ll03_hist_norm = correlations_ll03_hist / sum(
#    correlations_ll03_hist)
expl_correlation_hist_norm = expl_correlation_hist / sum(expl_correlation_hist)
expl_correlation_direct_hist_norm = expl_correlation_direct_hist / sum(
    expl_correlation_direct_hist)
expl_correlation_intermediate_hist_norm = \
    expl_correlation_intermediate_hist / sum(expl_correlation_intermediate_hist)
expl_correlation_AXB_hist_norm = expl_correlation_AXB_hist / sum(
    expl_correlation_AXB_hist)

# pylab.plot(binedges[:-1], correlations_ll03_hist,
#            label='all corr > 0.3')
pylab.plot(binedges[:-1], expl_correlation_hist, label='explained')
pylab.plot(binedges[:-1], expl_correlation_direct_hist, label='direct')
pylab.plot(binedges[:-1], expl_correlation_intermediate_hist,
           label='intermediate')
pylab.plot(binedges[:-1], expl_correlation_AXB_hist, label='AX, XB corr > 0.2')
pylab.xlabel('Correlation')
pylab.ylabel('Count')
pylab.legend(loc='upper left')
pylab.savefig('correlation_comparisons.png', dpi=500)
pylab.show()

print("x_is_intermediate_count_hist[:10]")
print(x_is_intermediate_count_hist[:10])
print("x_is_downstream_count_hist[:10]")
print(x_is_downstream_count_hist[:10])

pylab.figure(2)
pylab.bar(x_count_binedges[1:], x_is_downstream_count_hist,
          align='edge',
          log=True,
          label='Shared targets')
pylab.bar(x_count_binedges[1:], x_is_intermediate_count_hist,
          align='edge',
          log=True,
          label='Intermediates')
pylab.xlabel('Count of X')
pylab.ylabel('log(occurences of count)')
pylab.legend(loc='upper right')
pylab.title('log of occurences of X count (excluding shared regulator)')
pylab.savefig('X_count.png', dpi=500)
pylab.show()
x_count_df = x_count_df.append(other=x_count_df_list, ignore_index=True)
x_count_df.sort_values(inplace=True, by='count', ascending=False)
x_count_df.to_csv('X_counts.csv', index=None)
