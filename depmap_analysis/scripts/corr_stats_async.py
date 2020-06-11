from ctypes import c_wchar_p
from time import time
from datetime import datetime
from multiprocessing import Pool, cpu_count, Array
import logging

import numpy as np

logger = logging.getLogger(__name__)

global_results = []
global_results_pairs = []
global_vars = {}
list_of_genes = []


def _list_chunk_gen(lst, size, shuffle=False):
    """Given list, generate chunks <= size
    If shuffle is True, randomize input list before creating gnereator
    """
    if shuffle:
        np.random.shuffle(lst)
    n = max(1, size)
    return (lst[k:k+n] for k in range(0, len(lst), n))


def success_callback(res):
    global_results.append(res)


def success_callback_pairs(res):
    global_results_pairs.append(res)


class GlobalVars(object):
    def __init__(self, df=None, z_cm=None, sampl=10):
        if df is not None:
            global_vars['df'] = df
        if sampl:
            global_vars['subset_size'] = sampl
        if z_cm is not None:
            global list_of_genes
            global_vars['z_cm'] = z_cm
            list_of_genes = Array(c_wchar_p,
                                  np.array(z_cm.columns.values),
                                  lock=False)

    @staticmethod
    def update_global_vars(**kwargs):
        for varkey, obj in kwargs.items():
            global_vars[varkey] = obj
            if varkey == 'z_cm':
                global list_of_genes
                list_of_genes = Array(c_wchar_p,
                                      np.array(obj.columns.values),
                                      lock=False)

    @staticmethod
    def get_global_var_names():
        return set(global_vars.keys())

    @staticmethod
    def assert_global_vars(varnames):
        """

        varname : set(str)
            Set of names of variables to checkif they exists

        Returns
        -------
        bool
            True if variables in varnames are initialized
        """
        return all([global_vars.get(k, None) is not None for k in varnames])

    @staticmethod
    def assert_vars():
        """Same as assert_global_vars but with the shared array as well"""
        df_exists = global_vars.get('df', False) is not False
        z_cm_exists = global_vars.get('z_cm', False) is not False
        ssize_exists = global_vars.get('subset_size', False) is not False
        shared_ar_exists = bool(len(list_of_genes[:]))
        return df_exists and z_cm_exists and \
            ssize_exists and shared_ar_exists


# ToDo: make one work submitting function as a wrapper and provide the inner
#  loop as argument
def get_pairs_mp(ab_corr_pairs, max_proc=cpu_count(), max_pairs=10000):
    logger.info("Stratifying correlations by interaction type")
    logger.info(
        f'Starting workers for pairs at '
        f'{datetime.now().strftime("%H:%M:%S")} '
        f'with about {len(ab_corr_pairs)} pairs to check'
    )
    tstart = time()
    max_proc = min(cpu_count(), max_proc) if max_proc is not None else \
        cpu_count()
    if max_proc < 1:
        logger.warning('Max processes is set to < 1, resetting to 1')
        max_proc = 1

    if max_pairs and len(ab_corr_pairs) > max_pairs:
        corr_pairs = np.random.choice(ab_corr_pairs, size=max_pairs,
                                      replace=False)
    else:
        corr_pairs = ab_corr_pairs

    # Loop workers
    with Pool(max_proc) as pool:
        # Split up number of pairs
        size = len(corr_pairs) // max_proc + 1 if max_proc > 1 else 1
        lst_gen = _list_chunk_gen(lst=list(corr_pairs),
                                  size=size,
                                  shuffle=True)
        for corr_pairs in lst_gen:
            pool.apply_async(
                func=get_pairs,
                args=(corr_pairs, ),
                callback=success_callback_pairs
            )
        logger.info('Done submitting work to pool of workers')
        pool.close()
        logger.info('Pool is closed')
        pool.join()
        logger.info('Pool is joined')

    logger.info(f'Execution time: {time() - tstart} seconds')
    logger.info(f'Done at {datetime.now().strftime("%H:%M:%S")}')

    # Assemble results
    logger.info('Assembling results')
    results_pairs = set()
    for s in global_results_pairs:
        results_pairs.update(s)
    assert len(results_pairs) <= max_pairs, len(results_pairs)
    return results_pairs


def get_pairs(corr_pairs):
    # Get global args
    expl_df = global_vars['df']

    # Pairs where a-x-b AND a-b explanation exists
    pairs_axb_direct = set()

    # Pairs where a-x-b AND NOT a-b explanation exists
    pairs_axb_only = set()

    # all a-x-b "pathway" explanations, should be union of the above two
    pairs_any_axb = set()

    for s, o in corr_pairs:
        # Make sure we don't try to explain self-correlations
        if s == o:
            continue
        # Get all interaction types associated with given subject s and
        # object o
        int_types = set(expl_df['expl type'][(expl_df['agA'] == s) &
                                             (expl_df['agB'] == o)].values)
        # Check intersection of types
        axb_types = {'a-x-b', 'b-x-a', 'shared target'}.intersection(int_types)
        if axb_types:
            pairs_any_axb.add((s, o))
            if 'direct' in int_types:
                # Direct and pathway
                pairs_axb_direct.add((s, o))
            else:
                # Pathway and NOT direct
                pairs_axb_only.add((s, o))
        else:
            # Skip direct only
            continue

    # The union should be all pairs where a-x-b explanations exist
    ab_axb_union = pairs_axb_direct.union(pairs_axb_only)
    assert ab_axb_union == pairs_any_axb
    return pairs_axb_only


def get_corr_stats_mp(so_pairs, max_proc=cpu_count()):
    logger.info(
        f'Starting workers at {datetime.now().strftime("%H:%M:%S")} '
        f'with about {len(so_pairs)} pairs to check')
    tstart = time()

    max_proc = min(cpu_count(), max_proc)
    if max_proc < 1:
        logger.warning('Max processes is set to < 1, resetting to 1')
        max_proc = 1

    with Pool(max_proc) as pool:
        # Split up so_pairs in equal chunks
        size = len(so_pairs) // max_proc + 1 if max_proc > 1 else 1
        lst_gen = _list_chunk_gen(lst=list(so_pairs),
                                  size=size,
                                  shuffle=True)
        for pairs in lst_gen:
            pool.apply_async(
                func=get_corr_stats,
                args=(pairs, ),
                callback=success_callback,
            )
        logger.info('Done submitting work to pool of workers')
        pool.close()
        logger.info('Pool is closed')
        pool.join()
        logger.info('Pool is joined')
    logger.info(f'Execution time: {time() - tstart} seconds')
    logger.info(f'Done at {datetime.now().strftime("%H:%M:%S")}')

    logger.info('Assembling results')
    results = [[], [], [], [], []]
    for done_res in global_results:
        # Var name: all_x_corrs; Dict key: 'all_axb_corrs'
        results[0].extend(done_res['all_axb_corrs'])
        # Var name: avg_x_corrs; Dict key: axb_avg_corrs
        results[1].extend(done_res['axb_avg_corrs'])
        # Var name: top_x_corrs; Dict key: top_axb_corrs
        results[2].extend(done_res['top_axb_corrs'])
        # Var name: all_azb_corrs; Dict key: all_azb_corrs
        results[3].extend(done_res['all_azb_corrs'])
        # Var name: azb_avg_corrs; Dict key: azb_avg_corrs
        results[4].extend(done_res['azb_avg_corrs'])

    return results


def get_corr_stats(so_pairs):
    global list_of_genes
    df = global_vars['df']
    z_corr = global_vars['z_cm']
    subset_size = global_vars['subset_size']
    chunk_size = max(len(list_of_genes[:]) // subset_size, 1)

    all_axb_corrs = []
    top_axb_corrs = []
    azb_avg_corrs = []
    all_azb_corrs = []
    axb_avg_corrs = []
    for subj, obj in so_pairs:
        ab_avg_corrs = []
        path_rows = df[(df['agA'] == subj) &
                       (df['agB'] == obj) &
                       ((df['expl type'] == 'a-x-b') |
                        (df['expl type'] == 'b-x-a') |
                       (df['expl type'] == 'shared target'))]

        # Make sure we don't double-count Xs such that X is shared target
        # and also a pathway; also, don't include X if X = subj or obj
        x_set = set()
        for ix, path_row in path_rows.iterrows():
            x_set.update([x for x in path_row['expl data']
                          if x not in (subj, obj)])

        for x in x_set:
            if x in z_corr.columns:
                ax_corr = z_corr.loc[subj, x]
                xb_corr = z_corr.loc[x, obj]
                all_axb_corrs += [ax_corr, xb_corr]
                avg_corr = 0.5 * abs(ax_corr) + 0.5 * abs(xb_corr)
                ab_avg_corrs.append(avg_corr)
                axb_avg_corrs.append(avg_corr)
            else:
                # if warn < 3:
                #     warn += 1
                #     logger.warning('%s does not exist in the crispr or '
                #                    'rnai correlation matrices.' % x)
                # else:
                #     logger.warning('Muting warnings...')
                continue

        # Get a random subset of the possible correlation z scores
        for z in np.random.choice(list_of_genes[:], chunk_size, False):
            if z == subj or z == obj:
                continue
            az_corr = z_corr.loc[z, subj]
            bz_corr = z_corr.loc[z, obj]
            if np.isnan(az_corr) or np.isnan(bz_corr):
                continue  # Is there a more efficient way of doing this?
            all_azb_corrs.extend([az_corr, bz_corr])
            azb_avg_corrs.append(0.5 * abs(az_corr) + 0.5 * abs(bz_corr))

        # if warn:
        #     logger.warning('%d missing X genes out of %d in correlation '
        #                    'matrices' % (warn, len(x_list)))
        if len(ab_avg_corrs) > 0:
            max_magn_avg = max(ab_avg_corrs)
            top_axb_corrs.append((subj, obj, max_magn_avg))
            #ab_to_axb_avg.append((subj, obj, comb_zsc, max_magn_avg))

    return {'all_axb_corrs': all_axb_corrs,
            'axb_avg_corrs': axb_avg_corrs,
            'top_axb_corrs': top_axb_corrs,
            'all_azb_corrs': all_azb_corrs,
            'azb_avg_corrs': azb_avg_corrs}
