import numpy as np
from time import time
from datetime import datetime
from multiprocessing import Pool, cpu_count

import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

global_results = []
global_vars = {}


def _list_chunk_gen(lst, size):
    """Given list, generate chunks <= size"""
    n = max(1, size)
    return (lst[k:k+n] for k in range(0, len(lst), n))


def success_callback(res):
    global_results.append(res)


class GlobalVars(object):
    def __init__(self, df, z_cm):
        global_vars['df'] = df
        global_vars['z_cm'] = z_cm

    @staticmethod
    def assert_vars():
        df_exists = global_vars.get('df', False) is not False
        z_cm_exists = global_vars.get('z_cm', False) is not False
        return df_exists and z_cm_exists


def get_corr_stats_mp(so_pairs):
    logger.info(
        f'Starting workers at {datetime.now().strftime("%H:%M:%S")} '
        f'with about {len(so_pairs)} pairs to check')
    tstart = time()

    promises = []
    with Pool() as pool:
        # Split up so_pairs in equal chunks
        lst_gen = _list_chunk_gen(lst=list(so_pairs),
                                  size=len(so_pairs) // cpu_count() + 1)
        for n, pairs in enumerate(lst_gen):
            # async_res = pool.apply_async(
            pool.apply_async(
                func=get_corr_stats,
                args=(pairs, ),
                callback=success_callback,
            )
            max_n = n
            # promises.append(async_res)
        logger.info('Done submitting work to pool of workers')
        pool.close()
        logger.info('Pool is closed')
        pool.join()
        logger.info('Pool is joined')
    logger.info(f'Execution time: {time() - tstart} seconds')
    logger.info(f'Done at {datetime.now().strftime("%H:%M:%S")}')
    # results = [async_res.get()]

    # retries = 0
    # while len(global_results) < max_n:
    #     sleep(0.5)
    #     retries += 1
    #     if retries == retries_await:
    #         logger.warning('Timed out waiting for results')
    #         break
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
    df = global_vars['df']
    z_corr = global_vars['z_cm']

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

        # Get the set of possible correlation z scores
        for z in set(z_corr.columns):
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

    return {'all_axb_corrs': all_axb_corrs, 'axb_avg_corrs': axb_avg_corrs,
            'top_axb_corrs': top_axb_corrs, 'all_azb_corrs': all_azb_corrs,
            'azb_avg_corrs': azb_avg_corrs}
