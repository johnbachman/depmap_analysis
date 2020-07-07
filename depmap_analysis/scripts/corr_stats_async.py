from time import time
from ctypes import c_wchar_p
from datetime import datetime
from collections import Counter
from multiprocessing import Pool, cpu_count, Array, current_process
import logging
import random

import numpy as np
from pybel.dsl.node_classes import CentralDogma

from indra.util.multiprocessing_traceback import WrapException
from indra.databases.hgnc_client import get_current_hgnc_id, get_uniprot_id,\
    uniprot_ids, get_hgnc_name

logger = logging.getLogger(__name__)

uniprot_ids_reverse = {v: k for k, v in uniprot_ids.items()}

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


def error_callback(err):
    logger.error(f'An error occurred in process {current_process().pid}')
    logger.exception(err)


class GlobalVars(object):
    def __init__(self, df=None, z_cm=None, reactome=None, sampl=10):
        if df is not None:
            global_vars['df'] = df
        if sampl:
            global_vars['subset_size'] = sampl
        if reactome:
            global_vars['reactome'] = reactome
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
            Set of names of variables to check if they exists

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
        f'with {len(ab_corr_pairs)} pairs to check'
    )
    tstart = time()
    max_proc = min(cpu_count(), max_proc) if max_proc is not None else \
        cpu_count()
    if max_proc < 1:
        logger.warning('Max processes is set to < 1, resetting to 1')
        max_proc = 1

    if max_pairs and len(ab_corr_pairs) > max_pairs:
        logger.info(f'Down sampling ab_corr_pairs to {max_pairs}')
        corr_pairs = random.sample(
            ab_corr_pairs, max_pairs
        )
    else:
        corr_pairs = ab_corr_pairs

    # Loop workers
    with Pool(max_proc) as pool:
        # Split up number of pairs
        size = len(corr_pairs) // max_proc + 1 if max_proc > 1 else 1
        lst_gen = _list_chunk_gen(lst=list(corr_pairs),
                                  size=size,
                                  shuffle=True)
        for chunk_of_pairs in lst_gen:
            pool.apply_async(
                func=get_pairs,
                args=(chunk_of_pairs, ),
                callback=success_callback_pairs,
                error_callback=error_callback
            )
        logger.info('Done submitting work to pool of workers')
        pool.close()
        logger.info('Pool is closed')
        pool.join()
        logger.info('Pool is joined')

    logger.info(f'Execution time: {time() - tstart} seconds')
    logger.info(f'Done at {datetime.now().strftime("%H:%M:%S")}')

    # Assemble results
    logger.info(f'Assembling {len(global_results_pairs)} results')
    results_pairs = set()
    for s in global_results_pairs:
        results_pairs.update(s)
    assert len(results_pairs) <= len(ab_corr_pairs)
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
                error_callback=error_callback
            )
        logger.info('Done submitting work to pool of workers')
        pool.close()
        logger.info('Pool is closed')
        pool.join()
        logger.info('Pool is joined')
    logger.info(f'Execution time: {time() - tstart} seconds')
    logger.info(f'Done at {datetime.now().strftime("%H:%M:%S")}')

    logger.info(f'Assembling {len(global_results)} results')
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
    try:
        global list_of_genes
        df = global_vars['df']
        z_corr = global_vars['z_cm']
        reactome = global_vars['reactome']
        subset_size = global_vars['subset_size']
        chunk_size = max(len(list_of_genes[:]) // subset_size, 1)

        all_axb_corrs = []
        top_axb_corrs = []
        axb_avg_corrs = []

        azb_avg_corrs = []
        all_azb_corrs = []

        reactome_avg_corrs = []
        all_reactome_corrs = []

        # reset counters
        counter = Counter({'r_skip': 0,
                           'z_skip': 0,
                           'x_skip': 0})

        for subj, obj in so_pairs:
            # Get x values
            (avg_x_corrs_per_ab, axb_corrs), x_len = \
                get_interm_corr_stats_x(subj, obj, z_corr, df)
            all_axb_corrs += axb_corrs
            if len(avg_x_corrs_per_ab) > 0:
                max_magn_avg = max(avg_x_corrs_per_ab)
                axb_avg_corrs += avg_x_corrs_per_ab
                top_axb_corrs.append((subj, obj, max_magn_avg))
            counter['x_skip'] += x_len - len(avg_x_corrs_per_ab)

            # Get z values
            z_iter = np.random.choice(list_of_genes[:], chunk_size, False)
            avg_z_corrs_per_ab, azb_corrs = \
                get_interm_corr_stats_z(subj, obj, z_iter, z_corr)
            azb_avg_corrs += avg_z_corrs_per_ab
            all_azb_corrs += azb_corrs
            counter['z_skip'] += len(z_iter) - len(avg_z_corrs_per_ab)

            # Get reactome values
            (avg_reactome_corrs_per_ab, reactome_corrs), r_len = \
                get_interm_corr_stats_reactome(subj, obj, reactome, z_corr)
            reactome_avg_corrs += avg_reactome_corrs_per_ab
            counter['r_skip'] += r_len - len(avg_reactome_corrs_per_ab)

            ## OLD ##

            # path_rows = df[(df['agA'] == subj) &
            #                (df['agB'] == obj) &
            #                ((df['expl type'] == 'a-x-b') |
            #                 (df['expl type'] == 'b-x-a') |
            #                (df['expl type'] == 'shared target'))]
            #
            # # Make sure we don't double-count Xs such that X is shared target
            # # and also a pathway; also, don't include X if X = subj or obj
            # x_set = set()
            # for ix, path_row in path_rows.iterrows():
            #     x_set.update([x for x in path_row['expl data']
            #                   if x not in (subj, obj)])
            #
            # for x in x_set:
            #     x_name = x.name if isinstance(x, CentralDogma) else x
            #     if x_name in z_corr.columns:
            #         ax_corr = z_corr.loc[subj, x_name]
            #         xb_corr = z_corr.loc[x_name, obj]
            #         all_axb_corrs += [ax_corr, xb_corr]
            #         avg_corr = 0.5 * abs(ax_corr) + 0.5 * abs(xb_corr)
            #         avg_x_corrs_per_ab.append(avg_corr)
            #     else:
            #         continue
            #
            # # Get a random subset of the possible correlation z scores
            # for z in np.random.choice(list_of_genes[:], chunk_size, False):
            #     try:
            #         if z == subj or z == obj:
            #             counter.update('z_skip')
            #             continue
            #         az_corr = z_corr.loc[z, subj]
            #         bz_corr = z_corr.loc[z, obj]
            #         if np.isnan(az_corr) or np.isnan(bz_corr):
            #             # Is there a more efficient way of doing this?
            #             logger.info(
            #                 f'NaN correlations for '
            #                 f'subj-z ({str(subj)}-{str(z)}) or '
            #                 f'obj-z ({str(obj)}-{str(z)})'
            #             )
            #             continue
            #         all_azb_corrs.extend([az_corr, bz_corr])
            #         azb_avg_corrs.append(0.5*abs(az_corr) + 0.5*abs(bz_corr))
            #     except KeyError as err:
            #         raise KeyError(
            #             f'KeyError was raised trying to sample background '
            #             f'correlation distribution with subject {str(subj)}'
            #             f'({subj.__class__}), object {str(obj)} '
            #             f'({obj.__class__}) and z {z} ({z.__class__})'
            #         ) from err

        # TODO Add reactome output to returned dict
        assert all(len(cl) for cl in [all_axb_corrs, axb_avg_corrs,
                                      top_axb_corrs, all_azb_corrs,
                                      azb_avg_corrs]),\
            logger.error(f'No results for process {current_process().pid}. '
                         f'Stats:\nall_axb_corrs: {len(all_axb_corrs)}, '
                         f'axb_avg_corrs: {len(axb_avg_corrs)}, '
                         f'top_axb_corrs: {len(top_axb_corrs)}, '
                         f'all_azb_corrs: {len(all_azb_corrs)}, '
                         f'azb_avg_corrs: {len(azb_avg_corrs)}')

        for k, v in counter.items():
            if v > 0:
                logger.info(f'Skipped {k} {v} times')

        return {'all_axb_corrs': all_axb_corrs,
                'axb_avg_corrs': axb_avg_corrs,
                'top_axb_corrs': top_axb_corrs,
                'all_azb_corrs': all_azb_corrs,
                'azb_avg_corrs': azb_avg_corrs}
    except Exception:
        raise WrapException()


def get_interm_corr_stats_x(subj, obj, z_corr, df):
    path_rows = df[(df['agA'] == subj) &
                   (df['agB'] == obj) &
                   ((df['expl type'] == 'a-x-b') |
                    (df['expl type'] == 'b-x-a') |
                    (df['expl type'] == 'shared target'))]
    x_set = set()
    for ix, path_row in path_rows.iterrows():
        x_set.update([x for x in path_row['expl data']
                      if x not in (subj, obj)])
    return _get_interm_corr_stats(subj, obj, x_set, z_corr), len(x_set)


def get_interm_corr_stats_z(subj, obj, z_set, z_corr):
    return _get_interm_corr_stats(subj, obj, z_set, z_corr)


def get_interm_corr_stats_reactome(subj, obj, reactome, z_corr):
    pathways_by_gene, genes_by_pathway, _ = reactome
    subj_up = _hgncsym2up(subj)
    if subj_up is None:
        return [], []
    obj_up = _hgncsym2up(obj)
    if obj_up is None:
        return [], []

    paths = set(pathways_by_gene[subj_up]) & set(pathways_by_gene[obj_up])
    gene_set = set()
    for rp in paths:
        gene_set.update([g for g in genes_by_pathway[rp]])
    hgnc_gene_set = [_up2hgncsym(up) for up in gene_set]
    return _get_interm_corr_stats(subj, obj, hgnc_gene_set, z_corr),\
           len(hgnc_gene_set)


def _get_interm_corr_stats(a, b, y_set, z_corr):
    # Get a list of the maximum ax-bx average per pair
    avg_y_corrs_per_ab = []

    # Get all ax and bc correlations
    all_ayb_corrs = []

    for y in y_set:
        try:
            # Skip self correlations and non-existing names
            if y is None or y == a or y == b or y not in z_corr.columns:
                continue
            ay_corr = z_corr.loc[y, a]
            by_corr = z_corr.loc[y, b]
            if np.isnan(ay_corr) or np.isnan(by_corr):
                # Is there a more efficient way of doing this?
                logger.info(
                    f'NaN correlations for subj-y ({str(a)}-{str(y)}) or '
                    f'obj-y ({str(b)}-{str(y)})'
                )
                continue

            all_ayb_corrs += [ay_corr, by_corr]
            avg_y_corrs_per_ab.append(0.5*(abs(ay_corr) + abs(by_corr)))

        except KeyError as ke:
            raise KeyError(
                f'KeyError was raised trying to sample '
                f'correlation distribution with subject {str(a)}'
                f'({a.__class__}), object {str(b)} '
                f'({b.__class__}) and intermediate {str(y)} ({y.__class__})'
            ) from ke
    return avg_y_corrs_per_ab, all_ayb_corrs


def _hgncsym2up(hgnc_symb):
    hgnc_id = get_current_hgnc_id(hgnc_symb)
    if isinstance(hgnc_id, list):
        ix = 0
        upid = None
        while upid is None:
            try:
                upid = get_uniprot_id(hgnc_id[ix])
            except IndexError:
                break
            ix += 1
    else:
        upid = get_uniprot_id(hgnc_id)
    return upid


def _up2hgncsym(up_id):
    return get_hgnc_name(uniprot_ids_reverse.get(up_id))
