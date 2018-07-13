import csv
import argparse as ap
import logging
import pandas as pd
import numpy as np
import depmap_network_functions as dnf
from indra.tools.assemble_corpus import dump_statements
from indra.tools.assemble_corpus import load_statements as ac_load_stmts
from time import time
import pdb

# Temp fix because indra doesn't import when script is run from terminal
import sys
sys.path.append('~/repos/indra/')

logger = logging.getLogger('SlowClap')

# 1. no geneset -> use corr from full DepMap data, no filtering needed
# 2. geneset, not strict -> interaction has to contain at least one gene from
#    the loaded list
# 3. geneset loaded, strict -> each correlation can only contain genes from
#    the loaded data set


def _read_gene_set_file(gf, data):
    gset = []
    with open(gf, 'rt') as f:
        for g in f.readlines():
            gn = g.upper().strip()
            if gn in data:
                gset.append(gn)
    return gset


def _is_float(n):
    if type(n) is np.float64 or type(n) is float:
        return True
    else:
        return False


def main(args):

    # Open data file
    data = pd.read_csv(args.ceres_file, index_col=0, header=0)
    data = data.T

    if args.geneset_file:
        # Read gene set to look at
        gene_filter_list = _read_gene_set_file(gf=args.geneset_file, data=data)

    # 1. no loaded gene list OR 2. loaded gene list but not strict -> data.corr
    if not args.geneset_file or (args.geneset_file and not args.strict):
        # Calculate the full correlations, or load from cached
        if args.recalc:
            corr = data.corr()
            corr.to_hdf('correlations.h5', 'correlations')
        else:
            corr = pd.read_hdf(args.corr_file, 'correlations')
        # No gene set file, leave 'corr' intact and unstack
        if not args.geneset_file:
            fcorr_list = corr.unstack()

        # Gene set file present: filter and unstack
        elif args.geneset_file and not args.strict:
            fcorr_list = corr[gene_filter_list].unstack()

    # 3. Strict: both genes in interaction must be from loaded set
    elif args.geneset_file and args.strict:
        fcorr_list = data[gene_filter_list].corr().unstack()

    # Remove self correlation, correlations below ll, sort on magnitude,
    # leave correlation intact
    flarge_corr = fcorr_list[fcorr_list != 1.0]
    flarge_corr = flarge_corr[flarge_corr.abs() > args.ll]
    if args.ul < 1.0:
        flarge_corr = flarge_corr[flarge_corr.abs() < args.ul]
    fsort_corrs = flarge_corr[
        flarge_corr.abs().sort_values(ascending=False).index]

    # Find out if the HGNC pairs are connected and if they are how
    all_hgnc_ids = set()
    uniq_pairs = set()
    with open(args.outbasename+'_all.csv', 'w', newline='') as csvf:
        wrtr = csv.writer(csvf, delimiter=',')
        for pair in fsort_corrs.items():
            (id1, id2), correlation = pair
            if frozenset([id1, id2, correlation]) not in uniq_pairs:
                uniq_pairs.add(frozenset([id1, id2, correlation]))
                all_hgnc_ids.update([id1, id2])
                wrtr.writerow([id1, id2, correlation])

    # Get statements from file or from database that contain any gene from
    # provided list as set
    if args.statements_in:  # Get statments from file
        stmts_all = set(ac_load_stmts(args.statements_in))
    else:  # Use api to get statements. NOT the same as querying for each ID
        if args.geneset_file:
            stmts_all = dnf.dbc_load_statements(gene_filter_list)
        else:
            # if there is no gene set file, restrict to gene ids in
            # correlation data
            stmts_all = dnf.dbc_load_statements(list(all_hgnc_ids))

    # Dump statements to pickle file if output name has been given
    if args.statements_out:
        dump_statements(stmts=stmts_all, fname=args.statements_out)

    # Get nested dicts from statements
    nested_dict_statements = dnf.nested_dict_gen(stmts_all)

    # Loop through the unique pairs
    dir_conn_pairs = []
    dir_neg_conn_pairs = []
    unexplained = []
    npairs = len(uniq_pairs)

    f_con = open(args.outbasename + '_connections_text.txt', 'w')

    f_neg_c = open(args.outbasename + '_neg_conn_text.txt', 'w')

    logger.info('Looking for connections between %i pairs' % npairs)
    for pair in uniq_pairs:
        pl = list(pair)
        for li in pl:
            if _is_float(li):
                correlation = li
                break
        pl.remove(correlation)
        id1, id2 = pl

        forward_fail = False
        backward_fail = False

        # nested_dict_statements.get(id1).get(id2) raises AttributeError
        # if nested_dict_statements.get(id1) returns {}

        # Checks subj=id1, obj=id2
        if nested_dict_statements.get(id1) and \
                nested_dict_statements.get(id1).get(id2):
            stmts = nested_dict_statements[id1][id2]
            logger.info('Found connection between %s and %s' % (id1, id2))
            dir_conn_pairs.append((id1, id2, correlation, stmts))
            f_con.write(dnf.str_output(subj=id1, obj=id2, corr=correlation,
                                       stmts=stmts, ignore_str='parent'))

            if correlation < 0:
                dir_neg_conn_pairs.append((id1, id2, correlation, stmts))
                f_neg_c.write(dnf.str_output(subj=id1, obj=id2,
                                             corr=correlation,
                                             stmts=stmts,
                                             ignore_str='parent'))
        else:
            forward_fail = True

        # Checks subj=id2, obj=id1
        if nested_dict_statements.get(id2) and \
                nested_dict_statements.get(id2).get(id1):
            stmts = nested_dict_statements[id2][id1]
            logger.info('Found connection between %s and %s' % (id2, id1))
            dir_conn_pairs.append((id2, id1, correlation, stmts))
            f_con.write(dnf.str_output(subj=id2, obj=id1, corr=correlation,
                                       stmts=stmts, ignore_str='parent'))

            if correlation < 0:
                dir_neg_conn_pairs.append((id2, id1, correlation, stmts))
                f_neg_c.write(dnf.str_output(subj=id2, obj=id1,
                                             corr=correlation,
                                             stmts=stmts,
                                             ignore_str='parent'))

        else:
            backward_fail = True

        # If both failed, count as unexaplined
        if forward_fail and backward_fail:
            unexplained.append([id1, id2, correlation])

    with open(args.outbasename+'_connections.csv', 'w', newline='') as csvf:
        wrtr = csv.writer(csvf, delimiter=',')
        wrtr.writerows(dir_conn_pairs)

    with open(args.outbasename+'_neg_conn.csv', 'w', newline='') as csvf:
        wrtr = csv.writer(csvf, delimiter=',')
        wrtr.writerows(dir_neg_conn_pairs)

    with open(args.outbasename+'_unexplained.csv', 'w', newline='') as csvf:
        wrtr = csv.writer(csvf, delimiter=',')
        wrtr.writerows(unexplained)

    f_con.close()
    f_neg_c.close()


if __name__ == '__main__':
    parser = ap.ArgumentParser()
    parser.add_argument('-c', '--corr-file', required=True,
                        help='Precalculated correlations in h5 format')
    parser.add_argument('-f', '--ceres-file', required=True,
                        help='Ceres correlation score in csv format')
    parser.add_argument('-g', '--geneset-file',
                        help='Filter to interactions with gene set data file.')
    parser.add_argument('-o', '--outbasename', default=str(int(time())),
                        help='Base name for outfiles. Default: UTC timestamp')
    parser.add_argument('-r', '--recalc', action='store_true',
                        help='With \'-r\', recalculate full gene-gene '
                             'correlation of data set.')
    parser.add_argument('-s', '--strict', action='store_true',
                        help='With \'-s\', the correlations are restricted to '
                             'only be between loaded gene set. If no gene set '
                             'is loaded, this option has no effect.')
    parser.add_argument('-sti', '--statements-in', help='Loads a pickle file '
                                                        'to use instead of '
                                                        'quering a database.')
    parser.add_argument('-sto', '--statements-out', help='Saves the used '
                                                         'statements read from '
                                                         'the database')
    parser.add_argument('-ll', type=float, default=0.5,
                        help='Lower limit CERES correlation score filter.')
    parser.add_argument('-ul', type=float, default=1.0,
                        help='Upper limit CERES correlation score filter.')

    a = parser.parse_args()
    main(a)
