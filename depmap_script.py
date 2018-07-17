import csv
import logging
import argparse as ap
from time import time
import numpy as np
import pandas as pd
import networkx as nx
from indra.tools import assemble_corpus as ac
import depmap_network_functions as dnf

# Temp fix because indra doesn't import when script is run from terminal
import sys
sys.path.append('~/repos/indra/')

logger = logging.getLogger('depmap_script')

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


def get_correlations(ceres_file, geneset_file, corr_file, strict, outbasename,
                     recalc=False, lower_limit=0.3, upper_limit=1.0):
    # Open data file
    logger.info("Reading CERES data from %s" % ceres_file)
    data = pd.read_csv(ceres_file, index_col=0, header=0)
    data = data.T

    if geneset_file:
        # Read gene set to look at
        gene_filter_list = _read_gene_set_file(gf=geneset_file, data=data)

    # 1. no loaded gene list OR 2. loaded gene list but not strict -> data.corr
    if not geneset_file or (geneset_file and not strict):
        # Calculate the full correlations, or load from cached
        if recalc:
            logger.info("Calculating correlations (may take a long time)")
            corr = data.corr()
            corr.to_hdf('correlations.h5', 'correlations')
        else:
            logger.info("Loading correlations from %s" % corr_file)
            corr = pd.read_hdf(corr_file, 'correlations')
        # No gene set file, leave 'corr' intact and unstack
        if not geneset_file:
            fcorr_list = corr.unstack()

        # Gene set file present: filter and unstack
        elif geneset_file and not strict:
            fcorr_list = corr[gene_filter_list].unstack()

    # 3. Strict: both genes in interaction must be from loaded set
    elif geneset_file and strict:
        fcorr_list = data[gene_filter_list].corr().unstack()

    # Remove self correlation, correlations below ll, sort on magnitude,
    # leave correlation intact
    logger.info("Removing self correlations")
    flarge_corr = fcorr_list[fcorr_list != 1.0] # Self correlations
    logger.info("Filtering correlations to %.1f < C < %.1f" %
                (lower_limit, upper_limit))
    flarge_corr = flarge_corr[flarge_corr.abs() > lower_limit]
    if upper_limit < 1.0:
        flarge_corr = flarge_corr[flarge_corr.abs() < upper_limit]
    # Sort by absolute value
    logger.info("Sorting correlations by absolute value")
    fsort_corrs = flarge_corr[
        flarge_corr.abs().sort_values(ascending=False).index]

    # Compile set of correlations to be explained without duplicates A-B, B-A
    logger.info("Compiling deduplicated set of correlations")
    all_hgnc_ids = set()
    uniq_pairs = []
    with open(outbasename+'_all.csv', 'w', newline='') as csvf:
        wrtr = csv.writer(csvf, delimiter=',')
        for pair in fsort_corrs.items():
            (id1, id2), correlation = pair
            if (id2, id1, correlation) not in uniq_pairs:
                uniq_pairs.append((id1, id2, correlation))
                all_hgnc_ids.update([id1, id2])
                wrtr.writerow([id1, id2, correlation])
    return uniq_pairs, all_hgnc_ids, fsort_corrs


def main(args):
    uniq_pairs, all_hgnc_ids, fsort_corrs = \
            get_correlations(args.ceres_file, args.geneset_file,
                             args.corr_file, args.strict,
                             args.outbasename, args.recalc, args.ll, args.ul)

    # Get statements from file or from database that contain any gene from
    # provided list as set
    if args.statements_in:  # Get statments from file
        stmts_all = set(ac.load_statements(args.statements_in))
    else:  # Use api to get statements. NOT the same as querying for each ID
        if args.geneset_file:
            stmts_all = dnf.dbc_load_statements(gene_filter_list)
        else:
            # if there is no gene set file, restrict to gene ids in
            # correlation data
            stmts_all = dnf.dbc_load_statements(list(all_hgnc_ids))

    # Dump statements to pickle file if output name has been given
    if args.statements_out:
        ac.dump_statements(stmts=stmts_all, fname=args.statements_out)

    # Get nested dicts from statements
    nested_dict_statements = dnf.nested_dict_gen(stmts_all)

    # Loop through the unique pairs
    dir_conn_pairs = []
    dir_neg_conn_pairs = []
    unexplained = []
    npairs = len(uniq_pairs)

    f_con = open(args.outbasename + '_connections_latex.tex', 'w')

    f_neg_c = open(args.outbasename + '_neg_conn_latex.tex', 'w')

    logger.info('Looking for connections between %i pairs' % npairs)
    for pair in uniq_pairs:
        pl = list(pair)
        for li in pl:
            if _is_float(li):
                correlation = li
                fmt_corr = '{0:.04}'.format(correlation)
                break
        pl.remove(correlation)
        id1, id2 = pl

        forward_fail = False
        backward_fail = False

        if (nested_dict_statements.get(id1) and
                nested_dict_statements.get(id1).get(id2)) or \
                (nested_dict_statements.get(id2) and
                 nested_dict_statements.get(id2).get(id1)):
            new_pair = r'\section{{{}, {}: {}}}'.format(id1, id2, fmt_corr) \
                 +'\n'+ \
                 r'See correlation plot \href{{' \
                 r'https://depmap.org/portal/interactive/?xDataset=Avana' \
                 r'&xFeature={}&yDataset=Avana&yFeature={}&colorDataset=' \
                 r'lineage&colorFeature=all&filterDataset=context' \
                 r'&filterFeature=&regressionLine=false&statisticsTable=false' \
                 r'&associationTable=true&plotOnly=false}}{{here}}'.format(
                     id1, id2) + '\n\n'
            f_con.write(new_pair)
            if correlation < 0:
                f_neg_c.write(new_pair)

        # nested_dict_statements.get(id1).get(id2) raises AttributeError
        # if nested_dict_statements.get(id1) returns {}

        ev_fltr = 0

        # Checks subj=id1, obj=id2
        if nested_dict_statements.get(id1) and \
                nested_dict_statements.get(id1).get(id2):
            stmts = nested_dict_statements[id1][id2]
            logger.info('Found connection between %s and %s' % (id1, id2))
            dir_conn_pairs.append((id1, id2, correlation, stmts))
            output = dnf.latex_output(subj=id1, obj=id2, corr=correlation,
                                      ev_len_fltr=ev_fltr, stmts=stmts,
                                      ignore_str='parent')
            f_con.write(output)

            if correlation < 0:
                dir_neg_conn_pairs.append((id1, id2, correlation, stmts))
                f_neg_c.write(output)
        else:
            forward_fail = True

        # Checks subj=id2, obj=id1
        if nested_dict_statements.get(id2) and \
                nested_dict_statements.get(id2).get(id1):
            stmts = nested_dict_statements[id2][id1]
            logger.info('Found connection between %s and %s' % (id2, id1))
            dir_conn_pairs.append((id2, id1, correlation, stmts))
            output = dnf.latex_output(subj=id2, obj=id1, corr=correlation,
                                      ev_len_fltr=ev_fltr, stmts=stmts,
                                      ignore_str='parent')
            f_con.write(output)

            if correlation < 0:
                dir_neg_conn_pairs.append((id2, id1, correlation, stmts))
                f_neg_c.write(output)

        else:
            backward_fail = True

        # If both failed, count as unexplained
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
    parser.add_argument('-ll', type=float, default=0.3,
                        help='Lower limit CERES correlation score filter.')
    parser.add_argument('-ul', type=float, default=1.0,
                        help='Upper limit CERES correlation score filter.')
    a = parser.parse_args()

    with open('dep_map_script_log{}'.format(str(int(time()))), 'w') as f:
        f.write('Command line options used:\n')
        f.write('Correlation file: {}'.format(a.corr_file))
        f.write('Ceres file: {}'.format(a.ceres_file))
        f.write('Geneset file: {}'.format(a.geneset_file))
        f.write('Ignore correlations below: {}'.format(a.ll))
        f.write('Ignore correlations above: {}'.format(a.ul))
        f.write('Output basename: {}'.format(a.outbasename))
        f.write('Recalculate correlations: {}'.format(a.recalc))
        f.write('Strict: {}'.format(a.strict))
        f.write('Pickled statement file: {}'.format(a.statements_in))
        f.write('Output loaded statements to: {}'.format(a.statements_out))

    main(a)
