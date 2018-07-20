import csv
import logging
import numpy as np
import pandas as pd
from time import time
import networkx as nx
import argparse as ap
import itertools as itt
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


def _is_float(n):
    if type(n) is np.float64 or type(n) is float:
        return True
    else:
        return False


def _corr_web_latex(id1, id2, fmtcorr):
    web_text = r'\section{{{}, {}: {}}}'.format(id1, id2, fmtcorr) + '\n' + \
                 r'See correlation plot \href{{' \
                 r'https://depmap.org/portal/interactive/?xDataset=Avana' \
                 r'&xFeature={}&yDataset=Avana&yFeature={}&colorDataset=' \
                 r'lineage&colorFeature=all&filterDataset=context' \
                 r'&filterFeature=&regressionLine=false&statisticsTable=false' \
                 r'&associationTable=true&plotOnly=false}}{{here}}'.format(
                     id1, id2) + '\n\n'

    return web_text


def main(args):

    # Evidence list must have length above this value to be part of the
    # text output
    ev_fltr = 0

    # Prepare data
    gene_filter_list, uniq_pairs, all_hgnc_ids, fsort_corrs = \
        dnf.get_correlations(args.ceres_file, args.geneset_file,
                             args.corr_file, args.strict, args.outbasename,
                             args.recalc, args.ll, args.ul)

    # Get statements from file or from database that contain any gene from
    # provided list as set
    if args.statements_in:  # Get statments from file
        stmts_all = set(ac.load_statements(args.statements_in))
    else:  # Use api to get statements. _NOT_ the same as querying for each ID
        if args.geneset_file:
            stmts_all = dnf.dbc_load_statements(gene_filter_list)
        else:
            # if there is no gene set file, restrict to gene ids in
            # correlation data
            stmts_all = dnf.dbc_load_statements(list(all_hgnc_ids))

    # Dump statements to pickle file if output name has been given
    if args.statements_out:
        logger.info('Dumping read statements')
        ac.dump_statements(stmts=stmts_all, fname=args.statements_out)

    # Get nested dicts from statements
    nested_dict_statements = dnf.nested_dict_gen(stmts_all)

    # Get undirected graph from nested dict
    undir_nx_graph = dnf.nx_undirected_graph_from_nested_dict(
        nest_d=nested_dict_statements)
    node_set = set(undir_nx_graph.nodes())

    # Loop through the unique pairs
    dir_conn_pairs = []  # Save pairs that are directly connected
    two_step_conn_pairs = []  # Extend to pairs connected by two direct edges
    dir_neg_conn_pairs = []  # Directly connected pairs with correlation < 0
    unexplained = []  # Unexplained correlations
    npairs = len(uniq_pairs)

    # Open files to write text/latex output
    with open(args.outbasename + '_connections_latex.tex', 'w') as f_con, \
            open(args.outbasename + '_neg_conn_latex.tex', 'w') as f_neg_c:

        logger.info('Looking for connections between %i pairs' % npairs)
        for pair in uniq_pairs:
            correlation, fmt_corr = None, None
            pl = list(pair)
            for li in pl:
                if _is_float(li):
                    correlation = li
                    fmt_corr = '{0:.04}'.format(correlation)
                    break
            pl.remove(correlation)
            id1, id2 = pl

            if (nested_dict_statements.get(id1) and
                    nested_dict_statements.get(id1).get(id2)) or \
                    (nested_dict_statements.get(id2) and
                     nested_dict_statements.get(id2).get(id1)):
                new_pair = _corr_web_latex(id1, id2, fmt_corr)
                f_con.write(new_pair)

                if correlation < 0:
                    f_neg_c.write(new_pair)

            # nested_dict_statements.get(id1).get(id2) raises AttributeError
            # if nested_dict_statements.get(id1) returns {}

            found = []

            for subj, obj in itt.permutations((id1, id2), r=2):
                if nested_dict_statements.get(subj) and \
                        nested_dict_statements.get(subj).get(obj):
                    stmts = nested_dict_statements[subj][obj]
                    logger.info('Found direct connection between %s and %s' % (
                        subj, obj))
                    found.append(True)
                    stmt_tuple = (subj, obj, correlation, stmts)
                    dir_conn_pairs.append(stmt_tuple)

                    # ToDo Should two_step_conn_pairs contain connections of
                    # length == 2 or <= 2? ==2 for now.
                    two_step_conn_pairs.append((subj, obj, correlation,
                                                0, 'direct'))
                    output = dnf.latex_output(subj=subj,
                                              obj=obj,
                                              corr=correlation,
                                              ev_len_fltr=ev_fltr,
                                              stmts=stmts,
                                              ignore_str='parent')
                    f_con.write(output)
                    if correlation < 0:
                        dir_neg_conn_pairs.append(stmt_tuple)
                        f_neg_c.write(output)

                else:
                    found.append(False)

                # any(found) is True if at least one connection was found and
                # therefore "not any" is only True when no connection was found
                if not any(found):
                    # Check if connection A-X-B exists and save the set of
                    # middle nodes.
                    #
                    # This can be built to a recursive search for n middle
                    # nodes by inserting an 'else' at the bottom where the
                    # neighbors of the neighbors are searched.
                    if subj in node_set and obj in node_set:
                        neighbor_node_connections = \
                            set(undir_nx_graph[subj].keys()) & \
                            set(undir_nx_graph[obj].keys())
                        if neighbor_node_connections:
                                logger.info('Found connection of length 2 '
                                            'between %s and %s' % (subj, obj))
                                two_step_conn_pairs.\
                                    append((subj, obj, correlation,
                                            len(neighbor_node_connections),
                                            neighbor_node_connections))

                        # no connection A-B or A-X-B found
                        else:
                            unexplained.append([id1, id2, correlation])

    with open(args.outbasename+'_conn_correlations.csv', 'w', newline='')\
            as csvf:
        wrtr = csv.writer(csvf, delimiter=',')
        wrtr.writerows(dir_conn_pairs)

    with open(args.outbasename + '_two_conn_correlations.csv', 'w', newline='')\
            as csvf:
        wrtr = csv.writer(csvf, delimiter=',')
        wrtr.writerows(two_step_conn_pairs)

    with open(args.outbasename+'_neg_conn_correlations.csv', 'w', newline='') \
            as csvf:
        wrtr = csv.writer(csvf, delimiter=',')
        wrtr.writerows(dir_neg_conn_pairs)

    with open(args.outbasename+'_unexpl_correlations.csv', 'w', newline='') \
            as csvf:
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
