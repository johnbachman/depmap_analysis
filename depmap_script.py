import csv
import json
import logging
import pandas as pd
import pickle as pkl
from time import time
import networkx as nx
import argparse as ap
import itertools as itt
from numpy import float64
from collections import defaultdict
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


def nest_dict():
    """Returns a nested dictionary of arbitrary depth

    Returns
    -------
    defaultdict(nest_dict)
    """
    return defaultdict(nest_dict)


def _entry_exist(nest_dict, outer_key, inner_key):
    if nest_dict.get(outer_key) and nest_dict.get(outer_key).get(inner_key):
        return True
    else:
        return False


def _dump_it_to_pickle(fname, pyobj):
    with open(fname, 'wb') as po:
        pkl.dump(obj=pyobj, file=po)


def _dump_it_to_json(fname, pyobj):
    with open(fname, 'w') as json_out:
        json.dump(pyobj, json_out)


def _dump_it_to_csv(fname, pyobj, separator=','):
    with open(fname, 'w', newline='') as csvf:
        wrtr = csv.writer(csvf, delimiter=separator)
        wrtr.writerows(pyobj)


def _is_float(n):
    if type(n) is float64 or type(n) is float:
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

    # Prepare data (we need uniq_pairs to look for explainable interactions)
    gene_filter_list, uniq_pairs, all_hgnc_ids, fsort_corrs = \
        dnf.get_correlations(args.ceres_file, args.geneset_file,
                             args.corr_file, args.strict, args.outbasename,
                             args.unique_depmap_pairs, args.recalc, args.ll,
                             args.ul)

    # Get statements from file or from database that contain any gene from
    # provided list as set unless you're already loading a pre-calculated
    # nested dict and/or precalculated directed graph.

    if not (args.light_weight_stmts or args.nested_dict_in):
        if args.statements_in:  # Get statments from file
            stmts_all = set(ac.load_statements(args.statements_in))
        # Use api to get statements._NOT_ the same as querying for each ID
        else:
            if args.geneset_file:
                stmts_all = dnf.dbc_load_statements(gene_filter_list)
            else:
                # if there is no gene set file, restrict to gene ids in
                # correlation data
                stmts_all = dnf.dbc_load_statements(list(all_hgnc_ids))

        # Dump statements to pickle file if output name has been given
        if args.statements_out:
            logger.info('Dumping read raw statements')
            ac.dump_statements(stmts=stmts_all, fname=args.statements_out)

    # Get nested dicts from statements
    if args.light_weight_stmts:
        hash_df = pd.read_csv(args.light_weight_stmts, delimiter='\t')
        nested_dict_statements = dnf.nested_hash_dict_from_pd_dataframe(hash_df)
    elif args.nested_dict_in:
        with open(args.nested_dict_in, 'rb') as rpkl:
            nested_dict_statements = pkl.load(rpkl)
    else:
        nested_dict_statements = dnf.dedupl_nested_dict_gen(stmts_all)
        if args.nested_dict_out:
            with open(args.nested_dict_out, 'wb') as wpkl:
                pkl.dump(obj=nested_dict_statements, file=wpkl)

    # Get undirected graph from nested dict
    undir_nx_graph = dnf.nx_undirected_graph_from_nested_dict(
        nest_d=nested_dict_statements)
    undir_node_set = set(undir_nx_graph.nodes())

    # Get directed simple graph
    if args.directed_graph_in:
        with open(args.directed_graph_in, 'rb') as rpkl:
            nx_dir_graph = pkl.load(rpkl)
    else:
        nx_dir_graph = dnf.nx_directed_graph_from_nested_dict_3layer(
            nest_d=nested_dict_statements)
        # Save as pickle file
        if args.directed_graph_out:
            with open(args.directed_graph_out, 'wb') as pklout:
                pkl.dump(obj=nx_dir_graph, file=pklout)
    dir_node_set = set(nx_dir_graph.nodes())

    # Loop through the unique pairs
    dir_expl_count, im_expl_count = 0, 0
    explained_pairs = []  # Saves all explanations
    explained_neg_pairs = []  # Saves all explanations with correlation < 0
    unexplained = []  # Unexplained correlations
    npairs = len(uniq_pairs)

    # The explained nested dict: (1st key = subj, 2nd key = obj, 3rd key =
    # connection type)
    #
    # d[subj][obj] = {direct: [stmts/stmt hashes],
    #                 x_is_intermediary: [X],
    #                 x_is_downstream: [X],
    #                 x_is_upstream: [X]}
    #
    # Then in javascript you can for example do:
    # if SUBJ_is_subj_dict.obj.direct.length <-- should return zero if []
    #
    # Used to get: directed graph
    # 1. all nodes of directed graph -> 1st dropdown
    # 2. dir -> undir graph -> jsons to check all corr neighbors -> 2nd dropdown
    # 3. jsons to check if connection is direct or intermediary

    explained_nested_dict = nest_dict()

    # Open files to write text/latex output
    # with open(args.outbasename + '_connections_latex.tex', 'w') as f_con, \
    #         open(args.outbasename + '_neg_conn_latex.tex', 'w') as f_neg_c:

    logger.info('Looking for connections between %i pairs' % npairs)
    for pair in uniq_pairs:
        # Store bool(s) for found connection (either A-B or A-X-B)
        found = set()
        im_found = False  # Flag bool set for intermediate connections

        correlation, fmt_corr = None, None
        pl = list(pair)
        for li in pl:
            if _is_float(li):
                correlation = li
                fmt_corr = '{0:.04}'.format(correlation)
                break
        pl.remove(correlation)
        id1, id2 = pl

        # nested_dict_statements.get(id1).get(id2) raises AttributeError
        # if nested_dict_statements.get(id1) returns {}
        for subj, obj in itt.permutations((id1, id2), r=2):
            if _entry_exist(nested_dict_statements, subj, obj):
                # Get the statements
                stmts = nested_dict_statements[subj][obj]

                # Put in the explained nested dict
                explained_nested_dict[subj][obj]['direct'] = stmts

                logger.info('Found direct connection between %s and %s' % (
                    subj, obj))
                direct = True
                dir_expl_count += 1
                found.add(True)
                stmt_tuple = (subj, obj, correlation, 'direct', [])
                explained_pairs.append(stmt_tuple)

                if correlation < 0:
                    explained_neg_pairs.append(stmt_tuple)
                    # f_neg_c.write(output)

            # Checking 1. "pathway": A -> X -> B and B -> X -> A
            if subj in dir_node_set and obj in dir_node_set:
                dir_path_nodes = list(set(nx_dir_graph.succ[subj]) &
                                      set(nx_dir_graph.pred[obj]))
                if dir_path_nodes:
                    found.add(True)
                    im_found = True
                    logger.info('Found directed path of length 2 '
                                'between %s and %s' % (subj, obj))
                    explained_nested_dict[subj][obj]['x_is_intermediary']\
                        = dir_path_nodes
                    stmt_tuple = (subj, obj, correlation, 'pathway',
                                  dir_path_nodes)
                    explained_pairs.append(stmt_tuple)
                    if correlation < 0:
                        explained_neg_pairs.append(stmt_tuple)
                else:
                    found.add(False)

            else:
                found.add(False)

        if id1 in dir_node_set and id2 in dir_node_set:
            # Checking 2: share target/coregulator A -> X <- B
            downstream_share = list(set(nx_dir_graph.succ[id1]) &
                                    set(nx_dir_graph.succ[id2]))
            # Checking 3: No correlator A <- X -> B
            upstream_share = list(set(nx_dir_graph.pred[id1]) &
                                  set(nx_dir_graph.pred[id2]))
            if downstream_share:
                found.add(True)
                im_found = True
                stmt_tuple = (id1, id2, correlation, 'shared_target',
                              downstream_share)
                logger.info('Found downstream share: %s and %s share %i '
                            'targets' % (id1, id2, len(downstream_share)))
                explained_nested_dict[id1][id2]['x_is_downstream'] = \
                    downstream_share
                explained_nested_dict[id2][id1]['x_is_downstream'] = \
                    downstream_share
                explained_pairs.append(stmt_tuple)
                if correlation < 0:
                    explained_neg_pairs.append(stmt_tuple)

            if upstream_share:
                found.add(True)
                im_found = True
                stmt_tuple = (id1, id2, correlation, 'shared_upstream',
                              upstream_share)
                logger.info('Found upstream share: %s and %s are both directly '
                            'downstream of %i nodes' %
                            (id1, id2, len(upstream_share)))
                explained_nested_dict[id1][id2]['x_is_upstream'] = \
                    upstream_share
                explained_nested_dict[id2][id1]['x_is_upstream'] = \
                    upstream_share
                explained_pairs.append(stmt_tuple)
                if correlation < 0:
                    explained_neg_pairs.append(stmt_tuple)

            if not downstream_share and not upstream_share:
                found.add(False)
        else:
            found.add(False)

        if im_found:
            im_expl_count += 1

        # Make sure the connection types we didn't find are empty lists
        if any(found):
            for s, o in itt.permutations((id1, id2), r=2):
                # Direct
                if not _entry_exist(explained_nested_dict[s], o, 'direct'):
                    explained_nested_dict[s][o]['direct'] = []
                # x_is_intermediary
                if not _entry_exist(explained_nested_dict[s], o,
                                    'x_is_intermediary'):
                    explained_nested_dict[s][o]['x_is_intermediary'] = []
                # x_is_upstream
                if not _entry_exist(explained_nested_dict[s], o,
                                    'x_is_upstream'):
                    explained_nested_dict[s][o]['x_is_upstream'] = []
                # x_is_downstream
                if not _entry_exist(explained_nested_dict[s], o,
                                    'x_is_downstream'):
                    explained_nested_dict[s][o]['x_is_downstream'] = []

        # any(found) is True if at least one connection was found and
        # therefore "not any" is only True when no connection was found
        if not any(found):
            unexplained.append((id1, id2, correlation))
            if args.verbosity > 1:
                logger.info('No explainable path found between %s and '
                            '%s.' % (id1, id2))

    logger.info('Total unexplained: %i. Total explained: %i (%i direct and %i '
                'mediated by an intermediate node). Total number of pairs '
                'checked: %i' %
                (len(unexplained), len(explained_pairs), dir_expl_count,
                 im_expl_count, npairs))

    # Here create directed graph from explained nested dict
    nx_expl_dir_graph = dnf.nx_directed_graph_from_nested_dict_3layer(
        nest_d=explained_nested_dict)

    # 'explained_nodes' are used to produce first drop down
    explained_nodes = nx_expl_dir_graph.nodes()
    logger.info('Dumping json "explainable_ids.json" for first dropdown.')
    _dump_it_to_json('explainable_ids.json', explained_nodes)

    # Get undir graph and save each neighbor lookup as json for 2nd dropdown
    nx_expl_undir_graph = nx_expl_dir_graph.to_undirected()
    dnf.nx_undir_to_neighbor_lookup_json(expl_undir_graph=nx_expl_undir_graph)

    _dump_it_to_pickle(fname=args.outbasename+'_explained_nest_dict.pkl',
                       pyobj=explained_nested_dict)

    _dump_it_to_csv(fname=args.outbasename+'_expl_correlations.csv',
                    pyobj=explained_pairs)
    _dump_it_to_csv(fname=args.outbasename+'_expl_neg_correlations.csv',
                    pyobj=explained_neg_pairs)
    _dump_it_to_csv(fname=args.outbasename+'_unexpl_correlations.csv',
                    pyobj=unexplained)


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
        'only be between loaded gene set. If no gene set is loaded, '
        'this option has no effect.')
    parser.add_argument('-v', '--verbosity', action="count",
                        help='increase output verbosity (e.g., -vv is more '
                             'than -v)')
    parser.add_argument('-dgi', '--directed-graph-in', help='Load a'
        'precalculated directed graph of indra subjec/object network.')
    parser.add_argument('-dgo', '--directed-graph-out', help='Save the '
        'calculated directed graph of the indra statement network.')
    parser.add_argument('-ndi', '--nested-dict-in', help='Load precalculated '
        'nested dict of statements of the form d[subj][obj] = [stmts/stmt '
        'hashes].')
    parser.add_argument('-ndo', '--nested-dict-out', help='Save the '
        'calculated nested dict of statements')
    parser.add_argument('-lw', '--light-weight-stmts', help='A lightweight '
        'file with (hgnc_id1, hgnc_id2, stmt_type, stmt_hash) as columns.')
    parser.add_argument('-sti', '--statements-in', help='Loads a pickle file '
        'to use instead of quering a database. This option is unused if a '
        'lightweight file [-lw] is provided.')
    parser.add_argument('-sto', '--statements-out', help='Saves the used '
        'statements read from the database')
    parser.add_argument('-up', '--unique-depmap-pairs', help='Uses a '
        'previously saved file to read the unique pairs to read over. Useful '
        'if you are running the script on the full data with no filters.')
    parser.add_argument('-ll', type=float, default=0.3,
                        help='Lower limit CERES correlation score filter.')
    parser.add_argument('-ul', type=float, default=1.0,
                        help='Upper limit CERES correlation score filter.')
    a = parser.parse_args()

    with open('dep_map_script_log{}.log'.format(str(int(time()))), 'w',
              newline='\n') as f:
        f.write('Command line option - value\n')
        for arg in vars(a):
            f.write('{} - {}\n'.format(arg, getattr(a, arg)))
    main(a)
