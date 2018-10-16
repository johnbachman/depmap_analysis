import os
import csv
import json
import logging
import pandas as pd
import pickle as pkl
import networkx as nx
import argparse as ap
import itertools as itt
from numpy import float64
from time import time, strftime
from collections import defaultdict
from indra.tools import assemble_corpus as ac
import depmap_network_functions as dnf
from depmap_network_functions import create_nested_dict as nest_dict
# There are pickled files using "nest_dict" in their preserve import settings
# and we can therefore not use another name when using those files

# todo remove pdb import before merging PR
import pdb

logger = logging.getLogger('depmap_script')

# 1. no geneset -> use corr from full DepMap data, no filtering needed
# 2. geneset, not strict -> interaction has to contain at least one gene from
#    the loaded list
# 3. geneset loaded, strict -> each correlation can only contain genes from
#    the loaded data set


def _dump_it_to_pickle(fname, pyobj):
    with open(fname, 'wb') as po:
        pkl.dump(obj=pyobj, file=po)


def _dump_it_to_json(fname, pyobj):
    with open(fname, 'w') as json_out:
        json.dump(pyobj, json_out)


def _dump_it_to_csv(fname, pyobj, separator=',', header=None):
    if header:
        with open(fname, 'w') as fo:
            fo.write(','.join(header)+'\n')
    with open(fname, 'a', newline='') as csvf:
        wrtr = csv.writer(csvf, delimiter=separator)
        wrtr.writerows(pyobj)


def _dump_nest_dict_to_csv(fname, nested_dict, separator=',', header=None):
    if header:
        with open(fname, 'w') as fo:
            fo.write(separator.join(header)+'\n')
    with open(fname, 'a') as fo:
        for ok in nested_dict:
            for ik in nested_dict[ok]:
                cd = nested_dict[ok][ik]['correlations']
                fo.write('%s,%s,%f,%f\n' % (ok, ik, cd['crispr'], cd['rnai']))


def _pickle_open(file_path_to_pickle):
    with open(file_path_to_pickle, 'rb') as pi:
        return pkl.load(file=pi)


def _json_open(file_path_to_json):
    with open(file_path_to_json, 'r') as jo:
        return json.load(fp=jo)


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
    if not args.unique_depmap_crispr_pairs:
        unique_pairs_fpath = input('\nNo file of unique pair combinations for '
                                   'CRISPR data is provided.\nPlease provide '
                                   'one here or press [Enter] to proceed with '
                                   'the calculation of\nunique pairs (may take '
                                   'a while on full data set).\n')
        if unique_pairs_fpath:
            args.unique_depmap_crispr_pairs = unique_pairs_fpath
    if not args.unique_depmap_rnai_pairs:
        unique_pairs_fpath = input('\nNo file of unique pair combinations for '
                                   'RNAi data is provided.\nPlease provide '
                                   'one here or press [Enter] to proceed with '
                                   'the calculation of\nunique pairs (may take '
                                   'a while on full data set).\n')
        if unique_pairs_fpath:
            args.unique_depmap_rnai_pairs = unique_pairs_fpath

    # Check if belief dict is provided
    if not args.belief_score_dict and not args.nested_dict_in:
        logger.error('belief dict must be provided through the `-b ('
                     '--belief-score-dict)` argument if no nested dict '
                     'of statements with belief score is provided through the '
                     '`-ndi (--nested-dict-in)` argument.')
        raise FileNotFoundError

    # Prepare data (we need uniq_pairs to look for explainable interactions)
    filter_settings = {'margin': args.margin,
                       'filter_type': args.filter_type}

    args_dict = dnf.create_nested_dict()
    args_dict['crispr']['data'] = args.crispr_data_file
    args_dict['crispr']['corr'] = args.crispr_corr_file
    args_dict['crispr']['filter_gene_set'] = (args.geneset_file if
                                              args.geneset_file else [])
    args_dict['crispr']['ll'] = args.crispr_corr_range[0]
    args_dict['crispr']['ul'] = (args.crispr_corr_range[1] if len(
        args.crispr_corr_range) == 2 else 1.0)
    args_dict['crispr']['outbasename'] = args.outbasename+'_crispr'
    args_dict['crispr']['sigma'] = args.crispr_mean_sigma[0] if \
        args.crispr_mean_sigma else None
    args_dict['crispr']['mean'] = args.crispr_mean_sigma[1] if \
        args.crispr_mean_sigma else None
    args_dict['crispr']['strict'] = args.strict

    args_dict['rnai']['data'] = args.rnai_data_file
    args_dict['rnai']['corr'] = args.rnai_corr_file
    args_dict['rnai']['filter_gene_set'] = (args.geneset_file if
                                            args.geneset_file else [])
    args_dict['rnai']['ll'] = args.rnai_corr_range[0]
    args_dict['rnai']['ul'] = (args.rnai_corr_range[1] if len(
        args.rnai_corr_range) == 2 else 1.0)
    args_dict['rnai']['outbasename'] = args.outbasename+'_rnai'
    args_dict['rnai']['sigma'] = args.rnai_mean_sigma[0] if \
        args.rnai_mean_sigma else None
    args_dict['rnai']['mean'] = args.rnai_mean_sigma[1] if \
        args.rnai_mean_sigma else None
    args_dict['rnai']['strict'] = args.strict

    master_corr_dict, all_hgnc_ids, stats_dict = dnf.get_combined_correlations(
        dict_of_data_sets=args_dict, filter_settings=filter_settings)

    dnf._dump_master_corr_dict_to_pairs_in_csv(
        fname=args.outbasename+'_merged_corr_pairs.csv',
        nest_dict=master_corr_dict)

    if args.geneset_file:
        gene_filter_list = set(dnf._read_gene_set_file(
                gf=args_dict['crispr']['filter_gene_set'],
                data=args_dict['crispr']['data'])) & \
            set(dnf._read_gene_set_file(
                gf=args_dict['rnai']['filter_gene_set'],
                data=args_dict['crispr']['data']))
    else:
        gene_filter_list = None

    # Get dict of {hash: belief score}
    belief_dict = None  # ToDo use api to query belief scores if not loaded
    if args.belief_score_dict:
        if args.belief_score_dict.endswith('.json'):
            belief_dict = _json_open(args.belief_score_dict)
        elif args.belief_score_dict.endswith('.pkl'):
            belief_dict = _pickle_open(args.belief_score_dict)

    # LOADING INDRA STATEMENTS
    # Get statements from file or from database that contain any gene from
    # provided list as set unless you're already loading a pre-calculated
    # nested dict and/or precalculated directed graph.

    if not (args.light_weight_stmts or args.nested_dict_in):
        if args.statements_in:  # Get statments from file
            stmts_all = set(ac.load_statements(args.statements_in))
        # Use api to get statements. _NOT_ the same as querying for each ID
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
        nested_dict_statements = _pickle_open(args.nested_dict_in)
    else:
        if belief_dict is None:
            logger.error('belief dict must be provided through the `-b ('
                         '--belief-score-dict)` argument if no nested dict '
                         'of statements is provided through the `-ndi ('
                         '--nested-dict-in)` argument.')
            raise FileNotFoundError
        else:
            nested_dict_statements = dnf.dedupl_nested_dict_gen(stmts_all,
                                                            belief_dict)
        if args.nested_dict_out:
            _dump_it_to_pickle(fname=args.nested_dict_out,
                               pyobj=nested_dict_statements)

    # Get directed simple graph
    if args.directed_graph_in:
        with open(args.directed_graph_in, 'rb') as rpkl:
            nx_dir_graph = pkl.load(rpkl)
    else:
        # Create directed graph from statement dict
        nx_dir_graph = dnf.nx_directed_graph_from_nested_dict_2layer(
            nest_d=nested_dict_statements, belief_dict=belief_dict)
        # Save as pickle file
        if args.directed_graph_out:
            _dump_it_to_pickle(fname=args.directed_graph_out,
                               pyobj=nx_dir_graph)
    dir_node_set = set(nx_dir_graph.nodes)

    # LOOP THROUGH THE UNIQUE CORRELATION PAIRS, MATCH WITH INDRA NETWORK
    dir_expl_count, im_expl_count = 0, 0
    explained_pairs = []  # Saves all explanations
    explained_neg_pairs = []  # Saves all explanations with correlation < 0
    unexplained = []  # Unexplained correlations
    # npairs = dnf.rawincount(args.outbasename+'_merged_corr_pairs.csv')
    npairs = len(master_corr_dict)

    # The explained nested dict: (1st key = subj, 2nd key = obj, 3rd key =
    # connection type or correlation).
    #
    # directed: any A->B or B->A
    # undirected: any of complex, selfmodification, parent
    # x_is_intermediary: A->X->B or B->X->A
    # x_is_downstream: A->X<-B
    # x_is_upstream: A<-X->B
    #
    # d[subj][obj] = {correlation: {gene_set1: corr, gene_set2: corr, ...},
    #                 directed: [(stmt/stmt hash, belief score)],
    #                 undirected: [(stmt/stmt hash, belief score)],
    #                 x_is_intermediary: [(X, belief rank)],
    #                 x_is_downstream: [(X, belief rank)],
    #                 x_is_upstream: [(X, belief rank)]}
    #
    # Then in javascript you can for example do:
    # if SUBJ_is_subj_dict.obj.direct.length <-- should return zero if []
    #
    # Used to get: directed graph
    # 1. all nodes of directed graph -> 1st dropdown
    # 2. dir -> undir graph -> jsons to check all corr neighbors -> 2nd dropdown
    # 3. jsons to check if connection is direct or intermediary

    explained_nested_dict = dnf.create_nested_dict()

    # Open files to write text/latex output
    # with open(args.outbasename + '_connections_latex.tex', 'w') as f_con, \
    #         open(args.outbasename + '_neg_conn_latex.tex', 'w') as f_neg_c:

    logger.info('Looking for connections between %i pairs (length of dict)' %
                len(master_corr_dict))

    skipped = 0

    for outer_id, do in master_corr_dict.items():
        for inner_id, corr_dict in do.items():
            if len(corr_dict.keys()) == 0:
                skipped += 1
                if args.verbosity:
                    logger.info('Skipped outer_id=%s and inner_id=%s' %
                            (outer_id, inner_id))
                continue
            avg_corrs = []
            for set_name in corr_dict:
                avg_corrs.append(corr_dict[set_name])

            avg_corr = sum(avg_corrs)/len(avg_corrs)
            id1, id2 = outer_id, inner_id

            # Store bool(s) for found connection (either A-B or A-X-B)
            found = set()
            im_found = False  # Flag bool set for intermediate connections

            for subj, obj in itt.permutations((id1, id2), r=2):
                if dnf._entry_exist_dict(nested_dict_statements, subj, obj):
                    # Get the statements
                    stmts = nested_dict_statements[subj][obj]

                    # check if directed, put in the explained nested dict
                    dir_stmts, undir_stmts = dnf.get_directed(stmts)
                    explained_nested_dict[subj][obj]['directed'] = dir_stmts
                    explained_nested_dict[subj][obj]['undirected'] = undir_stmts

                    if args.verbosity:
                        logger.info('Found direct connection between %s and '
                                    '%s' % (subj, obj))
                    dir_expl_count += 1
                    found.add(True)
                    stmt_tuple = (subj, obj, corr_dict['crispr'],
                                  corr_dict['rnai'], 'direct', [])
                    explained_pairs.append(stmt_tuple)

                    if avg_corr < 0:
                        explained_neg_pairs.append(stmt_tuple)
                        # f_neg_c.write(output)

                # Checking 1. "pathway": A -> X -> B and B -> X -> A
                if subj in dir_node_set and obj in dir_node_set:
                    dir_path_nodes = list(set(nx_dir_graph.succ[subj]) &
                                          set(nx_dir_graph.pred[obj]))
                    if dir_path_nodes:
                        found.add(True)
                        im_found = True
                        if args.verbosity:
                            logger.info('Found directed path of length 2 '
                                        'between %s and %s' % (subj, obj))

                        dir_path_nodes_wb = dnf.rank_nodes(
                            node_list=dir_path_nodes,
                            nested_dict_stmts=nested_dict_statements,
                            gene_a=subj,
                            gene_b=obj,
                            x_type='x_is_intermediary')

                        explained_nested_dict[subj][obj]['x_is_intermediary']\
                            = dir_path_nodes_wb
                        stmt_tuple = (subj, obj, corr_dict['crispr'],
                                      corr_dict['rnai'], 'pathway',
                                      dir_path_nodes_wb)
                        explained_pairs.append(stmt_tuple)
                        if avg_corr < 0:
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
                    downstream_share_wb = dnf.rank_nodes(
                        node_list=downstream_share,
                        nested_dict_stmts=nested_dict_statements,
                        gene_a=id1,
                        gene_b=id2,
                        x_type='x_is_downstream')
                    stmt_tuple = (id1, id2, corr_dict['crispr'], 
                                  corr_dict['rnai'], 'shared_target',
                                  downstream_share_wb)
                    if args.verbosity:
                        logger.info('Found downstream share: %s and %s share '
                                    '%i targets' %
                                    (id1, id2, len(downstream_share)))
                    explained_nested_dict[id1][id2]['x_is_downstream'] = \
                        downstream_share_wb
                    explained_nested_dict[id2][id1]['x_is_downstream'] = \
                        downstream_share_wb
                    explained_pairs.append(stmt_tuple)
                    if avg_corr < 0:
                        explained_neg_pairs.append(stmt_tuple)

                if upstream_share:
                    found.add(True)
                    im_found = True
                    upstream_share_wb = dnf.rank_nodes(
                        node_list=upstream_share,
                        nested_dict_stmts=nested_dict_statements,
                        gene_a=id1,
                        gene_b=id2,
                        x_type='x_is_upstream')
                    stmt_tuple = (id1, id2, corr_dict['crispr'], 
                                  corr_dict['rnai'], 'shared_upstream',
                                  upstream_share_wb)
                    if args.verbosity:
                        logger.info('Found upstream share: %s and %s are both '
                                    'directly downstream of %i nodes' %
                                    (id1, id2, len(upstream_share)))
                    explained_nested_dict[id1][id2]['x_is_upstream'] = \
                        upstream_share_wb
                    explained_nested_dict[id2][id1]['x_is_upstream'] = \
                        upstream_share_wb
                    explained_pairs.append(stmt_tuple)
                    if avg_corr < 0:
                        explained_neg_pairs.append(stmt_tuple)

                if not downstream_share and not upstream_share:
                    found.add(False)
            else:
                found.add(False)

            # Count intermediate connections found
            if im_found:
                im_expl_count += 1

            # Make sure the connection types we didn't find are empty lists.
            # Also add correlation so it can be queried for at the same time
            # as the items for the second drop down.
            if any(found):
                for s, o in itt.permutations((id1, id2), r=2):
                    # Correlation
                    explained_nested_dict[s][o]['correlations'] = corr_dict
                    # Directed
                    if not dnf._entry_exist_dict(explained_nested_dict[s], o,
                                            'directed'):
                        explained_nested_dict[s][o]['directed'] = []
                    # Undirected
                    if not dnf._entry_exist_dict(explained_nested_dict[s], o,
                                            'undirected'):
                        explained_nested_dict[s][o]['undirected'] = []
                    # x_is_intermediary
                    if not dnf._entry_exist_dict(explained_nested_dict[s], o,
                                        'x_is_intermediary'):
                        explained_nested_dict[s][o]['x_is_intermediary'] = []
                    # x_is_upstream
                    if not dnf._entry_exist_dict(explained_nested_dict[s], o,
                                        'x_is_upstream'):
                        explained_nested_dict[s][o]['x_is_upstream'] = []
                    # x_is_downstream
                    if not dnf._entry_exist_dict(explained_nested_dict[s], o,
                                        'x_is_downstream'):
                        explained_nested_dict[s][o]['x_is_downstream'] = []

            # any(found) is True if at least one connection was found and
            # therefore "not any" is only True when no connection was found
            if not any(found):
                unexplained.append((id1, id2, corr_dict['crispr'],
                                    corr_dict['rnai']))
                if args.verbosity and args.verbosity > 1:
                    logger.info('No explainable path found between %s and '
                                '%s.' % (id1, id2))
    long_string = ''
    long_string += '-' * 63 + '\n'
    long_string += 'Summary:' + '\n'
    long_string += '> Skipped %i empty doublets in corr dict\n' % skipped
    long_string += '> Total unexplained: %i' % len(unexplained) + '\n'
    long_string += '> Total explained: %i,' % len(explained_pairs) + '\n'
    long_string += '> with %i direct and %i mediated by an intermediate ' \
                   'node.' % (dir_expl_count, im_expl_count) + '\n'
    long_string += '> Total number of pairs checked: %i' % npairs + '\n\n'
    long_string += 'Statistics:' + '\n\n'
    long_string += '  RNAi data ' + '\n'
    long_string += '  ----------' + '\n'
    long_string += '> mean: %f\n' % stats_dict['rnai']['mean']
    long_string += '> SD: %f\n' % stats_dict['rnai']['sigma']
    long_string += '  CRISPR data ' + '\n'
    long_string += '  ------------' + '\n'
    long_string += '> mean: %f\n' % stats_dict['crispr']['mean']
    long_string += '> SD: %f\n' % stats_dict['crispr']['sigma']
    long_string += '' + '\n'
    long_string += '-' * 63 + '\n\n'

    logger.info('\n' + long_string)

    # Here create directed graph from explained nested dict
    nx_expl_dir_graph = dnf.nx_directed_graph_from_nested_dict_3layer(
        nest_d=explained_nested_dict)

    # 'explained_nodes' are used to produce first drop down
    explained_nodes = list(nx_expl_dir_graph.nodes)
    logger.info('Dumping json "explainable_ids.json" for first dropdown.')
    _dump_it_to_json('explainable_ids.json', explained_nodes)

    # Get undir graph and save each neighbor lookup as json for 2nd dropdown
    nx_expl_undir_graph = nx_expl_dir_graph.to_undirected()
    dnf.nx_undir_to_neighbor_lookup_json(expl_undir_graph=nx_expl_undir_graph)

    _dump_nest_dict_to_csv(fname=args.outbasename+'_explained_pairs.csv',
                           nested_dict=explained_nested_dict,
                           header=['gene1', 'gene2',
                                   'crispr_corr', 'rnai_corr'])
    _dump_it_to_pickle(fname=args.outbasename+'_explained_nest_dict.pkl',
                       pyobj=explained_nested_dict)
    headers = ['subj', 'obj', 'crispr_corr', 'rnai_corr', 'type', 'X']
    _dump_it_to_csv(fname=args.outbasename+'_expl_correlations.csv',
                    pyobj=explained_pairs, header=headers)
    _dump_it_to_csv(fname=args.outbasename+'_expl_neg_correlations.csv',
                    pyobj=explained_neg_pairs, header=headers)
    _dump_it_to_csv(fname=args.outbasename+'_unexpl_correlations.csv',
                    pyobj=unexplained, header=headers[:-2])
    with open(args.outbasename+'_script_summary.txt', 'w') as fo:
        fo.write(long_string)
    return 0


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        description='Script to analyze and try to explain gene knockout data '
                    'from depmap.org. Minimum Working Example for running '
                    'script:    python depmap_script.py -cf <crispr gene '
                    'data csv file> -rf <rnai gene data csv file> '
                    '-o <output file name> Other good options are:    '
                    '-cc/-rc: precalculated correlation matrices in hdf '
                    'format    -ndi: nested dictionary of INDRA statements of '
                    'the format  `d[gene][gene] = [stmts/stmt hashes]`  OR    '
                    '-lw: a csv file with  `gene,gene,stmt type,stmt hash`  as '
                    'columns.')

    either_of = parser.add_argument_group('One of the following two arguments')
    required_args = parser.add_argument_group('Required Arguments')
    required_args.add_argument('-cf', '--crispr-data-file', required=True,
                               help='CRISPR gene dependency data in csv format')
    required_args.add_argument('-rf', '--rnai-data-file', required=True,
                               help='RNAi gene dependency data in csv format')
    either_of.add_argument('-b', '--belief-score-dict', help='Load a dict with '
        'stmt hash: belief score to be incorporated in the explainable '
        'network dict.')
    either_of.add_argument('-ndi', '--nested-dict-in',
                           help='Load precalculated nested dict of'
                                'statements of the form '
                                'd[subj][obj] = [(stmt/stmt hash, belief '
                                'score)].')
    parser.add_argument('-cc', '--crispr-corr-file',
                        help='Precalculated CRISPR correlations in h5 format')
    parser.add_argument('-rc', '--rnai-corr-file',
                        help='Precalculated RNAi correlations in h5 format')
    parser.add_argument('-g', '--geneset-file',
                        help='Filter to interactions with gene set data file.')
    parser.add_argument('--margin', type=float, default=1.0,
                        help='How large diff in terms of standard deviations '
                             'to accept between data sets when filtering for '
                             'correlations during merge. Default is 1 SD.')
    parser.add_argument('--filter-type', default='sigma-diff', type=str,
                        help='Type of filtering. Currently only supports '
                             '"sigma-diff"')
    parser.add_argument('-o', '--outbasename', default=str(int(time())),
                        help='Base name for outfiles. Default: UTC timestamp')
    parser.add_argument('-rec', '--recalc-crispr', action='store_true',
                        help='With \'-r\', recalculate full gene-gene '
                             'correlations of CRISPR data set.')
    parser.add_argument('-rer', '--recalc-rnai', action='store_true',
                        help='With \'-r\', recalculate full gene-gene '
                             'correlations of RNAi data set.')
    parser.add_argument('-s', '--strict', action='store_true',
                        help='With \'-s\', the correlations are restricted to '
        'only be between loaded gene set. If no gene set is loaded, '
        'this option has no effect.')
    parser.add_argument('-v', '--verbosity', action='count',
                        help='increase output verbosity (-vv is more '
                             'than -v)')
    parser.add_argument('-dgi', '--directed-graph-in', help='Load a'
        'precalculated directed graph of indra subject/object network.')
    parser.add_argument('-dgo', '--directed-graph-out', help='Save the '
        'calculated directed graph of the indra statement network.')
    parser.add_argument('-ndo', '--nested-dict-out', help='Save the '
        'calculated nested dict of statements')
    parser.add_argument('-lw', '--light-weight-stmts', help='A lightweight '
        'file with (hgnc_id1, hgnc_id2, stmt_type, stmt_hash) as columns.')
    parser.add_argument('-sti', '--statements-in', help='Loads a pickle file '
        'to use instead of quering a database. This option is unused if a '
        'lightweight file [-lw] is provided.')
    parser.add_argument('-sto', '--statements-out', help='Saves the used '
        'statements read from the database')
    parser.add_argument('-cup', '--unique-depmap-crispr-pairs', help='Uses a '
        'previously saved file to read the unique pairs to read over. Useful '
        'if you are running the script on the full data with no filters.')
    parser.add_argument('-rup', '--unique-depmap-rnai-pairs', help='Uses a '
        'previously saved file to read the unique pairs to read over. Useful '
        'if you are running the script on the full data with no filters.')
    parser.add_argument('-lls', action='store_true', help='Use exactly one '
        'standard deviation of full set as lower limit.')
    parser.add_argument('-crange', '--crispr-corr-range', default=[0.3],
                        type=float, nargs="+",
                        help='LOWER_LIM UPPER_LIM\nTwo decimal numbers '
                             'denoting the range of correlations to consider '
                             'in the crispr data.')
    parser.add_argument('-rrange', '--rnai-corr-range', default=[0.2],
                        type=float, nargs="+",
                        help='LOWER_LIM UPPER_LIM\nTwo decimal numbers '
                             'denoting the range of correlations to consider '
                             'in the rnai data.')
    parser.add_argument('-cstats', '--crispr-mean-sigma', type=float, nargs=2,
                        help='-rstats <mean> <stdev> | Provide a value of the '
                             'mean and standard deviation for the crispr data '
                             'instead of calculating it from the full data '
                             'set.')
    parser.add_argument('-rstats', '--rnai-mean-sigma', type=float, nargs=2,
                        help='-rstats <mean> <stdev> | Provide a value of the '
                             'mean and standard deviation for the rnai data '
                             'instead of calculating it from the full data '
                             'set.')
    a = parser.parse_args()

    with open(a.outbasename+'dep_map_script_log{}.log'.format(
            str(int(time()))), 'w', newline='\n') as f:
        f.write('Created on {}\n'.format(strftime('%Y %b %d, %H:%M:%S')))
        f.write('Execution path: {}\n\n'.format(os.getcwd()))
        f.write('Command line option : value\n---------------------------\n')
        for arg in vars(a):
            f.write('{} : {}\n'.format(arg, getattr(a, arg)))
    done = main(a)
    if done == 0 or done is None:
        logger.info('Script finished without errors')
