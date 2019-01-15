"""Script to analyze and try to explain gene knockout data from
depmap.org. Minimum Working Example for running script:

`python depmap_script.py -cf <crispr gene data csv file> -rf <rnai
gene data csv file> -o <output file name>`


Other important options are:

- -cc/-rc: precalculated correlation matrices in hdf format
- -ndi: nested dictionary of INDRA statements of the format
  `d[gene][gene] = [stmts/stmt hashes]` OR -lw: a csv file with `gene,gene,
  stmt type,stmt hash` as columns.
"""

import os
import sys
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
from random import choice as rnd_choice
from indra.tools import assemble_corpus as ac
import depmap_network_functions as dnf
from depmap_network_functions import create_nested_dict as nest_dict
# There are pickled files using "nest_dict" in their preserved import settings
# and we can therefore not use another name when using those files until we
# create new pickle files

logger = logging.getLogger('DepMap Script')

# 1. no geneset -> use corr from full DepMap data, no filtering needed
# 2. geneset, not strict -> interaction has to contain at least one gene from
#    the loaded list
# 3. geneset loaded, strict -> each correlation can only contain genes from
#    the loaded data set

global any_expl, any_expl_not_sr, common_parent, tuple_dir_expl_count, \
    both_dir_expl_count, tuple_im_expl_count, both_im_dir_expl_count, \
    tuple_im_st_expl_count, tuple_im_sr_expl_count, \
    tuple_sr_expl_only_count, explained_pairs, unexplained, \
    explained_nested_dict, id1, id2, nested_dict_statements, dataset_dict, \
    dir_node_set, nx_dir_graph, uninteresting_set, uninteresting


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
                cd = nested_dict[ok][ik]['meta_data']
                fo.write('%s,%s,%s,%s\n' %
                         (ok, ik, str(cd['crispr']), str(cd['rnai'])))


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


def _rnd_pair_gen(rgs):
    return rnd_choice(rgs), rnd_choice(rgs)


def _parse_cell_filter(cl_file, id2depmapid_pickle=None, namespace='CCLE_Name'):
    logger.info('Parsing cell lines...')
    cell_lines = []
    with open(cl_file, 'r') as fi:
        first_line = fi.readline()
        if ',' in first_line:
            col_separator = ','
        else:
            col_separator = None
    with open(cl_file, 'r') as fi:
        if id2depmapid_pickle:
            id2depmapid_dict = _pickle_open(id2depmapid_pickle)
            assert namespace in id2depmapid_dict.keys()
            for n, cl in enumerate(fi.readlines()):
                cl_name = cl.split(sep=col_separator)[0]
                if n == 0:
                    continue
                else:
                    try:
                        cell_lines.append(id2depmapid_dict[namespace][cl_name])
                    except KeyError:
                        logger.warning('Could not find mapping for %s' %
                                       cl_name)
                        continue
        else:
            for n, cl in enumerate(fi.readlines()):
                cl_name = cl.split(sep=col_separator)[0]
                if n == 0:
                    continue
                else:
                    cell_lines.append(cl_name)

    return cell_lines


def _parse_uninteresting_genes(gene_set_file, check_column):
    logger.info('Parsing uninteresting genes assuming column %s can be '
                'mapped to boolean' % check_column)
    gene_df = pd.read_csv(gene_set_file, header=0)
    if check_column in gene_df.columns.values:
        try:
            genes = set(gene_df['gene'][gene_df[check_column].apply(
                bool)].values)
        except KeyError:
            try:
                genes = set(gene_df['genes'][gene_df[check_column].apply(
                    bool)].values)
            except KeyError as exception:
                logger.warning(exception)
                sys.exit('No column with name "gene" or "genes" found!')
    else:
        sys.exit('Cannot find column %s!' % check_column)
    logger.info('Loaded %i uninteresting genes' % len(genes))
    return genes


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


def _arg_dict(args_struct):
    args_dict = dnf.create_nested_dict()

    if not args_struct.crispr_data_file and not args_struct.rnai_data_file \
            and not args_struct.brca_dependencies \
            and not args_struct.sampling_gene_file:
        sys.exit('Must provide at least one data set or a gene dependency '
                 'dataset or a list of genes for random sampling!')

    # CRISPR
    if args_struct.crispr_data_file:
        args_dict['crispr']['data'] = args_struct.crispr_data_file
        args_dict['crispr']['corr'] = args_struct.crispr_corr_file
        args_dict['crispr']['ll'] = max(args_struct.crispr_corr_range[0], 0.0)
        args_dict['crispr']['ul'] = (args_struct.crispr_corr_range[1]
            if len(args_struct.crispr_corr_range) == 2 else None)
        args_dict['crispr']['max_pairs'] = args_struct.max_pairs
        args_dict['crispr']['mean'] = args_struct.crispr_mean_sigma[0] if \
            args_struct.crispr_mean_sigma else None
        args_dict['crispr']['sigma'] = args_struct.crispr_mean_sigma[1] if \
            args_struct.crispr_mean_sigma else None

    # RNAi
    if args_struct.rnai_data_file:
        args_dict['rnai']['data'] = args_struct.rnai_data_file
        args_dict['rnai']['corr'] = args_struct.rnai_corr_file
        args_dict['rnai']['ll'] = max(args_struct.rnai_corr_range[0], 0.0)
        args_dict['rnai']['ul'] = (args_struct.rnai_corr_range[1]
            if len(args_struct.rnai_corr_range) == 2 else None)
        args_dict['rnai']['max_pairs'] = args_struct.max_pairs
        args_dict['rnai']['mean'] = args_struct.rnai_mean_sigma[0] if \
            args_struct.rnai_mean_sigma else None
        args_dict['rnai']['sigma'] = args_struct.rnai_mean_sigma[1] if \
            args_struct.rnai_mean_sigma else None

    # RANDOM SAMPLING
    if args_struct.sampling_gene_file:
        if not args_struct.max_pairs:
            sys.exit('Must specify a maximum number of pairs for random '
                     'sampling')
        args_dict['sampling_gene_file'] = args_struct.sampling_gene_file

    args_dict.default_factory = None
    return args_dict


def loop_body(args):

    global any_expl, any_expl_not_sr, common_parent, tuple_dir_expl_count, \
        both_dir_expl_count, tuple_im_expl_count, both_im_dir_expl_count, \
        tuple_im_st_expl_count, tuple_im_sr_expl_count, \
        tuple_sr_expl_only_count, explained_pairs, unexplained, \
        explained_nested_dict, id1, id2, nested_dict_statements, \
        dataset_dict, dir_node_set, nx_dir_graph, uninteresting_set, \
        uninteresting

    # Store bool(s) for found connection (either A-B or A-X-B)
    found = set()  # Flag anythin found
    dir_found = False  # Flag direct/complex connection
    im_found = False  # Flag intermediate connections
    sr_found = False  # Flag shared regulator connection
    non_sr_found = False  # Flag any non shared regulator connection

    for subj, obj in itt.permutations((id1, id2), r=2):
        if dnf._entry_exist_dict(nested_dict_statements, subj, obj):
            both_dir_expl_count += 1

            # Get the statements
            stmts = nested_dict_statements[subj][obj]

            # check if directed, put in the explained nested dict
            dir_stmts, undir_stmts = dnf.get_directed(stmts)
            explained_nested_dict[subj][obj]['directed'] = dir_stmts
            explained_nested_dict[subj][obj]['undirected'] = undir_stmts

            if args.verbosity:
                logger.info('Found direct connection between %s and '
                            '%s' % (subj, obj))
            found.add(True)
            dir_found = True
            non_sr_found = True
            stmt_tuple = (subj, obj, 'direct', [], dataset_dict)
            explained_pairs.append(stmt_tuple)

        # Checking 1. "pathway": A -> X -> B and B -> X -> A
        if subj in dir_node_set and obj in dir_node_set:
            dir_path_nodes = list(set(nx_dir_graph.succ[subj]) &
                                  set(nx_dir_graph.pred[obj]))
            if dir_path_nodes:
                found.add(True)
                im_found = True
                non_sr_found = True
                both_im_dir_expl_count += 1
                if args.verbosity:
                    logger.info('Found directed path of length 2 '
                                'between %s and %s' % (subj, obj))

                dir_path_nodes_wb = dnf.rank_nodes(
                    node_list=dir_path_nodes,
                    nested_dict_stmts=nested_dict_statements,
                    gene_a=subj,
                    gene_b=obj,
                    x_type='x_is_intermediary')

                explained_nested_dict[subj][obj]['x_is_intermediary'] \
                    = dir_path_nodes_wb
                stmt_tuple = (subj, obj, 'pathway',
                              dir_path_nodes_wb, dataset_dict)
                explained_pairs.append(stmt_tuple)

            else:
                found.add(False)

        else:
            found.add(False)

    # Check common parent
    if dnf.has_common_parent(id1=id1, id2=id2):
        common_parent += 1
        non_sr_found = True
        found.add(True)
        parents = list(dnf.common_parent(id1=id1, id2=id2))
        explained_nested_dict[id1][id2]['common_parents'] = parents
        explained_nested_dict[id2][id1]['common_parents'] = parents
    else:
        found.add(False)

    # Check if both A and B are in list of "uninteresting genes"
    if id1 in uninteresting_set and id2 in uninteresting_set:
        explained_nested_dict[id1][id2]['uninteresting'] = True
        explained_nested_dict[id2][id1]['uninteresting'] = True
        uninteresting += 1
        found.add(True)
    else:
        explained_nested_dict[id1][id2]['uninteresting'] = False
        explained_nested_dict[id2][id1]['uninteresting'] = False
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
            non_sr_found = True
            tuple_im_st_expl_count += 1
            downstream_share_wb = dnf.rank_nodes(
                node_list=downstream_share,
                nested_dict_stmts=nested_dict_statements,
                gene_a=id1,
                gene_b=id2,
                x_type='x_is_downstream')
            stmt_tuple = (id1, id2, 'shared_target',
                          downstream_share_wb, dataset_dict)
            if args.verbosity:
                logger.info('Found downstream share: %s and %s share '
                            '%i targets' %
                            (id1, id2, len(downstream_share)))
            explained_nested_dict[id1][id2]['x_is_downstream'] = \
                downstream_share_wb
            explained_nested_dict[id2][id1]['x_is_downstream'] = \
                downstream_share_wb
            explained_pairs.append(stmt_tuple)

        if upstream_share:
            found.add(True)
            im_found = True
            sr_found = True
            tuple_im_sr_expl_count += 1
            upstream_share_wb = dnf.rank_nodes(
                node_list=upstream_share,
                nested_dict_stmts=nested_dict_statements,
                gene_a=id1,
                gene_b=id2,
                x_type='x_is_upstream')
            stmt_tuple = (id1, id2, 'shared_upstream',
                          upstream_share_wb, dataset_dict)
            if args.verbosity:
                logger.info('Found upstream share: %s and %s are both '
                            'directly downstream of %i nodes' %
                            (id1, id2, len(upstream_share)))
            explained_nested_dict[id1][id2]['x_is_upstream'] = \
                upstream_share_wb
            explained_nested_dict[id2][id1]['x_is_upstream'] = \
                upstream_share_wb
            explained_pairs.append(stmt_tuple)

        if not downstream_share and not upstream_share:
            found.add(False)
    else:
        found.add(False)

    # Make sure the connection types we didn't find are empty lists.
    # Also add correlation so it can be queried for at the same time
    # as the items for the second drop down.
    if any(found):
        # Any explanation found
        any_expl += 1

        # Count A-B or B-A connections found per set(A,B)
        if dir_found:
            tuple_dir_expl_count += 1

        # Count A-X-B connections found per set(A,B)
        if im_found:
            tuple_im_expl_count += 1

        # Count non shared regulators found
        if non_sr_found:
            any_expl_not_sr += 1

        # Count only shared regulators found
        if sr_found and not non_sr_found:
            tuple_sr_expl_only_count += 1

        for s, o in itt.permutations((id1, id2), r=2):
            # Correlation/meta data
            explained_nested_dict[s][o]['meta_data'] = dataset_dict
            # common_parents
            if not dnf._entry_exist_dict(explained_nested_dict[s], o,
                                         'common_parents'):
                explained_nested_dict[s][o]['common_parents'] = []
            # directed
            if not dnf._entry_exist_dict(explained_nested_dict[s], o,
                                         'directed'):
                explained_nested_dict[s][o]['directed'] = []
            # undirected
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

    # any(found) == True if at least one connection was found
    # not any(found) == True is only True when no connection was found
    if not any(found):
        unexplained.append((id1, id2, dataset_dict))
        if args.verbosity and args.verbosity > 1:
            logger.info('No explainable path found between %s and '
                        '%s.' % (id1, id2))


def main(args):

    global any_expl, any_expl_not_sr, common_parent, tuple_dir_expl_count, \
        both_dir_expl_count, tuple_im_expl_count, both_im_dir_expl_count, \
        tuple_im_st_expl_count, tuple_im_sr_expl_count, \
        tuple_sr_expl_only_count, explained_pairs, unexplained, \
        explained_nested_dict, id1, id2, nested_dict_statements, dataset_dict, \
        avg_corr, dir_node_set, nx_dir_graph, uninteresting_set, uninteresting

    if args.cell_line_filter and not len(args.cell_line_filter) > 2:
        logger.info('Filtering to provided cell lines in correlation '
                    'calculations.')
        cell_lines = _parse_cell_filter(*args.cell_line_filter)
        assert len(cell_lines) > 0
    elif args.cell_line_filter and len(args.cell_line_filter) > 2:
        sys.exit('Argument --cell-line-filter only takes one or two arguments')
    # No cell line dictionary and rnai data and filtering is requested
    elif args.cell_line_filter and len(args.cell_line_filter) == 1 and \
            args.rnai_data_file:
        sys.exit('Need a translation dictionary if RNAi data is provided and '
                 'filter is requested')
    else:
        # Should be empty only when --cell-line-filter is not provided
        logger.info('No cell line filter provided. Using all cell lines in '
                    'correlation calculations.')
        cell_lines = []

    # Parse "uninteresting genes"
    if args.uninteresting and len(args.uninteresting) == 2:
        uninteresting_set = _parse_uninteresting_genes(
            gene_set_file=args.uninteresting[0],
            check_column=args.uninteresting[1])
        logger.info('Loading "uninteresting pairs."')
    elif args.uninteresting and len(args.uninteresting) != 2:
        sys.exit('Argument --uninteresting takes exactly two arguments: '
                 '--uninteresting <file> <column name>')

    # Check if belief dict is provided
    if not args.belief_score_dict and not args.nested_dict_in:
        logger.error('Belief dict must be provided through the `-b ('
                     '--belief-score-dict)` argument if no nested dict '
                     'of statements with belief score is provided through the '
                     '`-ndi (--nested-dict-in)` argument.')
        raise FileNotFoundError

    # Get dict of {hash: belief score}
    belief_dict = None  # ToDo use api to query belief scores if not loaded
    if args.belief_score_dict:
        if args.belief_score_dict.endswith('.json'):
            belief_dict = _json_open(args.belief_score_dict)
        elif args.belief_score_dict.endswith('.pkl'):
            belief_dict = _pickle_open(args.belief_score_dict)

    args_dict = _arg_dict(args)
    npairs = 0

    filter_settings = {'gene_set_filter': args.gene_set_filter,
                       'strict': args.strict,
                       'cell_line_filter': cell_lines,
                       'cell_line_translation_dict': _pickle_open(
                                                     args.cell_line_filter[1])
                       if args.cell_line_filter and len(args.cell_line_filter)
                       == 2 else None,
                       'margin': args.margin,
                       'filter_type': (args.filter_type
                                       if args.filter_type in ['sigma-diff',
                                                               'corr-corr-corr']
                                       else None)
                       }

    output_settings = {'dump_unique_pairs': args.dump_unique_pairs,
                       'outbasename': args.outbasename}

    # Parse CRISPR and/or RNAi data
    if args_dict.get('crispr') or args_dict.get('rnai'):
        if not filter_settings['filter_type'] and \
            args.crispr_data_file and \
                args.rnai_data_file:
            logger.info('No merge filter set. Output will be intersection of '
                        'the two data sets.')

        master_corr_dict, all_hgnc_ids, stats_dict = \
            dnf.get_combined_correlations(dict_of_data_sets=args_dict,
                                          filter_settings=filter_settings,
                                          output_settings=output_settings)

        # Count pairs in merged correlation dict and dum it
        npairs = dnf._dump_master_corr_dict_to_pairs_in_csv(
            fname=args.outbasename+'_merged_corr_pairs.csv',
            nest_dict=master_corr_dict)

        if args.gene_set_filter:
            gene_filter_list = None
            if args_dict.get('crispr') and not args_dict.get('rnai'):
                gene_filter_list = dnf._read_gene_set_file(
                    gf=filter_settings['gene_set_filter'],
                    data=pd.read_csv(args_dict['crispr']['data'],
                                         index_col=0, header=0))
            elif args_dict.get('rnai') and not args_dict.get('crispr'):
                gene_filter_list = dnf._read_gene_set_file(
                        gf=filter_settings['gene_set_filter'],
                        data=pd.read_csv(args_dict['rnai']['data'],
                                         index_col=0, header=0))
            elif args_dict.get('crispr') and args_dict.get('rnai'):
                gene_filter_list = \
                    set(dnf._read_gene_set_file(
                        gf=filter_settings['gene_set_filter'],
                        data=pd.read_csv(args_dict['crispr']['data'],
                                         index_col=0, header=0))) & \
                    set(dnf._read_gene_set_file(
                        gf=filter_settings['gene_set_filter'],
                        data=pd.read_csv(args_dict['rnai']['data'],
                                         index_col=0, header=0)))
            assert gene_filter_list is not None

        else:
            gene_filter_list = None
    else:
        stats_dict = None

    # LOADING INDRA STATEMENTS
    # Get statements from file or from database that contain any gene from
    # provided list as set unless you're already loading a pre-calculated
    # nested dict and/or precalculated directed graph.

    if not (args.light_weight_stmts or args.nested_dict_in):
        if args.statements_in:  # Get statments from file
            stmts_all = set(ac.load_statements(args.statements_in))
        # Use api to get statements. _NOT_ the same as querying for each ID
        else:
            if args.gene_set_filter:
                stmts_all = dnf.dbc_load_statements(gene_filter_list)
            else:
                # if there is no gene set file, restrict to gene ids in
                # input data
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
    any_expl = 0  # Count if any explanation per (A,B) correlation found
    # Count any explanation per (A,B) found, excluding shared regulator
    any_expl_not_sr = 0
    common_parent = 0  # Count if common parent found per set(A,B)
    uninteresting = 0  # Count pairs removed because they were "uninteresting"
    tuple_dir_expl_count = 0  # Count A-B/B-A as one per set(A,B)
    both_dir_expl_count = 0  # Count A-B and B-A separately per set(A,B)
    tuple_im_expl_count = 0  # Count any A->X->B,B->X->A as one per set(A,B)
    both_im_dir_expl_count = 0  # Count A->X->B,B->X->A separately per set(A,B)
    tuple_im_st_expl_count = 0  # Count if shared target found per set(A,B)
    tuple_im_sr_expl_count = 0  # Count if shared regulator found per set(A,B)
    tuple_sr_expl_only_count = 0  # Count if only shared regulator found
    explained_pairs = []  # Saves all explanations
    unexplained = []  # Unexplained correlations
    skipped = 0

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
    #                 common_parents: [list of parents]
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

    # Using the following loop structure for counter variables:
    # a = 2
    # def for_loop_body():
    #     global a
    #     a += 1
    # # Then loop like:
    # if dict:
    #     for pairs in dict:
    #         for_loop_body(args)
    # elif random:
    #     for random pair:
    #         for_loop_body(args)

    explained_nested_dict = dnf.create_nested_dict()

    # Loop rnai and/or crispr only
    if args_dict.get('rnai') or args_dict.get('crispr') and \
            not args.brca_dependencies:
        logger.info('Gene pairs generated from DepMap knockout screening data '
                    'sets')
        logger.info('Looking for connections between %i pairs' % (
            npairs if npairs > 0 else args.max_pairs)
        )
        for outer_id, do in master_corr_dict.items():
            for inner_id, dataset_dict in do.items():
                if len(dataset_dict.keys()) == 0:
                    skipped += 1
                    if args.verbosity:
                        logger.info('Skipped outer_id=%s and inner_id=%s' %
                                (outer_id, inner_id))
                    continue

                id1, id2 = outer_id, inner_id
                loop_body(args)

    # Loop rnai and/or crispr AND BRCA cell line dependencies
    elif args_dict.get('rnai') or args_dict.get('crispr') and \
            args.brca_dependencies:
        logger.info('Gene pairs generated from combined knockout screens. '
                    'Output data will incluide BRCA cell line dependency\n'
                    'data as well as correlation data from knockout screens.')
        logger.info('Looking for connections between %i pairs' % (
            npairs if npairs > 0 else args.max_pairs)
        )

        # Load BRCA dependency data
        brca_data_set = pd.read_csv(args.brca_dependencies, header=0)
        depend_in_breast_genes = brca_data_set.drop(
            axis=1, labels=['Url Label', 'Type'])[brca_data_set['Type'] ==
                                                  'gene']
        genes = set(depend_in_breast_genes['Gene/Compound'].values)

        for outer_id, do in master_corr_dict.items():
            for inner_id, knockout_dict in do.items():
                if len(knockout_dict.keys()) == 0:
                    skipped += 1
                    if args.verbosity:
                        logger.info('Skipped outer_id=%s and inner_id=%s' %
                                (outer_id, inner_id))
                    continue

                id1, id2 = outer_id, inner_id
                dataset_dict = {}
                gene1_data = []
                gene2_data = []

                # Get BRCA dep data
                if id1 in genes:
                    for row in depend_in_breast_genes[
                        depend_in_breast_genes[
                            'Gene/Compound'] == id1].iterrows():
                        gene1_data.append((row[1]['Dataset'],
                                           row[1]['T-Statistic'],
                                           row[1]['P-Value']))
                if id2 in genes:
                    for row in depend_in_breast_genes[
                        depend_in_breast_genes[
                            'Gene/Compound'] == id2].iterrows():
                        gene2_data.append((row[1]['Dataset'],
                                           row[1]['T-Statistic'],
                                           row[1]['P-Value']))

                dataset_dict[id1] = gene1_data
                dataset_dict[id2] = gene2_data

                dataset_dict['crispr'] = (knockout_dict['crispr'] if
                                          knockout_dict.get('crispr')
                                          else None),
                dataset_dict['rnai'] = (knockout_dict['rnai'] if
                                        knockout_dict.get('rnai')
                                        else None)

                if id1 not in genes and id2 not in genes:
                    dataset_dict = knockout_dict

                # Run loop body
                loop_body(args)

    # loop brca dependency ONLY
    elif args.brca_dependencies and not \
            (args_dict.get('rnai') or args_dict.get('crispr')):
        logger.info('Gene pairs generated from BRCA gene enrichment data only.')
        brca_data_set = pd.read_csv(args.brca_dependencies, header=0)
        depend_in_breast_genes = brca_data_set.drop(
            axis=1, labels=['Url Label', 'Type'])[brca_data_set['Type'] ==
                                                  'gene']
        genes = set(depend_in_breast_genes['Gene/Compound'].values)
        npairs = len(list(itt.combinations(genes, 2)))
        logger.info('Looking for connections between %i pairs' % (
            npairs if npairs > 0 else args.max_pairs)
        )
        for id1, id2 in itt.combinations(genes, 2):
            gene1_data = []
            gene2_data = []
            # For each non-diagonal pair in file, insert in dataset_dict:
            # geneA, geneB,
            # dataset for A, dataset for B,
            # T-stat for A, T-stat for B,
            # P-value for A, P-value
            for row in depend_in_breast_genes[
                 depend_in_breast_genes['Gene/Compound'] == id1].iterrows():
                gene1_data.append((row[1]['Dataset'],
                                   row[1]['T-Statistic'],
                                   row[1]['P-Value']))

            for row in depend_in_breast_genes[
                 depend_in_breast_genes['Gene/Compound'] == id2].iterrows():
                gene2_data.append((row[1]['Dataset'],
                                   row[1]['T-Statistic'],
                                   row[1]['P-Value']))
            # dataset_dict = {id1:
            #                 [(dataset1, T-stat1, P-value1),
            #                  (dataset2, T-stat2, P-value2)],
            #                 id2:
            #                  [(..., ...)],
            #                  ...}
            dataset_dict = {id1: gene1_data, id2: gene2_data}
            loop_body(args)

    # loop random pairs from data set
    elif args_dict.get('sampling_gene_file'):
        logger.info('Gene pairs generated at random from %s' %
                    args_dict['sampling_gene_file'])
        with open(args_dict['sampling_gene_file'], 'r') as fi:
            rnd_gene_set = [l.strip() for l in fi.readlines()]

        npairs = args.max_pairs
        dataset_dict = None
        logger.info('Looking for connections between %i pairs' % (
            npairs if npairs > 0 else args.max_pairs)
        )
        for _ in range(npairs):
            id1, id2 = _rnd_pair_gen(rnd_gene_set)
            assert not isinstance(id1, list)
            loop_body(args)

    long_string = ''
    long_string += '-' * 63 + '\n'
    long_string += 'Summary for matching INDRA network to correlation pairs:'\
                   + '\n\n'
    long_string += '> Total number of correlation pairs checked: %i' % npairs\
                   + '\n'
    if args.verbosity:
        long_string += '> Skipped %i empty doublets in corr dict\n' % skipped

    long_string += '> Total correlations unexplained: %i' % len(unexplained)\
                   + '\n'
    long_string += '> Total correlations explained: %i' % any_expl + '\n'
    long_string += '> Total correlations explained, excluding shared ' \
                   'regulator: %i' % any_expl_not_sr + '\n'
    long_string += '>    %i correlations have an explanation involving a ' \
                   'common parent' % common_parent + '\n'
    if args.uninteresting:
        long_string += '>    %i gene pairs were considered explained as part ' \
                       'of the "uninterested genes"' % uninteresting + '\n'
    long_string += '>    %i correlations have an explanation involving a ' \
                   'direct connection' % tuple_dir_expl_count + \
                   '\n'
    long_string += '>    %i direct connections found (count A-B and B-A ' \
                   'separately, including complexes)' % both_dir_expl_count + \
                   '\n'
    long_string += '>    %i correlations have an explanation ' \
                   'involving an intermediate node (A-X-B).' \
                   % tuple_im_expl_count + '\n'
    long_string += '>    %i A->X->B or B->X->A connections found (one count ' \
                   'per direction)' % both_im_dir_expl_count + '\n'
    long_string += '>    %i correlations have an explanation involving a ' \
                   'shared target (A->X<-B)' % tuple_im_st_expl_count + '\n'
    long_string += '>    %i correlations have an explanation involving a ' \
                   'shared regulator (A<-X->B)' % tuple_im_sr_expl_count + '\n'
    long_string += '>    %i correlations have shared regulator as only ' \
                   'explanation' % tuple_sr_expl_only_count + '\n\n'

    long_string += 'Statistics of input data:' + '\n\n'
    if stats_dict and stats_dict.get('rnai'):
        long_string += '  RNAi data ' + '\n'
        long_string += '  ----------' + '\n'
        long_string += '> mean: %f\n' % stats_dict['rnai']['mean']
        long_string += '> SD: %f\n' % stats_dict['rnai']['sigma']
        long_string += '> lower bound: %.3f*SD = %.4f\n' % (
            args_dict['rnai']['ll'],
            args_dict['rnai']['ll']*stats_dict['rnai']['sigma']
        )
        if args_dict['rnai']['ul']:
            long_string += '> upper bound: %.3f*SD = %.4f\n\n' % (
                args_dict['rnai']['ul'],
                args_dict['rnai']['ul'] * stats_dict['rnai']['sigma']
            )
    if stats_dict and stats_dict.get('crispr'):
        long_string += '  CRISPR data ' + '\n'
        long_string += '  ------------' + '\n'
        long_string += '> mean: %f\n' % stats_dict['crispr']['mean']
        long_string += '> SD: %f\n' % stats_dict['crispr']['sigma']
        long_string += '> lower bound: %.3f*SD = %.4f\n' % (
            args_dict['crispr']['ll'],
            args_dict['crispr']['ll']*stats_dict['crispr']['sigma']
        )
        if args_dict['crispr']['ul']:
            long_string += '> upper bound: %.3f*SD = %.4f\n\n' % (
                args_dict['crispr']['ul'],
                args_dict['crispr']['ul'] * stats_dict['crispr']['sigma']
            )
    long_string += '-' * 63 + '\n\n'

    logger.info('\n' + long_string)

    # Here create directed graph from explained nested dict
    nx_expl_dir_graph = dnf.nx_directed_graph_from_nested_dict_3layer(
        nest_d=explained_nested_dict)

    if not args.no_web_files:
        # 'explained_nodes' are used to produce first drop down
        explained_nodes = list(nx_expl_dir_graph.nodes)
        logger.info('Dumping json "explainable_ids.json" for first dropdown.')
        _dump_it_to_json(args.outbasename+'explainable_ids.json',
                         explained_nodes)

        # Get undir graph and save each neighbor lookup as json for 2nd dropdown
        nx_expl_undir_graph = nx_expl_dir_graph.to_undirected()
        dnf.nx_undir_to_neighbor_lookup_json(
            expl_undir_graph=nx_expl_undir_graph, outbasename=args.outbasename)

        _dump_nest_dict_to_csv(fname=args.outbasename+'_explained_pairs.csv',
                               nested_dict=explained_nested_dict,
                               header=['gene1', 'gene2', 'meta_data'])

    _dump_it_to_pickle(fname=args.outbasename+'_explained_nest_dict.pkl',
                       pyobj=explained_nested_dict)
    headers = ['subj', 'obj', 'type', 'X', 'meta_data']
    _dump_it_to_csv(fname=args.outbasename+'_expl_correlations.csv',
                    pyobj=explained_pairs, header=headers)
    _dump_it_to_csv(fname=args.outbasename+'_unexpl_correlations.csv',
                    pyobj=unexplained, header=headers[:-2])
    with open(args.outbasename+'_script_summary.txt', 'w') as fo:
        fo.write(long_string)
    return 0


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        description=
        """Script to analyze and try to explain gene knockout data from 
        depmap.org. Minimum Working Example for running script:
        
        
        `python depmap_script.py -cf <crispr gene data csv file> -rf <rnai 
        gene data csv file> -o <output file name>`
        
        Other important options are:
        
        - -cc/-rc: precalculated correlation matrices in hdf format
        - Pre-saved statement structures: -ndi: nested dictionary of INDRA 
          statements of the format `d[gene][gene] = [stmts/stmt hashes]` OR 
          -lw: a csv file with `gene,gene,stmt type,stmt hash` as columns.
        """
    )

    either_of = parser.add_argument_group('One of the following arguments')
    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument('-cf', '--crispr-data-file',
                               help='CRISPR gene dependency data in csv format')
    required_args.add_argument('-rf', '--rnai-data-file',
                               help='RNAi gene dependency data in csv format')
    required_args.add_argument('-brca', '--brca-dependencies',
                               help='`Dependencies enriched in breast` from '
        'https://depmap.org/portal/context/breast')
    required_args.add_argument('-sampl', '--sampling-gene-file',
                               help='A file containing hgnc symbols that will '
        'have pairs sampled from it at random. The pairs are then used to '
        'generate statistics over the explanation ratio for random pairs.')
    either_of.add_argument('-b', '--belief-score-dict', help='Load a dict with '
        'stmt hash: belief score to be incorporated in the explainable '
        'network dict.')
    either_of.add_argument('-ndi', '--nested-dict-in', help='Load a nested '
        'dict of statements of the form  d[subj][obj] = '
        '[(stmt/stmt hash, belief score)].')
    parser.add_argument('-cc', '--crispr-corr-file',
                        help='Precalculated CRISPR correlations in h5 format')
    parser.add_argument('-rc', '--rnai-corr-file',
                        help='Precalculated RNAi correlations in h5 format')
    parser.add_argument('--uninteresting', type=str, nargs='+',
                        help='--uninteresting <filepath> <column name> | Load '
        'a gene set file. The genes in this set will be considered '
        '"uninteresting" when looking for explanations for a pair. If both '
        'genes in the pair are members of the "uniunteresting set", the pair '
        'will be considered explained.')
    parser.add_argument('-gsf', '--gene-set-filter', type=str, help='Load a '
        'file with a gene set to filter to interactions with the gene set.')
    parser.add_argument('-clf', '--cell-line-filter', type=str, nargs='+',
        help='If one argument is provided, the script assumes it is a text '
        'file with cell line identifiers in the DepMap ID format. If two '
        'arguments are provided, the script assumes the first argument is a '
        'text file with cell line identifiers of non DepMap ID format and the '
        'second argument is a dictionary with mappings from the provided cell '
        'line ID format to DepMap ID.')
    parser.add_argument('--margin', type=float, default=1.0, help='How large '
        'diff in terms of standard deviations to accept between data sets '
        'when filtering for correlations during merge. Default is 1 SD.')
    parser.add_argument('--filter-type', default='None', type=str,
                        help='Type of filtering. Options are: `sigma-diff` - '
        'The difference in the distances from the mean measured in number of '
        'standard deviations must be smaller than given by --margin. '
        '`corr-corr-corr` - The product of the scaled correlations* must be '
        'greater than given by --margin. `None` - No filter is applied when '
        'merging the data sets. The resulting correlation dictionary will '
        'simply be the intersection of the provided data sets. | *Scaled '
        'Correlation = (corr-mean)/SD')
    parser.add_argument('-o', '--outbasename', default=str(int(time())),
                        help='Base name for outfiles. Default: UTC timestamp.')
    parser.add_argument('-rec', '--recalc-crispr', action='store_true',
                        help='With \'-r\', recalculate full gene-gene '
                             'correlations of CRISPR data set.')
    parser.add_argument('-rer', '--recalc-rnai', action='store_true',
                        help='With \'-r\', recalculate full gene-gene '
                             'correlations of RNAi data set.')
    parser.add_argument('-s', '--strict', action='store_true', help='With '
        '\'-s\', the correlations are restricted to only be between loaded '
        'gene set. If no gene set is loaded, this option has no effect.')
    parser.add_argument('--dump-unique-pairs', action='store_true', help='With '
        '\'--dump-unique-pairs\', you can save the unique '
        'gene-gene-correlation pairs for each data set. Default: False (no '
        'output)')
    parser.add_argument('-v', '--verbosity', action='count', help='increase '
        'output verbosity (-vv is more than -v)')
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
        'statements to a pickle file.')
    parser.add_argument('-cup', '--unique-depmap-crispr-pairs', help='Uses a '
        'previously saved file to read the unique pairs to read over. Useful '
        'if you are running the script on the full data with no filters.')
    parser.add_argument('-rup', '--unique-depmap-rnai-pairs', help='Uses a '
        'previously saved file to read the unique pairs to read over. Useful '
        'if you are running the script on the full data with no filters.')
    parser.add_argument('--max-pairs', type=int, default=None, help='Limit the '
        'maximum number of gene-gene pairs to explain. If used, the pairs '
        'used will be sampled at random.')
    parser.add_argument('-crange', '--crispr-corr-range', default=[1.0],
                        type=float, nargs="+", help='-crange LOWER_LIM ('
        'UPPER_LIM) | One or two decimal numbers denoting the range of '
        'correlations to consider in terms of number of SD in the CRISPR data.')
    parser.add_argument('-rrange', '--rnai-corr-range', default=[1.0],
                        type=float, nargs="+", help='-rrange LOWER_LIM ('
        'UPPER_LIM) | One or two decimal numbers denoting the range of '
        'correlations to consider in terms of number of SD in the RNAi data.')
    parser.add_argument('-cstats', '--crispr-mean-sigma', type=float, nargs=2,
                        help='-rstats <mean> <stdev> | Provide a value of the '
        'mean and standard deviation for the CRISPR data instead of '
        'calculating it from the full data set.')
    parser.add_argument('-rstats', '--rnai-mean-sigma', type=float, nargs=2,
                        help='-rstats <mean> <stdev> | Provide a value of the '
        'mean and standard deviation for the RNAi data instead of calculating '
        'it from the full data set.')
    parser.add_argument('-noweb', '--no-web-files', action='store_true',
                        help='With the flag active, no output files aimed for '
        'web services are produced')
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
