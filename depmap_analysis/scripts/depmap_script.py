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
import json
import logging
import pickle as pkl
import argparse as ap
import itertools as itt
from datetime import datetime
from random import choice as rnd_choice

import pandas as pd
from numpy import float64
from indra.tools import assemble_corpus as ac
from indra_db.util.dump_sif import load_pickle_from_s3

import depmap_analysis.util.io_functions as io
import depmap_analysis.network_functions.net_functions as nf
import depmap_analysis.network_functions.famplex_functions as ff
from depmap_analysis.network_functions import depmap_network_functions as dnf

logger = logging.getLogger('DepMap Script')

KEY_ERR_TOL = 100

# 1. no geneset -> use corr from full DepMap data, no filtering needed
# 2. geneset, not strict -> interaction has to contain at least one gene from
#    the loaded list
# 3. geneset loaded, strict -> each correlation can only contain genes from
#    the loaded data set


def _dump_nest_dict_to_csv(fname, nested_dict, separator=',', header=None,
                           excl_sr=True):
    if header:
        with open(fname, 'w') as fo, \
                open(fname.split('.')[0]+'_sr_only.csv', 'w') as fosro:
            fo.write(separator.join(header)+'\n')
            fosro.write(separator.join(header)+'\n')
    with open(fname, 'a') as fo, \
            open(fname.split('.')[0]+'_sr_only.csv', 'a') as fosro:
        skip = 0
        for ok in nested_dict:
            for ik in nested_dict[ok]:
                cd = nested_dict[ok][ik]['meta_data']
                if excl_sr and nested_dict[ok][ik]['sr_only']:
                    skip += 1
                    if skip < 10:
                        logger.info('Skipping sr only: %s and %s' % (ik, ok))
                    elif skip == 10:
                        logger.info('Muting skip messages...')
                    fosro.write('%s,%s,%s,%s\n' %
                        (ok, ik,
                         str(cd['crispr']) if cd and json.dumps(cd.get(
                             'crispr'))
                         else '0',
                         str(cd['rnai']) if cd and json.dumps(cd.get('rnai'))
                         else '0'))
                    continue
                fo.write('%s,%s,%s,%s\n' %
                         (ok, ik,
                          str(cd['crispr']) if cd and json.dumps(cd.get(
                              'crispr'))
                          else '0',
                          str(cd['rnai']) if cd and json.dumps(cd.get('rnai'))
                          else '0'))


def _is_float(n):
    if type(n) is float64 or type(n) is float:
        return True
    else:
        return False


def _rnd_pair_gen(rgs):
    return rnd_choice(rgs), rnd_choice(rgs)


def _parse_cell_filter(cl_file, id2depmapid_pkl=None, namespace='CCLE_Name'):
    logger.info('Parsing cell lines...')
    cell_lines = []
    with open(cl_file, 'r') as fi:
        first_line = fi.readline()
        if ',' in first_line:
            col_separator = ','
        else:
            col_separator = None
    with open(cl_file, 'r') as fi:
        if id2depmapid_pkl:
            id2depmapid_dict = io.pickle_open(id2depmapid_pkl)
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


def _parse_explained_genes(gene_set_file, check_column):
    logger.info('Parsing explaind genes assuming column %s can be '
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
    logger.info('Loaded %i explained genes' % len(genes))
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


def loop_body(args, hgnc_name2id=None, fplx_name2id=None):

    global any_expl, any_expl_not_sr, common_parent, ab_expl_count, \
        directed_im_expl_count, both_im_dir_expl_count, \
        any_axb_non_sr_expl_count, sr_expl_count, \
        shared_regulator_only_expl_count, explanations_of_pairs, unexplained, \
        explained_nested_dict, id1, id2, nested_dict_statements, \
        dataset_dict, dir_node_set, nx_dir_graph, explained_set, \
        part_of_explained, sr_explanations, any_expl_ign_sr

    # Store booleans for each found connection
    found = False  # Flag anything found
    directed = False  # Flag direct
    undirected = False  # Flag complex
    x_is_intermediary = False  # Flag intermediate connections
    x_is_downstream = False  # Flag shared target
    x_is_upstream = False  # Flag shared regulator connection
    part_of_explained_set = False  # Flag uniteresting
    has_common_parent = False  # Flag common parent

    for subj, obj in itt.permutations((id1, id2), r=2):
        if dnf.entry_exist_dict(nested_dict_statements, subj, obj):

            # Get the statements
            stmts = nested_dict_statements[subj][obj]

            # check if directed, put in the explained nested dict
            dir_stmts, undir_stmts = dnf.get_directed(stmts)
            directed = bool(dir_stmts)
            undirected = bool(undir_stmts)
            explained_nested_dict[subj][obj]['directed'] = dir_stmts
            explained_nested_dict[subj][obj]['undirected'] = undir_stmts

            if args.verbosity:
                logger.info('Found direct connection between %s and '
                            '%s' % (subj, obj))
            found = True
            stmt_tuple = (subj, obj, 'direct', [], json.dumps(dataset_dict))
            explanations_of_pairs.append(stmt_tuple)

        # Checking 1. "pathway": A -> X -> B and B -> X -> A
        if subj in dir_node_set and obj in dir_node_set:
            dir_path_nodes = list(set(nx_dir_graph.succ[subj]) &
                                  set(nx_dir_graph.pred[obj]))
            if dir_path_nodes and \
                    dnf.entry_exist_dict(nested_dict_statements, subj, obj):
                found = True
                x_is_intermediary = True
                if args.verbosity:
                    logger.info('Found directed path of length 2 '
                                'between %s and %s' % (subj, obj))

                dir_path_nodes_wb = nf.rank_nodes(
                    node_list=dir_path_nodes,
                    nested_dict_stmts=nested_dict_statements,
                    gene_a=subj,
                    gene_b=obj,
                    x_type='x_is_intermediary')

                explained_nested_dict[subj][obj]['x_is_intermediary'] \
                    = dir_path_nodes_wb
                stmt_tuple = (subj, obj, 'pathway',
                              dir_path_nodes_wb, json.dumps(dataset_dict))
                explanations_of_pairs.append(stmt_tuple)

    # Check common parent

    # Get ns:id
    id1_ns, id1_id = None, None
    id2_ns, id2_id = None, None

    # HGNC
    if hgnc_name2id:
        id1_ns, id1_id = (hgnc_name2id[id1], 'HGNC') if \
            hgnc_name2id.get(id1) else (None, None)
        id2_ns, id2_id = (hgnc_name2id[id2], 'HGNC') if \
            hgnc_name2id.get(id2) else (None, None)

    # Try FPLX if not HGNC
    if id1_id is None and fplx_name2id:
        id1_ns, id1_id = (fplx_name2id[id1], 'FPLX') if \
            fplx_name2id.get(id1) else (None, None)
    if id2_id is None and fplx_name2id:
        id2_ns, id2_id = (fplx_name2id[id2], 'FPLX') if \
            fplx_name2id.get(id2) else (None, None)

    # Last resort: GILDA
    if id1_id is None:
        id1_ns, id1_id = nf.ns_id_from_name(id1)
    if id2_id is None:
        id2_ns, id2_id = nf.ns_id_from_name(id2)

    if id1_id and id2_id and ff.has_common_parent(ns1=id1_ns, id1=id1_id,
                                                  ns2=id2_ns, id2=id2_id):
        has_common_parent = True
        found = True
        parents = list(ff.common_parent(id1=id1, id2=id2))
        explained_nested_dict[id1][id2]['common_parents'] = parents
        explained_nested_dict[id2][id1]['common_parents'] = parents
        stmt_tuple = (id1, id2, 'common_parents', parents, [])
        explanations_of_pairs.append(stmt_tuple)

    # Check if both A and B are in list of "uninteresting genes"
    if explained_set and \
            id1 in explained_set and id2 in explained_set:
        explained_nested_dict[id1][id2]['explained_set'] = True
        explained_nested_dict[id2][id1]['explained_set'] = True
        found = True
        part_of_explained_set = True
        stmt_tuple = (id1, id2, 'explained_set', [], [])
        explanations_of_pairs.append(stmt_tuple)

    if id1 in dir_node_set and id2 in dir_node_set:
        # Checking 2: shared target A -> X <- B
        downstream_share = list(set(nx_dir_graph.succ[id1]) &
                                set(nx_dir_graph.succ[id2]))
        # Checking 3: shared regulator A <- X -> B
        upstream_share = list(set(nx_dir_graph.pred[id1]) &
                              set(nx_dir_graph.pred[id2]))
        if downstream_share:
            found = True
            x_is_downstream = True
            downstream_share_wb = nf.rank_nodes(
                node_list=downstream_share,
                nested_dict_stmts=nested_dict_statements,
                gene_a=id1,
                gene_b=id2,
                x_type='x_is_downstream')
            stmt_tuple = (id1, id2, 'shared_target',
                          downstream_share_wb, json.dumps(dataset_dict))
            if args.verbosity:
                logger.info('Found downstream share: %s and %s share '
                            '%i targets' %
                            (id1, id2, len(downstream_share)))
            explained_nested_dict[id1][id2]['x_is_downstream'] = \
                downstream_share_wb
            explained_nested_dict[id2][id1]['x_is_downstream'] = \
                downstream_share_wb
            explanations_of_pairs.append(stmt_tuple)

        if upstream_share:
            found = True
            x_is_upstream = True
            upstream_share_wb = nf.rank_nodes(
                node_list=upstream_share,
                nested_dict_stmts=nested_dict_statements,
                gene_a=id1,
                gene_b=id2,
                x_type='x_is_upstream')
            stmt_tuple = (id1, id2, 'shared_upstream',
                          upstream_share_wb, json.dumps(dataset_dict))
            if args.verbosity:
                logger.info('Found upstream share: %s and %s are both '
                            'directly downstream of %i nodes' %
                            (id1, id2, len(upstream_share)))
            explained_nested_dict[id1][id2]['x_is_upstream'] = \
                upstream_share_wb
            explained_nested_dict[id2][id1]['x_is_upstream'] = \
                upstream_share_wb
            sr_explanations.append(stmt_tuple)

    # Make sure the connection types we didn't find are empty lists.
    # Also add correlation so it can be queried for at the same time
    # as the items for the second drop down.
    if found:
        # Any explanation found
        any_expl += 1

        # Count explanations with only non-shared regulators
        if not x_is_upstream and any([directed, undirected, x_is_intermediary,
                                     x_is_downstream, has_common_parent,
                                      part_of_explained_set]):
            any_expl_not_sr += 1

        # Count explanations with common parents
        if has_common_parent:
            common_parent += 1

        # Count explanations due to uninteresting genes
        if part_of_explained_set:
            part_of_explained += 1
        else:
            explained_nested_dict[id1][id2]['explained_set'] = False
            explained_nested_dict[id2][id1]['explained_set'] = False

        if any([directed, undirected, x_is_intermediary, x_is_downstream,
                has_common_parent, part_of_explained_set]):
            any_expl_ign_sr += 1

        # Count A-B or B-A connections found per set(A,B)
        if directed or undirected:
            ab_expl_count += 1

        # Count directed A-X-B connections found per set(A,B)
        if x_is_intermediary:
            directed_im_expl_count += 1

        # Count A-X-B, ignoring shared regulators
        if x_is_intermediary or x_is_downstream:
            any_axb_non_sr_expl_count += 1

        if x_is_upstream:
            sr_expl_count += 1

        # Count when shared regulator is the only explanation
        # NOTE: this should be equal to any_expl - any_expl_ign_sr
        if all([not directed, not undirected, not x_is_intermediary,
                not x_is_downstream, not has_common_parent,
                not part_of_explained_set]) and x_is_upstream:
            shared_regulator_only_expl_count += 1
            explained_nested_dict[id1][id2]['sr_only'] = True
            explained_nested_dict[id2][id1]['sr_only'] = True
        else:
            explained_nested_dict[id1][id2]['sr_only'] = False
            explained_nested_dict[id2][id1]['sr_only'] = False

        assert shared_regulator_only_expl_count == (any_expl - any_expl_ign_sr)

        for s, o in itt.permutations((id1, id2), r=2):
            # Correlation/meta data
            explained_nested_dict[s][o]['meta_data'] = dataset_dict
            # common_parents
            if not dnf.entry_exist_dict(explained_nested_dict[s], o,
                                         'common_parents'):
                explained_nested_dict[s][o]['common_parents'] = []
            # directed
            if not dnf.entry_exist_dict(explained_nested_dict[s], o,
                                         'directed'):
                explained_nested_dict[s][o]['directed'] = []
            # undirected
            if not dnf.entry_exist_dict(explained_nested_dict[s], o,
                                         'undirected'):
                explained_nested_dict[s][o]['undirected'] = []
            # x_is_intermediary
            if not dnf.entry_exist_dict(explained_nested_dict[s], o,
                                         'x_is_intermediary'):
                explained_nested_dict[s][o]['x_is_intermediary'] = []
            # x_is_upstream
            if not dnf.entry_exist_dict(explained_nested_dict[s], o,
                                         'x_is_upstream'):
                explained_nested_dict[s][o]['x_is_upstream'] = []
            # x_is_downstream
            if not dnf.entry_exist_dict(explained_nested_dict[s], o,
                                         'x_is_downstream'):
                explained_nested_dict[s][o]['x_is_downstream'] = []

    # found == True if at least one connection was found
    # not found == True is only True when no connection was found
    if not found:
        unexplained.append((id1, id2, json.dumps(dataset_dict)))
        if args.verbosity and args.verbosity > 1:
            logger.info('No explainable path found between %s and '
                        '%s.' % (id1, id2))


def main(args, hgnc_name2id=None, fplx_name2id=None):
    if hgnc_name2id or fplx_name2id:
        logger.info('Using loaded dict(s) for mapping')
    global any_expl, any_expl_not_sr, common_parent, ab_expl_count, \
        directed_im_expl_count, both_im_dir_expl_count, \
        any_axb_non_sr_expl_count, sr_expl_count, \
        shared_regulator_only_expl_count, explanations_of_pairs, unexplained, \
        explained_nested_dict, id1, id2, nested_dict_statements, dataset_dict, \
        avg_corr, dir_node_set, nx_dir_graph, explained_set, part_of_explained,\
        sr_explanations, any_expl_ign_sr

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

    # Parse "explained genes"
    if args.explained_set and len(args.explained_set) == 2:
        explained_set = _parse_explained_genes(
            gene_set_file=args.explained_set[0],
            check_column=args.explained_set[1])
        logger.info('Loading "explained pairs."')
    elif args.explained_set and len(args.explained_set) != 2:
        sys.exit('Argument --explained-set takes exactly two arguments: '
                 '--explained-set <file> <column name>')

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
        if args.belief_score_dict == 's3':
            belief_key = 'belief_dict.pkl'
            belief_dict = load_pickle_from_s3(key=belief_key)
        elif args.belief_score_dict.endswith('.json'):
            belief_dict = io.json_open(args.belief_score_dict)
        elif args.belief_score_dict.endswith('.pkl'):
            belief_dict = io.pickle_open(args.belief_score_dict)

    args_dict = _arg_dict(args)
    npairs = 0

    gene_filter_list = []
    all_hgnc_ids = set()
    stmts_all = []
    master_corr_dict = dnf.create_nested_dict()

    filter_settings = {'gene_set_filter': args.gene_set_filter,
                       'strict': args.strict,
                       'cell_line_filter': cell_lines,
                       'cell_line_translation_dict': io.pickle_open(
                                                     args.cell_line_filter[1])
                       if args.cell_line_filter and len(args.cell_line_filter)
                       == 2 else None,
                       'margin': args.margin,
                       'filter_type': (args.filter_type
                                       if args.filter_type
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
        elif filter_settings.get('filter_type'):
            logger.info('Using filter type "%s"' %
                        filter_settings['filter_type'])
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
                gene_filter_list = io.read_gene_set_file(
                    gf=filter_settings['gene_set_filter'],
                    data=pd.read_csv(args_dict['crispr']['data'],
                                         index_col=0, header=0))
            elif args_dict.get('rnai') and not args_dict.get('crispr'):
                gene_filter_list = io.read_gene_set_file(
                        gf=filter_settings['gene_set_filter'],
                        data=pd.read_csv(args_dict['rnai']['data'],
                                         index_col=0, header=0))
            elif args_dict.get('crispr') and args_dict.get('rnai'):
                gene_filter_list = \
                    set(io.read_gene_set_file(
                        gf=filter_settings['gene_set_filter'],
                        data=pd.read_csv(args_dict['crispr']['data'],
                                         index_col=0, header=0))) & \
                    set(io.read_gene_set_file(
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

    if not args.light_weight_stmts and not args.nested_dict_in and \
            not args.sif_df_in:
        if args.statements_in:  # Get statments from file
            stmts_all = set(ac.load_statements(args.statements_in))
        # Use api to get statements. _NOT_ the same as querying for each ID
        else:
            if args.gene_set_filter and gene_filter_list:
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
        with open(args.light_weight_stmts, 'r') as csvf:
            first_line = csvf.readline()
            delim = '\t' if '\t' in first_line else ','

        hash_df = pd.read_csv(args.light_weight_stmts, delimiter=delim)
        nested_dict_statements = \
            dnf.nested_hash_dict_from_pd_dataframe(hash_df)
    elif args.nested_dict_in:
        nested_dict_statements = io.pickle_open(args.nested_dict_in)
    elif args.sif_df_in:
        nested_dict_statements = dnf.sif_dump_df_to_nest_d(args.sif_df_in,
                                                           belief_dict)
    else:
        nested_dict_statements = dnf.dedupl_nested_dict_gen(stmts_all,
                                                            belief_dict)
        if args.nested_dict_out:
            io.dump_it_to_pickle(fname=args.nested_dict_out,
                                 pyobj=nested_dict_statements)
    assert nested_dict_statements is not None

    # Get directed simple graph
    if args.directed_graph_in:
        with open(args.directed_graph_in, 'rb') as rpkl:
            nx_dir_graph = pkl.load(rpkl)
    else:
        # Create directed graph from statement dict
        nx_dir_graph = dnf.nested_stmt_dict_to_nx_digraph(
            nest_d=nested_dict_statements, belief_dict=belief_dict)
        # Save as pickle file
        if args.directed_graph_out:
            io.dump_it_to_pickle(fname=args.directed_graph_out,
                                 pyobj=nx_dir_graph)
    dir_node_set = set(nx_dir_graph.nodes)

    # LOOP THROUGH THE UNIQUE CORRELATION PAIRS, MATCH WITH INDRA NETWORK
    any_expl = 0  # Count if any explanation per (A,B) correlation found
    any_expl_not_sr = 0  # Count any explanation, exlcuding when shared
    # regulator is the only explanation
    any_expl_ign_sr = 0  # Count any explanation, ingoring shared regulator
    # explanations
    common_parent = 0  # Count if common parent found per set(A,B)
    part_of_explained = 0  # Count pairs part the "explained set"
    ab_expl_count = 0  # Count A-B/B-A as one per set(A,B)
    directed_im_expl_count = 0  # Count any A>X>B,B>X>A as one per set(A,B)
    any_axb_non_sr_expl_count = 0  # if shared target found per set(A,B)
    sr_expl_count = 0  # Count if shared regulator found per set(A,B)
    shared_regulator_only_expl_count = 0  # Count if only shared regulator
    explanations_of_pairs = []  # Saves all non shared regulator explanations
    sr_explanations = []  # Saves all shared regulator explanations
    unexplained = []  # Unexplained correlations
    skipped = 0

    # KeyError count
    key_errs = 0

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
    # 2. dir -> undir graph -> jsons to check all corr neighbors ->
    #    2nd dropdown
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
        try:
            for outer_id, do in master_corr_dict.items():
                for inner_id, dataset_dict in do.items():
                    if len(dataset_dict.keys()) == 0:
                        skipped += 1
                        if args.verbosity:
                            logger.info('Skipped outer_id=%s and inner_id=%s'
                                        % (outer_id, inner_id))
                        continue

                    id1, id2 = outer_id, inner_id
                    loop_body(args, hgnc_name2id=hgnc_name2id,
                              fplx_name2id=fplx_name2id)
        except KeyError as err:
            key_errs += 1
            if key_errs > KEY_ERR_TOL:
                raise KeyError('Too many KeyErrors registered. '
                               'Aborting...\n%s' % repr(err))

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

        try:
            for outer_id, do in master_corr_dict.items():
                for inner_id, knockout_dict in do.items():
                    if len(knockout_dict.keys()) == 0:
                        skipped += 1
                        if args.verbosity:
                            logger.info('Skipped outer_id=%s and inner_id=%s'
                                        % (outer_id, inner_id))
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
                    loop_body(args, hgnc_name2id=hgnc_name2id,
                              fplx_name2id=fplx_name2id)
        except KeyError as err:
            key_errs += 1
            if key_errs > KEY_ERR_TOL:
                raise KeyError('Too many KeyErrors registered. '
                               'Aborting...\n%s' % repr(err))

    # loop brca dependency ONLY
    elif args.brca_dependencies and not \
            (args_dict.get('rnai') or args_dict.get('crispr')):
        logger.info('Gene pairs generated from BRCA gene enrichment data '
                    'only.')
        brca_data_set = pd.read_csv(args.brca_dependencies, header=0)
        depend_in_breast_genes = brca_data_set.drop(
            axis=1, labels=['Url Label', 'Type'])[brca_data_set['Type'] ==
                                                  'gene']
        genes = set(depend_in_breast_genes['Gene/Compound'].values)
        npairs = len(list(itt.combinations(genes, 2)))
        logger.info('Looking for connections between %i pairs' % (
            npairs if npairs > 0 else args.max_pairs)
        )
        try:
            for id1, id2 in itt.combinations(genes, 2):
                gene1_data = []
                gene2_data = []
                # For each non-diagonal pair in file, insert in dataset_dict:
                # geneA, geneB,
                # dataset for A, dataset for B,
                # T-stat for A, T-stat for B,
                # P-value for A, P-value
                for row in depend_in_breast_genes[
                     depend_in_breast_genes[
                         'Gene/Compound'] == id1].iterrows():
                    gene1_data.append((row[1]['Dataset'],
                                       row[1]['T-Statistic'],
                                       row[1]['P-Value']))

                for row in depend_in_breast_genes[
                    depend_in_breast_genes[
                        'Gene/Compound'] == id2].iterrows():
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
                loop_body(args, hgnc_name2id=hgnc_name2id,
                          fplx_name2id=fplx_name2id)
        except KeyError as err:
            key_errs += 1
            if key_errs > KEY_ERR_TOL:
                raise KeyError('Too many KeyErrors registered. '
                               'Aborting...\n%s' % repr(err))

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
        try:
            for _ in range(npairs):
                id1, id2 = _rnd_pair_gen(rnd_gene_set)
                assert not isinstance(id1, list)
                loop_body(args, hgnc_name2id=hgnc_name2id,
                          fplx_name2id=fplx_name2id)
        except KeyError as err:
            key_errs += 1
            if key_errs > KEY_ERR_TOL:
                raise KeyError('Too many KeyErrors registered. '
                               'Aborting...\n%s' % repr(err))

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
    long_string += '> Total correlations explained, ignoring shared ' \
                   'regulator: %i' % any_expl_ign_sr + '\n'
    long_string += '> Total correlations explained, excluding shared ' \
                   'regulator (total - shared only): %i' % \
                   (any_expl - shared_regulator_only_expl_count) + '\n'
    long_string += '>    %i correlations have an explanation involving a ' \
                   'common parent' % common_parent + '\n'
    if args.explained_set:
        long_string += '>    %i gene pairs were considered explained as ' \
                       'part of the "explained set"' % part_of_explained + '\n'
    long_string += '>    %i explanations involving direct connection or ' \
                   'complex' % ab_expl_count + '\n'
    long_string += '>    %i correlations have a directed explanation ' \
                   'involving an intermediate node (A->X->B/A<-X<-B)' \
                   % directed_im_expl_count + '\n'
    long_string += '>    %i correlations have an explanation involving an ' \
                   'intermediate node excluding shared regulators' % \
                   any_axb_non_sr_expl_count + '\n'
    long_string += '>    %i correlations have an explanation involving a ' \
                   'shared regulator (A<-X->B)' % sr_expl_count + '\n'
    long_string += '>    %i correlations have shared regulator as only ' \
                   'explanation' % shared_regulator_only_expl_count + '\n\n'

    if stats_dict and (stats_dict.get('rnai') or stats_dict.get('crispr')):
        long_string += 'Statistics of input data:' + '\n\n'
    if stats_dict and stats_dict.get('rnai'):
        long_string += '  RNAi data ' + '\n'
        long_string += ' -----------' + '\n'
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
        long_string += ' -------------' + '\n'
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
    nx_expl_dir_graph = dnf.nested_stmt_explained_dict_nx_digraph(
        nest_d=explained_nested_dict)

    if not args.no_web_files:
        # 'explained_nodes' are used to produce first drop down
        explained_nodes = list(nx_expl_dir_graph.nodes)
        logger.info('Dumping json "explainable_ids.json" for first dropdown.')
        io.dump_it_to_json(args.outbasename + '_explainable_ids.json',
                           explained_nodes)

        # Get undir graph, save each neighbor lookup as json for 2nd dropdown
        nx_expl_undir_graph = nx_expl_dir_graph.to_undirected()
        dnf.nx_undir_to_neighbor_lookup_json(
            expl_undir_graph=nx_expl_undir_graph, outbasename=args.outbasename)

    # Easiest way to check if pairs are explained or not is to loop explained
    # dict. Skip shared regulators.
    _dump_nest_dict_to_csv(
        fname=args.outbasename+'_explained_correlations.csv',
        nested_dict=explained_nested_dict,
        header=['gene1', 'gene2', 'meta_data'],
        excl_sr=True)

    io.dump_it_to_pickle(fname=args.outbasename + '_explained_nest_dict.pkl',
                      pyobj=explained_nested_dict)
    headers = ['subj', 'obj', 'type', 'X', 'meta_data']
    io.dump_it_to_csv(fname=args.outbasename + '_explanations_of_pairs.csv',
                   pyobj=explanations_of_pairs, header=headers)
    io.dump_it_to_csv(fname=
                    args.outbasename+'_explanations_of_shared_regulators.csv',
                   pyobj=sr_explanations, header=headers)
    io.dump_it_to_csv(fname=args.outbasename + '_unexpl_correlations.csv',
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
        
        -cc/-rc: precalculated correlation matrices in hdf format of 
        crispr/rnai
        
        Either of
        -ndi: nested dictionary of INDRA statements of the format `d[gene][
        gene] = [stmts/stmt hashes]` OR 
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
    parser.add_argument('--explained-set', type=str, nargs='+',
                        help='--explained-set <filepath> <column name> | Load '
        'a gene set file. The genes in this set will be considered '
        '"explained" when looking for explanations for a pair. If both '
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
        'difference in z-score between data sets to accept when filtering for '
        'correlations during merge. Default is 1 (i.e. 1 SD).')
    parser.add_argument('--filter-type', default='None', type=str,
                        help='Type of filtering for merging correlations from '
        'multiple data sets. Options are: `z_score_diff` - The difference in '
        'z-score must be smaller than given by --margin. `z_score_product` - '
        'The product of the z-scores must be greater than given by --margin. '
        '`z_score_mean` - The mean of the z-scores must be smaller than given '
        'by --margin. `sign` - Only filter out correlations if they have '
        'opposite signs. `None` - No filter is applied when merging the data '
        'sets. The resulting correlation dictionary will simply be the '
        'intersection of the provided data sets.')
    parser.add_argument('-o', '--outbasename',
                        default=datetime.utcnow().strftime('%Y%m%d'),
                        help='Base name for outfiles. Default: UTC '
                             'timestamp.')
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
    parser.add_argument('--sif-df-in',
                        help='Use a sif dump dataframe for generating the '
                             'nested dict and the dir_graph')
    parser.add_argument('--hgnc-name2id',
                        help='Pickled dict mapping hgnc names to hgnc ids '
                             'used in the db dump input.')
    parser.add_argument('--fplx-name2id',
                        help='Pickled dict mapping fplx names to fplx ids '
                             'used in the db dump input.')
    a = parser.parse_args()

    ymd_date = datetime.utcnow().strftime('%Y%m%d')
    with open(a.outbasename+'dep_map_script_log_%s.log' % ymd_date,
              'w', newline='\n') as f:
        f.write('Created on %s\n' %
                datetime.utcnow().strftime('%Y %b %d, %H:%M:%S'))
        f.write('Execution path: {}\n\n'.format(os.getcwd()))
        f.write('Command line option : value\n---------------------------\n')
        for arg in vars(a):
            f.write('{} : {}\n'.format(arg, getattr(a, arg)))

    # Load mappings
    if a.hgnc_name2id:
        hgnc_nm2id_map = io.pickle_open(a.hgnc_name2id)
    else:
        hgnc_nm2id_map = None
    if a.fplx_name2id:
        fplx_nm2id_map = io.pickle_open(a.fplx_name2id)
    else:
        fplx_nm2id_map = None
    done = main(a, hgnc_name2id=hgnc_nm2id_map, fplx_name2id=fplx_nm2id_map)
    if done == 0 or done is None:
        logger.info('Script finished without errors')
