"""The DepMap script

The script matches observed correlations with different graph represenations
if the IndraNetwork

# Inputs:
#   1. Pre-processed correlation matrix with only (gene, gene, z-score)
#   2. nx.DiGraph (or nx.MultiDiGraph?) of IndraNetwork containing at least
#      agA/B: (name, ns, id), hash, type, belief, sign

# Processing:
#   0. If signed graph: match edge sign with correlation sign
#   1. Find direct links both ways
#   2. Find A-X-B links: both ways (2), common target, common regulator
#   3. Find famplex links (have common parent)

# Questions:
#   Q1. Where would the cutting down to specific SD ranges be done?
#   A: Probably outside match correlations, somewhere inside or after
#      preprocessing. Better to do it all at once for one

# Output:
#   An object of a new class that wraps a dataframe that can generate
#   different explanations statistics
"""
import logging
import argparse
import numpy as np
import pandas as pd
import networkx as nx
from pathlib import Path
from datetime import datetime
from depmap_analysis.util.io_functions import pickle_open, dump_it_to_pickle
from depmap_analysis.network_functions.net_functions import \
    SIGNS_TO_INT_SIGN, ns_id_from_name
from depmap_analysis.network_functions.famplex_functions import common_parent
from depmap_analysis.network_functions.depmap_network_functions import \
    corr_matrix_to_generator, same_sign, get_sign
from depmap_analysis.util.statistics import DepMapExplainer
from depmap_analysis.scripts.depmap_preprocessing import run_corr_merge


logger = logging.getLogger('DepMap Script')


def match_correlations(corr_z, indranet, **kwargs):
    """The main loop for matching correlations with INDRA explanations

    Parameters
    ----------
    corr_z : pd.DataFrame
        The pre-processed correlation matrix. No more processing of the
        matrix should have to be done here, i.e. it should already have
        filtered the correlations to the proper SD ranges and removed the
        genes that are not applicable for this explanation,
        self correlations should also have been removed.
    indranet : nx.DiGraph
        The graph representation of the indra network. Each edge should
        have an attribute named 'statements' containing a list of sources
        supporting that edge. If signed search, indranet is expected to be an
        nx.MultiDiGraph with edges keyes by (gene, gene, sign) tuples.

    Returns
    -------
    depmap_explainer
        An instance of the DepMapExplainer class containing the explanations
        for the correlations.
    """
    min_columns = ('agA', 'agB', 'z-score')
    id_columns = min_columns + ('agA_ns', 'agA_id', 'agB_ns', 'agB_id')
    # Map each expl type to a function that handles that explanation
    expl_types = {'a-b': expl_ab,
                  'b-a': expl_ba,
                  'common parent': find_cp,
                  'explained set': None,  # a priori explained
                  'a-x-b': expl_axb,
                  'b-x-a': expl_bxa,
                  'shared regulator': get_sr,
                  'shared target': get_st,
                  }
    bool_columns = ('not in graph', 'explained') + tuple(expl_types.keys())
    stats_columns = id_columns + bool_columns
    expl_columns = ('agA', 'agB', 'z-score', 'expl type', 'expl data')

    # The generator skips self correlations as it yields pairs from the
    # upper triangle of the square matrix
    corr_iter = corr_matrix_to_generator(corr_z)
    signed_search = kwargs.get('signed_search', False)
    ymd_date = kwargs['indra_datetime'] if kwargs.get('indra_datetime') else\
        datetime.now().strftime('%Y%m%d')
    explainer = DepMapExplainer(stats_columns=stats_columns,
                                expl_columns=expl_columns,
                                indra_network_date=ymd_date)
    stats_dict = {k: [] for k in stats_columns}
    expl_dict = {k: [] for k in expl_columns}

    for A, B, zsc in corr_iter:
        # Initialize current iteration stats
        stats = {k: False for k in bool_columns}

        # Append to stats_dict
        stats_dict['agA'].append(A)
        stats_dict['agB'].append(B)
        stats_dict['z-score'].append(zsc)

        # Skip if A or B not in graph
        if A not in indranet.nodes or B not in indranet.nodes:
            for k in set(stats_dict.keys()).difference(set(min_columns)):
                if k == 'not in graph':
                    # Flag not in graph
                    stats_dict[k].append(True)
                else:
                    # All columns are NaN's
                    stats_dict[k].append(np.nan)
            continue

        # Get ns:id
        a_ns, a_id, b_ns, b_id = get_ns_id(A, B, indranet)

        # Append to stats dict
        stats_dict['agA_ns'].append(a_ns)
        stats_dict['agB_ns'].append(b_ns)
        stats_dict['agA_id'].append(a_id)
        stats_dict['agB_id'].append(b_id)

        # If in expl set, skip other explanations
        if kwargs.get('explained_set'):
            if A in kwargs['explained_set'] and B in kwargs['eplained_Set']:
                # Set explained set = True
                stats_dict['explained set'].append(True)

                # Set overall explained = True
                stats_dict['explained'].append(True)

                # All other columns to False
                for k in set(bool_columns).difference(
                        {'explained set', 'explained'}):
                    stats_dict[k].append(False)

                # Set explanation type and data
                # Append to expl_dict
                expl_dict['agA'].append(A)
                expl_dict['agB'].append(B)
                expl_dict['z-score'].append(zsc)
                expl_dict['expl type'].append('explained set')
                expl_dict['expl data'].append(np.nan)

                # And skip the rest of explanations
                continue

        # Loop expl functions
        for expl_type, expl_func in expl_types.items():
            # Function signature: s, o, corr, net, signed, **kwargs
            # Function should return what will be kept in the 'expl_data'
            # column of the expl_df

            # Skip if 'explained set', which is caught above
            if expl_type == 'explained set':
                continue

            # Some functions reverses A, B hence the s, o assignment
            s, o, expl_data = expl_func(A, B, zsc, indranet, signed_search,
                                        **kwargs)
            if expl_data:
                expl_dict['agA'].append(s)
                expl_dict['agB'].append(o)
                expl_dict['z-score'].append(zsc)
                expl_dict['expl type'].append(expl_type)
                expl_dict['expl data'].append(expl_data)

            stats[expl_type] = bool(expl_data)

        # Set explained column
        stats['explained'] = any([b for b in stats.values()])

        # Add stats to stats_dict
        for expl_tp in stats:
            stats_dict[expl_tp].append(stats[expl_tp])

        # Assert that all columns are the same length
        if not all(len(ls) for ls in stats_dict.values()):
            raise IndexError('Unequal column lengths in stats_dict after '
                             'iteration')

    explainer.stats_df = explainer.stats_df.append(other=pd.DataFrame(
        data=stats_dict))
    explainer.expl_df = explainer.expl_df.append(other=pd.DataFrame(
        data=expl_dict))
    explainer.has_data = True
    return explainer


def find_cp(s, o, corr, net, signed, **kwargs):
    s_ns, s_id, o_ns, o_id = get_ns_id(s, o, net)
    if not s_id or not o_id:
        s_ns, s_id = ns_id_from_name(s)
        o_ns, o_id = ns_id_from_name(o)
    parents = list(common_parent(ns1=s_ns, id1=s_id, ns2=o_ns, id2=o_id))
    if parents:
        return s, o, parents
    else:
        return s, o, None


def expl_axb(s, o, corr, net, signed, **kwargs):
    if signed:
        pass  # Todo implement signed check
    else:
        x_nodes = set(net.succ[s]) & set(net.pred[o])
        if x_nodes:
            return s, o, list(x_nodes)
        else:
            return s, o, None


def expl_bxa(s, o, corr, net, signed, **kwargs):
    return expl_axb(o, s, corr, net, signed, **kwargs)


# Shared regulator: A<-X->B
def get_sr(s, o, corr, net, signed, **kwargs):
    if signed:
        pass  # Todo: implement for signed
    else:
        x_nodes = set(net.pred[s]) & set(net.pred[o])
        if x_nodes:
            return s, o, list(x_nodes)
        else:
            return s, o, None


# Shared target: A->X<-B
def get_st(s, o, corr, net, signed, **kwargs):
    if signed:
        pass  # Todo: implement for signed
    else:
        x_nodes = set(net.succ[s]) & set(net.succ[o])
        if x_nodes:
            return s, o, list(x_nodes)
        else:
            return s, o, None


def expl_ab(s, o, corr, net, signed, **kwargs):
    edge_dict = net.edges.get((s, o, kwargs['sign']), None) if signed else \
        net.edges.get((s, o), None)
    if edge_dict:
        return s, o, edge_dict.get('statements')
    return s, o, None


def expl_ba(s, o, corr, net, signed, **kwargs):
    return expl_ab(o, s, corr, net, signed, **kwargs)


def get_edge_statements(s, o, corr, net, signed, **kwargs):
    if signed:
        corr_sign = SIGNS_TO_INT_SIGN[get_sign(corr)]
        return net.edges.get((s, o, corr_sign), None)
    else:
        return net.edges.get((s, o))


def get_ns_id(subj, obj, net):
    """Get ns:id for both subj and obj

    Parameters
    ----------

    subj : str
        The subject node
    obj : str
        The source node
    net : nx.Graph
        A networkx graph object that at least contains node entries.

    Returns
    -------
    tuple
        A tuple with four entries:
        (subj namespace, subj id, obj namespace, obj id)
    """
    s_ns = net.nodes[subj]['ns'] if net.nodes.get(subj) else None
    s_id = net.nodes[subj]['id'] if net.nodes.get(subj) else None
    o_ns = net.nodes[obj]['ns'] if net.nodes.get(obj) else None
    o_id = net.nodes[obj]['id'] if net.nodes.get(obj) else None

    return s_ns, s_id, o_ns, o_id


if __name__ == '__main__':
    parser = argparse.ArgumentParser('DepMap Explainer main script')
    #   1a Load depmap data from scratch | load crispr/rnai raw corr | z-score
    corr_group = parser.add_mutually_exclusive_group(
        required=True)
    corr_group.add_argument(
        '--raw-data', nargs=2, type=str,
        help='--raw-data CRISPR RNAI | File paths to CRISPR raw data and RNAi '
             'raw data from the DepMap Portal. The CRISPR file name should '
             'match '
             '*gene_effect.csv. The RNAi file name should match '
             '*gene_dep_scores.csv'
    )
    corr_group.add_argument(
        '--raw-corr', nargs=2, type=str,
        help='--raw-corr CRISPR RNAI | File paths to raw correlation data ('
             'before z-score conversion) containing hdf compressed '
             'correlation data. These files contain the result of running '
             '`raw_df.corr()`'
    )
    corr_group.add_argument(
        '--z-score', type=str,
        help='The file path to the fully merged correlation matrix '
             'containing z-scores.')

    #   1b Load indranet
    parser.add_argument(
        '--indranet', type=str, required=True,
        help='The indra network to use for explanations. Should be either a '
             'DiGraph or signed DiGraph (a MultiDiGraph with max two edges '
             'per node pair, one edge per sign).'
    )

    #   2. Filter to SD range
    parser.add_argument('--sd-range', nargs='+', type=float, required=True,
                        help='SD range to filter to')
    #   3. Ignore list as file
    parser.add_argument(
        '--ignore-list', type=str,
        help='Provide a text file with one gene name per line to skip in the'
             'explanations')

    # 4 output
    parser.add_argument(
        '--outname', required=True, type=str,
        help='The output name (could contain a path as well) of the pickle '
             'dump of the explainer object')

    args = parser.parse_args()
    run_options = {}

    # Check options
    sd_l, sd_u = args.sd_range if len(args.sd_range) == 2 else\
        ((args.sd_range[0], None) if len(args.sd_range) == 1 else (None, None))

    outname = args.outname if args.outname.endswith('.pkl') else \
        args.outname + '.pkl'
    outpath = Path(outname)

    if not Path(args.indranet).is_file():
        raise FileNotFoundError(f'{args.indranet} does not exist')

    if not sd_l and not sd_u:
        raise ValueError('Must specify at least a lower bound for the SD '
                         'range')

    # Get the z-score corr matrix
    if args.z_score:
        z_corr = pd.read_hdf(args.z_score)
    else:
        z_sc_options = {
            'crispr_raw': args.raw_data[0],
            'rnai_raw': args.raw_data[1],
            'crispr_corr': args.raw_corr[0],
            'rnai_corr': args.raw_corr[1]
        }
        z_corr = run_corr_merge(**z_sc_options)

    # Get indranet
    net = pickle_open(args.indranet)
    run_options['indranet'] = net
    run_options['signed_search'] = net.is_multigraph()  # Todo check signed

    # 2. Filter to SD range
    if sd_l and sd_u:
        z_corr = z_corr[((z_corr > sd_l) & (z_corr < sd_u)) |
                        ((z_corr < -sd_l) & (z_corr > -sd_u))]
    elif sd_l and not sd_u:
        z_corr = z_corr[(z_corr > sd_l) | (z_corr < -sd_l)]
    run_options['corr_z'] = z_corr

    # 3. Ignore list as file
    ignore_file = args.ignore_list
    if ignore_file:
        with open(ignore_file, 'r') as f:
            run_options['explained_set'] = set(f.readlines())

    explanations = match_correlations(**run_options)

    # mkdir in case it  doesn't exist
    outpath.parent.mkdir(parents=True, exist_ok=True)
    dump_it_to_pickle(fname=outpath.as_posix(), pyobj=explanations)
