import logging
from decimal import Decimal
from collections import defaultdict

import requests
import subprocess

import requests
import numpy as np
import pandas as pd
from requests.exceptions import ConnectionError

from indra.config import CONFIG_DICT
from indra.ontology.bio import bio_ontology
from indra.belief import load_default_probs
from indra.assemblers.indranet import IndraNet
from indra.databases import get_identifiers_url
from indra_reading.readers import get_reader_classes

from indra.explanation.model_checker.model_checker import \
    signed_edges_to_signed_nodes
from depmap_analysis.util.io_functions import pickle_open
import depmap_analysis.network_functions.famplex_functions as fplx_fcns

logger = logging.getLogger(__name__)

bio_ontology.initialize()
np.seterr(all='raise')
NP_PRECISION = 10 ** -np.finfo(np.longfloat).precision  # Numpy precision
INT_PLUS = 0
INT_MINUS = 1
SIGN_TO_STANDARD = {INT_PLUS: '+', '+': '+', 'plus': '+',
                    '-': '-', 'minus': '-', INT_MINUS: '-'}
SIGNS_TO_INT_SIGN = {INT_PLUS: INT_PLUS, '+': INT_PLUS, 'plus': INT_PLUS,
                     '-': INT_MINUS, 'minus': INT_MINUS, 1: INT_MINUS,
                     None: None}
REVERSE_SIGN = {INT_PLUS: INT_MINUS, INT_MINUS: INT_PLUS,
                '+': '-', '-': '+',
                'plus': 'minus', 'minus': 'plus'}

READERS = set([rc.name.lower() for rc in get_reader_classes()])
READERS.update(['medscan', 'rlimsp'])


def _get_smallest_belief_prior():
    def_probs = load_default_probs()
    return min([v for v in def_probs['syst'].values() if v > 0])# +
               #[[v for v in def_probs['rand'].values() if v > 0]])


MIN_BELIEF = _get_smallest_belief_prior()


GRND_URI = None
try:
    GRND_URI = CONFIG_DICT['GILDA_URL']
except KeyError:
    logger.warning('Indra Grounding service not available. Add '
                   'INDRA_GROUNDING_SERVICE_URL to `indra/config.ini`')


def pinger(domain, timeout=1):
    """Returns True if host at domain is responding"""
    return subprocess.run(["ping", "-c", "1", '-w%d' % int(timeout),
                           domain]).returncode == 0


def gilda_pinger():
    """Return True if the gilda service is available"""
    try:
        logger.info('Trying to reach GILDA service again...')
        return requests.post(GRND_URI, json={'text': 'erk'}).status_code == 200
    except ConnectionError:
        return False


def _curated_func(ev_dict):
    """Return False if no source dict exists, or if all sources are
    readers, otherwise return True."""
    return False if not ev_dict else \
        (False if all(s.lower() in READERS for s in ev_dict) else True)



def _weight_from_belief(belief):
    """Map belief score 'belief' to weight. If the calculation goes below
    precision, return longfloat precision insted to avoid making the
    weight zero."""
    return np.max([NP_PRECISION, -np.log(belief, dtype=np.longfloat)])


def _weight_mapping(G, verbosity=0):
    """Mapping function for adding the weight of the flattened edges

    Parameters
    ----------
    G : IndraNet
        Incoming graph

    Returns
    -------
    G : IndraNet
        Graph with updated belief
    """
    for edge in G.edges:
        try:
            G.edges[edge]['weight'] = \
                _weight_from_belief(G.edges[edge]['belief'])
        except FloatingPointError as err:
            logger.warning('FloatingPointError from unexpected belief '
                           '%s. Resetting ag_belief to 10*np.longfloat '
                           'precision (%.0e)' %
                           (G.edges[edge]['belief'],
                            Decimal(NP_PRECISION * 10)))
            if verbosity == 1:
                logger.error('Error string: %s' % err)
            elif verbosity > 1:
                logger.error('Exception output follows:')
                logger.exception(err)
            G.edges[edge]['weight'] = NP_PRECISION
    return G


def sif_dump_df_merger(df, strat_ev_dict, belief_dict, set_weights=True,
                       verbosity=0):
    """Merge the sif dump df with the provided dictionaries

    Parameters
    ----------
    df : str|pd.DataFrame
        A dataframe, either as a file path to a pickle or a pandas
        DataFrame object.
    belief_dict : str|dict
        The file path to a pickled dict or a dict object keyed by statement
        hash containing the belief score for the corresponding statements.
        The hashes should correspond to the hashes in the loaded dataframe.
    strat_ev_dict : str|dict
        The file path to a pickled dict or a dict object keyed by statement
        hash containing the stratified evidence count per statement. The
        hashes should correspond to the hashes in the loaded dataframe.
    set_weights : bool
        If True, set the edge weights. Default: True.
    verbosity : int
        Output various extra messages if > 1.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame with new columns from the merge
    """
    sed = None

    if isinstance(df, str):
        merged_df = pickle_open(df)
    else:
        merged_df = df

    if 'hash' in merged_df.columns:
        merged_df.rename(columns={'hash': 'stmt_hash'}, inplace=True)

    if isinstance(belief_dict, str):
        belief_dict = pickle_open(belief_dict)
    elif isinstance(belief_dict, dict):
        belief_dict = belief_dict

    if isinstance(strat_ev_dict, str):
        sed = pickle_open(strat_ev_dict)
    elif isinstance(strat_ev_dict, dict):
        sed = strat_ev_dict

    # Extend df with these columns:
    #   belief score from provided dict
    #   stratified evidence count by source
    # Extend df with famplex rows
    # 'stmt_hash' must exist as column in the input dataframe for merge to work
    # Preserve all rows in merged_df, so do left join:
    # merged_df.merge(other, how='left', on='stmt_hash')

    hashes = []
    beliefs = []
    for k, v in belief_dict.items():
        hashes.append(k)
        beliefs.append(v)

    merged_df = merged_df.merge(
        right=pd.DataFrame(data={'stmt_hash': hashes, 'belief': beliefs}),
        how='left',
        on='stmt_hash'
    )
    # Check for missing hashes
    if merged_df['belief'].isna().sum() > 0:
        logger.warning('%d rows with missing belief score found' %
                       merged_df['belief'].isna().sum())
        if verbosity > 1:
            logger.info('Missing hashes in belief dict: %s' % list(
                merged_df['stmt_hash'][merged_df['belief'].isna() == True]))
        logger.info('Setting missing belief scores to 1/evidence count')

    hashes = []
    strat_dicts = []
    for k, v in sed.items():
        hashes.append(k)
        strat_dicts.append(v)

    merged_df = merged_df.merge(
        right=pd.DataFrame(data={'stmt_hash': hashes,
                                 'source_counts': strat_dicts}),
        how='left',
        on='stmt_hash'
    )
    # Check for missing hashes
    if merged_df['source_counts'].isna().sum() > 0:
        logger.warning('%d rows with missing evidence found' %
                       merged_df['source_counts'].isna().sum())
        if verbosity > 1:
            logger.info(
                'Missing hashes in stratified evidence dict: %s' %
                list(merged_df['stmt_hash'][
                         merged_df['source_counts'].isna() == True]))

    logger.info('Setting "curated" flag')
    # Map to boolean 'curated' for reader/non-reader
    merged_df['curated'] = merged_df['source_counts'].apply(func=_curated_func)

    if set_weights:
        logger.info('Setting edge weights')
        # Add weight: -log(belief) or 1/evidence count if no belief
        has_belief = (merged_df['belief'].isna() == False)
        has_no_belief = (merged_df['belief'].isna() == True)
        merged_df['weight'] = 0
        if has_belief.sum() > 0:
            merged_df.loc[has_belief, 'weight'] = merged_df['belief'].apply(
                func=_weight_from_belief)
        if has_no_belief.sum() > 0:
            merged_df.loc[has_no_belief, 'weight'] = merged_df['evidence_count'].apply(
                func=lambda ec: 1/np.longfloat(ec))
    else:
        logger.info('Skipping setting edge weight')

    return merged_df


def sif_dump_df_to_digraph(df, strat_ev_dict, belief_dict,
                           graph_type='digraph',
                           include_entity_hierarchies=True,
                           verbosity=0):
    """Return a NetworkX digraph from a pandas dataframe of a db dump

    Parameters
    ----------
    df : str|pd.DataFrame
        A dataframe, either as a file path to a pickle or a pandas
        DataFrame object.
    belief_dict : str|dict
        The file path to a pickled dict or a dict object keyed by statement
        hash containing the belief score for the corresponding statements.
        The hashes should correspond to the hashes in the loaded dataframe.
    strat_ev_dict : str|dict
        The file path to a pickled dict or a dict object keyed by statement
        hash containing the stratified evidence count per statement. The
        hashes should correspond to the hashes in the loaded dataframe.
    graph_type : str
        Return type for the returned graph. Currently supports:
            - 'digraph': IndraNet(nx.DiGraph) (Default)
            - 'multidigraph': IndraNet(nx.MultiDiGraph)
            - 'signed': SignedGraphModelChecker(ModelChecker)
    include_entity_hierarchies : bool
        Default: True
    verbosity: int
        Output various messages if > 0. For all messages, set to 4

    Returns
    -------
    indranet_graph : IndraNet(graph_type)
        The type is determined by the graph_type argument"""
    graph_options = ['digraph', 'multidigraph', 'signed']
    if graph_type not in graph_options:
        raise ValueError('Graph type %s not supported. Can only chose between'
                         ' %s' % (graph_type, graph_options))

    sif_df = sif_dump_df_merger(df, strat_ev_dict, belief_dict, verbosity)

    # Map ns:id to node name
    logger.info('Creating dictionary with mapping from (ns,id) to node name')
    ns_id_name_tups = set(
        zip(sif_df.agA_ns, sif_df.agA_id, sif_df.agA_name)).union(
        set(zip(sif_df.agB_ns, sif_df.agB_id, sif_df.agB_name)))
    ns_id_to_nodename = {(ns, _id): name for ns, _id, name in ns_id_name_tups}

    # Create graph from df
    if graph_type == 'multidigraph':
        indranet_graph = IndraNet.from_df(sif_df)
    elif graph_type is 'digraph':
        # Flatten
        indranet_graph = IndraNet.digraph_from_df(sif_df,
                                                  'complementary_belief',
                                                  _weight_mapping)
    elif graph_type == 'signed':
        signed_edge_graph = IndraNet.signed_from_df(
            df=sif_df,
            flattening_method='complementary_belief',
            weight_mapping=_weight_mapping
        )
        signed_node_graph = signed_edges_to_signed_nodes(
            graph=signed_edge_graph, copy_edge_data={'weight'})
        signed_edge_graph.graph['node_by_ns_id'] = ns_id_to_nodename
        signed_node_graph.graph['node_by_ns_id'] = ns_id_to_nodename
        return signed_edge_graph, signed_node_graph

    # Add hierarchy relations to graph (not applicable for signed graphs)
    if include_entity_hierarchies and graph_type in ('multidigraph',
                                                     'digraph'):
        logger.info('Fetching entity hierarchy relationsships')
        full_entity_list = fplx_fcns.get_all_entities()
        logger.info('Adding entity hierarchy manager as graph attribute')
        node_by_uri = {uri: _id for (ns, _id, uri) in full_entity_list}
        added_pairs = set()  # Save (A, B, URI)
        logger.info('Building entity relations to be added to data frame')
        entities = 0
        for ns, _id, uri in full_entity_list:
            node = _id
            # Get name in case it's different than id
            if ns_id_to_nodename.get((ns, _id), None):
                node = ns_id_to_nodename[(ns, _id)]
            else:
                ns_id_to_nodename[(ns, _id)] = node

            # Add famplex edge
            for pns, pid in bio_ontology.get_parents(ns, _id):
                puri = get_identifiers_url(pns, pid)
                pnode = pid
                if ns_id_to_nodename.get((pns, pid), None):
                    pnode = ns_id_to_nodename[(pns, pid)]
                else:
                    ns_id_to_nodename[(pns, pid)] = pnode
                # Check if edge already exists
                if (node, pnode, puri) not in added_pairs:
                    entities += 1
                    # Belief and evidence are conditional
                    added_pairs.add((node, pnode, puri))  # A, B, uri of B
                    ed = {'agA_name': node, 'agA_ns': ns, 'agA_id': _id,
                          'agB_name': pnode, 'agB_ns': pns, 'agB_id': pid,
                          'stmt_type': 'fplx', 'evidence_count': 1,
                          'source_counts': {'fplx': 1}, 'stmt_hash': puri,
                          'belief': 1.0, 'weight': NP_PRECISION,
                          'curated': True}
                    # Add non-existing nodes
                    if ed['agA_name'] not in indranet_graph.nodes:
                        indranet_graph.add_node(ed['agA_name'],
                                                ns=ed['agA_ns'],
                                                id=ed['agA_id'])
                    if ed['agB_name'] not in indranet_graph.nodes:
                        indranet_graph.add_node(ed['agB_name'],
                                                ns=ed['agB_ns'],
                                                id=ed['agB_id'])
                    # Add edges
                    ed.pop('agA_id')
                    ed.pop('agA_ns')
                    ed.pop('agB_id')
                    ed.pop('agB_ns')
                    if indranet_graph.is_multigraph():
                        # MultiDiGraph
                        indranet_graph.add_edge(ed['agA_name'],
                                                ed['agB_name'],
                                                **ed)
                    else:
                        # DiGraph
                        u = ed.pop('agA_name')
                        v = ed.pop('agB_name')

                        # Check edge
                        if indranet_graph.has_edge(u, v):
                            indranet_graph.edges[(u, v)]['statements'].append(
                                ed)
                        else:
                            indranet_graph.add_edge(u,
                                                    v,
                                                    belief=1.0,
                                                    weight=1.0,
                                                    statements=[ed])

        logger.info('Loaded %d entity relations into dataframe' % entities)
        indranet_graph.graph['node_by_uri'] = node_by_uri
    indranet_graph.graph['node_by_ns_id'] = ns_id_to_nodename
    return indranet_graph


def rank_nodes(node_list, nested_dict_stmts, gene_a, gene_b, x_type):
    """Returns a list of tuples of nodes and their rank score

    The provided node list should contain the set of nodes that connects subj
    and obj through an intermediate node found in nested_dict_stmts.

    nested_dict_stmts

        d[subj][obj] = [stmts/stmt hashes]

    node_list : list[nodes]
    nested_dict_stmts : defaultdict(dict)
        Nested dict of statements: nest_d[subj][obj]
    gene_a : str
        Name of node A in an A-X-B connection
    gene_b : str
        Name of node B in an A-X-B connection
    x_type : str
        One of 'x_is_intermediary', 'x_is_downstream' or 'x_is_upstream'

    -------
    Returns
    dir_path_nodes_wb : list[(node, rank)]
        A list of node, rank tuples.
    """
    def _tuple_rank(ax_stmts, xb_stmts):
        def _body(t):
            assert len(t) == 2 or len(t) == 3
            bel = MIN_BELIEF
            if len(t) == 2:
                tp, hs = t
            elif len(t) == 3:
                tp, hs, bel = t
            else:
                raise IndexError('Tuple must have len(t) == 2,3 Tuple: %s' %
                                 repr(t))
            return tp, hs, bel
        ax_score_list = []
        xb_score_list = []
        for tup in ax_stmts:
            typ, hsh_a, belief = _body(tup)
            ax_score_list.append(belief)
        for tup in xb_stmts:
            typ, hsh_b, belief = _body(tup)
            xb_score_list.append(belief)
        return ax_score_list, xb_score_list

    def _dict_rank(ax_stmts, xb_stmts):
        ax_score_list = []
        xb_score_list = []
        for sd in ax_stmts:
            ax_score_list.append(float(sd.get('belief', MIN_BELIEF)))
        for sd in xb_stmts:
            xb_score_list.append(float(sd.get('belief', MIN_BELIEF)))
        return ax_score_list, xb_score_list

    def _calc_rank(nest_dict_stmts, subj_ax, obj_ax, subj_xb, obj_xb):
        ax_stmts = nest_dict_stmts[subj_ax][obj_ax]
        xb_stmts = nest_dict_stmts[subj_xb][obj_xb]
        hsh_a, hsh_b = None, None

        # The statment with the highest belief score should
        # represent the edge (potentially multiple stmts per edge)

        if isinstance(ax_stmts[0], tuple):
            ax_score_list, xb_score_list = _tuple_rank(ax_stmts, xb_stmts)
        elif isinstance(ax_stmts[0], (dict, defaultdict)):
            ax_score_list, xb_score_list = _dict_rank(ax_stmts, xb_stmts)

        # Rank by multiplying the best two belief scores for each edge
        rank = max(ax_score_list) * max(xb_score_list)

        # No belief score should be zero, thus rank should never be zero
        try:
            assert rank != 0
        except AssertionError:
            logger.warning('Combined rank == 0 for hashes %s and %s, implying '
                           'belief score is 0 for at least one of the '
                           'following statements: ' % (hsh_a, hsh_b))
        return rank

    dir_path_nodes_wb = []

    if x_type is 'x_is_intermediary':  # A->X->B or A<-X<-B
        for gene_x in node_list:
            x_rank = _calc_rank(nest_dict_stmts=nested_dict_stmts,
                                subj_ax=gene_a, obj_ax=gene_x,
                                subj_xb=gene_x, obj_xb=gene_b)
            dir_path_nodes_wb.append((gene_x, x_rank))

    elif x_type is 'x_is_downstream':  # A->X<-B
        for gene_x in node_list:
            x_rank = _calc_rank(nest_dict_stmts=nested_dict_stmts,
                                subj_ax=gene_a, obj_ax=gene_x,
                                subj_xb=gene_b, obj_xb=gene_x)
            dir_path_nodes_wb.append((gene_x, x_rank))
    elif x_type is 'x_is_upstream':  # A<-X->B

        for gene_x in node_list:
            x_rank = _calc_rank(nest_dict_stmts=nested_dict_stmts,
                                subj_ax=gene_x, obj_ax=gene_a,
                                subj_xb=gene_x, obj_xb=gene_b)
            dir_path_nodes_wb.append((gene_x, x_rank))

    return dir_path_nodes_wb


def ag_belief_score(belief_list):
    """Each item in `belief_list` should be a float"""
    # Aggregate belief score: 1-prod(1-belief_i)
    try:
        ag_belief = np.longfloat(1.0) - np.prod(np.fromiter(map(
            lambda belief: np.longfloat(1.0) - belief, belief_list),
            dtype=np.longfloat)
        )
    except FloatingPointError as err:
        logger.warning('%s: Resetting ag_belief to 10*np.longfloat '
                       'precision (%.0e)' %
                       (err, Decimal(NP_PRECISION * 10)))
        ag_belief = NP_PRECISION * 10

    return ag_belief


def ns_id_from_name(name, gilda_retry=False):
    """Query the groudning service for the most likely ns:id pair for name"""
    global GILDA_TIMEOUT
    if gilda_retry and GILDA_TIMEOUT and gilda_pinger():
        logger.info('GILDA is responding again!')
        GILDA_TIMEOUT = False

    if GRND_URI and not GILDA_TIMEOUT:
        try:
            res = requests.post(GRND_URI, json={'text': name})
            if res.status_code == 200:
                rj = res.json()[0]
                return rj['term']['db'], rj['term']['id']
            else:
                logger.warning('Grounding service responded with code %d, '
                               'check your query format and URL' %
                               res.status_code)
        except IndexError:
            logger.info('No grounding exists for %s' % name)
        except ConnectionError as err2:
            logger.warning('GILDA has timed out, ignoring future requests')
            GILDA_TIMEOUT = True
    else:
        if GILDA_TIMEOUT:
            logger.warning('Indra Grounding service not available.')
        else:
            logger.warning('Indra Grounding service not available. Add '
                           'GILDA_URL to `indra/config.ini`')
    return None, None
