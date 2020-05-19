import logging
from decimal import Decimal

import requests
import subprocess
import numpy as np
import pandas as pd
import networkx as nx
from collections import deque
from requests.exceptions import ConnectionError
from networkx.classes.reportviews import OutMultiEdgeView, OutEdgeView, \
    NodeView
from indra.config import CONFIG_DICT
from indra.preassembler import hierarchy_manager as hm
from indra.assemblers.indranet import IndraNet
from indra.explanation.pathfinding_util import signed_edges_to_signed_nodes, \
    signed_nodes_to_signed_edge
from depmap_analysis.util.io_functions import pickle_open
import depmap_analysis.network_functions.famplex_functions as fplx_fcns

logger = logging.getLogger(__name__)

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


GILDA_TIMEOUT = False


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
    sed = None
    readers = {'medscan', 'rlimsp', 'trips', 'reach', 'sparser', 'isi'}

    def _curated_func(ev_dict):
        return False if not ev_dict else \
            (False if all(s.lower() in readers for s in ev_dict) else True)

    def _weight_from_belief(belief):
        """Map belief score 'belief' to weight. If the calculation goes below
        precision, return longfloat precision insted to avoid making the
        weight zero."""
        return np.max([NP_PRECISION, -np.log(belief, dtype=np.longfloat)])

    def _weight_mapping(G):
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
            G.edges[edge]['weight'] = \
                _weight_from_belief(G.edges[edge]['belief'])
        return G

    if isinstance(df, str):
        sif_df = pickle_open(df)
    else:
        sif_df = df

    if 'hash' in sif_df.columns:
        sif_df.rename(columns={'hash': 'stmt_hash'}, inplace=True)

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
    # Preserve all rows in sif_df, so do left join:
    # sif_df.merge(other, how='left', on='stmt_hash')

    hashes = []
    beliefs = []
    for k, v in belief_dict.items():
        hashes.append(k)
        beliefs.append(v)

    sif_df = sif_df.merge(
        right=pd.DataFrame(data={'stmt_hash': hashes, 'belief': beliefs}),
        how='left',
        on='stmt_hash'
    )
    # Check for missing hashes
    if sif_df['belief'].isna().sum() > 0:
        logger.warning('%d rows with missing belief score found' %
                       sif_df['belief'].isna().sum())
        if verbosity > 1:
            logger.info('Missing hashes in belief dict: %s' % list(
                sif_df['stmt_hash'][sif_df['belief'].isna() == True]))
        logger.info('Setting missing belief scores to 1/evidence count')

    hashes = []
    strat_dicts = []
    for k, v in sed.items():
        hashes.append(k)
        strat_dicts.append(v)

    sif_df = sif_df.merge(
        right=pd.DataFrame(data={'stmt_hash': hashes,
                                 'source_counts': strat_dicts}),
        how='left',
        on='stmt_hash'
    )
    # Check for missing hashes
    if sif_df['source_counts'].isna().sum() > 0:
        logger.warning('%d rows with missing evidence found' %
                       sif_df['source_counts'].isna().sum())
        if verbosity > 1:
            logger.info('Missing hashes in stratified evidence dict: %s' %
                        list(sif_df['stmt_hash'][sif_df['source_counts'].isna()
                                            == True]))
    # Map ns:id to node name
    logger.info('Creating dictionary with mapping from (ns,id) to node name')
    ns_id_name_tups = set(
        zip(sif_df.agA_ns, sif_df.agA_id, sif_df.agA_name)).union(
        set(zip(sif_df.agB_ns, sif_df.agB_id, sif_df.agB_name)))
    ns_id_to_nodename = {(ns, _id): name for ns, _id, name in ns_id_name_tups}

    logger.info('Setting "curated" flag')
    # Map to boolean 'curated' for reader/non-reader
    sif_df['curated'] = sif_df['source_counts'].apply(func=_curated_func)

    logger.info('Setting edge weights')
    # Add weight: -log(belief) or 1/evidence count if no belief
    has_belief = (sif_df['belief'].isna() == False)
    has_no_belief = (sif_df['belief'].isna() == True)
    sif_df['weight'] = 0
    if has_belief.sum() > 0:
        sif_df.loc[has_belief, 'weight'] = sif_df['belief'].apply(
            func=_weight_from_belief)
    if has_no_belief.sum() > 0:
        sif_df.loc[has_no_belief, 'weight'] = sif_df['evidence_count'].apply(
            func=lambda ec: 1/np.longfloat(ec))

    # Create graph from df
    if graph_type == 'multidigraph':
        indranet_graph = IndraNet.from_df(sif_df)
    elif graph_type is 'digraph':
        # Flatten
        indranet_graph = IndraNet.digraph_from_df(sif_df,
                                                  'complementary_belief',
                                                  _weight_mapping)
    elif graph_type == 'signed':
        signed_edge_graph = IndraNet.signed_from_df(sif_df,
            flattening_method='complementary_belief',
            weight_mapping=_weight_mapping)
        signed_node_graph = signed_edges_to_signed_nodes(
            graph=signed_edge_graph, copy_edge_data={'weight'})
        signed_edge_graph.graph['node_by_ns_id'] = ns_id_to_nodename
        signed_node_graph.graph['node_by_ns_id'] = ns_id_to_nodename
        return signed_edge_graph, signed_node_graph

    # Add hierarchy relations to graph (not applicable for signed graphs)
    if include_entity_hierarchies and graph_type != 'signed':
        logger.info('Fetching entity hierarchy relationsships')
        full_entity_list = fplx_fcns.get_all_entities()
        ehm = hm.hierarchies['entity']
        ehm.initialize()
        logger.info('Adding entity hierarchy manager as graph attribute')
        node_by_uri = {}
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
            node_by_uri[uri] = node

            # Add famplex edge
            for puri in ehm.get_parents(uri):
                pns, pid = ehm.ns_id_from_uri(puri)
                pnode = pid
                if ns_id_to_nodename.get((pns, pid), None):
                    pnode = ns_id_to_nodename[(pns, pid)]
                else:
                    ns_id_to_nodename[(pns, pid)] = pnode
                node_by_uri[puri] = pnode
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
        HGNC name of gene A in an A-X-B connection
    gene_b : str
        HGNC name of gene B in an A-X-B connection
    x_type : str
        One of 'x_is_intermediary', 'x_is_downstream' or 'x_is_upstream'

    -------
    Returns
    dir_path_nodes_wb : list[(node, rank)]
        A list of node, rank tuples.
    """

    def _calc_rank(nest_dict_stmts, subj_ax, obj_ax, subj_xb, obj_xb):
        ax_stmts = nest_dict_stmts[subj_ax][obj_ax]
        xb_stmts = nest_dict_stmts[subj_xb][obj_xb]
        ax_score_list = []
        xb_score_list = []
        hsh_a, hsh_b = None, None

        # The statment with the highest belief score should
        # represent the edge (potentially multiple stmts per edge)
        # To get latest belief score: see indra_db.belief
        for tup in ax_stmts:
            assert len(tup) == 2 or len(tup) == 3
            belief = 1
            if len(tup) == 2:
                typ, hsh_a = tup
            elif len(tup) == 3:
                typ, hsh_a, belief = tup
            ax_score_list.append(belief)
        for tup in xb_stmts:
            assert len(tup) == 2 or len(tup) == 3
            belief = 1
            if len(tup) == 2:
                typ, hsh_b = tup
            elif len(tup) == 3:
                typ, hsh_b, belief = tup
            xb_score_list.append(belief)

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


def get_sorted_neighbors(G, node, reverse, g_edges):
    # better sorted key
    """Sort by aggregated belief per edge"""
    neighbors = G.predecessors(node) if reverse else G.successors(node)
    # Check signed node
    if isinstance(node, tuple):
        if reverse:
            return sorted(
                neighbors,
                key=lambda n:
                    g_edges[signed_nodes_to_signed_edge(n, node)]['belief'],
                reverse=True
            )
        else:
            return sorted(
                neighbors,
                key=lambda n:
                    g_edges[signed_nodes_to_signed_edge(node, n)]['belief'],
                reverse=True)

    else:
        if reverse:
            return sorted(neighbors,
                          key=lambda n: g_edges[(n, node)]['belief'],
                          reverse=True)
        else:
            return sorted(neighbors,
                          key=lambda n: g_edges[(node, n)]['belief'],
                          reverse=True)


# Implementation inspired by networkx's
# networkx.algorithms.traversal.breadth_first_search::generic_bfs_edges
def bfs_search(g, source_node, g_nodes=None, g_edges=None, reverse=False,
               depth_limit=2, path_limit=None, max_per_node=5,
               node_filter=None, node_blacklist=None, terminal_ns=None,
               sign=None, **kwargs):
    """Do breadth first search from a given node and yield paths

    Parameters
    ----------
    g : nx.Digraph
        An nx.DiGraph to search in. Can also be a signed node graph.
    source_node : node
        Node in the graph to start from.
    g_nodes : nx.classes.reportviews.nodesNodeView
        The nodes property to look up nodes from. Set this if the node
        attribute 'ns' needs to be looked up from another graph object than
        the one provided as `g`. Default: g.nodes
    g_edges : nx.classes.reportviews.OutMultiEdgeView|OutEdgeView
        The edges property to look up edges and their data from. Set this if
        the edge beliefs needs to be looked up from another grapth object
        than `g`. Default: d.edges
    reverse : bool
        If True go upstream from source, otherwise go downstream. Default:
        False.
    depth_limit : int
        Stop when all paths with this many edges have been found. Default: 2.
    path_limit : int
        The maximum number of paths to return. Default: no limit.
    max_per_node : int
        The maximum number of paths to yield per parent node. If 1 is
        chosen, the search only goes down to the leaf node of its first
        encountered branch. Default: 5
    node_filter : list[str]
        The allowed namespaces (node attribute 'ns') for the nodes in the
        path
    node_blacklist : set[node]
        A set of nodes to ignore. Default: None.
    terminal_ns : list[str]
        Force a path to terminate when any of the namespaces in this list
        are encountered.
    sign : int
        If set, defines the search to be a signed search. Default: None.

    Yields
    ------
    path : tuple(node)
        Paths in the bfs search starting from `source`.
    """
    int_plus = 0
    int_minus = 1
    g_nodes = g.nodes if g_nodes is None else g_nodes
    g_edges = g.edges if g_edges is None else g_edges
    if not isinstance(g_nodes, NodeView):
        raise ValueError('Provided object for g_nodes is not a valid '
                         'NodeView object')
    if not isinstance(g_edges, (OutEdgeView, OutMultiEdgeView)):
        raise ValueError('Provided object for g_edges is not a valid '
                         'OutEdgeView or OutMultiEdgeView object')

    queue = deque([(source_node,)])
    visited = ({source_node}).union(node_blacklist) \
        if node_blacklist else {source_node}
    yielded_paths = 0
    while queue:
        cur_path = queue.popleft()
        last_node = cur_path[-1]
        node_name = last_node[0] if isinstance(last_node, tuple) else \
            last_node

        # if last node is in terminal_ns, continue to next path
        if terminal_ns and g_nodes[node_name]['ns'].lower() in terminal_ns:
            # Check correct leaf sign for signed search
            continue

        sorted_neighbors = get_sorted_neighbors(G=g, node=last_node,
                                                reverse=reverse,
                                                g_edges=g_edges)

        yielded_neighbors = 0
        # for neighb in neighbors:
        for neighb in sorted_neighbors:
            neig_name = neighb[0] if isinstance(neighb, tuple) else neighb

            # Check cycles
            if sign is not None:
                # Avoid signed paths ending up on the opposite sign of the
                # same node
                if (neig_name, int_minus) in cur_path or \
                        (neig_name, int_plus) in cur_path:
                    continue
            elif neighb in visited:
                continue

            # Check namespace
            if node_filter and len(node_filter) > 0:
                if g_nodes[neig_name]['ns'].lower() not in node_filter:
                    continue

            # Add to visited nodes and create new path
            visited.add(neighb)
            new_path = cur_path + (neighb,)

            # Check yield and break conditions
            if len(new_path) > depth_limit + 1:
                continue
            else:
                # Yield newest path and recieve new ignore values

                # Signed search yield
                if sign is not None:
                    if reverse:
                        # Upstream signed search should not end in negative
                        # node
                        if new_path[-1][1] == int_minus:
                            ign_vals = None
                            pass
                        else:
                            ign_vals = yield new_path
                            yielded_paths += 1
                            yielded_neighbors += 1

                    else:
                        # Downstream signed search has to end on node with
                        # requested sign
                        if new_path[-1][1] != sign:
                            ign_vals = None
                            pass
                        else:
                            ign_vals = yield new_path
                            yielded_paths += 1
                            yielded_neighbors += 1

                # Unsigned search
                else:
                    ign_vals = yield new_path
                    yielded_paths += 1
                    yielded_neighbors += 1

                # If new ignore nodes are recieved, update set
                if ign_vals is not None:
                    ign_nodes, ign_edges = ign_vals
                    visited.update(ign_nodes)

                # Check max paths reached, no need to add to queue
                if path_limit and yielded_paths >= path_limit:
                    break

            # Append yielded path
            queue.append(new_path)

            # Check if we've visited enough neighbors
            # Todo: add all neighbors to 'visited' and add all skipped
            #  paths to queue? Currently only yielded paths are
            #  investigated deeper
            if max_per_node and yielded_neighbors >= max_per_node:
                break

        # Check path limit again to catch the inner break for path_limit
        if path_limit and yielded_paths >= path_limit:
            break
