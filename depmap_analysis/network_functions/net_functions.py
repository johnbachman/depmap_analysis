import logging
import requests
import numpy as np
import pandas as pd
import networkx as nx
from decimal import Decimal
import networkx.algorithms.simple_paths as simple_paths

from indra.config import CONFIG_DICT
from indra.preassembler import hierarchy_manager as hm
from indra.assemblers.indranet import IndraNet
from indra.explanation.model_checker import SignedGraphModelChecker

from depmap_analysis.util.io_functions import pickle_open
import depmap_analysis.network_functions.famplex_functions as fplx_fcns

logger = logging.getLogger(__name__)

np.seterr(all='raise')
NP_PRECISION = 10 ** -np.finfo(np.longfloat).precision  # Numpy precision


GRND_URI = None
try:
    GRND_URI = CONFIG_DICT['INDRA_GROUNDING_SERVICE_URL']
except KeyError:
    logger.warning('Indra Grounding service not available. Add '
                   'INDRA_GROUNDING_SERVICE_URL to `indra/config.ini`')


def sif_dump_df_to_nx_digraph(df, strat_ev_dict, belief_dict,
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
            - 'signed': IndraNet(nx.MultiDiGraph)
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
    if graph_type is 'multidigraph':
        indranet_graph = IndraNet.from_df(sif_df)
    elif graph_type is 'digraph':
        # Flatten
        indranet_graph = IndraNet.digraph_from_df(sif_df,
                                                  'complementary_belief',
                                                  _weight_mapping)
    elif graph_type is 'signed':
        model = IndraNet.signed_from_df(sif_df,
            flattening_method='complementary_belief',
            weight_mapping=_weight_mapping)
        indranet_graph =\
            SignedGraphModelChecker(model=model).get_graph()

    # Add hierarchy relations to graph
    if include_entity_hierarchies and graph_type is not 'signed':
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


def ns_id_from_name(name):
    """Query the groudninfg service for the most likely ns:id pair for name"""
    if GRND_URI:
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
    else:
        logger.warning('Indra Grounding service not available. Add '
                       'INDRA_GROUNDING_SERVICE_URL to `indra/config.ini`')
    return None, None


# Copy from networkx.algorithms.simple_paths
# Added ignore_nodes and ignore_edges arguments
def shortest_simple_paths(G, source, target, weight=None, ignore_nodes=None,
                          ignore_edges=None):
    """Generate all simple paths in the graph G from source to target,
       starting from shortest ones.

    A simple path is a path with no repeated nodes.

    If a weighted shortest path search is to be used, no negative weights
    are allawed.

    Parameters
    ----------
    G : NetworkX graph

    source : node
       Starting node for path

    target : node
       Ending node for path

    weight : string
        Name of the edge attribute to be used as a weight. If None all
        edges are considered to have unit weight. Default value None.

    ignore_nodes : container of nodes
       nodes to ignore, optional

    ignore_edges : container of edges
       edges to ignore, optional

    Returns
    -------
    path_generator: generator
       A generator that produces lists of simple paths, in order from
       shortest to longest.

    Raises
    ------
    NetworkXNoPath
       If no path exists between source and target.

    NetworkXError
       If source or target nodes are not in the input graph.

    NetworkXNotImplemented
       If the input graph is a Multi[Di]Graph.

    Examples
    --------

    >>> G = nx.cycle_graph(7)
    >>> paths = list(nx.shortest_simple_paths(G, 0, 3))
    >>> print(paths)
    [[0, 1, 2, 3], [0, 6, 5, 4, 3]]

    You can use this function to efficiently compute the k shortest/best
    paths between two nodes.

    >>> from itertools import islice
    >>> def k_shortest_paths(G, source, target, k, weight=None):
    ...     return list(islice(nx.shortest_simple_paths(G, source, target,
    ...         weight=weight), k))
    >>> for path in k_shortest_paths(G, 0, 3, 2):
    ...     print(path)
    [0, 1, 2, 3]
    [0, 6, 5, 4, 3]

    Notes
    -----
    This procedure is based on algorithm by Jin Y. Yen [1]_.  Finding
    the first $K$ paths requires $O(KN^3)$ operations.

    See Also
    --------
    all_shortest_paths
    shortest_path
    all_simple_paths

    References
    ----------
    .. [1] Jin Y. Yen, "Finding the K Shortest Loopless Paths in a
       Network", Management Science, Vol. 17, No. 11, Theory Series
       (Jul., 1971), pp. 712-716.

    """
    if source not in G:
        raise nx.NodeNotFound('source node %s not in graph' % source)

    if target not in G:
        raise nx.NodeNotFound('target node %s not in graph' % target)

    if weight is None:
        length_func = len
        shortest_path_func = simple_paths._bidirectional_shortest_path
    else:
        def length_func(path):
            return sum(G.adj[u][v][weight] for (u, v) in zip(path, path[1:]))
        shortest_path_func = simple_paths._bidirectional_dijkstra

    culled_ignored_nodes = set() if ignore_nodes is None else set(ignore_nodes)
    culled_ignored_edges = set() if ignore_edges is None else set(ignore_edges)
    listA = list()
    listB = simple_paths.PathBuffer()
    prev_path = None
    while True:
        cur_ignore_nodes = culled_ignored_nodes.copy()
        cur_ignore_edges = culled_ignored_edges.copy()
        if not prev_path:
            length, path = shortest_path_func(G, source, target, weight=weight,
                                              ignore_nodes=cur_ignore_nodes,
                                              ignore_edges=cur_ignore_edges)
            listB.push(length, path)
        else:
            for i in range(1, len(prev_path)):
                root = prev_path[:i]
                root_length = length_func(root)
                for path in listA:
                    if path[:i] == root:
                        cur_ignore_edges.add((path[i - 1], path[i]))
                try:
                    length, spur = shortest_path_func(
                        G, root[-1], target, ignore_nodes=cur_ignore_nodes,
                        ignore_edges=cur_ignore_edges, weight=weight)
                    path = root[:-1] + spur
                    listB.push(root_length + length, path)
                except nx.NetworkXNoPath:
                    pass
                cur_ignore_nodes.add(root[-1])
        if listB:
            path = listB.pop()
            rcvd_ignore_values = yield path
            if rcvd_ignore_values is not None:
                culled_ignored_nodes = culled_ignored_nodes.union(
                    rcvd_ignore_values[0])
                culled_ignored_edges = culled_ignored_edges.union(
                    rcvd_ignore_values[1])
            listA.append(path)
            prev_path = path
        else:
            break
