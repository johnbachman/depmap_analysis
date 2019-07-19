import logging
import requests
import numpy as np
import networkx as nx
from decimal import Decimal
import networkx.algorithms.simple_paths as simple_paths

from indra.config import CONFIG_DICT
from indra.preassembler import hierarchy_manager as hm

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


def sif_dump_df_to_nx_digraph(df, belief_dict=None, strat_ev_dict=None,
                              multi=False, include_entity_hierarchies=True,
                              verbosity=0):
    """Return a NetworkX digraph from a pandas dataframe of a db dump

    Parameters
    ----------
    df : str|pandas.DataFrame
        A dataframe, either as a file path to a pickle or a pandas
        DataFrame object
    belief_dict : str
        The file path to a belief dict that is keyed by statement hashes
        corresponding to the statement hashes loaded in df
    strat_ev_dict : str
        The file path to a dict keyed by statement hashes containing the
        stratified evidence count per statement
    multi : bool
        Default: False; Return an nx.MultiDiGraph if True, otherwise
        return an nx.DiGraph
    include_entity_hierarchies : bool
        Default: True
    verbosity: int
        Output various messages if > 0. For all messages, set to 4

    Returns
    -------
    nx_graph : nx.DiGraph or nx.MultiDiGraph
        By default an nx.DiGraph is returned. By setting multi=True,
        an nx.MultiDiGraph is returned instead."""
    bsd = None
    sed = None
    ns_id_to_nodename = {}
    readers = {'medscan', 'rlimsp', 'trips', 'reach', 'sparser', 'isi'}

    if isinstance(df, str):
        sif_df = pickle_open(df)
    else:
        sif_df = df

    if isinstance(belief_dict, str):
        bsd = pickle_open(belief_dict)
    elif isinstance(belief_dict, dict):
        bsd = belief_dict
    else:
        logger.warning('No belief dict provided, weights will be set to '
                       '1/evidence count')

    if isinstance(strat_ev_dict, str):
        sed = pickle_open(strat_ev_dict)
    elif isinstance(strat_ev_dict, dict):
        sed = strat_ev_dict
    else:
        logger.info('No stratified evidence dict provided')

    # Add as nodes:
    #   'agA_name', 'agB_name'
    # Columns to be added as node attributes:
    #   'agA_ns', 'agA_id', 'agB_ns', 'agB_id'
    # Columns to be added as edge attributes
    #   'stmt_type', 'evidence_count', 'hash'
    # Add from external source:
    #   belief score from provided dict
    #   stratified evidence count by source
    #   famplex edges using entity hierarchies

    if multi:
        nx_graph = nx.MultiDiGraph()
    else:
        nx_graph = nx.DiGraph()
    index = 0
    skipped = 0
    for index, row in sif_df.iterrows():
        if row['agA_name'] is None or row['agB_name'] is None:
            skipped += 1
            if verbosity > 3:
                logger.warning('Skip: %d None agent found: %s' %
                               (skipped, repr(row)))
            elif verbosity > 1:
                logger.warning('Skip: %d None agent found: %s,%s (%d)' %
                               (skipped, row['agA_name'],
                                row['agB_name'],
                                row['hash']))
            continue
        # Add non-existing nodes
        if row['agA_name'] not in nx_graph.nodes:
            nx_graph.add_node(row['agA_name'],
                            ns=row['agA_ns'], id=row['agA_id'])
            ns_id_to_nodename[(row['agA_ns'], row['agA_id'])] = row['agA_name']
        if row['agB_name'] not in nx_graph.nodes:
            nx_graph.add_node(row['agB_name'],
                            ns=row['agB_ns'], id=row['agB_id'])
            ns_id_to_nodename[(row['agB_ns'], row['agB_id'])] = row['agB_name']
        # Add edges
        if bsd:
            try:
                if bsd[row['hash']] == 1 and verbosity:
                    logger.info('Resetting weight from belief score to 1.0 '
                                'for %s' % str(row['hash']))
                b_s = bsd[row['hash']]
                weight = -np.log(max(b_s - 1e-7, 1e-7))
                bs = b_s
            except KeyError:
                if verbosity > 3:
                    logger.warning('Index: %d; Skipped: %d; Hash: %s is '
                                   'missing from belief dict %s' %
                                   (index, skipped, str(row['hash']),
                                    repr(row)))
                elif verbosity:
                    logger.warning('Index: %d; Skipped: %d; Hash: %s is '
                                   'missing from belief dict' %
                                   (index, skipped, str(row['hash'])))
                weight = 1/row['evidence_count']
                bs = None
        else:
            weight = 1 / row['evidence_count']
            bs = None

        if sed:
            try:
                evidence = sed[row['hash']]
            except KeyError:
                if verbosity > 3:
                    logger.warning('Index %d; Skipped %d; Hash: %s is missing '
                                   'from stratified evidence count dict %s' %
                                   (index, skipped, str(row['hash']),
                                    repr(row)))
                elif verbosity:
                    logger.warning('Index %d; Skipped %d; Hash: %s is missing '
                                   'from stratified evidence count dict' %
                                   (index, skipped, str(row['hash'])))
                evidence = {}
        else:
            evidence = {}
        curated = False if not evidence else (False if all(src.lower() in
            readers for src in evidence) else True)
        ed = {'u_for_edge': row['agA_name'],
              'v_for_edge': row['agB_name'],
              'weight': weight,
              'stmt_type': row['stmt_type'],
              'stmt_hash': row['hash'],
              'evidence_count': row['evidence_count'],
              'evidence': evidence,
              'curated': curated,
              'bs': bs}

        if multi:
            nx_graph.add_edge(**ed)
        else:
            if ed.pop('u_for_edge', None):
                ed.pop('v_for_edge', None)
            # Add non-existing edges, append to existing
            if (row['agA_name'], row['agB_name']) in nx_graph.edges:
                nx_graph.edges[(row['agA_name'],
                                row['agB_name'])]['stmt_list'].append(ed)
            else:
                nx_graph.add_edge(row['agA_name'], row['agB_name'],
                                  stmt_list=[ed])
    if skipped:
        logger.warning('Skipped %d edges with None as node' % skipped)
    logger.info('Loaded %d statements into %sDiGraph' %
                (index-skipped, 'Multi' if multi else ''))

    if not multi:
        logger.info('Aggregating belief score for DiGraph edge weights')
        for e in nx_graph.edges:
            # weight = -log(bs)
            ag_belief = ag_belief_score(
                [s['bs'] for s in nx_graph.edges[e]['stmt_list']]
            )
            nx_graph.edges[e]['bs'] = ag_belief
            nx_graph.edges[e]['weight'] = -np.log(ag_belief)

    if include_entity_hierarchies:
        def _fplx_edge_in_list(mg, edge, check_uri, g):
            if mg:
                ec = 0
                es = g.edges.get((*edge, ec), None)
                while es:
                    if es['stmt_type'] == 'fplx' and \
                            es['stmt_hash'] == check_uri:
                        return True
                    else:
                        ec += 1
                        es = g.edges.get((*edge, ec), None)
            else:
                if g.edges.get(edge):
                    for es in g.edges.get(edge).get('stmt_list'):
                        if es['stmt_type'] == 'fplx' and \
                                es['stmt_hash'] == check_uri:
                            return True
                        else:
                            continue
                else:
                    return False

        logger.info('Fetching entity hierarchy relationsships')
        full_entity_list = fplx_fcns.get_all_entities()
        ehm = hm.hierarchies['entity']
        ehm.initialize()
        logger.info('Adding entity hierarchy manager as graph attribute')
        nx_graph.graph['entity_hierarchy_manager'] = ehm
        node_by_uri = {}
        logger.info('Adding entity relations as edges in graph')
        entities = 0
        for ns, _id, uri in full_entity_list:
            node = _id
            # Get name in case it's different than id
            if ns_id_to_nodename.get((ns, _id), None):
                node = ns_id_to_nodename[(ns, _id)]

            if node not in nx_graph.nodes:
                nx_graph.add_node(node, ns=ns, id=_id)
            node_by_uri[uri] = node

            # Add famplex edge
            for puri in ehm.get_parents(uri):
                pns, pid = ehm.ns_id_from_uri(puri)
                pnode = pid
                if ns_id_to_nodename.get((pns, pid), None):
                    pnode = ns_id_to_nodename[(pns, pid)]
                node_by_uri[puri] = pnode
                if pnode not in nx_graph.nodes:
                    nx_graph.add_node(pnode, ns=pns, id=pid)
                # Check if edge already exists
                if not _fplx_edge_in_list(multi, (node, pnode), puri,
                                          nx_graph):
                    entities += 1
                    ed = {'u_for_edge': node,
                          'v_for_edge': pnode,
                          'weight': 1.0,
                          'stmt_type': 'fplx',
                          'stmt_hash': puri,
                          'evidence_count': 1,
                          'evidence': {'fplx': 1},
                          'curated': True,
                          'bs': 1.0}
                    if multi:
                        nx_graph.add_edge(**ed)
                    else:
                        if ed.pop('u_for_edge', None):
                            ed.pop('v_for_edge', None)
                        if (node, pnode) in nx_graph.edges:
                            nx_graph.edges[(node,
                                            pnode)]['stmt_list'].append(ed)
                        else:
                            # The fplx edge is the only edge, add custom
                            # aggregate bs, weight
                            nx_graph.add_edge(node, pnode, stmt_list=[ed],
                                              bs=1.0, weight=1.0)

        logger.info('Loaded %d entity relations into %sDiGraph' %
                    (entities, 'Multi' if multi else ''))
        nx_graph.graph['node_by_uri'] = node_by_uri
        nx_graph.graph['node_by_ns_id'] = ns_id_to_nodename
    return nx_graph


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
            bs = 1
            if len(tup) == 2:
                typ, hsh_a = tup
            elif len(tup) == 3:
                typ, hsh_a, bs = tup
            ax_score_list.append(bs)
        for tup in xb_stmts:
            assert len(tup) == 2 or len(tup) == 3
            bs = 1
            if len(tup) == 2:
                typ, hsh_b = tup
            elif len(tup) == 3:
                typ, hsh_b, bs = tup
            xb_score_list.append(bs)

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
    # Aggregate belief score: 1-prod(1-bs_i)
    try:
        ag_belief = np.longfloat(1.0) - np.prod(np.fromiter(map(
            lambda bs: np.longfloat(1.0) - bs, belief_list),
            dtype=np.longfloat)
        )
    except FloatingPointError as err:
        logger.warning('%s: Resetting ag_belief to 10*np.longfloat '
                       'precision (%.0e)' %
                       (err, Decimal(NP_PRECISION * 10)))
        ag_belief = NP_PRECISION * 10

    return ag_belief


def _ns_id_from_name(name):
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
    ...     return list(islice(nx.shortest_simple_paths(G, source, target, weight=weight), k))
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
