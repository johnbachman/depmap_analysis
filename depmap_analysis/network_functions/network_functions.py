import logging
import numpy as np
import networkx as nx
from decimal import Decimal

from indra.preassembler import hierarchy_manager as hm

from depmap_analysis.util.io_functions import pickle_open
import depmap_analysis.network_functions.famplex_functions as fplx_fcns


np.seterr(all='raise')
NP_PRECISION = 10 ** -np.finfo(np.longfloat).precision  # Numpy precision


logger = logging.getLogger('INDRA Network Search')


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
    np_prec = 10 ** -np.finfo(np.longfloat).precision  # Numpy precision
    ns_id_to_nodename = {}
    readers = {'medscan', 'rlimsp', 'trips', 'reach', 'sparser', 'isi'}
    if isinstance(df, str):
        sif_df = pickle_open(df)
    else:
        sif_df = df
    if belief_dict:
        bsd = pickle_open(belief_dict)
    else:
        logger.warning('No belief dict provided, weights will be set to '
                       '1/evidence count')
    if strat_ev_dict:
        sed = pickle_open(strat_ev_dict)
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
                bs = np_prec*10
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
        def _fplx_edge_in_list(multi, edge, check_uri, nx_graph):
            if multi:
                e = 0
                es = nx_graph.edges.get((*edge, e), None)
                while es:
                    if es['stmt_type'] == 'fplx' and \
                            es['stmt_hash'] == check_uri:
                        return True
                    else:
                        e += 1
                        es = nx_graph.edges.get((*edge, e), None)
            else:
                if nx_graph.edges.get(edge):
                    for es in nx_graph.edges.get(edge).get('stmt_list'):
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
        for ns, id, uri in full_entity_list:
            node = id
            # Get name in case it's different than id
            if ns_id_to_nodename.get((ns, id), None):
                node = ns_id_to_nodename[(ns, id)]

            if node not in nx_graph.nodes:
                nx_graph.add_node(node, ns=ns, id=id)
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
