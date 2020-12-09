"""Explainer and helper functions for depmap_script2.py When adding new
explanation functions, please also add them to the mapping at the end"""
import logging
import networkx as nx
from typing import Set, Union, Tuple, List, Optional

import pandas as pd
from networkx import DiGraph, MultiDiGraph
from pybel.dsl import CentralDogma

from depmap_analysis.network_functions.famplex_functions import common_parent
from depmap_analysis.network_functions.net_functions import gilda_normalization, \
    INT_PLUS, INT_MINUS


__all__ = ['explained', 'expl_ab', 'expl_ba', 'expl_axb', 'expl_bxa',
           'find_cp', 'get_sd', 'get_sr', 'get_st', 'get_ns_id_pybel_node',
           'get_ns_id', 'normalize_corr_names', 'expl_functions',
           'funcname_to_colname']
logger = logging.getLogger(__name__)


def explained(s, o, corr, net, _type, **kwargs):
    # This function is used for a priori explained relationships
    return s, o, 'explained_set'


def find_cp(s, o, corr, net, _type, **kwargs):
    if _type == 'pybel':
        s_name = kwargs['s_name']
        s_ns, s_id = get_ns_id_pybel_node(s_name, s)
        o_name = kwargs['o_name']
        o_ns, o_id = get_ns_id_pybel_node(o_name, o)
    else:
        s_name = s
        o_name = o
        s_ns, s_id, o_ns, o_id = get_ns_id(s, o, net)

    if not s_id:
        s_ns, s_id, s_norm_name = gilda_normalization(s_name)
    if not o_id:
        o_ns, o_id, o_norm_name = gilda_normalization(o_name)

    if s_id and o_id:
        # Possible kwargs:
        #   - immediate_only : bool
        #         Determines if all or just the immediate parents should be
        #         returned. Default: False.
        #   - is_a_part_of : iterable
        #         If provided, the parents must be in this set of ids. The
        #         set is assumed to be valid ontology labels (see
        #         ontology.label()).
        parents = list(common_parent(
            ns1=s_ns, id1=s_id, ns2=o_ns, id2=o_id,
            immediate_only=kwargs.get('immediate_only', False),
            is_a_part_of=kwargs.get('is_a_part_of')
        ))
        if parents:
            # if kwargs.get('ns_set'):
            #     parents = {(ns, _id) for ns, _id in parents if ns.lower() in
            #                kwargs['ns_set']} or None
            return s, o, parents

    return s, o, None


def expl_axb(s, o, corr, net, _type, **kwargs):
    x_set = set(net.succ[s]) & set(net.pred[o])
    if _type in {'signed', 'pybel'}:
        x_nodes = _get_signed_interm(s, o, corr, net, x_set)
    else:
        x_nodes = x_set

    # Filter ns
    if kwargs.get('ns_set'):
        x_nodes = {x for x in x_nodes if
                   net.nodes[x]['ns'].lower() in kwargs['ns_set']} or None

    if x_nodes:
        return s, o, list(x_nodes)
    else:
        return s, o, None


def expl_bxa(s, o, corr, net, _type, **kwargs):
    if _type == 'pybel':
        s_name = kwargs.pop('s_name')
        o_name = kwargs.pop('o_name')
        options = {'o_name': s_name, 's_name': o_name}
    else:
        options = {}
    return expl_axb(o, s, corr, net, _type, **kwargs, **options)


# Shared regulator: A<-X->B
def get_sr(s, o, corr, net, _type, **kwargs):
    # Filter ns
    if kwargs.get('ns_set'):
        ns_filt_args = (net, kwargs['ns_set'])
        s_pred = set(_node_ns_filter(net.pred[s], *ns_filt_args))
        o_pred = set(_node_ns_filter(net.pred[o], *ns_filt_args))
        x_set = _node_ns_filter(s_pred & o_pred, *ns_filt_args)
        x_set_union = _node_ns_filter(s_pred | o_pred, *ns_filt_args)
    else:
        s_pred = set(net.pred[s])
        o_pred = set(net.pred[o])
        x_set = s_pred & o_pred
        x_set_union = s_pred | o_pred

    # Sort out sign
    if _type in {'signed', 'pybel'}:
        x_nodes = _get_signed_shared_regulators(s, o, corr, net, x_set, False)
        x_nodes_union = _get_signed_shared_regulators(s, o, corr, net,
                                                      x_set_union, True)
    else:
        x_nodes = x_set
        x_nodes_union = x_set_union

    # Return if there is anything in union, s_pred or o_pred
    if x_nodes_union or s_pred or o_pred:
        return s, o, (list(s_pred), list(o_pred),
                      list(x_nodes or []), list(x_nodes_union or []))
    else:
        return s, o, None


# Shared target: A->X<-B
def get_st(s, o, corr, net, _type, **kwargs):
    # Filter ns
    if kwargs.get('ns_set'):
        ns_filt_args = (net, kwargs['ns_set'])
        s_succ = set(_node_ns_filter(net.succ[s], *ns_filt_args))
        o_succ = set(_node_ns_filter(net.succ[o], *ns_filt_args))
        x_set = _node_ns_filter(s_succ & o_succ, *ns_filt_args)
        x_set_union = _node_ns_filter(s_succ | o_succ, *ns_filt_args)
    else:
        s_succ = set(net.succ[s])
        o_succ = set(net.succ[o])
        x_set = s_succ & o_succ
        x_set_union = s_succ | o_succ

    # Sort out sign
    if _type in {'signed', 'pybel'}:
        x_nodes = _get_signed_shared_targets(s, o, corr, net, x_set, False)
        x_nodes_union = _get_signed_shared_targets(s, o, corr, net,
                                                   x_set_union, True)
    else:
        x_nodes = x_set
        x_nodes_union = x_set_union

    # Return if there is any node in union, s_succ or o_succ
    if x_nodes_union or s_succ or o_succ:
        return s, o, (list(s_succ), list(o_succ),
                      list(x_nodes or []), list(x_nodes_union or []))
    else:
        return s, o, None


def get_sd(s, o, corr, net, _type, **kwargs):
    # Get next-nearest-neighborhood for subject
    args = (net, _type in {'signed', 'pybel'}, kwargs.get('ns_set'))
    s_x_set = _get_nnn_set(s, *args)
    o_x_set = _get_nnn_set(o, *args)

    # Get intersection and union of each nodes' 1st & 2nd layer neighbors
    x_set = s_x_set & o_x_set
    x_set_union = s_x_set | o_x_set

    if _type in {'signed', 'pybel'}:
        x_nodes = _get_signed_deep_interm(s, o, corr, net, x_set, False)
        x_nodes_union = _get_signed_deep_interm(s, o, corr, net, x_set_union,
                                                True)
    else:
        x_nodes = x_set
        x_nodes_union = x_set_union

    if x_nodes_union or s_x_set or o_x_set:
        s_x_list = set()
        o_x_list = set()
        if _type in {'signed', 'pybel'}:
            for xy in s_x_set:
                s_x_list.update(xy)
            for xy in o_x_set:
                o_x_list.update(xy)
        else:
            s_x_list = s_x_set
            o_x_list = o_x_set

        return s, o, (list(s_x_list or []), list(o_x_list or []),
                      list(x_nodes or []), list(x_nodes_union or []))
    else:
        return s, o, None


def expl_ab(s, o, corr, net, _type, **kwargs):
    edge_dict = get_edge_statements(s, o, corr, net, _type, **kwargs)
    if edge_dict:
        return s, o, edge_dict.get('stmt_hash') if _type == 'pybel' else \
            edge_dict.get('statements')
    return s, o, None


def expl_ba(s, o, corr, net, _type, **kwargs):
    if _type == 'pybel':
        s_name = kwargs.pop('s_name')
        o_name = kwargs.pop('o_name')
        options = {'o_name': s_name, 's_name': o_name}
    else:
        options = {}
    return expl_ab(o, s, corr, net, _type, **kwargs, **options)


def get_edge_statements(s, o, corr, net, _type, **kwargs):
    if _type in {'signed', 'pybel'}:
        int_sign = INT_PLUS if corr >= 0 else INT_MINUS
        return net.edges.get((s, o, int_sign), None)
    else:
        return net.edges.get((s, o))


def _get_signed_interm(s, o, corr, sign_edge_net, x_set):
    # Used for a->x->b and b->x->a relations
    # Make sure we have the right sign type
    int_sign = INT_PLUS if corr >= 0 else INT_MINUS

    # ax and xb sign need to match correlation sign
    x_approved = set()
    for x in x_set:
        ax_plus = (s, x, INT_PLUS) in sign_edge_net.edges
        ax_minus = (s, x, INT_MINUS) in sign_edge_net.edges
        xb_plus = (x, o, INT_PLUS) in sign_edge_net.edges
        xb_minus = (x, o, INT_MINUS) in sign_edge_net.edges

        if int_sign == INT_PLUS:
            if ax_plus and xb_plus or ax_minus and xb_minus:
                x_approved.add(x)
        if int_sign == INT_MINUS:
            if ax_plus and xb_minus or ax_minus and xb_plus:
                x_approved.add(x)
    return x_approved


def _get_signed_shared_regulators(s: str, o: str, corr: float,
                                  sign_edge_net: nx.MultiDiGraph,
                                  x_set: Set, union: bool) -> Set[str]:
    x_approved = set()

    for x in x_set:
        xs_plus = (x, s, INT_PLUS) in sign_edge_net.edges
        xo_plus = (x, o, INT_PLUS) in sign_edge_net.edges
        xs_minus = (x, s, INT_MINUS) in sign_edge_net.edges
        xo_minus = (x, o, INT_MINUS) in sign_edge_net.edges

        if union:
            if any([xs_plus, xo_plus, xs_minus, xo_minus]):
                x_approved.add(x)
            else:
                pass
        else:
            if corr > 0:
                if xs_plus and xo_plus or xs_minus and xo_minus:
                    x_approved.add(x)
            else:
                if xs_plus and xo_minus or xs_minus and xo_plus:
                    x_approved.add(x)

    return x_approved


def _get_signed_shared_targets(s: str, o: str, corr: float,
                               sign_edge_net: nx.MultiDiGraph,
                               x_set: Set, union: bool) -> Set[str]:
    x_approved = set()

    for x in x_set:
        sx_plus = (s, x, INT_PLUS) in sign_edge_net.edges
        ox_plus = (o, x, INT_PLUS) in sign_edge_net.edges
        sx_minus = (s, x, INT_MINUS) in sign_edge_net.edges
        ox_minus = (o, x, INT_MINUS) in sign_edge_net.edges

        if union:
            if any([sx_plus, ox_plus, sx_minus, ox_minus]):
                x_approved.add(x)
            else:
                pass
        else:
            if corr > 0:
                if sx_plus and ox_plus or sx_minus and ox_minus:
                    x_approved.add(x)
            else:
                if sx_plus and ox_minus or sx_minus and ox_plus:
                    x_approved.add(x)

    return x_approved


def _get_signed_deep_interm(
        s: str, o: str, corr: float, sign_edge_net: nx.MultiDiGraph,
        xy_set: Set[Tuple[str, str]], union: bool) -> Set[str]:
    # Make sure we have the right sign type
    path_sign = INT_PLUS if corr >= 0 else INT_MINUS

    # a-x-y and b-x-y need to both match path sign
    x_approved = set()
    for x, y in xy_set:
        sx_plus = (s, x, INT_PLUS) in sign_edge_net.edges
        ox_plus = (o, x, INT_PLUS) in sign_edge_net.edges
        xy_plus = (x, y, INT_PLUS) in sign_edge_net.edges
        sx_minus = (s, x, INT_MINUS) in sign_edge_net.edges
        ox_minus = (o, x, INT_MINUS) in sign_edge_net.edges
        xy_minus = (x, y, INT_MINUS) in sign_edge_net.edges

        # Match args for _approve_signed_paths
        args = (sx_minus, sx_plus, ox_minus, ox_plus, xy_minus, xy_plus,
                path_sign, union)

        # Add nodes that form paths with the correct sign
        if _approve_signed_paths(*args):
            x_approved.update({x, y})
    return x_approved


def _approve_signed_paths(sxm: bool, sxp: bool, oxm: bool, oxp: bool,
                          xym: bool, xyp: bool, sign: int, union: bool) \
        -> bool:
    def _asp(n1: bool, n2: bool, p1: bool, p2: bool, s: int) -> bool:
        # Approve Signed Path
        if s == INT_PLUS:
            return p1 and p2 or n1 and n2
        else:
            return p1 and n2 or n1 and p2

    # Match args for _asp
    sargs = (sxm, xym, sxp, xyp, sign)
    oargs = (oxm, xym, oxp, xyp, sign)
    if union:
        return _asp(*sargs) or _asp(*oargs)
    else:
        return _asp(*sargs) and _asp(*oargs)


def _get_nnn_set(
        n: str, g: nx.MultiDiGraph, signed: bool,
        ns_set: Optional[Union[Set[str], List[str], Tuple[str]]] = None
) -> Set[Union[str, Tuple[str, str]]]:
    n_x_set = set()
    for x in g.succ[n]:
        if ns_set and g.nodes[x]['ns'].lower() not in ns_set:
            # Skip if x is not in allowed name space
            continue

        # If signed, add edges instead and match sign in helper
        if signed:
            for y in g.succ[x]:
                if ns_set and g.nodes[y]['ns'].lower() not in ns_set:
                    # Skip if y is not in allowed name space
                    continue
                n_x_set.add((x, y))
        # Just add nodes for unsigned
        else:
            n_x_set.add(x)
            if ns_set:
                n_x_set.update({
                    y for y in g.succ[x] if g.nodes[y]['ns'].lower() in ns_set
                })
            else:
                n_x_set.update(g.succ[x])
    return n_x_set


def _node_ns_filter(node_list: Union[Set[str], List[str]],
                    net: Union[DiGraph, MultiDiGraph],
                    allowed_ns: Union[Set[str], List[str], Tuple[str]]) \
        -> Set[str]:
    return {x for x in node_list if net.nodes[x]['ns'].lower()
            in allowed_ns}


def get_ns_id(subj, obj, net):
    """Get ns:id for both subj and obj

    Note: should *NOT* be used with PyBEL nodes

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


def get_ns_id_pybel_node(hgnc_sym, node):
    """

    Parameters
    ----------
    hgnc_sym : str
        Name to match
    node : CentralDogma|tuple
        PyBEL node or tuple of PyBEL nodes

    Returns
    -------
    tuple
        Tuple of ns, id for node
    """
    # If tuple of nodes, recursive call until match is found
    if isinstance(node, tuple):
        for n in node:
            ns, _id = get_ns_id_pybel_node(hgnc_sym, n)
            if ns is not None:
                return ns, _id
        logger.warning('None of the names in the tuple matched the HGNC '
                       'symbol')
        return None, None
    # If PyBEL node, check name match, return if match, else None tuple
    elif isinstance(node, CentralDogma):
        if node.name == hgnc_sym:
            try:
                return node.namespace, node.identifier
            except AttributeError:
                return None, None
    # Not recognized
    else:
        logger.warning(f'Type {node.__class__} not recognized')
        return None, None


def normalize_corr_names(corr_m: pd.DataFrame,
                         graph: Union[DiGraph, MultiDiGraph],
                         ns: str = None) -> pd.DataFrame:
    # todo:
    #  1. Move this function, together with get_ns_id,
    #     get_ns_id_pybel_node, normalize_entitites to net_functions
    #  2. Provide ns and id to the correlation matrix here too (requires
    #     overhaul of depmap script)
    #  3. Add support for pybel
    """

    Parameters
    ----------
    corr_m : pd.DataFrame
        A square pandas dataframe representing a correlation matrix. It is
        assumed that columns and indices are identical.
    graph : Union[DiGraph, MultiDiGraph]
        A graph to look in to see if the names are there
    ns : str
        The assumed namespace of the names in corr_m

    Returns
    -------
    pd.DataFrame
    """
    def _get_ns_id(n: str, g: Union[DiGraph, MultiDiGraph]) -> Tuple[str, str]:
        return g.nodes[n]['ns'], g.nodes[n]['id']

    col_names = corr_m.columns.values
    normalized_names = []
    for name in col_names:
        if name in graph.nodes:
            normalized_names.append(name)
        else:
            ns, _id, nn = gilda_normalization(name)
            if nn:
                normalized_names.append(nn)
            else:
                normalized_names.append(name)

    # Reset the normalized names
    corr_m.columns = normalized_names
    corr_m.index = normalized_names

    return corr_m


# Add new function to the tuple
expl_func_list = (explained, expl_ab, expl_ba, expl_axb, expl_bxa, find_cp,
                  get_sr, get_st, get_sd)
# Map the name of the function to a more human friendly column name
funcname_to_colname = {
    'explained': 'explained set',
    'expl_ab': 'a-b',
    'expl_ba': 'b-a',
    'expl_axb': 'a-x-b',
    'expl_bxa': 'b-x-a',
    'find_cp': 'common parent',
    'get_sr': 'shared regulator',
    'get_st': 'shared target',
    'get_sd': 'shared downstream'
}
expl_functions = {f.__name__: f for f in expl_func_list}
assert len(expl_func_list) == len(funcname_to_colname)
