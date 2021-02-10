"""
Explainer and helper functions for depmap_script2.py

Explanation functions must have the signature:
s, o, cor, net, _type, **kwargs

When adding new explanation functions, please also add them to the mapping
at the end
"""
import inspect
import logging
import networkx as nx
from typing import Set, Union, Tuple, List, Optional, Dict
from itertools import product

import pandas as pd
from networkx import DiGraph, MultiDiGraph
from pybel.dsl import CentralDogma

from indra.databases.hgnc_client import get_current_hgnc_id, get_uniprot_id
from depmap_analysis.util.io_functions import dump_it_to_pickle
from depmap_analysis.network_functions.famplex_functions import common_parent
from depmap_analysis.network_functions.net_functions import \
    gilda_normalization, INT_PLUS, INT_MINUS


__all__ = ['get_ns_id_pybel_node', 'get_ns_id', 'normalize_corr_names',
           'expl_functions', 'funcname_to_colname', 'apriori', 'axb_colname',
           'bxa_colname', 'ab_colname', 'ba_colname', 'st_colname',
           'sr_colname', 'sd_colname', 'cp_colname', 'react_colname',
           'react_funcname']

logger = logging.getLogger(__name__)


class FunctionRegistrationError(Exception):
    """Raise when a function does not adhere to the explainer function rules"""


def apriori_explained(s: str, o: str, corr: float,
                      net: Union[DiGraph, MultiDiGraph], _type: str, **kwargs)\
        -> Tuple[str, str, bool, Union[str, None]]:
    """A mock function that is used for a-priori explained pairs

    Parameters
    ----------
    s: str
        Subject node
    o: str
        Object node
    corr: float
        Correlation, either as [-1.0, 1.0] or z-score
    net: Union[DiGraph, MultiDiGraph]
        The indra graph used to explain the correlation between s and o
    _type: str
        The graph type used

    Returns
    -------
    Tuple[str, str, Union[str, None]]
        If explained, the tuple (s, o, explanation) is returned,
        where `explanation` is a string mapped from the subject or object.
    """
    expl_mapping: Union[Dict[str, str], None] = kwargs.get('expl_mapping')
    if expl_mapping is None:
        return s, o, False, None

    # Get a-priori explanations for s and o
    why_s = expl_mapping.get(s)
    why_o = expl_mapping.get(o)

    if why_s or why_o:
        explanation = f'{s}: {why_s}, {o}: {why_o}'
        return s, o, True, explanation
    else:
        return s, o, False, None


def common_reactome_paths(s: str, o: str, corr: float,
                          net: Union[DiGraph, MultiDiGraph],
                          _type: str, reactome_dict: Dict[str, List[str]],
                          **kwargs) \
        -> Tuple[str, str, bool, Union[None, List[str]]]:
    """Explain pair by matching common reactome pathways

    The pair is explained if they have any common reactome pathways

    Parameters
    ----------
    s: str
        Subject node
    o: str
        Object node
    corr: float
        Correlation, either as [-1.0, 1.0] or z-score
    net: Union[DiGraph, MultiDiGraph]
        The indra graph used to explain the correlation between s and o
    _type: str
        The graph type used
    reactome_dict: Dict[str, List[str]]
        Dict mapping gene names to reactome pathway identifiers. The dict is
        keyed by gene UP IDs, so s and o must be translated to UP.

    Returns
    -------

    """
    s_up = _get_upid_from_hgnc_symbol(s)
    if s_up is None:
        return s, o, False, None

    o_up = _get_upid_from_hgnc_symbol(o)
    if o_up is None:
        return s, o, False, None

    common_reactome = set(reactome_dict.get(s_up, [])) & \
        set(reactome_dict.get(o_up, []))
    return s, o, bool(common_reactome), common_reactome or None


def find_cp(s: str, o: str, corr: float, net: Union[DiGraph, MultiDiGraph],
            _type: str, **kwargs) -> Tuple[str, str, bool, Union[None,
                                                                 List[str]]]:
    """Explain pair by looking for ontological parents

    The pair is explained if the two entities have common ontological parents

    Parameters
    ----------
    s: str
        Subject node
    o: str
        Object node
    corr: float
        Correlation, either as [-1.0, 1.0] or z-score
    net: Union[DiGraph, MultiDiGraph]
        The indra graph used to explain the correlation between s and o
    _type: str
        The graph type used

    Returns
    -------
    Tuple[str, str, Union[None, List[str]]]
        A tuple of s, o and, if any, the list of common parents
    """
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
            return s, o, True, parents

    return s, o, False, None


def expl_axb(s: str, o: str, corr: float, net: Union[DiGraph, MultiDiGraph],
             _type: str, **kwargs) -> Tuple[str, str, bool, Union[None,
                                                                  List[str]]]:
    """Explain pair by looking for intermediate nodes connecting a to b

    The pair is considered explained if there is at least one node x
    connecting s with o in a directed sense: s -> x -> o

    Parameters
    ----------
    s: str
        Subject node
    o: str
        Object node
    corr: float
        Correlation, either as [-1.0, 1.0] or z-score
    net: Union[DiGraph, MultiDiGraph]
        The indra graph used to explain the correlation between s and o
    _type: str
        The graph type used

    Returns
    -------
    Tuple[str, str, Union[None, List[str]]]
        A tuple of s, o and a list of the nodes connecting s and o (if any)
    """
    s_succ = set(net.succ[s])
    o_pred = set(net.pred[o])
    # Filter ns
    if kwargs.get('ns_set'):
        ns_filt_args = (net, kwargs['ns_set'])
        s_succ = set(_node_ns_filter(s_succ, *ns_filt_args))
        o_pred = set(_node_ns_filter(o_pred, *ns_filt_args))
    # Filter sources
    if kwargs.get('src_set'):
        # Use reverse=False for downstream
        # net, reverse, allowed_src
        s_succ = _src_filter(s, s_succ, net, False, kwargs['src_set'])
        o_pred = _src_filter(o, o_pred, net, True, kwargs['src_set'])
    # Get intersection
    x_set = s_succ & o_pred

    # Sort out sign
    if _type in {'signed', 'pybel'}:
        x_nodes = _get_signed_interm(s, o, corr, net, x_set)
    else:
        x_nodes = x_set

    if x_nodes:
        return s, o, True, list(x_nodes)
    else:
        return s, o, False, None


def expl_bxa(s: str, o: str, corr: float, net: Union[DiGraph, MultiDiGraph],
             _type: str, **kwargs) -> Tuple[str, str, bool, Union[None,
                                                                  List[str]]]:
    """Reversal of expl_axb

    Parameters
    ----------
    s: str
        Subject node
    o: str
        Object node
    corr: float
        Correlation, either as [-1.0, 1.0] or z-score
    net: Union[DiGraph, MultiDiGraph]
        The indra graph used to explain the correlation between s and o
    _type: str
        The graph type used

    Returns
    -------
    Tuple[str, str, Union[None, List[str]]]
        A tuple of o, s and a list of the nodes connecting o and s (if any)
    """
    if _type == 'pybel':
        s_name = kwargs.pop('s_name')
        o_name = kwargs.pop('o_name')
        options = {'o_name': s_name, 's_name': o_name}
    else:
        options = {}
    u, v, expl, data = expl_axb(o, s, corr, net, _type, **kwargs, **options)
    return u, v, expl, data


# Shared regulator: A<-X->B
def get_sr(s: str, o: str, corr: float, net: Union[DiGraph, MultiDiGraph],
           _type: str, **kwargs) -> \
        Tuple[str, str, bool,
              Union[None, Tuple[List[str], List[str], List[str], List[str]]]]:
    """Explain pair by finding common upstream nodes

    The pair is explained if there is at least one common upstream node.
    However, the function returns data even if there are no common upstream
    nodes and there are predecessors to either of s or o.

    Parameters
    ----------
    s: str
        Subject node
    o: str
        Object node
    corr: float
        Correlation, either as [-1.0, 1.0] or z-score
    net: Union[DiGraph, MultiDiGraph]
        The indra graph used to explain the correlation between s and o
    _type: str
        The graph type used

    Returns
    -------
    Tuple[str, str, Union[None, Tuple[List[str], List[str], List[str],
    List[str]]]]
        A tuple of s, o and a tuple of the predecessors of s, o and their
        intersection and union
    """
    # Filter ns
    if kwargs.get('ns_set'):
        ns_filt_args = (net, kwargs['ns_set'])
        s_pred = set(_node_ns_filter(net.pred[s], *ns_filt_args))
        o_pred = set(_node_ns_filter(net.pred[o], *ns_filt_args))
    else:
        s_pred = set(net.pred[s])
        o_pred = set(net.pred[o])

    # Filter sources
    if kwargs.get('src_set'):
        # Use reverse=False for downstream
        # net, reverse, allowed_src
        src_args = (net, True, kwargs['src_set'])
        s_pred = _src_filter(s, s_pred, *src_args)
        o_pred = _src_filter(o, o_pred, *src_args)

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

    # Return x_nodes if strict, else if there is anything in union or s_pred
    # or o_pred
    strict = kwargs.get('strict_intermediates', False)
    if (strict and x_nodes) or (not strict and (s_pred or o_pred)):
        expl = True
    else:
        expl = False
    return s, o, expl, (list(s_pred), list(o_pred), list(x_nodes or []),
                        list(x_nodes_union or []))


# Shared target: A->X<-B
def get_st(s: str, o: str, corr: float, net: Union[DiGraph, MultiDiGraph],
           _type: str, **kwargs) -> \
        Tuple[str, str, bool,
              Union[None, Tuple[List[str], List[str], List[str], List[str]]]]:
    """Explain pair by finding common downstream nodes

    The pair is explained if there is at least one common downstream node.
    However, the function returns data even if there are no common downstream
    nodes and there are successors to either of s or o.

    Parameters
    ----------
    s: str
        Subject node
    o: str
        Object node
    corr: float
        Correlation, either as [-1.0, 1.0] or z-score
    net: Union[DiGraph, MultiDiGraph]
        The indra graph used to explain the correlation between s and o
    _type: str
        The graph type used

    Returns
    -------
    Tuple[str, str, Union[None, Tuple[List[str], List[str], List[str],
    List[str]]]]
        A tuple of s, o and a tuple of the successors of s, o and their
        intersection and union
    """
    s_succ = set(net.succ[s])
    o_succ = set(net.succ[o])
    # Filter ns
    if kwargs.get('ns_set'):
        ns_filt_args = (net, kwargs['ns_set'])
        s_succ = set(_node_ns_filter(s_succ, *ns_filt_args))
        o_succ = set(_node_ns_filter(o_succ, *ns_filt_args))
    # Filter sources
    if kwargs.get('src_set'):
        # Use reverse=False for downstream
        # net, reverse, allowed_src
        src_args = (net, False, kwargs['src_set'])
        s_succ = _src_filter(s, s_succ, *src_args)
        o_succ = _src_filter(o, o_succ, *src_args)
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

    # Return x_nodes if strict, else if there is anything in union or s_succ
    # or o_succ
    strict = kwargs.get('strict_intermediates', False)
    if (strict and x_nodes) or (not strict and (s_succ or o_succ)):
        expl = True
    else:
        expl = False

    return s, o, expl, (list(s_succ), list(o_succ), list(x_nodes or []),
                        list(x_nodes_union or []))


def get_sd(s: str, o: str, corr: float, net: Union[DiGraph, MultiDiGraph],
           _type: str, **kwargs) -> \
        Tuple[str, str, bool,
              Union[None, Tuple[List[str], List[str], List[str], List[str]]]]:
    """Explain pair by finding common downstream nodes two edges from s and o

    The pair is explained if there is at least one common node two edges
    away from s and o. However, the function returns data even if there
    are no common downstream nodes and there are successors to either of s
    or o.

    Parameters
    ----------
    s: str
        Subject node
    o: str
        Object node
    corr: float
        Correlation, either as [-1.0, 1.0] or z-score
    net: Union[DiGraph, MultiDiGraph]
        The indra graph used to explain the correlation between s and o
    _type: str
        The graph type used

    Returns
    -------
    Tuple[str, str, Union[None, Tuple[List[str], List[str], List[str],
    List[str]]]]
        A tuple of s, o and a tuple of the two-edge successors of s,
        o and their intersection and union
    """
    # Get nodes two edges away for subject
    args = (net, _type in {'signed', 'pybel'}, kwargs.get('ns_set'),
            kwargs.get('src_set'))
    s_y_set = _get_nnn_set(s, *args)
    o_y_set = _get_nnn_set(o, *args)

    # Get intersection and union of each nodes' 1st & 2nd layer neighbors
    y_set = s_y_set & o_y_set
    y_set_union = s_y_set | o_y_set

    if _type in {'signed', 'pybel'}:
        y_nodes = _get_signed_deep_interm(s, o, corr, net, y_set, False)
        y_nodes_union = _get_signed_deep_interm(s, o, corr, net, y_set_union,
                                                True)
    else:
        y_nodes = y_set
        y_nodes_union = y_set_union

    if y_nodes_union or s_y_set or o_y_set:
        s_y_list = set()
        o_y_list = set()
        if _type in {'signed', 'pybel'}:
            for _, sy in s_y_set:
                s_y_list.add(sy)
            for _, oy in o_y_set:
                o_y_list.update(oy)
        else:
            s_y_list = s_y_set
            o_y_list = o_y_set

        return s, o, True, (list(s_y_list or []), list(o_y_list or []),
                            list(y_nodes or []), list(y_nodes_union or []))
    else:
        return s, o, False, None


def expl_ab(s: str, o: str, corr: float, net: Union[DiGraph, MultiDiGraph],
            _type: str, **kwargs) -> Tuple[str, str, bool, Union[None,
                                                                 Tuple[List]]]:
    """Explain pair by checking for an edge between s and o

    The pair is explained if there exists and edge between s and o. The edge
    meta data is returned if the edge exists

    Parameters
    ----------
    s: str
        Subject node
    o: str
        Object node
    corr: float
        Correlation, either as [-1.0, 1.0] or z-score
    net: Union[DiGraph, MultiDiGraph]
        The indra graph used to explain the correlation between s and o
    _type: str
        The graph type used

    Returns
    -------
    Tuple[str, str, Union[None, Tuple[List]]]
        A tuple of s, o and, if the edge s-o exists, the edge meta data
    """
    edge_dict = _get_edge_statements(s, o, corr, net, _type, **kwargs)
    if edge_dict:
        return s, o, True, edge_dict.get('stmt_hash') if _type == 'pybel' \
            else edge_dict.get('statements')
    return s, o, False, None


def expl_ba(s: str, o: str, corr: float, net: Union[DiGraph, MultiDiGraph],
            _type: str, **kwargs) -> Tuple[str, str, bool, Union[None,
                                                                 Tuple[List]]]:
    """Reversal of expl_ab

    The pair is explained if there exists and edge between o and s. The edge
    meta data is returned if the edge exists

    Parameters
    ----------
    s: str
        Subject node
    o: str
        Object node
    corr: float
        Correlation, either as [-1.0, 1.0] or z-score
    net: Union[DiGraph, MultiDiGraph]
        The indra graph used to explain the correlation between s and o
    _type: str
        The graph type used

    Returns
    -------
    Tuple[str, str, Union[None, Tuple[List]]]
        A tuple of o, s and, if the edge o-s exists, the edge meta data
    """
    if _type == 'pybel':
        s_name = kwargs.pop('s_name')
        o_name = kwargs.pop('o_name')
        options = {'o_name': s_name, 's_name': o_name}
    else:
        options = {}
    u, v, expl, data = expl_ab(o, s, corr, net, _type, **kwargs, **options)
    return u, v, expl, data


def _get_edge_statements(s: str, o: str, corr: float,
                         net: Union[DiGraph, MultiDiGraph], _type: str,
                         **kwargs) -> \
        Dict[str, Union[
            str, int, float,
            List[Dict[str, Union[int, float, str, Dict[str, int]]]]
        ]]:
    """

    Parameters
    ----------
    s: str
        Subject node
    o: str
        Object node
    corr: float
        Correlation, either as [-1.0, 1.0] or z-score
    net: Union[DiGraph, MultiDiGraph]
        The indra graph used to explain the correlation between s and o
    _type: str
        The graph type used

    Returns
    -------
    Dict[str, Union[str, int, float,
                    List[Dict[str, Union[int, float, str, Dict[str, int]]]]
    ]]
    """
    if _type in {'signed', 'pybel'}:
        int_sign = INT_PLUS if corr >= 0 else INT_MINUS
        return net.edges.get((s, o, int_sign), None)
    else:
        return net.edges.get((s, o))


def _get_signed_interm(s: str, o: str, corr: float,
                       sign_edge_net: MultiDiGraph, x_set: Set[str]) -> \
        Set[str]:
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
    # Used for a<-x->b type relationships
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
    # Used for a->x<-b type relationships
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
    # Used for a->()->x<-()<-b type relationships
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

        # Add node that form paths with the correct sign
        if _approve_signed_paths(*args):
            x_approved.add(y)
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


def _get_nnn_set(n: str,
                 g: nx.MultiDiGraph,
                 signed: bool,
                 ns_set: Optional[Set[str]] = None,
                 src_set: Optional[Set[str]] = None) \
        -> Set[Union[str, Tuple[str, str]]]:
    # Filter node ns at all levels, only filter edge stmt sources at first edge
    # Todo make signed check here instead of late in the caller
    # For signed, have to keep pairs
    n_x_set = set()
    for x in g.succ[n]:
        if ns_set and g.nodes[x]['ns'].lower() not in ns_set:
            # Skip if x is not in allowed name space
            continue
        if src_set and not _src_in_edge(g.edges[(n, x)]['statements'],
                                        src_set):
            # Skip if there are no sources from the allowed set in edge
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
            if ns_set:
                n_x_set.update({
                    y for y in g.succ[x] if g.nodes[y]['ns'].lower() in ns_set
                })
            else:
                n_x_set.update(g.succ[x])
    return n_x_set


def _get_upid_from_hgnc_symbol(hgnc_gene: str) -> Union[str, None]:
    hgnc_id = get_current_hgnc_id(hgnc_gene)
    if isinstance(hgnc_id, list):
        ix = 0
        while True:
            try:
                up_id = get_uniprot_id(hgnc_id[ix])
            except IndexError:
                up_id = None
                break
            if up_id is None:
                ix += 1
    else:
        up_id = get_uniprot_id(hgnc_id)
    return up_id


def _node_ns_filter(node_list: Union[Set[str], List[str]],
                    net: Union[DiGraph, MultiDiGraph],
                    allowed_ns: Union[Set[str], List[str], Tuple[str]]) \
        -> Set[str]:
    return {x for x in node_list if net.nodes[x]['ns'].lower()
            in allowed_ns}


def _src_filter(start_node: str, neighbor_nodes: Set[str],
                net: Union[DiGraph, MultiDiGraph], reverse: bool,
                allowed_src: Set[str]) -> Set[str]:
    # Filter out nodes from 'nodes' if they don't have any sources from the
    # allowed sources

    # Make a sorted list of the neighbors
    node_list = sorted(neighbor_nodes)

    # Create an edge iterator with correct order
    edge_iter = \
        product(node_list, [start_node]) if reverse else \
        product([start_node], node_list)

    # Check which edges have the allowed sources
    filtered_nodes = set()
    for n, edge in zip(node_list, edge_iter):
        stmt_list = net.edges[edge]['statements']
        if _src_in_edge(stmt_list, allowed_src):
            filtered_nodes.add(n)
    return filtered_nodes


def _src_in_edge(
        stmt_list: List[Dict[str, Union[int, float, str, Dict[str, int]]]],
        allowed_src: Set[str]
) -> bool:
    """Assumes a list of stmt meta data as dicts"""
    for stmt_dict in stmt_list:
        # Catch empty source counts and check if any of the wanted sources
        if isinstance(stmt_dict['source_counts'], dict) \
                and any([s.lower() in allowed_src for s in
                         stmt_dict['source_counts']]):
            return True
    return False


def get_ns_id(subj: str, obj: str, net: Union[DiGraph, MultiDiGraph]) -> \
        Tuple[str, str, str, str]:
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
                         ns: str = None,
                         name_mapping: Dict[str, str] = None,
                         dump_mapping: bool = False,
                         dump_name: str = None) -> pd.DataFrame:
    # todo:
    #  - Move this function, together with get_ns_id,
    #    get_ns_id_pybel_node, normalize_entitites to net_functions
    #  - Provide ns and id to the correlation matrix here too (requires
    #    overhaul of depmap script)
    #  - Add support for pybel
    #  - If ns is provided loop through results and get the data for the
    #    matching name space, otherwise use the function used on Agent
    #    normalization
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
    name_mapping : Optional[Dict[str, str]]
        If provided try to map names from this dict
    dump_mapping : bool
        If True, save the mapping to pickle, Default: False.
    dump_name : Optional[str]
        The file path to save the mapping at

    Returns
    -------
    pd.DataFrame
    """
    def _get_ns_id(n: str, g: Union[DiGraph, MultiDiGraph]) -> Tuple[str, str]:
        return g.nodes[n]['ns'], g.nodes[n]['id']

    if name_mapping is None and dump_mapping and not dump_name:
        raise ValueError('Must provide file path with variable `dump_name` '
                         'if name mapping is dumped')

    col_names = corr_m.columns.values
    normalized_names = []
    mapping = {}
    for name in col_names:
        # If mapping is provided and name is in the mapping
        if name_mapping and name in name_mapping:
            normalized_names.append(name_mapping[name])
        # Otherwise use gilda
        else:
            if name in graph.nodes:
                normalized_names.append(name)
                # If we want to save the mapping
                if dump_mapping and name_mapping is None:
                    mapping[name] = name
            else:
                ns, _id, nn = gilda_normalization(name)
                if nn:
                    normalized_names.append(nn)
                    # If we want to save the mapping
                    if dump_mapping and name_mapping is None:
                        mapping[name] = nn
                else:
                    normalized_names.append(name)
                    # If we want to save the mapping
                    if dump_mapping and name_mapping is None:
                        mapping[name] = name

    # Reset the normalized names
    corr_m.columns = normalized_names
    corr_m.index = normalized_names

    if dump_mapping and mapping:
        dump_it_to_pickle(fname=dump_name, pyobj=mapping, overwrite=True)

    return corr_m


# Add new function to the tuple
expl_func_list = (apriori_explained, expl_ab, expl_ba, expl_axb, expl_bxa,
                  find_cp, get_sr, get_st, get_sd, common_reactome_paths)

# Map the name of the function to a more human friendly column name
funcname_to_colname = {
    'apriori_explained': 'apriori_explained',
    'expl_ab': 'a_b',
    'expl_ba': 'b_a',
    'expl_axb': 'a_x_b',
    'expl_bxa': 'b_x_a',
    'find_cp': 'common_parent',
    'get_sr': 'shared_regulator',
    'get_st': 'shared_target',
    'get_sd': 'shared_downstream',
    'common_reactome_paths': 'reactome_paths'
}

# Set reactome funcname
react_funcname = 'common_reactome_paths'

# Map function name to function
expl_functions = {f.__name__: f for f in expl_func_list}

# Set colnames to variables
apriori = funcname_to_colname['apriori_explained']
axb_colname = funcname_to_colname['expl_axb']
bxa_colname = funcname_to_colname['expl_bxa']
ab_colname = funcname_to_colname['expl_ab']
ba_colname = funcname_to_colname['expl_ba']
st_colname = funcname_to_colname['get_st']
sr_colname = funcname_to_colname['get_sr']
sd_colname = funcname_to_colname['get_sd']
cp_colname = funcname_to_colname['find_cp']
react_colname = funcname_to_colname['common_reactome_paths']

# Check that functions added to expl_func_list also exist in name to func map
try:
    assert len(expl_func_list) == len(funcname_to_colname)
except AssertionError:
    raise FunctionRegistrationError(
        'Missing function(s) in either expl_func_list or funcname_to_colname'
    )

# Check that function names are matched in column name mapping
try:
    assert set(expl_functions.keys()) == set(funcname_to_colname.keys())
except AssertionError:
    raise FunctionRegistrationError(
        'Function name mismatch between explanation functions tuple and '
        'function name to column name mapping'
    )

# Check that all functions have the fixed arg structure

for func in expl_func_list:
    try:
        assert \
            {'_type', 'corr', 'net', 'o', 's'}.issubset(
                set(inspect.signature(func).parameters.keys()))
    except AssertionError:
        raise FunctionRegistrationError(
            f'Function "{func.__name__}" does not have the required minimum '
            f'signature for its arguments. The required signature is '
            f'(s, o, corr, net, _type, **kwargs)'
        )

# Check that all func return str, str, bool as first three values

