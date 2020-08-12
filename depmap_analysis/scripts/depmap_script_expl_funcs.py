"""Explainer and helper functions for depmap_script2.py"""
import logging
from pybel.dsl import CentralDogma

from depmap_analysis.network_functions.famplex_functions import common_parent
from depmap_analysis.network_functions.net_functions import ns_id_from_name, \
    INT_PLUS, INT_MINUS


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
        s_ns, s_id, o_ns, o_id = get_ns_id(s, o, net)

    if not s_id:
        s_ns, s_id = ns_id_from_name(s_name) if _type == 'pybel' else \
            ns_id_from_name(s)
    if not o_id:
        o_ns, o_id = ns_id_from_name(o_name) if _type == 'pybel' else \
            ns_id_from_name(o)

    if s_id and o_id:
        parents = list(common_parent(ns1=s_ns, id1=s_id, ns2=o_ns, id2=o_id))
        if parents:
            if kwargs.get('ns_set'):
                parents = {(ns, _id) for ns, _id in parents if ns.lower() in
                           kwargs['ns_set']} or None
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
    x_set = set(net.pred[s]) & set(net.pred[o])

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


# Shared target: A->X<-B
def get_st(s, o, corr, net, _type, **kwargs):
    x_set = set(net.succ[s]) & set(net.succ[o])

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


def get_sd(s, o, corr, net, _type, **kwargs):
    # Get next-nearest-neighborhood for subject
    s_x_set = set()
    for x in net.succ[s]:
        # If signed, add edges instead and match sign in helper
        if _type in {'signed', 'pybel'}:
            for y in net.succ[x]:
                s_x_set.add((x, y))
        # Just add nodes
        else:
            s_x_set.add(x)
            s_x_set.update(net.succ[x])
    # Get next-nearest-neighborhood for object
    o_x_set = set()
    for x in net.succ[o]:
        # If signed, add edges instead and match sign in helper
        if _type in {'signed', 'pybel'}:
            for y in net.succ[x]:
                o_x_set.add((x, y))
        else:
            o_x_set.add(x)
            o_x_set.update(net.succ[x])

    # Get intersection of each nodes' 1st & 2nd layer neighbors
    x_set = s_x_set & o_x_set

    if _type in {'signed', 'pybel'}:
        x_nodes = _get_signed_deep_interm(s, o, corr, net, x_set)
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


def _get_signed_deep_interm(s, o, corr, sign_edge_net, xy_set):
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

        # Add nodes that form paths with the correct sign
        if path_sign == INT_PLUS:
            if (sx_plus and xy_plus or sx_minus and xy_minus) and \
                    (ox_plus and xy_plus or ox_minus and xy_minus):
                x_approved.update({x, y})
        else:
            if (sx_plus and xy_minus or sx_minus and xy_plus) and \
                    (ox_plus and xy_minus or ox_minus and xy_plus):
                x_approved.update({x, y})

    return x_approved


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