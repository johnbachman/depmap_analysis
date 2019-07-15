from indra.preassembler import hierarchy_manager as hm


def get_all_entities(eh=hm.hierarchies['entity']):
    """Get a list of all entities included in HierarchyManager['entity']

    Parameters
    ----------
    eh : HierarchyManager object
        A HierarchyManager object initialized to an entities HierarchyManager.

    Returns
    -------
    entity_list : list
        A list of namespace, id, uri_id tuples
    """
    ent_list = []
    eh.initialize()
    entities_keyed_by_parent = eh._children
    for puri, children_uris in entities_keyed_by_parent.items():
        pns, pid = eh.ns_id_from_uri(puri)
        ent_list.append((pns, pid, puri))
        ns_id_child_set = set([(*eh.ns_id_from_uri(p), p)
                               for p in entities_keyed_by_parent[puri]])
        for cns, cid, curi in ns_id_child_set:
            ent_list.append((cns, cid, curi))

    return ent_list


def find_parent(ho=hm.hierarchies['entity'], ns='HGNC',
                id_=None, type_='all'):
    """A wrapper function for he.get_parents to make the functionality more
    clear.

    Parameters
    ----------
    ho : HierarchyManager object
        A HierarchyManager object. Default: entity hierarchy object
    ns : str
        namespace id. Default: HGNC
    id_ : str
        id to check parents for. Default: None
    type_ : str
        'all': (Default) return all parents irrespective of level;
        'immediate': return only the immediate parents;
        'top': return only the highest level parents

    Returns
    -------
    set
        set of parents of database id in namespace ns
    """
    return ho.get_parents(ho.get_uri(ns, id_), type_)


def common_parent(ho=hm.hierarchies['entity'], ns1='HGNC',
                  id1=None, ns2='HGNC', id2=None, type_='all'):
    """Returns the set of common parents.

    Parameters
    ----------
    ho : HierarchyManager object
        A HierarchyManager object. Default: entity hierarchy object
    ns1 : str
        namespace id. Default: HGNC
    id1 : str
        First id to check parents for. Default: None
    ns2 : str
        namespace id. Default: HGNC
    id2 : str
        Second id to check parents for. Default: None
    type_ : str
        'all': (Default) return all parents irrespective of level;
        'immediate': return only the immediate parents;
        'top': return only the highest level parents

    Returns
    -------
    set
        set of common parents in uri format
    """
    return find_parent(ho, ns1, id1, type_) & find_parent(ho, ns2, id2, type_)


def has_common_parent(ho=hm.hierarchies['entity'], ns1='HGNC', id1=None,
                      ns2='HGNC', id2=None, type='all'):

    """Returns True if id1 and id2 has at least one common parent.

    Parameters
    ----------
    ho : HierarchyManager object
        A HierarchyManager object. Default: entity hierarchy object
    ns1 : str
        namespace id. Default: HGNC
    id1 : str
        First id to check parents for. Default: None
    ns2 : str
        namespace id. Default: HGNC
    id2 : str
        Second id to check parents for. Default: None
    type : str
        'all': return all parents irrespective of level;
        'immediate': return only the immediate parents;
        'top': return only the highest level parents

    Returns
    -------
    bool
        True if hgnc1 and hgnc2 has one or more common parents.
    """
    return bool(common_parent(ho, ns1, id1, ns2, id2, type))


