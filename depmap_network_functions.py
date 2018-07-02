from indra.preassembler import hierarchy_manager as hm
from indra.sources.indra_db_rest import client_api as capi


def find_parent(ho=hm.hierarchies['entity'], ns='HGNC',
                id=None, type='all'):
    """A wrapper function for he.get_parents to make the functionilty more
    clear.

    ho : HierarchyManager object
        A HierarchyManager object. Default: entity hierarchy object
    ns : str
        namespace id. Default: HGNC
    id : str
        id to check parents for. Default: None
    type : str
        'all': (Default) return all parents irrespective of level;
        'immediate': return only the immediate parents;
        'top': return only the highest level parents

    Returns
    -------
    set
        set of parents of database id in namespace ns
    """
    return ho.get_parents(ho.get_uri(ns, id), type)


def common_parent(ho=hm.hierarchies['entity'], ns1='HGNC',
                  id1=None, ns2='HGNC', id2=None, type='all'):
    """Returns the set of common parents.

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
        'all': (Default) return all parents irrespective of level;
        'immediate': return only the immediate parents;
        'top': return only the highest level parents

    Returns
    -------
    set
        set of common parents in uri(?) format  # ToDo Format name is uri?
    """
    return find_parent(ho, ns1, id1, type) & find_parent(ho, ns2, id2, type)


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


def direct_relation(id1, id2, on_limit='sample'):
    """Queries INDRA DB for Statements linking two genes.

    Parameters
    ----------
    id1/id2 : str
        Strings of the two ids to check a direct relation between.
        Default: None
    on_limit : str
        There are four options for handling the a query that is to large:
        `sample` - (default) take a sample of statements from the result,
        `truncate` - simply return the first 10,000 statements of the result,
        `error` - raise an error if the query is too large, or
        `persist` - perform as many queries as needed to get all the statements.
        Note that this last option generally takes much much longer to execute

    Returns
    -------
    stmts : list[:py:class:`indra.statements.Statement`]
        A list of INDRA Statement instances.
    """
    stmts = capi.get_statements(subject=id1, object=id2, on_limit=on_limit)
    stmts + capi.get_statements(subject=id2, object=id1, on_limit=on_limit)
    return stmts


def relation_type(indra_stmt):
    """
    Return the relation type in an INDRA statement as a string.

    Parameters
    ----------
    indra_stmt : :py:class:`indra.statements.Statement`

    Returns
    -------
    relation type : str
        A string containing an INDRA relation type
    """
    return indra_stmt.to_json()['type']


def relation_types(stmts):
    """Returns INDRA Statement types associated with a list of Statements.

    Parameters
    ----------
    stmts : list[:py:class:`indra.statements.Statement`]
        A list of INDRA Statement instances

    Returns
    -------
    types : list[INDRA statement types]
        A list of strings containing the INDRA statement types
    """
    types = []
    for stmt in stmts:
        types.append(relation_type(stmt))
    return types


def has_direct_relation(id1, id2):
    """Indicates whether two genes are linked by Statements in the INDRA DB.

    Parameters
    ----------
    id1/id2 : str
        HGNC names for the two genes.

    Returns
    -------
    bool
        True if the HGNC ids has a direct relation found in the
        indra.sources.indra_db_rest.client_api databases.
    """
    return bool(direct_relation(id1, id2))


def are_connected(id1, id2):
    """Indicates whether two genes have a connection.

    Parameters
    ----------
    id1/i2 : str
        HGNC id

    Returns
    -------
    bool
        True if the two HGNC ids either have a common parent or if they have a
        direct found in the indra.sources.indra_db_rest.client_api databases.
    """
    return has_common_parent(ns1='HGCN', id1=id1, ns2='HGCN', id2=id2) or \
           has_direct_relation(id1=id1, id2=id2)


def connection_type(id1, id2):
    """Indicates the connection type linking two genes.

    Parameters
    ----------
    id1/i2 : str
        HGNC id

    Returns
    -------
    ctype : str
        Returns one of four options:
        `None` - No connection.
        `direct` - id1 and id2 are directly connected in the
            indra.sources.indra_db_rest.client_api databases, but are *not*
            connected by family or complex.
        `parent` - id1 and id2 are connected by family or complex, but are
            *not* directly connected in the
            indra.sources.indra_db_rest.client_api databases.
        `both` - id1 and id2 are connected both directly and by family/complex.
    """

    ctype = None

    if has_direct_relation(id1=id1, id2=id2) and \
       not has_common_parent(id1=id1, id2=id2):
        ctype = 'direct'

    elif has_common_parent(id1=id1, id2=id2) and \
         not has_direct_relation(id1=id1, id2=id2):
        ctype = 'parent'

    elif has_common_parent(id1=id1, id2=id2) and \
         has_direct_relation(id1=id1, id2=id2):
        ctype = 'both'

    return ctype
