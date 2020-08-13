import logging

from indra.ontology.bio import bio_ontology
from indra.databases import get_identifiers_url

logger = logging.getLogger(__name__)


def get_all_entities(ontology=None):
    """Get a list of all entities included in HierarchyManager['entity']

    Parameters
    ----------
    ontology : IndraOntology object
        An IndraOntology object. Default: INDRA BioOntology

    Returns
    -------
    entity_list : list
        A list of namespace, id, uri_id tuples
    """
    ontology = bio_ontology if not ontology else ontology
    ent_list = []
    bio_ontology.initialize()
    for node in ontology.nodes:
        db_ns, db_id = ontology.get_ns_id(node)
        if db_ns in {'FPLX', 'HGNC'}:
            ent_list.append((db_ns, db_id, get_identifiers_url(db_ns, db_id)))
    return ent_list


def find_parent(id_, ns='HGNC', ontology=None, immediate_only=False,
                is_a_part_of=None):
    """A wrapper function for IndraOntology.get_parents()

    Parameters
    ----------
    id_ : str
        id to check parents for
    ns : str
        namespace id. Default: HGNC
    ontology : IndraOntology object
        An IndraOntology object. Default: INDRA BioOntology
    immediate_only : bool
        Determines if all or just the immediate parents should be returned
    is_a_part_of : iterable
        If provided, the parents must be in this set of ids

    Returns
    -------
    set[tuple]
        set of parents as (namespace, identifier) tuples
    """
    ontology = bio_ontology if not ontology else ontology

    if immediate_only:
        parents = {p for p in ontology.child_rel(ns, id_, {'isa', 'partof'})}
    else:
        parents = set(ontology.get_parents(ns, id_))

    if is_a_part_of:
        parents = {p for p in parents if p[1] in is_a_part_of}

    return parents


def common_parent(id1, id2, ns1='HGNC', ns2='HGNC', ontology=None,
                  immediate_only=False, is_a_part_of=None):
    """Returns the set of common parents.

    Parameters
    ----------
    id1 : str
        First id to check parents for
    id2 : str
        Second id to check parents for
    ns1 : str
        namespace id. Default: HGNC
    ns2 : str
        namespace id. Default: HGNC
    ontology : IndraOntology object
        An IndraOntology object. Default: INDRA BioOntology
    immediate_only : bool
        Determines if all or just the immediate parents should be returned.
        Default: False, i.e. all parents.
    is_a_part_of : iterable
        If provided, the parents must be in this set of ids

    Returns
    -------
    set
        set of common parents as (namespace, identifier) tuples
    """
    ontology = bio_ontology if not ontology else ontology
    return find_parent(ontology, ns1, id1, immediate_only, is_a_part_of) & \
        find_parent(ontology, ns2, id2, immediate_only, is_a_part_of)


def has_common_parent(id1, id2, ns1='HGNC', ns2='HGNC',
                      ontology=None, immediate_only=False):
    """Returns True if id1 and id2 has at least one common parent.

    Parameters
    ----------
    id1 : str
        First id to check parents for
    id2 : str
        Second id to check parents for
    ns1 : str
        namespace id. Default: HGNC
    ns2 : str
        namespace id. Default: HGNC
    ontology : IndraOntology object
        An IndraOntology object. Default: INDRA BioOntology
    immediate_only : bool
        Determines if all or just the immediate parents should be returned.
        Default: False, i.e. all parents.

    Returns
    -------
    bool
        True if id1 and id2 has one or more common parents.
    """
    ontology = bio_ontology if not ontology else ontology
    return bool(common_parent(ontology, ns1, id1, ns2, id2,
                              immediate_only=immediate_only))
