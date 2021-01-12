import logging

from typing import Tuple, Union
from indra.ontology.bio import bio_ontology
from indra.ontology.ontology_graph import IndraOntology
from indra.databases import get_identifiers_url

logger = logging.getLogger(__name__)


def get_all_entities(ontology=None):
    """Get a list of all entities included in an IndraOntology

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
    ontology.initialize()
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
        If provided, the parents must be in this set of ids. The set is
        assumed to contain valid ontology labels (see ontology.label()).

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
        If provided, the parents must be in this set of ids. The set is
        assumed to be valid ontology labels (see ontology.label()).

    Returns
    -------
    set
        set of common parents as (namespace, identifier) tuples
    """
    ontology = bio_ontology if not ontology else ontology
    return find_parent(ns=ns1, id_=id1, ontology=ontology,
                       immediate_only=immediate_only,
                       is_a_part_of=is_a_part_of) & \
        find_parent(ns=ns2, id_=id2, ontology=ontology,
                    immediate_only=immediate_only, is_a_part_of=is_a_part_of)


def has_common_parent(id1, id2, ns1='HGNC', ns2='HGNC',
                      ontology=None, immediate_only=False, is_a_part_of=None):
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
    is_a_part_of : iterable
        If provided, the parents must be in this set of ids. The set is
        assumed to be valid ontology labels (see ontology.label()).

    Returns
    -------
    bool
        True if id1 and id2 has one or more common parents.
    """
    ontology = bio_ontology if not ontology else ontology
    return bool(common_parent(id1, id2, ns1=ns1, ns2=ns2, ontology=ontology,
                              immediate_only=immediate_only,
                              is_a_part_of=is_a_part_of))


def ns_id_xref(from_ns: str, from_id: str, to_ns: str,
               ontology: IndraOntology = bio_ontology) \
        -> Union[Tuple[str, str], None]:
    """Get the id in another namespace given an ns-id pair

    Parameters
    ----------
    from_ns : str
        The namespace to translate from
    from_id : str
        The id in ns from_ns
    to_ns : str
        The namespace to find an id in
    ontology : IndraOntology
        The ontology to look in. Default: BioOnto

    Returns
    -------
    Union[Tuple[str, str], None]
        If found, a tuple of (ns, id), otherwise None.
    """
    return ontology.map_to(ns1=from_ns, id1=from_id, ns2=to_ns)


def ns_id_to_name(ns: str, _id: str, ontology: IndraOntology = bio_ontology) \
        -> Union[str, None]:
    """

    Parameters
    ----------
    ns
    _id
    ontology

    Returns
    -------
    Union[str, None]
        If found, the name is returned. Otherwise None is returned
    """
    return ontology.get_name(ns=ns, id=_id)
