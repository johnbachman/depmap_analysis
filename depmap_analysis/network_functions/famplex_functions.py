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


def find_parent(ontology=None, ns='HGNC',
                id_=None, type_='all'):
    """A wrapper function for IndraOntology.get_parents()

    Parameters
    ----------
    ontology : IndraOntology object
        An IndraOntology object. Default: INDRA BioOntology
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
    set[tuple]
        set of parents as (namespace, identifier) tuples
    """
    ontology = bio_ontology if not ontology else ontology
    return set(ontology.get_parents(ns, id_))


def common_parent(ontology=None, ns1='HGNC',
                  id1=None, ns2='HGNC', id2=None, type_='all'):
    """Returns the set of common parents.

    Parameters
    ----------
    ontology : IndraOntology object
        An IndraOntology object. Default: INDRA BioOntology
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
        set of common parents as (namespace, identifier) tuples
    """
    ontology = bio_ontology if not ontology else ontology
    return find_parent(ontology, ns1, id1, type_) & \
        find_parent(ontology, ns2, id2, type_)


def has_common_parent(ontology=None, ns1='HGNC', id1=None,
                      ns2='HGNC', id2=None, type='all'):

    """Returns True if id1 and id2 has at least one common parent.

    Parameters
    ----------
    ontology : IndraOntology object
        An IndraOntology object. Default: INDRA BioOntology
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
    ontology = bio_ontology if not ontology else ontology
    return bool(common_parent(ontology, ns1, id1, ns2, id2, type))


