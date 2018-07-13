from indra.preassembler import hierarchy_manager as hm
from indra.preassembler import Preassembler as pa
from indra.tools import assemble_corpus as ac
from indra.sources.indra_db_rest import client_api as capi
from indra.sources.indra_db_rest.client_api import IndraDBRestError
from collections import defaultdict
from math import ceil, log10
import itertools as itt
import logging
from indra.db import client as dbc
from indra.db import util as dbu
from sqlalchemy.exc import StatementError
import pdb
db_prim = dbu.get_primary_db()
dnf_logger = logging.getLogger('DepMapFunctionsLogger')


def agent_name_set(stmt):
    """Returns the list of agent names in a statement.

    stmt : :py:class:`indra.statements.Statement`

    Returns
    -------
    ags : list[agent names]

    """
    ags = []
    try:
        ags.update(list(map(lambda ag: ag.name, stmt.agent_list())))
    except AttributeError:
        for ag in stmt.agent_list():
            if ag is None:
                pass
            else:
                ags.append(ag.name)
    return ags


def nested_dict_gen(stmts):
    """Generates a nested dict of the form dict[key1][key2] = [statement list]
    from INDRA statements.

    stmts :  list[:py:class:`indra.statements.Statement`]
        List or set of INDRA statements to find connections in

    Returns
    -------
    stmts_dict : collections.defaultdict
         dict of the form dict[subj][obj] = list[stmts]
    """

    nested_stmt_dicts = defaultdict(dict)

    count = 0
    for st in stmts:
        count += 1
        # NOTE: If statement is complex, it migth have more than two agents
        # and the agents won't be distinguishable as subject,object

        agent_list = agent_name_set(stmt=st)

        # It takes two agents to tango
        if len(agent_list) > 1:

            # Is complex or selfmodification
            if st.to_json()['type'].lower in ['complex', 'selfmodification']:
                for agent, other_agent in itt.permutations(agent_list, r=2):
                    try:
                        nested_stmt_dicts[agent][other_agent].append(st)
                    except KeyError:  # If pair does not exist yet
                        nested_stmt_dicts[agent][other_agent] = [st]

            # Non-complex interaction
            else:
                subj = agent_list[0]
                obj = agent_list[1]

                try:
                    nested_stmt_dicts[subj][obj].append(st)
                except KeyError:
                    nested_stmt_dicts[subj][obj] = [st]

            # Check common parent (same familiy or complex)
            for agent, other_agent in itt.permutations(agent_list, r=2):
                if has_common_parent(id1=agent, id2=other_agent):
                    try:
                        if 'parent' not in \
                                nested_stmt_dicts[agent][other_agent]:
                            nested_stmt_dicts[agent][other_agent].append(
                                'parent')
                    except KeyError:
                        nested_stmt_dicts[agent][other_agent] = ['parent']

        # Ignore when we only have one agent
        else:
            continue

    dnf_logger.info('Created nested dict of length %i from %i statements.' %
                    (len(nested_stmt_dicts), len(stmts)))
    return nested_stmt_dicts


def deduplicate_stmt_list(stmts, ignore_str):
    """Takes a list of statements list[stmts] and runs
    indra.preassembler.Preassembler.combine_duplicate_stmts() while also
    taking care of non-statements

    stmts : list[:py:class:`indra.statements.Statement`]

    Returns
    -------
    stmts : list[:py:class:`indra.statements.Statement`]
         List of preassembled statments possibly including a non-statements
    """
    # subjects should be the outer keys and objects should be the inner

    if ignore_str in stmts:
        only_stmt_list = [s for s in stmts if type(s) is not str]
        stmts = pa_filter_unique_evidence(only_stmt_list)
        stmts += [ignore_str]
    else:
        stmts = pa_filter_unique_evidence(stmts)
    return stmts


def pa_filter_unique_evidence(stmts):

    # Ground statemtens:
    grounded_stmts = ac.map_grounding(stmts)

    # Use curated site information to standardize modification sites in stmts
    ms_stmts = ac.map_sequence(grounded_stmts)

    # Compiles together raw statements to one statement per type
    opa_stmts = ac.run_preassembly(ms_stmts, return_toplevel=False)
    return opa_stmts


def _old_str_output(subj, obj, corr, stmts, ignore_str='parent'):

    # Build up a string that shows explanations for each connection
    output = 'subj: %s; obj: %s; corr: %f \n' % (subj, obj, corr)

    if ignore_str in stmts:
        pure_list = [s for s in stmts if type(s) is not str]
        types = relation_types(pure_list)
        types += [ignore_str] * stmts.count(ignore_str)
    else:
        types = relation_types(stmts)

    cp_stmts = stmts.copy()
    dedupl_stmts = deduplicate_stmt_list(cp_stmts, ignore_str)

    types_set = set(types)
    types_sstmt = []
    for tp in types_set:
        for st in dedupl_stmts:
            if type(st) is not str:
                if st.to_json()['type'] == tp:
                    types_sstmt.append((tp, str(st)))
            elif type(st) is str and tp is ignore_str:
                types_sstmt.append((tp, str(st)))

    for tp, str_stmt in types_sstmt:
        if tp is not ignore_str:
            output += '- - - - - - - - - - - - - - - - - - - - - - - - - - - -'
            output += '\nInstances found of statement %s: %i\n' % \
                      (str_stmt, types.count(tp))
        for stmt in stmts:
            if type(stmt) is str and str(stmt) == ignore_str:
                output += '%s and %s are in the same complex or family\n' % \
                          (subj, obj)
            elif type(stmt) is not str and stmt.to_json()['type'] == tp:
                output += 'Evidence for uuid %s: ' % stmt.uuid
                ev = stmt.evidence[0].text
                output += ('N/A' if ev is None else ev)+'\n'
            else:
                continue

    # Add separator between each connection
    output += '\n\n#### #### #### #### #### ####\n'
    return output


def str_output(subj, obj, corr, stmts, ignore_str='parent'):

    # Build up a string that shows explanations for each connection
    output = 'subj: %s; obj: %s; corr: %f \n' % (subj, obj, corr)

    cp_stmts = stmts.copy()
    pa_stmts = deduplicate_stmt_list(stmts=cp_stmts, ignore_str=ignore_str)

    for stmt in pa_stmts:
        output += '- - - - - - - - - - - - - - - - - - - - - - - - - - - -'
        if type(stmt) is str and str(stmt) == ignore_str:
            output += '%s and %s are in the same complex or family\n' % \
                      (subj, obj)
        else:
            # Remove duplicate evidence text
            ev_text_set = set(['N/A' if ev.text is None else ev.text for ev in
                               stmt.evidence])
            ev_text_list = list(ev_text_set)
            if 'N/A' in ev_text_list:
                ev_text_list.remove('N/A')

            output += '\nInstances found of statement %s: %i; supports ' \
                      'count: %i; Supported by count: %i\n' % \
                      (str(stmt), len(ev_text_list), len(stmt.supports),
                       len(stmt.supported_by))
            for count, ev_text in enumerate(ev_text_list):
                try:
                    output += 'Evidence %i: ' % (count+1) + ev_text + '\n'
                except TypeError:
                    pdb.set_trace()

    # Add separator between each connection
    output += '\n\n#### #### #### #### #### ####\n'
    return output


def dbc_load_statements(hgnc_ids):
    """Load statements where hgnc id is subject or object from indra.db.client

    Parameters
    ----------
    hgnc_ids : iterable
        An iterable containing HGNC ids

    Returns
    -------
    stmts : set{:py:class:`indra.statements.Statement`}
        A set of all retrieved INDRA statemetents containing HGNC id
    """
    stmts = set()
    counter = 0
    n_hgnc_ids = len(hgnc_ids)
    try:
        for hgnc_id in hgnc_ids:
            stmts.update(dbc.get_statements_by_gene_role_type(agent_id=hgnc_id,
                                                              db=db_prim,
                                                              preassembled=
                                                              False,
                                                              fix_refs=False))
            counter += 1
            if counter % max(10, 10 ** ceil(log10(n_hgnc_ids)) // 100) == 0:
                dnf_logger.info(' : : : Finished %i queries out of %i '
                                ': : :' % (counter, n_hgnc_ids))

    except KeyboardInterrupt as e:
        db_prim.session.rollback()
        raise e
    except StatementError as e:
        db_prim.session.rollback()
        raise e
    return stmts


def find_parent(ho=hm.hierarchies['entity'], ns='HGNC',
                id=None, type='all'):
    """A wrapper function for he.get_parents to make the functionality more
    clear.

    Parameters
    ----------
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


def direct_relation(id1, id2, long_stmts=set()):
    """Returns a list of INDRA statements

    Parameters
    ----------
    id1/id2 : str
        Strings of the two ids to check a direct relation between.
    long_stmts : set[:py:class:`indra.statements.Statement`]
        (Optional) List or set of INDRA statements to find connections in

    Returns
    -------
    stmts : list[:py:class:`indra.statements.Statement`]
        List of INDRA statements that directly relate id1 and id2
    """
    if not long_stmts:
        stmts = direct_relation_from_api(id1=id1, id2=id2)
    else:
        stmts = direct_relation_from_stmts(id1=id1, id2=id2,
                                           stmts_in=long_stmts)
    return stmts


def direct_relation_from_api(id1, id2, on_limit='sample'):
    """Queries INDRA DB for Statements linking two genes and returns a list
    containing the matching statements.

    Parameters
    ----------
    id1/id2 : str
        Strings of the two ids to check a direct relation between.
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
    try:
        stmts = capi.get_statements(subject=id1, object=id2, on_limit=on_limit)
        stmts + capi.get_statements(subject=id2, object=id1, on_limit=on_limit)
    except IndraDBRestError:
        stmts = capi.get_statements(subject=id1 + '@TEXT', object=id2 + '@TEXT',
                                    on_limit=on_limit)
        stmts + capi.get_statements(subject=id2 + '@TEXT', object=id1 + '@TEXT',
                                    on_limit=on_limit)
    return stmts


def direct_relation_from_stmts(id1, id2, stmts_in):
    """Returns a list of INDRA statements that connect id1 and id2 queried
    from a provided list of statements,

    Parameters
    ----------
    id1/id2 : str
        Strings of the two ids to check a direct relation between.
    stmts_in : set[:py:class:`indra.statements.Statement`]
        List of INDRA statements to find connections in.

    Returns
    -------
    stmts_out : list[:py:class:`indra.statements.Statement`]
        List of INDRA statements that directly relate id1 and id2
    """
    target_ag = {id1, id2}
    stmts_out = []
    for stms in stmts_in:
        s_agents = agent_name_set(stms)
        if target_ag.issubset(s_agents):
            stmts_out.append(stms)
    return stmts_out


def relation_type(indra_stmt):
    """Return the statement type in an INDRA statement as a string.

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
    """Returns the corresponding list of INDRA Statement types associated
    with a list of Statements.

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


def has_direct_relation(id1, id2, long_stmts=set()):
    """Indicates whether two genes are linked by Statements in the INDRA data
    bases.

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
    return bool(direct_relation(id1, id2, long_stmts=long_stmts))


def are_connected(id1, id2, long_stmts=set()):
    """Indicates whether two genes have a connection either through a direct
    relation or a through a common parent.

    Parameters
    ----------
    id1/i2 : str
        HGNC id

    Returns
    -------
    bool
        True if the two HGNC ids either have a common parent or if they have a
        direct relation found in the indra.sources.indra_db_rest.client_api
        databases.
    """
    return has_common_parent(ns1='HGCN', id1=id1, ns2='HGCN', id2=id2) or \
        has_direct_relation(id1=id1, id2=id2, long_stmts=long_stmts)


def connection_types(id1, id2, long_stmts=set()):
    """Returns a list of the connection types linking two genes.

    Parameters
    ----------
    id1/i2 : str
        HGNC id

    Returns
    -------
    ctypes : list[type]
        Returns a list of connection types.
        `[]` - empty list if no connections.
        Type is any of:
        `INDRA statement` - Any INDRA statement type
        `parent` - id1 and id2 are connected through common parent(s)
    """

    ctypes = relation_types(direct_relation(id1=id1, id2=id2,
                                            long_stmts=long_stmts))
    if has_common_parent(id1=id1, id2=id2):
        ctypes += ['parent']
    return ctypes
