import logging
import numpy as np
import pandas as pd
import networkx as nx
import itertools as itt
from math import ceil, log10
from collections import defaultdict
from sqlalchemy.exc import StatementError
from pandas.core.series import Series as pd_Series_class
from pandas.core.frame import DataFrame as pd_DataFrame_class
from indra.db import util as dbu
from indra.db import client as dbc
from indra.tools import assemble_corpus as ac
from indra.preassembler import Preassembler as pa
from indra.preassembler import hierarchy_manager as hm
from indra.sources.indra_db_rest import client_api as capi
from indra.sources.indra_db_rest.client_api import IndraDBRestError

db_prim = dbu.get_primary_db()
dnf_logger = logging.getLogger('DepMapFunctionsLogger')


def nx_graph_from_corr_pd_series(corr_sr, source='id1', target='id2',
                                 edge_attr='correlation', use_abs_corr=False):
    """Return a graph from a pandas sereis containing correlaton between gene A
    and gene B, using the correlation as edge weight

    corr_sr : pandas Series or DataFrame
        Pandas Series/DataFrame containing A, B corr

    source : str
        which column to identify as source node (output is still undirected)

    target : str
        which column to identify as target nodes (output is still undirected)

    edge_attr : int or str
        Column to use for edge attributes
    absolute : Bool
        Use absolute value as edge weight. Otherwise magnitude is used.

    Returns
    -------
    corr_weight_graph : nx.Graph
        An undirected, weighted, networkx graph
    """

    # check if corr_sr is series or dataframe
    if type(corr_sr) == pd_Series_class:
        dnf_logger.info('Converting Pandas Series to Pandas DataFrame')
        corr_df = pd.DataFrame(corr_sr).reset_index()
    else:
        corr_df = corr_sr

    corr_df = corr_df.rename(
        columns={'level_0': source, 'level_1': target, 0: edge_attr})

    if use_abs_corr:
        dnf_logger.info('Using absolute correlation values')
        corr_df = corr_df.apply(lambda c: c.abs() if np.issubdtype(
            c.dtype, np.number) else c)

    if type(edge_attr) is list or type(edge_attr) is bool:
        dnf_logger.warning('More than one attribute might be added to edges. '
                           'Resulting networkx graph might not be usable as '
                           'simple weighted graph.')
    dnf_logger.info('Creating weighted undirected graph from network data')
    corr_weight_graph = nx.from_pandas_dataframe(df=corr_df,
                                                 source=source,
                                                 target=target,
                                                 edge_attr=edge_attr)
    return corr_weight_graph


def nx_graph_from_corr_tuple_list(corr_list, use_abs_corr=False):
    """Return a graph from a list of edges, using the correlation as weight

    corr_list : list or iterator
        Edge tuples

    absolute : Bool
        Use absolute value as edge weight. Otherwise magnitude is used.

    Returns
    -------
    corr_weight_graph : nx.Graph
        An undirected, weighted, networkx graph
    """
    corr_weight_graph = nx.Graph()

    if use_abs_corr:
        dnf_logger.info('Using absolute correlation values')
        corr_list = map(lambda t: (t[0], t[1], abs(t[2])), corr_list)

    dnf_logger.info('Converting tuples to an edge bunch')
    edge_bunch = map(lambda t: (t[0], t[1], {'weight': t[2]}), corr_list)

    dnf_logger.info('Creating weighted undirected graph from network data')
    corr_weight_graph.add_edges_from(ebunch=edge_bunch)

    return corr_weight_graph


def agent_name_set(stmt):
    """Return the list of agent names in a statement.

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


def _uniq_evidence_count(stmt):
    """Count the number of evidences listed.

    stmt : indra statement

    Returns
    """

    return


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
    """Wrapper function for chaining preassembly statements meant to reduce
    the number of statements.

    stmts : list[:py:class:`indra.statements.Statement`]

    Returns
    -------
    stmts : list[:py:class:`indra.statements.Statement`]
        List of preassembled indra statements
    """

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
    """Formats information about statements and returns a string.

    subj : str
        indra statement subject
    obj : str
        indra statment object
    corr : float
        Correlation between subject and object
    stmts : list[:py:class:`indra.statements.Statement`]
        List of indra statements
    ignore_str : str
        String to ignore if it appears in the list of indra statements

    Returns
    -------
    output : str
        string with information about the statements that connect subject and
        object formatted for printing or for writing to a text file.
    """

    output = ''

    # Build up a string that shows explanations for each connection
    output = 'subj: %s; obj: %s; corr: %f \n' % (subj, obj, corr) + \
             'https://depmap.org/portal/interactive/?xDataset=Avana&xFeature' \
             '={}&yDataset=Avana&yFeature={}&colorDataset=lineage' \
             '&colorFeature=all&filterDataset=context&filterFeature=' \
             '&regressionLine=false&statisticsTable=false&associationTable=' \
             'true&plotOnly=false\n'.format(subj, obj)

    cp_stmts = stmts.copy()
    pa_stmts = deduplicate_stmt_list(stmts=cp_stmts, ignore_str=ignore_str)

    for stmt in pa_stmts:
        output += '- - - - - - - - - - - - - - - - - - - - - - - - - - - -\n'
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
                output += 'Evidence %i: ' % (count+1) + ev_text + '\n'

    # Add separator between each connection
    output += '\n\n#### #### #### #### #### ####\n'
    return output


def latex_output(subj, obj, corr, stmts, ev_len_fltr, ignore_str='parent'):
    """Compiles information about statements and returns a LaTeX
    formatted string that can be written to a file.

    subj : str
        indra statement subject
    obj : str
        indra statment object
    corr : float
        Correlation between subject and object
    stmts : list[:py:class:`indra.statements.Statement`]
        List of indra statements
    ignore_str : str
        String to ignore if it appears in the list of indra statements

    Returns
    -------
    output : str
        string with information about the statements that connect subject and
        object formatted for printing or for writing to a text file.
    """

    output = ''

    # Build up a string that shows explanations for each connection
    # This string is put in the script instead
    # output = r'\section{{{}, {}: {}}}'.format(subj, obj, corr) + '\n' + \
    #          r'See correlation plot \href{{' \
    #          r'https://depmap.org/portal/interactive/?xDataset=Avana' \
    #          '&xFeature' \
    #          '={}&yDataset=Avana&yFeature={}&colorDataset=lineage' \
    #          '&colorFeature=all&filterDataset=context&filterFeature=' \
    #          '&regressionLine=false&statisticsTable=false&associationTable=' \
    #          'true&plotOnly=false}}{{here}}'.format(subj, obj) + '\n\n'

    cp_stmts = stmts.copy()
    pa_stmts = deduplicate_stmt_list(stmts=cp_stmts, ignore_str=ignore_str)

    # HERE: insert subsection A->B
    output += r'\subsection{{{A} $\rightarrow$ {B}}}'.format(A=subj, B=obj)+'\n'

    # Sort stmts by evidence length
    # stmts_dict = dict()
    # ev_lens = []
    # for stmt in pa_stmts:
    #     stmts_dict[str(stmt)] = stmt
    #     ev_text_list = list(set(['N/A' if ev.text is None else ev.text for ev in
    #                        stmt.evidence]))
    #     # pdb.set_trace()
    #     if 'N/A' in ev_text_list:
    #         ev_text_list.remove('N/A')
    #     # Save tuple (len, str(stmt))
    #     ev_lens.append((len(ev_text_list), str(stmt)))
    # ev_lens.sort(key=lambda tup: tup[0], reverse=True)

    # HERE: itemize per statement type in stmt for loop
    output += r'\begin{itemize}'+'\n'
    for stmt in pa_stmts:
    # for lene, st_key in ev_lens:
    #     stmt = stmts_dict[st_key]
        if type(stmt) is str and str(stmt) == ignore_str:
                output += r'\item {s} and {o} are in the same complex ' \
                          r'or family'.format(s=subj, o=obj) + '\n'
        else:
            # Remove duplicate evidence text
            ev_text_set = set(['N/A' if ev.text is None else ev.text for ev in
                               stmt.evidence])
            ev_text_list = list(ev_text_set)
            if 'N/A' in ev_text_list:
                ev_text_list.remove('N/A')

            # assert lene == len(ev_text_list)

            if len(ev_text_list) >= ev_len_fltr:

                output += r'\item {nstmts} instances found of statement '\
                          r'{stmt}; supports count: {supc}; Supported by ' \
                          r'count: {supbc}'.format(stmt=str(stmt),
                                                   nstmts=len(ev_text_list),
                                                   supc=len(stmt.supports),
                                                   supbc=len(
                                                       stmt.supported_by))+'\n'

                # There are statements with zero length evidence lists
                if len(ev_text_list) > 0:
                    # HERE enumerate evidence text
                    output += r'\begin{enumerate}'+'\n'

                    max_ev = 25  # Dont ouput more than 25 evidences
                    for count, ev_text in enumerate(ev_text_list):
                        output += r'\item Evidence: \texttt{' + ev_text + \
                                  r'}' + '\n'
                        if count+1 == max_ev:
                            break

                    output += r'\end{enumerate}'+'\n'
            else:
                output += 'Evidence count below threshold of {}.\n'\
                    .format(ev_len_fltr)

    output += r'\end{itemize}' + '\n'

    # Don't forget to escape latex characters in 'ev_text'
    output = output.replace('_', '\_')\
        .replace('%', '\%')\
        .replace('&', '\&')\
        .replace('^', '\^')\
        .replace('~', '\~')

    return output.encode('ascii', 'ignore').decode('ascii')


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
            if counter % max([10, 10 ** ceil(log10(n_hgnc_ids)) // 100]) == 0:
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
