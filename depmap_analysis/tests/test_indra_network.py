import logging
import unittest
import numpy as np
from random import random as rnd
from collections import defaultdict

import indra_db.tests.util as tu
from indra_db.util.dump_sif import load_db_content, make_ev_strata, \
    make_dataframe, NS_LIST

from depmap_analysis.network_functions.net_functions import \
    sif_dump_df_to_digraph
from depmap_analysis.network_functions.indra_network import IndraNetwork

logger = logging.getLogger('IndraNetworkSearch test')

# Get db
db = tu.get_db_with_views(1000)

# Get stratified evidence and belief scores
sed = make_ev_strata(pkl_filename=None, db=db)

# Get dataframe
df = make_dataframe(reconvert=True,
                    db_content=load_db_content(reload=True,
                                               ns_list=NS_LIST,
                                               pkl_filename=None,
                                               db=db),
                    pkl_filename=None)

# Create fake belief dict
bsd = {}
for n, h in df['stmt_hash'].iteritems():
    bsd[h] = rnd()

# Add custom row to df that can be checked later
test_edge = ('GENE_A', 'GENE_B')
test_medge = (*test_edge, 0)
test_node = test_edge[0]
test_hash = 1234567890
test_type = 'TestStatement'
test_row = {
    'agA_ns': 'TEST', 'agA_id': '1234', 'agA_name': test_edge[0],
    'agB_ns': 'TEST', 'agB_id': '2345', 'agB_name': test_edge[1],
    'stmt_type': test_type, 'evidence_count': 1, 'stmt_hash': test_hash
}
test_source = 'pc11'
test_evidence = {test_source: 1}
test_belief = 0.987654321
df = df.append(test_row,
    ignore_index=True)
sed[test_hash] = test_evidence
bsd[test_hash] = test_belief


class TestNetwork(unittest.TestCase):
    def setUp(self):
        self.df = df
        self.indra_network = IndraNetwork(
            indra_dir_graph=sif_dump_df_to_digraph(
                df=self.df, belief_dict=bsd, strat_ev_dict=sed, multi=False,
                include_entity_hierarchies=True),
            indra_multi_dir_graph=sif_dump_df_to_digraph(
                df=self.df, belief_dict=bsd, strat_ev_dict=sed, multi=True,
                include_entity_hierarchies=True)
        )
        self.indra_network.verbose = 2

    def test_network_search(self):
        query = {
            'source': test_edge[0],
            'target': test_edge[1],
            'stmt_filter': [],
            'edge_hash_blacklist': [],
            'node_filter': ['test'],
            'node_blacklist': [],
            'path_length': False,
            'sign': 'no_sign',
            'weighted': False,
            'bsco': 0.0,
            'curated_db_only': False,
            'fplx_expand': False,
            'k_shortest': 1,
            'two_way': True
        }

        result = self.indra_network.handle_query(**query)
        logger.info('Got result: %s' % result)
        assert result['timeout'] is False
        assert isinstance(result['paths_by_node_count'], (dict, defaultdict))
        assert 2 in result['paths_by_node_count']['forward']
        assert not result['paths_by_node_count']['backward']
        assert len(result['paths_by_node_count']['forward'][2]) == 1
        assert isinstance(result['paths_by_node_count']['forward'][2][0],
                          (dict, defaultdict))
        path_dict = result['paths_by_node_count']['forward'][2][0]
        assert path_dict['path'] == list(test_edge)
        assert isinstance(path_dict['cost'], str)
        assert isinstance(path_dict['sort_key'], str)
        stmts = path_dict['stmts']
        assert isinstance(stmts, list)
        assert isinstance(stmts[0][test_type], list)
        assert isinstance(stmts[0][test_type][0], dict)
        assert stmts[0]['subj'], stmts[0]['obj'] == test_edge
        stmt_dict = stmts[0][test_type][0]
        assert isinstance(stmt_dict['weight'], (np.longfloat, float))
        assert stmt_dict['stmt_type'] == test_row['stmt_type']
        assert stmt_dict['stmt_hash'] == str(test_row['stmt_hash'])
        assert stmt_dict['evidence_count'] == test_row['evidence_count']
        assert isinstance(stmt_dict['source_counts'], dict)
        assert stmt_dict['source_counts'] == test_evidence
        assert stmt_dict['source_counts'][test_source] == \
               test_evidence[test_source]
        assert isinstance(stmt_dict['curated'], bool)
        assert stmt_dict['curated'] is True
        assert stmt_dict['belief'] == test_belief

    def test_query_handling(self):
        result = self.indra_network.handle_query(test=True)
        assert {'paths_by_node_count', 'common_targets', 'common_parents',
                'timeout'} == set(result.keys())
        assert isinstance(result['paths_by_node_count'], dict)
        assert not result['paths_by_node_count'].get('forward', False)
        assert not result['paths_by_node_count'].get('backward', False)
        assert isinstance(result['common_targets'], list)
        assert isinstance(result['common_parents'], dict)
        assert result['timeout'] is False

    def test_dir_edge_structure(self):
        # Get an edge from test DB
        e = None
        for e in self.indra_network.dir_edges:
            if e != test_edge:
                break

        # Check basic edge
        assert isinstance(e, tuple)
        assert len(e) == 2

        # Check edge dict
        edge_dict = self.indra_network.dir_edges[e]
        edge_dict_test = self.indra_network.dir_edges[test_edge]
        assert isinstance(edge_dict, dict)
        assert isinstance(edge_dict_test, dict)
        assert isinstance(edge_dict['belief'], (np.longfloat, float))
        assert isinstance(edge_dict_test['belief'], (np.longfloat, float))
        assert isinstance(edge_dict['weight'], np.longfloat)
        assert isinstance(edge_dict_test['weight'], np.longfloat)

        # Check stmt meta data list
        stmt_list = edge_dict['statements']
        test_stmt_list = edge_dict_test['statements']
        assert isinstance(stmt_list, list)
        assert isinstance(test_stmt_list, list)
        assert isinstance(stmt_list[0], dict)
        assert isinstance(test_stmt_list[0], dict)

        # Check stmt meta data
        assert isinstance(stmt_list[0]['weight'], (float, np.longfloat))
        assert isinstance(test_stmt_list[0]['weight'], (float, np.longfloat))

        assert isinstance(stmt_list[0]['stmt_type'], str)
        assert test_stmt_list[0]['stmt_type'] == 'TestStatement'

        assert isinstance(stmt_list[0]['stmt_hash'], int)
        assert test_stmt_list[0]['stmt_hash'] == 1234567890

        assert isinstance(stmt_list[0]['evidence_count'], int)
        assert test_stmt_list[0]['evidence_count'] == 1

        assert isinstance(stmt_list[0]['source_counts'], dict)
        assert isinstance(test_stmt_list[0]['source_counts'], dict)
        assert len(test_stmt_list[0]['source_counts']) == 1
        assert test_source in test_stmt_list[0]['source_counts']
        assert test_stmt_list[0]['source_counts'][test_source] == 1

        assert isinstance(stmt_list[0]['curated'], bool)
        assert test_stmt_list[0]['curated'] is True

        assert isinstance(stmt_list[0]['belief'], (float, np.longfloat))
        assert isinstance(test_stmt_list[0]['belief'], (float, np.longfloat))
        assert test_stmt_list[0]['belief'] == 0.987654321

    def test_multi_dir_edge_structure(self):
        # Get an edge from test DB
        e = None
        for e in self.indra_network.mdg_edges:
            if e != test_medge:
                break

        # Check basic edge
        assert isinstance(e, tuple)
        assert len(e) == 3

        # Check edge dict
        edge_dict = self.indra_network.mdg_edges[e]
        edge_dict_test = self.indra_network.mdg_edges[test_medge]
        assert isinstance(edge_dict, dict)
        assert isinstance(edge_dict_test, dict)

        assert isinstance(edge_dict['belief'], (np.longfloat, float))
        assert isinstance(edge_dict_test['belief'], (np.longfloat, float))
        assert edge_dict_test['belief'] == 0.987654321

        assert isinstance(edge_dict['weight'], (np.longfloat, float))
        assert isinstance(edge_dict_test['weight'], (np.longfloat, float))

        assert isinstance(edge_dict['stmt_type'], str)
        assert edge_dict_test['stmt_type'] == 'TestStatement'

        assert isinstance(edge_dict['stmt_hash'], int)
        assert edge_dict_test['stmt_hash'] == 1234567890

        assert isinstance(edge_dict['evidence_count'], int)
        assert edge_dict_test['evidence_count'] == 1

        assert isinstance(edge_dict['source_counts'], dict)
        assert isinstance(edge_dict_test['source_counts'], dict)
        assert len(edge_dict_test['source_counts']) == 1
        assert test_source in edge_dict_test['source_counts']
        assert edge_dict_test['source_counts'][test_source] == 1

        assert isinstance(edge_dict['curated'], bool)
        assert edge_dict_test['curated'] is True

    def test_nodes(self):
        # Get a db node
        node = None
        for node in self.indra_network.nodes:
            if node != test_node:
                break

        # Check nodes
        node_dict = self.indra_network.nodes[node]
        test_node_dict = self.indra_network.nodes[test_node]
        assert isinstance(node_dict, dict)
        assert isinstance(test_node_dict, dict)

        assert isinstance(node_dict['ns'], str)
        assert test_node_dict['ns'] == 'TEST'

        assert isinstance(node_dict['id'], str)
        assert test_node_dict['id'] == '1234'
