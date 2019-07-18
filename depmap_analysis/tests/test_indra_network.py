import unittest
import numpy as np
from random import random as rnd

import indra_db.tests.util as tu
from indra_db.util.dump_sif import load_db_content, make_ev_strata, \
    make_dataframe, NS_LIST

from depmap_analysis.network_functions.network_functions import \
    sif_dump_df_to_nx_digraph
from depmap_analysis.network_functions.indra_network import IndraNetwork

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
for n, h in df['hash'].iteritems():
    bsd[h] = rnd()

# Add custom row to df that can be checked later
test_edge = ('GENE_A', 'GENE_B')
test_node = test_edge[0]
df = df.append({
        'agA_ns': 'TEST', 'agA_id': '1234', 'agA_name': test_edge[0],
        'agB_ns': 'TEST', 'agB_id': '2345', 'agB_name': test_edge[1],
        'stmt_type': 'TestStatement', 'evidence_count': 1, 'hash': 1234567890
    },
    ignore_index=True)
sed[1234567890] = {'tester': 1}
bsd[1234567890] = 0.987654321


class TestNetwork(unittest.TestCase):
    def setUp(self):
        self.df = df
        self.indra_network = IndraNetwork(
            indra_dir_graph=sif_dump_df_to_nx_digraph(
                df=self.df, belief_dict=bsd, strat_ev_dict=sed, multi=False,
                include_entity_hierarchies=True),
            indra_multi_dir_graph=sif_dump_df_to_nx_digraph(
                df=self.df, belief_dict=bsd, strat_ev_dict=sed, multi=True,
                include_entity_hierarchies=True)
        )

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
        assert isinstance(edge_dict['bs'], (np.longfloat, np.float))
        assert isinstance(edge_dict_test['bs'], (np.longfloat, np.float))
        assert isinstance(edge_dict['weight'], np.longfloat)
        assert isinstance(edge_dict_test['weight'], np.longfloat)

        # Check stmt meta data list
        stmt_list = edge_dict['stmt_list']
        test_stmt_list = edge_dict_test['stmt_list']
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

        assert isinstance(stmt_list[0]['evidence'], dict)
        assert isinstance(test_stmt_list[0]['evidence'], dict)
        assert len(test_stmt_list[0]['evidence']) == 1
        assert 'tester' in test_stmt_list[0]['evidence']
        assert test_stmt_list[0]['evidence']['tester'] == 1

        assert isinstance(stmt_list[0]['curated'], bool)
        assert test_stmt_list[0]['curated'] is True

        assert isinstance(stmt_list[0]['bs'], (float, np.longfloat))
        assert isinstance(test_stmt_list[0]['bs'], (float, np.longfloat))
        assert test_stmt_list[0]['bs'] == 0.987654321

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
