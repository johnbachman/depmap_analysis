from networkx import DiGraph
from depmap_analysis.scripts.depmap_script_expl_funcs import *

mek_erk_edge = {'statements': [{'stmt_hash': 32842005799216216,
                                'stmt_type': 'Complex',
                                'evidence_count': 1,
                                'belief': 0.65,
                                'source_counts': {'sparser': 1},
                                'residue': None,
                                'weight': 0.4307829160924542,
                                'curated': False,
                                'position': None,
                                'english': 'MEK binds ERK.'}],
                'belief': 0.65,
                'weight': 0.4307829160924542}
map2k1_erk_edge = {'statements': [{'stmt_hash': -1337529928737095,
                                   'stmt_type': 'Phosphorylation',
                                   'evidence_count': 7,
                                   'belief': 0.9994999999961258,
                                   'source_counts': {'reach': 7},
                                   'residue': None,
                                   'weight': 0.0005001250455584151,
                                   'curated': False,
                                   'position': None,
                                   'english': 'MAP2K1 phosphorylates ERK.'}],
                   'belief': 0.9994999999961258,
                   'weight': 0.0005001250455584151}
mapk1_mek_edge = {'statements': [{'stmt_hash': 1194824756615807,
                                  'stmt_type': 'Complex',
                                  'evidence_count': 1,
                                  'belief': 0.65,
                                  'source_counts': {'sparser': 1},
                                  'residue': None,
                                  'weight': 0.4307829160924542,
                                  'curated': False,
                                  'position': None,
                                  'english': 'MAPK1 binds MEK.'}],
                  'belief': 0.65,
                  'weight': 0.4307829160924542}


def test_parent_connections():
    dg = DiGraph()

    dg.add_node('MAP2K1', ns='HGNC', id='6840')
    dg.add_node('MEK', ns='FPLX', id='MEK')
    dg.add_node('MAPK1', ns='HGNC', id='6871')
    dg.add_node('ERK', ns='FPLX', id='ERK')

    # Add ns-id mapping
    dg.graph['node_by_ns_id'] = {('HGNC', '6840'): 'MAP2K1',
                                 ('HGNC', '6871'): 'MAPK1',
                                 ('FPLX', 'MEK'): 'MEK',
                                 ('FPLX', 'ERK'): 'ERK'}

    # Add mek-erk
    dg.add_edge('MEK', 'ERK', **mek_erk_edge)
    # Add MAP2K1-ERK
    dg.add_edge('MAP2K1', 'ERK', **map2k1_erk_edge)
    # Add MAPK1-MEK
    dg.add_edge('MAPK1', 'MEK', **mapk1_mek_edge)

    par_conn_func = expl_functions['parent_connections']

    a, b, explained, result = par_conn_func(s='MAP2K1', o='MAPK1', net=dg,
                                            corr=3.1415, _type='unsigned')

    # Check that results were found
    assert explained
    assert result is not None

    # Check keys
    assert set(result.keys()) == {'MAP2K1_op', 'MAPK1_sp', 'sp_op'}

    # Check mek-erk
    assert result['sp_op']['MEK_ERK'] == mek_erk_edge
    # Check map2k1-erk
    assert result['MAP2K1_op']['MAP2K1_ERK'] == map2k1_erk_edge
    # Check mapk1-mek
    assert result['MAPK1_sp']['MAPK1_MEK'] == mapk1_mek_edge
