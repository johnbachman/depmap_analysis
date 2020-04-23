import networkx as nx
from indra.explanation.model_checker import signed_edges_to_signed_nodes
from depmap_analysis.network_functions.net_functions import bfs_search, \
    INT_PLUS, INT_MINUS

# Set up
dg = nx.DiGraph()
all_ns = ['a', 'b', 'c', 'd', 'z']
edges = [('Z1', 'A1'), ('A1', 'B1'), ('A2', 'B1'), ('A3', 'B2'),
         ('A4', 'B2'), ('B1', 'C1'), ('B2', 'C1'), ('B3', 'C1'),
         ('C1', 'D1')]

# Ensures alphabetical order
edge_beliefs = {('Z1', 'A1'): 1-0.2,
                ('A1', 'B1'): 1-0.2,
                ('A2', 'B1'): 1-0.3,
                ('A3', 'B2'): 1-0.5,
                ('A4', 'B2'): 1-0.6,
                ('B1', 'C1'): 1-0.2,
                ('B2', 'C1'): 1-0.3,
                ('B3', 'C1'): 1-0.4,
                ('C1', 'D1'): 1-0.2}
dg.add_edges_from(edges)

# Add belief
for e in dg.edges:
    dg.edges[e]['belief'] = edge_beliefs[e]

# Add namespaces
nodes1, nodes2 = list(zip(*edges))
nodes = set(nodes1).union(nodes2)
for node in nodes:
    ns = node[0]
    _id = node[1]
    dg.nodes[node]['ns'] = ns
    dg.nodes[node]['id'] = _id


def test_bfs():
    # Test basic part of algorithm
    assert len([p for p in
                bfs_search(dg, 'C1', depth_limit=1, reverse=True)]) == 3
    assert len([p for p in
                bfs_search(dg, 'C1', depth_limit=2, reverse=True)]) == 7
    assert len([p for p in
                bfs_search(dg, 'C1', depth_limit=2, reverse=True,
                           path_limit=4)]) == 4

    # Test ns allowance list
    ans = ['c', 'b']
    assert len([p for p in bfs_search(dg, 'C1', depth_limit=2, reverse=True,
                                      node_filter=ans)]) == 3
    assert all(len(p) < 3 for p in
               bfs_search(dg, 'C1', depth_limit=2, reverse=True,
                          node_filter=ans))

    # Test longer paths
    assert len([p for p in bfs_search(dg, 'D1', depth_limit=5,
                                      reverse=True)]) == 9

    # Test node blacklist
    assert len([p for p in bfs_search(dg, 'D1', depth_limit=5, reverse=True,
                                      node_blacklist={'Z1'})]) == 8

    # Test max per node option
    # Should get 4 paths with max_per_node=1
    expected_paths = {('D1', 'C1'), ('D1', 'C1', 'B1'),
                      ('D1', 'C1', 'B1', 'A1'),
                      ('D1', 'C1', 'B1', 'A1', 'Z1')}
    paths = [p for p in bfs_search(g=dg, source='D1', depth_limit=5,
                                   reverse=True, max_per_node=1,
                                   node_filter=all_ns)]
    assert len(paths) == 6
    assert set(paths) == expected_paths

    # Test terminal NS
    # Terminate on 'b'
    expected_paths = {('D1', 'C1'), ('D1', 'C1', 'B1'), ('D1', 'C1', 'B2'),
                      ('D1', 'C1', 'B3')}
    paths = [p for p in bfs_search(g=dg, source='D1', depth_limit=5,
                                   reverse=True, terminal_ns=['b'],
                                   node_filter=all_ns)]
    assert len(paths) == 4
    assert set(paths) == expected_paths
    # Terminate on 'a'
    expected_paths = {('D1', 'C1'), ('D1', 'C1', 'B1'), ('D1', 'C1', 'B2'),
                      ('D1', 'C1', 'B3'), ('D1', 'C1', 'B1', 'A1'),
                      ('D1', 'C1', 'B1', 'A2'), ('D1', 'C1', 'B2', 'A3'),
                      ('D1', 'C1', 'B2', 'A4')}
    paths = [p for p in bfs_search(g=dg, source='D1', depth_limit=5,
                                   reverse=True, terminal_ns=['a'],
                                   node_filter=all_ns)]
    assert len(paths) == len(expected_paths)
    assert set(paths) == expected_paths


def test_signed_bfs():
    seg = nx.MultiDiGraph()
    signed_edges = [
        ('Z1', 'A1', INT_PLUS),  # 1
        ('Z1', 'A1', INT_MINUS),  # 2
        ('A1', 'B1', INT_PLUS),  # 3
        ('A2', 'B1', INT_MINUS),  # 4
        ('B1', 'C1', INT_PLUS),  # 5
        ('A3', 'B2', INT_PLUS),  # 6
        ('A4', 'B2', INT_MINUS),  # 7
        ('B2', 'C1', INT_PLUS),  # 8
        ('B2', 'C1', INT_MINUS),  # 9
        ('B3', 'C1', INT_MINUS),  # 10
        ('C1', 'D1', INT_PLUS),  # 11
        ('C1', 'D1', INT_MINUS),  # 12
    ]

    seg.add_edges_from(signed_edges)
    # ATTN!! seg.edges yields u, v, index while seg.edges() yields u, v
    for u, v, sign in seg.edges:
        seg.edges[(u,v,sign)]['sign'] = sign
        seg.edges[(u,v,sign)]['belief'] = random()

    sng = signed_edges_to_signed_nodes(graph=seg, prune_nodes=True,
                                       copy_edge_data=False)

    # D1 being upregulated: 12 paths
    paths = [p for p in bfs_search(
        g=sng, source=('D1', INT_PLUS), g_nodes=dg.nodes, reverse=True,
        depth_limit=5, node_filter=all_ns, sign=INT_PLUS)
    ]
    assert len(paths) == 12, len(paths)
