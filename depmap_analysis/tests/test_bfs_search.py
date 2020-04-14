import networkx as nx
from depmap_analysis.network_functions.net_functions import bfs_search


def test_bfs_search():
    dg = nx.DiGraph()
    edges = [('Z1', 'A1'), ('A1', 'B1'), ('A2', 'B1'), ('A3', 'B2'),
             ('A4', 'B2'), ('B1', 'C1'), ('B2', 'C1'), ('B3', 'C1'),
             ('C1', 'D1')]
    dg.add_edges_from(edges)

    nodes1, nodes2 = list(zip(*edges))
    nodes = set(nodes1).union(nodes2)
    for node in nodes:
        ns = node[0]
        _id = node[1]
        dg.nodes[node]['ns'] = ns
        dg.nodes[node]['id'] = _id

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
                                   node_filter=['a', 'b', 'c', 'd', 'z'])]
    assert len(paths) == 4
    assert set(paths) == expected_paths

    # Test terminal NS
    # Terminate on 'b'
    expected_paths = {('D1', 'C1'), ('D1', 'C1', 'B1'), ('D1', 'C1', 'B2'),
                      ('D1', 'C1', 'B3')}
    paths = [p for p in bfs_search(g=dg, source='D1', depth_limit=5,
                                   reverse=True, terminal_ns=['b'],
                                   node_filter=['a', 'b', 'c', 'd', 'z'])]
    assert len(paths) == 4
    assert set(paths) == expected_paths
    # Terminate on 'a'
    expected_paths = {('D1', 'C1'), ('D1', 'C1', 'B1'), ('D1', 'C1', 'B2'),
                      ('D1', 'C1', 'B3'), ('D1', 'C1', 'B1', 'A1'),
                      ('D1', 'C1', 'B1', 'A2'), ('D1', 'C1', 'B2', 'A3'),
                      ('D1', 'C1', 'B2', 'A4')}
    paths = [p for p in bfs_search(g=dg, source='D1', depth_limit=5,
                                   reverse=True, terminal_ns=['a'],
                                   node_filter=['a', 'b', 'c', 'd', 'z'])]
    assert len(paths) == len(expected_paths)
    assert set(paths) == expected_paths

