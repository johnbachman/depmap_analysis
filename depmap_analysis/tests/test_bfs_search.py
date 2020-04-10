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
