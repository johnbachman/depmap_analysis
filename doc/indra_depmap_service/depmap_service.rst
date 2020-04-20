INDRA Network Search
--------------------

The INDRA Network Search service is available at
https://network.indra.bio/

Available Service Endpoints
===========================

There are several endpoints served by the service, shown in the table below.
Each endpoint is described further in the sections below.

+----------------------+------------------+------------------------+
| Endpoint             | Available method | Purpose                |
+======================+==================+========================+
| `/` and              | GET              | Main page, path search |
| `/query`             |                  | with options:          |
|                      |                  |                        |
|                      |                  | - weighted search      |
|                      |                  | - signed search        |
|                      |                  | - restrict visited     |
|                      |                  |   nodes                |
+----------------------+------------------+------------------------+
| `/query/submit`      | POST             | Handles path queries   |
+----------------------+------------------+------------------------+
| `/multi_interactors` | POST             | Handles queries for    |
|                      |                  | immediate common up-   |
|                      |                  | or downstream          |
|                      |                  | interactors from a list|
|                      |                  | of targets or          |
|                      |                  | regulators             |
+----------------------+------------------+------------------------+
| `/node` and `/nodes` | POST             | Check if a node or a   |
|                      |                  | list of nodes is in the|
|                      |                  | network                |
+----------------------+------------------+------------------------+
| `/bfs_search`        | POST             | Handles breadth first  |
|                      |                  | queries from a source  |
|                      |                  | or target node.        |
+----------------------+------------------+------------------------+

Node Existence Lookup
.....................

The `/node` and `/nodes` endpoints serve the prupose of quickly looking up
if a particular node exists in the network. To use the `/node` endpoint you
have to POST a json to https://network.indra.bio/node::

    Method: POST
    JSON: {'node': '<name of node>'}

The return will be a json of the following format::

    {'node': '<name of node>', 'in_network': True/False}

To check a list of multiple nodes at the same time, the `/nodes` enpoint is
a better choice. This endpoint expects a list of node names in the incoming
json::

    Method: POST
    JSON: {'nodes': ['<name1>', '<name2>', ...]}

The resulting json will be a node name - boolean dictionary for each passed
node name in the initial json::

    {'<name1>': True/False, '<name2>': True/False, ...}


Direct Multi Interactors Lookup
...............................

the `/multi_interactors` endpoint allows for search of common upstream or
downstream interactors to the provided list of targets/regulators. If
regulators/targets are provided, the search will be for common
targets/regulators. The only required option is to provided a list of
regulators *or* targets. There are several options that are not required and
will defulat to different values. The table below describes the available
options.

+----------------+---------+----------+----------------+-------------------+
| Option         | Default | Required | Allowed values | Action            |
+================+=========+==========+================+===================+
|  regulators OR |         | Yes      | List of node   | Input for         |
|  targets       |         |          | names          | search            |
+----------------+---------+----------+----------------+-------------------+
|  allowed_ns    | any ns  | No       | *see below*    | Skip nodes with   |
|                |         |          |                | namespaces not in |
|                |         |          |                | the provided list |
+----------------+---------+----------+----------------+-------------------+
|  belief_cutoff |    0    | No       | 0 <= n <= 1    | Skip statements   |
|                |         |          |                | below threshold   |
+----------------+---------+----------+----------------+-------------------+
| skip_stmt_types|   [ ]   | No       | *see below*    | Skip statements   |
|                |         |          |                | of provided       |
|                |         |          |                | types             |
+----------------+---------+----------+----------------+-------------------+
| db_only        |  False  | No       | True/False     | Only allow        |
|                |         |          |                | statements from   |
|                |         |          |                | database sources  |
+----------------+---------+----------+----------------+-------------------+

The POST request to the endpoint should look like this::

    Method: POST
    JSON: {'regulators' OR 'targets': ['<name1>', '<name2>', ...],
           'allowed_ns': ['<name space1>', '<name space>', ...],
           'belief_cutoff': 0 <= float <=1,
           'skip_stmt_types': ['<statement type>', ...],
           'db_only': bool}


For `allowed_ns`, the following values are allowed:

- HGNC
- FPLX
- CHEBI
- PUBCHEM
- MIRBASE
- GO
- MESH
- HMDB

For `skip_stmt_types`, valid statement types need to be provided. To see a
full list of statement types, see
https://indra.readthedocs.io/en/latest/modules/statements.html#module-indra.statements.statements

By POSTing `"{'help': ''}"` to the endpoint, a small json is returned that
describes the options available.

Breadth First Search Endpoint
.............................
For doing breadth first searches the `/bfs_search` endpoint can be used. The
only required option is to provide a node name to start the search with.
Other options are the following:

+----------------+-----------+---------+----------------+-------------------+
| Option         | Default   |Required | Allowed values | Action            |
+================+===========+=========+================+===================+
|  source        |           |Yes      | Any node name  | Input for search  |
+----------------+-----------+---------+----------------+-------------------+
|  reverse       | False     |No       | True/False     | Reverse the search|
|                |           |         |                | to search upstream|
|                |           |         |                | instead of        |
|                |           |         |                | downstream        |
+----------------+-----------+---------+----------------+-------------------+
|  allowed_ns    | any ns    |No       | *see below*    | Skip nodes with   |
|                |           |         |                | namespaces not in |
|                |           |         |                | the provided list |
+----------------+-----------+---------+----------------+-------------------+
|  belief_cutoff |    0      |No       | 0 <= float <= 1| Skip statements   |
|                |           |         |                | below threshold   |
+----------------+-----------+---------+----------------+-------------------+
| skip_stmt_types|   [ ]     |No       | *see below*    | Skip statements   |
|                |           |         |                | of provided       |
|                |           |         |                | types             |
+----------------+-----------+---------+----------------+-------------------+
| depth_limit    |  2        |No       | 1 < int        | Only allow        |
|                |           |         |                | paths up to the   |
|                |           |         |                | provided length   |
+----------------+-----------+---------+----------------+-------------------+
| path_limit     | 100       |No       | 0 < int        | The maximum number|
|                |           |         |                | of paths to yield |
|                |           |         |                | from the generator|
+----------------+-----------+---------+----------------+-------------------+
| max_results    |  50       |No       | int            | The maximum number|
|                |           |         |                | of results to     |
|                |           |         |                | return            |
+----------------+-----------+---------+----------------+-------------------+
| node_blacklist |  []       |No       | str            | Node names to skip|
|                |           |         |                | in the search     |
+----------------+-----------+---------+----------------+-------------------+
| terminal_ns    | ['CHEBI', |No       | *see below*    | Node namespaces to|
|                | 'PUBCHEM']|         |                | end paths on      |
+----------------+-----------+---------+----------------+-------------------+
| max_per_node   |  5        |No       | True/False     | Maximum number of |
|                |           |         |                | paths from the    |
|                |           |         |                | same leaf parent  |
|                |           |         |                | node to yield     |
+----------------+-----------+---------+----------------+-------------------+
| sign           |  None     |No       | '+'/'-'        | Perform a signed  |
|                |           |         |                | search.           |
|                |           |         |                | *See below*       |
+----------------+-----------+---------+----------------+-------------------+


`'allowed_ns'` specifies which namespace a node in the resulting paths are
allowed to have. The following values are allowed:

- HGNC
- FPLX
- CHEBI
- PUBCHEM
- MIRBASE
- GO
- MESH
- HMDB

`'terminal_ns'` indicates a namespace that when encountered will not yield
further paths from it. For example, if one is only interested in paths with
a terminal node that is a chemical, 'PUBCHEM' and 'CHEBI' would be specified
here for example.

Signed Breadth First Search
^^^^^^^^^^^^^^^^^^^^^^^^^^^
By providing a sign with the keyword `'sign'`, the search is done over a
signed directed graph instead of a standard directed graph. The sign
determines if the result of the path is an up- or downregulation. To search
e.g. upregulations of the gene `'ACE2'`, the following JSON would have to be
sent::

    {'source': 'ACE2',
     'reverse': True,  # If reverse is True, search is upstream from 'source'
     'sign': '+'}

To search for any gene or gene family or protein complex that is upstream
of `'ACE2'`, the following JSON should be used::

    {'source': 'ACE2',
     'reverse': True,
     'allowed_ns': ['FPLX', 'HNGC'],
     'terminal_ns': []}  # 'terminal_ns' defaults to ['CHEBI', 'PUBCHEM']


Path Query Endpoint
...................
*Note:* This section describes the behavior of the POST endpoint
`/query/submit` handling path queries. For documentation of the web UI
endpoint `/query`, go to the web ui `documentation <./web_ui_introduction
.html>`_. The following list describes the available options for the POST
endpoint:

- **source** (str): the source node for the path.
- **target** (str): the target for the path.
- **stmt_filter** ([str]): a list of valid indra statement types or FamPlex
  child-parent connections (as 'fplx') *to exclude* in the path.
- **node_filter** ([str]): A list of node namespaces *to include* in the path
- **node_blacklist** ([str]): A list of node names to ignore. If a path
  contains a node in this list, the path will be discarded.
- **edge_hash_blacklist** ([str/int]): A list of statement hashes (as
  strings or ints) to ignore. If an edge statement hash is found in this
  list, it will be discarded from the assembled edge list.
- **cull_best_node** (int): A positive integer. Every x valid paths, cull the
  node with the highest (weighted): degree from the network. This increases
  the variety of paths found and reduces the impact of nodes with extremely
  high connectivity in the network.
- **path_length** (int|False): a positive integer stating the number of edges
  that should be in the returned path. If False, return paths with any number
  of edges.
- **sign** (None|str): If 'no_sign' or None, do regular unsigned graph search
  over the directed graph. If 'plus'/'minus', only paths with overall up/down
  regulation will be returned.
- **weighted** (Bool): If True, do a weighted path search. Weights in the
  network are assigned as -log(belief score).
- **bsco** (0 <= float <= 1.0): Belief Score Cut-Off, a positive decimal
  number <1.0 indicating at what belief score an edge statement should be
  ignored.
- **curated_db_only** (Bool): Filter results to only allow edges that are
  sourced from curated databases.
- **fplx_expand** (Bool): If True, when no path is found in the initial search,
  look for paths between the parents of the source and target.
- **k_shortest** (Bool|int): An integer stating the maximum number of directed
  paths to return in the result. The maximum allowed value is 50. If False,
  the maximum number of paths returned will be set to the maximum allowed
  value.
- **user_timeout** (float): A decimal specifying the number of seconds to use
  for timeout. If not provided, the default of 30 seconds is used.
- **two_way** (Bool): If True, search path both ways, i.e. search A->B and
  B->A.

The json would look like this::

    Method: POST
    {'source': <str>,
     'target': <str>,
     'stmt_filter': [str],
     'node_blacklist': [str],
     'edge_hash_blacklist': [str/int],
     'cull_best_node': [int],
     'path_length': <int|False>,
     'sign': <None|str>,
     'weighted': <Bool>,
     'bsco': <float>,
     'curated_db_only': <Bool>,
     'fplx_expand': <Bool>,
     'k_shortest': Bool|<int>,
     'user_timeout' : <float>,
     'two_way': <Bool>}


You can also read about the specific settings in the docstring of the
``IndraNetwork.handle_query`` method in
``depmap_analysis/indra_depmap_service/api.py``.


Running the Service Locally
===========================

To run the service locally, two things are needed:

1. Fetch the latest update to the branch
   `'master' <https://github.com/indralab/depmap_analysis/tree/master>`_
   of the depmap_analysis repository from one of the maintainers.
2. Download the latest network representations of the indra network
   (might require AWS S3 login):

   * ``indranet_dir_graph_latest.pkl``
   * ``indranet_sign_edge_graph_latest.pkl`` (optional)
   * ``indranet_sign_node_graph_latest.pkl`` (optional)

   The signed representations of the graph are only needed for signed path
   search.

Dependecies are Python 3.6+, but otherwise the same as for INDRA and
INDRA_DB. In the depmap_analysis.indra_depmap_service, run api.py with the
followung arguments::

  python -m api.py [-h] [--host HOST] [--port PORT] [--cache DG_GRAPH
  MDG_GRAPH|None SIGN_EDGE_GRAPH|None SIGN_NODE_GRAPH|None]

where ``HOST`` is the address to use (default is ``127.0.0.1``), ``PORT``
is the port to use (default is ``5000``) and ``DG_GRAPH``, ``MDG_GRAPH``,
``SIGN_EDGE_GRAPH`` and ``SING_NODE_GRAPH`` are pickled graphs representing
the INDRA knowledge network in DiGraph, MultiDiGraph and SignedGraph
representations, respectively. The ``--cache`` flag overrides the defaults
in the file so that any file can be provided. If default settings are used
for ``HOST`` and ``PORT``, a web ui is hosted on http://localhost:5000/query
and query submissions are done to http://localhost:5000/query/submit.
