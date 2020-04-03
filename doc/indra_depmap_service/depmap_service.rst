INDRA Network Search
--------------------

The INDRA Network Search service is available at
https://network.indra.bio/

Available Service Endpoints
===========================

There are several endpoints served by the service

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

To perform a path search, enter a source and a target in the appropriate text
boxes and apply the desired settings. You can read more about the specific
settings in the docstring of the ``IndraNetwork.handle_query`` method in
``depmap_analysis/indra_depmap_service/api.py``.

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
|                |         |          |                | name spaces not in|
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
