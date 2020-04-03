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

The `/node` and `/nodes` endpoint


Running the Service Locally
===========================

To run the service locally, two things are needed:

1. Fetch the latest update to the branch
   `'master' <https://github.com/indralab/depmap_analysis/tree/master>`_
   of the depmap_analysis repository from one of the maintainers.
2. Download the latest network representations of the indra network
   (might require AWS S3 login). The files needed are:

   * ``indranet_signed_nodes_graph_latest.pkl``
   * ``indranet_multi_digraph_latest.pkl``
   * ``indranet_dir_graph_latest.pkl``

Running the Service Locally
===========================

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
