INDRA Network Search
--------------------

Currently, during development, the service is only available to be run locally.

Setup
=====

To run the service locally, two things are needed:

1. Fetch the latest update to the branch
   `'master' <https://github.com/indralab/depmap_analysis/tree/master>`_
   of the depmap_analysis repository from one of the maintainers.
2. Download the latest network representations of the indra network
   (might require AWS S3 login). The files needed are:

   * ``nx_bs_fam_dir_graph_db_refresh_20190702.pkl``
   * ``nx_bs_fam_multi_digraph_db_refresh_20190702.pkl``

Running the Service Locally
===========================

Running the service locally requires Python 3.5+. Dependencies are the same
as for INDRA and INDRA_DB::

  python api.py [-h] [--host HOST] [--port PORT] [--cache DG_GRAPH MDG_GRAPH]

where ``HOST`` is the address to use (default is ``127.0.0.1``), ``PORT``
is the port to use (default is ``5000``) and ``DG_GRAPH`` and ``MDG_GRAPH`` are pickled NetworkX files representing
the INDRA knowledge network in DiGraph and MultiDiGraph representations, respectively. The ``--cache``
flag overrides the defaults in the file so that any file can be provided. If default settings are used for ``HOST``
and ``PORT``, a web ui is hosted on http://localhost:5000/query and query submissions are
done to http://localhost:5000/query/submit.

Searching
=========

To search, enter a source and a target in the appropriate text boxes and
apply the desired settings. You can read more about the specific settings in the docstring of the
``IndraNetwork.handle_query`` method in ``depmap_analysis/indra_depmap_service/api.py``
