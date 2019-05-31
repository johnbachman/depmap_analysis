INDRA DepMap Service
--------------------

Currently, while the being under development, the service is only available
to be run locally.

Setup
=====

To run the service locally, three things are needed:

1. Fetch the latest update to the branch
   `'web-overhaul' <https://github.com/kkaris/depmap_analysis/tree/web-overhaul>`_
   of the depmap_analysis repository from one of the maintainers.
2. Download the `latest network representations <https://s3.amazonaws.com/depmap-public/_cache/>`_
   of the indra network (might require AWS S3 login). The files needed are:

   * ``nx_bs_fam_multi_digraph_db_dump_20190417.pkl``
   * ``nx_bs_fam_dir_graph_db_dump_20190417.pkl``
   * (optional) ``test_mdg_network.pkl``
   * (optional) ``test_dir_network.pkl``
3. Make sure the constants ``INDRA_DG_CACHE`` and ``INDRA_MDG_CACHE`` in
   ``api.py`` corresond to the files downloaded in the previous step.

Running the Service
===================

Running the service locally requires Python 3.5+. Dependencies are the same
as for INDRA and INDRA_DB::

  python api.py [-h] [--host HOST] [--port PORT]

where ``HOST`` is the address to use (default is ``127.0.0.1``) and ``PORT``
is the port to use (default is ``5000``). If default settings are used, a
web ui is hosted on http://127.0.0.1:5000/query and query submissions are
done to http://127.0.0.1:5000/query/submit.

Searching
=========

To search, enter a source and a target in the appropriate text boxes and
apply the desired settings.
