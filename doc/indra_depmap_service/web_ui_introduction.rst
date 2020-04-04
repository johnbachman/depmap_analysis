=====================
Network Search Web UI
=====================
This documentation introduces the web user interface for the INDRA Network
Search.

.. figure:: ../static/images/indra_network_search_screenshot.png
  :align: center
  :figwidth: 100 %

  *The network search interface with no input or results.*

Search Options
--------------

Source and Target
~~~~~~~~~~~~~~~~~
The source and target are mandatory fields for the search. The source and
target are the nodes between which to find a path. Source and target does
not have to be of the allowed name spaces (see below). If no result is found
initially, greounding is done on the service backend to try to find an
alternative name for the provided node name.

Path Length
~~~~~~~~~~~
The path length to search for. Should be a positive integer. For the purpose
of this search interface, the path length is defined here as the number of
edges between the source and the target.

Max # Paths
~~~~~~~~~~~
The maximum number of results to return per category in the results. The
default and the maximum allowed is 50 results.

Statement Types to *Exclude*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is a multiselect dropdown which contains multiple statement type names
to exclude from the results. If an edge of a path only contains statement
types that are excluded, the whole path will be skipped from the result.

Node Name Spaces to *Exclude*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
abcd

Node Name Blacklist
~~~~~~~~~~~~~~~~~~~
abcd

Edge Hash Blacklist
~~~~~~~~~~~~~~~~~~~
abcd

Belief Score cut-off
~~~~~~~~~~~~~~~~~~~~
abcd

Cull Highest Degree Node
~~~~~~~~~~~~~~~~~~~~~~~~
abcd

Signed Search
~~~~~~~~~~~~~
abcd

Include Shared Regulators
~~~~~~~~~~~~~~~~~~~~~~~~~
abcd

Include Reverse Search
~~~~~~~~~~~~~~~~~~~~~~
abcd

Weighted Search
~~~~~~~~~~~~~~~
Uses a modified Djikstra's weighted search algorithm

Databases Only
~~~~~~~~~~~~~
abcd

Include Famplex Families and Complexes in Path Search
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
abcd

Expand search to FamPlex
~~~~~~~~~~~~~~~~~~~~~~~~
abcd

Timeout
~~~~~~~
abcd

Result Categories
-----------------
If there are not results for the specific section, that section card won't
show up.

Complexes and Families
~~~~~~~~~~~~~~~~~~~~~~

Common Targets
~~~~~~~~~~~~~~

Shared Regulators
~~~~~~~~~~~~~~~~~

N Edge Paths
~~~~~~~~~~~~


Download Results
----------------
You can download the search result json and the statements from the path
search (not the other searches) in a json format.
