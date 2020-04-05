=====================
Network Search Web UI
=====================
This documentation introduces the web user interface for the INDRA Network
Search Service.

.. figure:: ../static/images/indra_network_search_screenshot.png
  :align: center
  :figwidth: 100 %

  *The network search interface with no input or results.*


The Graphs Used
---------------
The multiple graphs used for the network search is assembled from a full
snapshot of the `INDRA DataBase <https://github.com/indralab/indra_db>`_ and
is updated regularly. Each statement that includes two or three agents are
assembled into the support for the edges for the graphs with one edge
possibly containing multiple statements. There are three graph types used:

1. DiGraph
2. signed edge DiGraph
3. signed node DiGraph

The **DiGraph** is used for unsigned causal search and for assembling the
statement data supporting the results of the search while the **signed node
graph** is used for signed causal search and the **signed edge graph** is
used for assembling the statement data supporting the signed node search
results.

The edges in the signed edge graphs only contain statements that have clear
up- or downreguations associated with them, which currently is
`IncreaseAmount` and `Activation` for upregulation, and `DecreaseAmount` and
`Inhibition` for downregulation.

The code assembling the graphs can be found `here <https://github
.com/indralab/depmap_analysis/blob/master/depmap_analysis/network_functions
/net_functions.py>`_ in the function `sif_dump_df_to_nx_digraph()`.

Search Options
--------------

Source and Target
~~~~~~~~~~~~~~~~~
The source and target are mandatory fields for the search. The source and
target are the nodes between which to find a path. Source and target does
not have to be of the allowed namespaces (see below). If no result is found
initially, grounding is done on the service backend to try to find an
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
to exclude from the results. If an edge of a path only contain statement
types that are excluded, the whole path will be skipped from the result.

Node Namespaces to *Include*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The namespaces included here are the ones that are allowed on any node
visited in the path search. The namespace of the source and target are
excluded from this restriction. A namespace in INDRA is the type of
identifier used to uniquely identify an entity. For example, a chemical can
be identified using a `CHEBI` identifier and would then be identified in the
`CHEBI` namespace.

Node Name Blacklist
~~~~~~~~~~~~~~~~~~~
Node names entered here are skipped in the path search. This is a good way
to avoid nodes of extremely high degree that overwhelmes the results and
effectively blocks out results including lower degree nodes. *See also Cull
Highest Degree Node below.*

Edge Hash Blacklist
~~~~~~~~~~~~~~~~~~~
To ignore a specific statement supporting an edge, the statement hash for
that statement can be added here.

Belief Score cut-off
~~~~~~~~~~~~~~~~~~~~
This option enables a belief score cut-off so that statements supporting an
egde has to have a belief score above this threshold. It is set to zero by
default. Read more about belief scores in the `belief module
<https://indra.readthedocs.io/en/latest/modules/belief/index.html>`_ of
INDRA.

Cull Highest Degree Node
~~~~~~~~~~~~~~~~~~~~~~~~
Entering an integer N here allows the path search to include the highest
degree node for the first N returned paths, after which it is added to the
**Node Name Blacklist**. This is repeated for the second highest degree node
for the following N paths, then for the third highest degree node and so
forth. *Note:* This option is currently only applied for unsigned path
searches.

Signed Search
~~~~~~~~~~~~~
To perform a signed search, click on the drop down menu that says "No sign"
and chose a sign. "+" means that all the returned paths are upregulations,
and "-" means that all the returned paths are downregulations. For the
purpose of signed search, only statements that imply a clear up- or
downregulation are considered. Currently this mean `IncreaseAmount` and
`Activation` for upregulation, and `DecreaseAmount` and `Inhibition` for
downregulation.

Include Shared Regulators
~~~~~~~~~~~~~~~~~~~~~~~~~
This checkbox adds results from a search of direct common shared regulators
of source and target. A direct shared regulator is defined as any node that
is exactly one edge upstream of *both* source and target.

Include Reverse Search
~~~~~~~~~~~~~~~~~~~~~~
With this option, the reverse search *from* target *to* source is done as
well as the original search from source to target. If the timeout is reached
(see below) before the reverse seach can start, the reverse search will
not return any paths. If the timeout is reached during the reverse search,
fewer paths than for the original search will be retured.

Weighted Search
~~~~~~~~~~~~~~~
When performing a weighted search, the cost along every path encountered is
calculated as the sum of the weights along the path. The paths are then
returned in ascending order of cost. The cost of a path is defined as the
sum of the weights of all the edges along the paths. The weigthed search
uses a slightly modified version of the Djikstra weighted search empolyed in
Networkx. *Note:* A weighted search is costly and usually takes longer than
a normal search. It is common that a very heavy weigthed search times out,
especially for a *signed* weighted search.

The code implemented for the weighted search is available `here <../.
./master/depmap_analysis/network_functions/net_functions.py>`_ in the
function `shortest_simple_paths()`.

Databases Only
~~~~~~~~~~~~~~
With this option, only statements that contain sources from curated
databases like PathwayCommons and Signor are allowed to support edges in the
returned paths.

Include Famplex Families and Complexes in Path Search
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This option allows for edges to be between a gene and its family or
beteween a gene and a complex formed by its encoded protein. For example: an
edge between `BRCA1` and its family `BRCA` would be allowed.

Expand search to FamPlex
~~~~~~~~~~~~~~~~~~~~~~~~
If a path search returns empty, this option will allow the path search to be
retried with parents if the source and/or target entities. For example, if a
search with `BRCA1` as source returns empty, the search would be retried
with the `BRCA` family as source instead.

Timeout
~~~~~~~
Setting a timeout allows to set a larger (or smaller) timeout than the
default 30 seconds timeout. The time since the path search was started is
checked after each path has been checked during the search. If the time
passed is larger than the allowed timeout, the search is interrupted and
returns as fast as possible. The timeout provided has to be a decimal number
smaller than or equeal to 120 seconds.

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
