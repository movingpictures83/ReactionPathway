# ReactionPathway
# Language: Python
# Input: Prefix (for multiple input files)
# Output: TXT (statistics)
# Tested with: PluMA 1.1, Python 3.6

PluMA plugin to perform statistical analysis of a network based on information regarding
"pathways" for the nodes.  This plugin is likely most useful for metabolomics data, or gene
expression data where corresponding pathway data is available.

The plugin accepts a prefix as its "inputfile" argument, and will then assume the following
four files:
1. prefix.csv: Network of nodes in CSV format, where rows and columns both represent nodes and
entry (i, j) corresponds to the edge weight from node i to node j.  A correlation network
works well here; since the statistics account for positive, zero and negative edge weights.
2. prefix.clusters.noa: A NOde Attribute (NOA) file that contains a simple table which
maps node names (these should be the same as in #1) to unique cluster identifiers.  Two
nodes that map to the same cluster ID are assumed to be in the same cluster.
3. prefix.pathways.noa: An NOA file with a simple table that maps node names (same as #1
and #2) to unique pathway identifiers.  Two nodes that map to the same pathway ID are
assumed to be on the same biological pathway (there can be various kinds of pathways,
depending on the input data type).
4. prefix.pathways.txt: A table that maps pathway identifiers (same as #3) to textual
descriptions of the pathway.

The plugin sends statistics to the screen, including:
* Percentage of metabolite pairs on the same pathway that also appear in the same cluster.
* Percentage of metabolite pairs on the same pathway that have positive edge weights.
* Average edge weight of metabolite pairs on the same pathway.
* Average edge weight of metabolite pairs on different pathways.
* Percentage of positive, negative and zero edges that contain metabolite pairs on the same pathway.

Finally, the plugin accepts an output file in plaintext (TXT) format, which will contain
the textual description of each pathway from #4 above, followed by the number of positive, zero
and negative edge weights between pairs of nodes on that pathway, the percentage of these
edges that are positive, and the average edge weight over all these edges.

