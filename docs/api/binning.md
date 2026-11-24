# Partition and validation

`metator.partition` partitions the network many times with the Louvain or
Leiden algorithm, and merges the contigs that always fall together into
clusters. `metator.clusters` holds the clustering rounds of a run, the
partition among them, as the map stores them. `metator.validation` measures
the clusters with miComplete, and partitions the contaminated ones again, and
`metator.figures` draws the plots of the final clusters.
`metator.markers` finds miComplete's marker genes on each contig of a map,
and gives miComplete's numbers of any group of contigs from them (see
[Marker genes of the contigs](../formats.md#marker-genes-of-the-contigs)).
See [Partition and validation](../guides/binning.md).

::: metator.partition

::: metator.clusters

::: metator.validation

::: metator.markers

::: metator.figures
