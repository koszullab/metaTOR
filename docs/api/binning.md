# Partition and validation

`metator.partition` partitions the network many times with the Louvain or
Leiden algorithm, and merges the contigs that always fall together into
clusters. `metator.clusters` holds the clustering rounds of a run, the
partition among them, as the map stores them. `metator.validation` measures
the clusters with miComplete, and partitions the contaminated ones again.
See [Partition and validation](../guides/binning.md).

::: metator.partition

::: metator.clusters

::: metator.validation
