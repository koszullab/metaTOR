# Contact maps

`metator.metacool` writes the contact map of a run, grouped by cluster, as a
`.metacool` file, from the map that the run started from,
with its clustering rounds and the normalization of its network. See
[The contact map of the clusters](../guides/metacool.md).

::: metator.metacool

`metator.inputs` reads the inputs of `metator network` and `metator
pipeline`, and makes the map of the contigs that a run starts from out of its
pairs files or BAM files (see [The inputs of a run](../cli.md#inputs)).

::: metator.inputs

`metator.dump` writes the tables of a run from the map of its clusters, as
`metator dump` does (see [`dump`](../cli.md#dump)).

::: metator.dump
