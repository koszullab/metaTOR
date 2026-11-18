# Alignment and network

`metator.align` aligns the metaHiC reads, with bowtie2 or bwa and hicstuff's
mapping modes, and writes the pairs of contigs each read pair joins.
`metator.network` builds the contig network from them, or from a `.metacool`
contact map, and normalises its edges. See
[Aligning reads and building the network](../guides/network.md).

::: metator.align

::: metator.network
