# Other tools

metaTOR 2.0 clusters contigs, and nothing else: the commands that did what
hicstuff, metacooler or instaGRAAL do were removed. This page lists each one,
with what replaces it. `MAP` is the contact map of a run, `metator.metacool`
in its output directory by default.

| Removed | Use instead |
|---|---|
| `metator qc` | `metacooler qc MAP --per-genome qc.tsv`: the contacts within contigs, within clusters (the genomes of the map) and between them, per cluster. For the uncut, loop and weird events of a library: hicstuff's `pipeline -f --plot`. |
| `metator contactmap` | `metacooler extract MAP --genomes NAME -o NAME.cool`: the contact map of a cluster, a `.cool` file, at the bin size of the map. A map of restriction fragments: hicstuff's `pipeline` on the FASTA file of the cluster. |
| The plots of `pipeline` and `validation`, the heatmap included | `metacooler show MAP -o map.png`: the contacts between the clusters (`--level genome`, the default), or of some of them (`--genomes`). `cluster_summary.txt` gives the completeness, contamination and size of each cluster. |
| `metator metacool` | `pipeline` and `validation` write the map of the clusters. By hand: `metacooler reassign MAP --contigs contig_data_final.txt --contigs-format Name,Final_cluster -o clusters.metacool`, then `metacooler annotate` for the columns of the tables. |
| `metator scaffold`, `pipeline --scaffold` | hicstuff's `pipeline -M graal` on the FASTA file of a cluster, then instaGRAAL on its matrix. |
| `metator pairs` | pairix's tools: sort the pairs by `chr1`, `chr2`, `pos1` and `pos2`, compress them with `bgzip`, then index them with `pairix -p pairs`. |

See [Migrating to 2.0](migration.md) for every removed option and output.
