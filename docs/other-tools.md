# Other tools

metaTOR 2.0 clusters contigs, and nothing else: the commands that did what
hicstuff, metacooler, instaGRAAL or pairix do were removed. This page lists
each one, with what replaces it. `MAP` is the contact map of a run,
`metator.metacool` in its output directory by default, `NAME` a cluster of
it, and `R1.fq.gz` and `R2.fq.gz` the Hi-C reads of the run.

| Removed | Use instead |
|---|---|
| `metator qc` | `metacooler qc MAP --per-genome qc.tsv`: the contacts within contigs, within clusters (the genomes of the map) and between them, per cluster. For the uncut, loop and weird events of a library: hicstuff's `pipeline -g assembly.fa -e ENZYME -f --plot R1.fq.gz R2.fq.gz`, whose `plots/` show them. |
| `metator contactmap` | `metacooler extract MAP --genomes NAME -o NAME.cool`: the contact map of a cluster, a `.cool` file, at the bin size of the map. A map of restriction fragments: hicstuff's `pipeline -g NAME.fa -e ENZYME R1.fq.gz R2.fq.gz` on the FASTA file of the cluster. |
| The plots of `pipeline` and `validation`, the heatmap included | `metacooler show MAP -o map.png`: the contacts between the clusters (`--level genome`, the default); between the contigs of some of them with `--genomes NAME,... --level contig`. `cluster_summary.txt` gives the completeness, redundancy and size of each cluster. |
| `metator metacool` | `pipeline` and `validation` write the map of the clusters. By hand, from the map of the contigs (`contigs.metacool`) and the tables of a run: `metacooler reassign contigs.metacool --contigs contig_data_final.txt --contigs-format Name,Final_cluster --genomes cluster_summary.txt --genomes-format 1 -o clusters.metacool`, then `metacooler annotate clusters.metacool --contigs contig_data_final.txt --contigs-format Name --columns contigs.GC_content,contigs.Hit,contigs.Shotgun_coverage,contigs.Restriction_site --inplace` for the columns of the contigs. |
| `metator scaffold`, `pipeline --scaffold` | hicstuff's `pipeline -g NAME.fa -e ENZYME -M graal -o NAME R1.fq.gz R2.fq.gz` on the FASTA file of a cluster, then instaGRAAL on its matrix (`NAME/abs_fragments_contacts_weighted.txt`, `fragments_list.txt` and `info_contigs.txt`). |
| `metator pairs` | hicstuff's `pipeline` writes its pairs sorted by `chr1`, `chr2`, `pos1` and `pos2`, and compressed with bgzip (`valid_idx.pairs.gz`, `valid_idx_filtered.pairs.gz` with `-f`): `pairix -p pairs` indexes them as they are. Another pairs file: sort it so, compress it with `bgzip`, then index it with `pairix -p pairs`. pip's `pypairix` installs neither the `pairix` nor the `bgzip` command; in Python, `pypairix.build_index(path, "pairs")` indexes a file. |

See [Migrating to 2.0](migration.md) for every removed option and output.
