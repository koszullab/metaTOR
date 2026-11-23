# Other tools

metaTOR 2.0 clusters contigs, and nothing else: it starts from a contact map
of the contigs, and the commands that did what hicstuff, metacooler,
instaGRAAL or pairix do were removed. This page lists each one, with what
replaces it. `MAP` is the contact map of a run: the map that it wrote in its
output directory (or the map that it started from, which it updates in place
with `--inplace`); `NAME` a cluster of it, and `R1.fq.gz` and `R2.fq.gz` the
Hi-C reads of the run.

| Removed | Use instead |
|---|---|
| Reads, BAM and pairs inputs (`--forward`, `--reverse`, `--start`, `--aligner`, `--aligner-mode`, `--min-quality`) | `hicstuff pipeline -M metacool` makes the map from the reads; `metacooler create pairs` from pairs files, or `metacooler create cool` from a `.cool` file; `metacooler merge` sums the maps of several libraries. See [From the reads to the map](#from-the-reads-to-the-map). |
| `metator qc` | `metacooler qc MAP --per-genome qc.tsv`: the contacts within contigs, within clusters (the genomes of the map) and between them, per cluster. For the uncut, loop and weird events of a library: hicstuff's `pipeline -g assembly.fa -e ENZYME -f --plot R1.fq.gz R2.fq.gz`, whose `plots/` show them. |
| `metator contactmap` | `metacooler extract MAP --genomes NAME -o NAME.cool`: the contact map of a cluster, a `.cool` file, at the bin size of the map. A map of restriction fragments: hicstuff's `pipeline -g NAME.fa -e ENZYME R1.fq.gz R2.fq.gz` on the FASTA file of the cluster. |
| The plots of `pipeline` and `validation`, the heatmap included | `metacooler show MAP -o map.png`: the contacts between the clusters (`--level genome`, the default); between the contigs of some of them with `--genomes NAME,... --level contig`. The genomes table of the map (`metacooler dump -t genomes -H MAP`) gives the completeness, redundancy and size of each cluster, as `cluster_summary.txt`. |
| `metator metacool` | `pipeline` and `validation` write the map of the clusters. By hand, from the map of the contigs (`contigs.metacool`) and the tables of a run of 1.x: `metacooler reassign contigs.metacool --contigs contig_data_final.txt --contigs-format metator --genomes bin_summary.txt --genomes-format 1 -o clusters.metacool` (`--contigs-format metator` reads the cluster of a contig in `Final_cluster`, then in 1.x's `Final_bin`; the summary of a 2.0 development build is `cluster_summary.txt`), then `metacooler annotate clusters.metacool --contigs contig_data_final.txt --contigs-format Name --columns contigs.GC_content,contigs.Hit,contigs.Shotgun_coverage,contigs.Restriction_site --inplace` for the columns of the contigs. |
| `metator scaffold`, `pipeline --scaffold` | hicstuff's `pipeline -g NAME.fa -e ENZYME -M graal -o NAME R1.fq.gz R2.fq.gz` on the FASTA file of a cluster, then instaGRAAL on its matrix (`NAME/abs_fragments_contacts_weighted.txt`, `fragments_list.txt` and `info_contigs.txt`). |
| `metator pairs` | hicstuff's `pipeline` writes its pairs sorted by `chr1`, `chr2`, `pos1` and `pos2`, and compressed with bgzip (`valid_idx.pairs.gz`, `valid_idx_filtered.pairs.gz` with `-f`): `pairix -p pairs` indexes them as they are. Another pairs file: sort it so, compress it with `bgzip`, then index it with `pairix -p pairs`. metaTOR reads no pairs file, and no longer comes with `pypairix`: `metacooler create pairs` makes the map from a pairs file, sorted or not. |

See [Migrating to 2.0](migration.md) for every removed option and output.

## From the reads to the map

hicstuff aligns the reads of a library on the assembly and makes its map, with
no clustering yet: one bin per restriction fragment, and every contig in the
genome `__unbinned__`. metacooler sums the maps of several libraries:

```bash
# 1. The map of each library
hicstuff pipeline -g assembly.fa -e DpnII -M metacool -o hicstuff_L1 -P L1 L1_R1.fq.gz L1_R2.fq.gz
hicstuff pipeline -g assembly.fa -e DpnII -M metacool -o hicstuff_L2 -P L2 L2_R1.fq.gz L2_R2.fq.gz

# 2. One map of both libraries, their contacts summed
metacooler merge hicstuff_L1/L1.metacool hicstuff_L2/L2.metacool -o sample.metacool

# 3. The clusters: binning/sample.metacool, the map of the clusters
metator pipeline sample.metacool -a assembly.fa -e DpnII -o binning
```

!!! note "`-M metacool` needs hicstuff 4.0.0"

    `hicstuff pipeline -M metacool`, and its mapping mode `5spm`, come with
    hicstuff 4.0.0: 3.3.0 and earlier have neither. With an earlier
    hicstuff, make the map from its pairs files with `metacooler create
    pairs`, as below.

hicstuff aligns the reads with bowtie2 by default; `-m 5spm` aligns them with
`bwa mem -5SPM`, close to metaTOR 1.x's `--aligner bwa`, which ran
`bwa mem -5SP`. Besides the map,
`<prefix>.metacool`, it writes the pairs of the library,
`<prefix>.valid_idx.pairs.gz`. The map holds the contact statistics of the
reads of each contig, which `metator mge` reads, at 1001 bp with metacooler
0.4.0 or later (see
[Clustering mobile genetic elements](guides/mge.md#running-it)).

From the pairs files of another aligner, `metacooler create pairs` makes the
map, one bin per contig or of a fixed size (`--bins 1000`), taking the lengths
of the contigs from the `#chromsize` lines of the pairs files, or from
`--fasta`:

```bash
metacooler create pairs --pairs L1.pairs.gz --pairs L2.pairs.gz --bins 1000 --no-ps \
    --assembly-name assembly -o sample.metacool
```

`--no-ps` leaves out the P(s) table that metacooler writes next to the map by
default, and `--assembly-name` names the assembly of the map, `unknown`
otherwise. A `.cool` file of the contigs becomes a map with
`metacooler create cool`. [Starting from a contact map](guides/metacool.md#starting-from-a-contact-map)
says what metaTOR reads from the map.
