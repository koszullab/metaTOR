# Other tools

metaTOR 2.0 only clusters contigs: the 1.x commands that duplicated hicstuff,
metacooler, instaGRAAL or pairix are gone. Below, `MAP` is the contact map a
run wrote, `NAME` one of its clusters, and `R1.fq.gz`, `R2.fq.gz` the Hi-C
reads.

| Removed | Use instead |
|---|---|
| Reads inputs, and `--forward`, `--reverse`, `--start`, `--aligner`, `--aligner-mode` | `hicstuff pipeline -M metacool` (from reads), `metacooler create pairs` (from pairs files) or `metacooler create cool` (from a `.cool` file) makes the map; `metacooler merge` sums libraries (see [below](#from-the-reads-to-the-map)). `network` and `pipeline` still take pairs or BAM files (see [The inputs of a run](cli.md#inputs)). |
| `metator qc` | `metacooler qc MAP --per-genome qc.tsv`: contacts within contigs, within clusters and between them. For uncut, loop and weird events: `hicstuff pipeline -g assembly.fa -e ENZYME -f --plot R1.fq.gz R2.fq.gz`. |
| `metator contactmap` | `metacooler extract MAP --genomes NAME -o NAME.cool`. For restriction fragments: `hicstuff pipeline -g NAME.fa -e ENZYME R1.fq.gz R2.fq.gz` on the cluster's FASTA file. |
| The heatmap of the clusters and the GC and coverage plots (`network_heatmap.png`, `MAGs_*_distribution.png`) | `metacooler show MAP -o map.png` (`--genomes NAME,... --level contig` for contigs). `metacooler dump -t contigs MAP` gives each contig's GC content, hits, coverage and cluster. The quality plots remain, in `plots/`. |
| `metator metacool` | `pipeline` and `validation` write the map of the clusters. From a 1.x run's tables and a map of the contigs: `metacooler reassign contigs.metacool --contigs contig_data_final.txt --contigs-format metator --genomes bin_summary.txt --genomes-format 1 -o clusters.metacool`, then `metacooler annotate clusters.metacool --contigs contig_data_final.txt --contigs-format Name --columns contigs.GC_content,contigs.Hit,contigs.Shotgun_coverage,contigs.Restriction_site --inplace`. |
| `metator scaffold`, `pipeline --scaffold` | `hicstuff pipeline -g NAME.fa -e ENZYME -M graal -o NAME R1.fq.gz R2.fq.gz` on the cluster's FASTA file, then instaGRAAL on its matrix (`NAME/abs_fragments_contacts_weighted.txt`, `fragments_list.txt`, `info_contigs.txt`). |
| `metator pairs` | hicstuff's pairs (`valid_idx.pairs.gz`) are sorted and bgzipped, ready for `pairix -p pairs`. Other pairs files: sort, `bgzip`, then `pairix -p pairs`. metaTOR no longer needs `pypairix`. |

[Migrating to 2.0](migration.md) lists every removed option and output.

## From the reads to the map

hicstuff makes the map of each library, one bin per restriction fragment,
every contig in the genome `__unbinned__`. metacooler sums them:

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

    hicstuff 3.3.0 and earlier have neither `-M metacool` nor the mapping
    mode `5spm`: make the map from their pairs files with
    `metacooler create pairs` (below).

- hicstuff aligns with bowtie2 by default; `-m 5spm` runs `bwa mem -5SPM`,
  close to 1.x's `--aligner bwa` (`bwa mem -5SP`).
- Besides the map, `<prefix>.metacool`, hicstuff writes the pairs,
  `<prefix>.valid_idx.pairs.gz`.
- The map holds the contact statistics that `metator mge` needs (metacooler
  0.4.0 or later; see [Clustering mobile genetic elements](guides/mge.md#running-it)).

### From pairs files

`metacooler create pairs` makes the map from the pairs files of any aligner,
one bin per contig or fixed bins (`--bins 1000`), with contig lengths from
the `#chromsize` lines or `--fasta`:

```bash
metacooler create pairs --pairs L1.pairs.gz --pairs L2.pairs.gz --no-ps \
    --assembly-name assembly -o sample.metacool
```

- One bin per contig without `--bins` needs metacooler 0.5.0; with 0.4.0,
  give `--bins contigs.bed` (see
  [From the reads to the map](guides/metacool.md#from-the-reads-to-the-map)).
- `--no-ps` skips the P(s) table; `--assembly-name` names the assembly
  (`unknown` otherwise).
- `metacooler create cool` makes a map from a `.cool` file of the contigs.

[Starting from a contact map](guides/metacool.md#starting-from-a-contact-map)
says what metaTOR reads from the map.

### Pairs or BAM files as the run's input

`metator network` and `metator pipeline` also make the map themselves, one
bin per contig (see [The inputs of a run](cli.md#inputs)):

```bash
# The pairs files of hicstuff: binning/sample.metacool, made from them, becomes the map of the clusters
metator pipeline hicstuff_L1/L1.valid_idx.pairs.gz hicstuff_L2/L2.valid_idx.pairs.gz -p sample \
    -a assembly.fa -e DpnII -o binning
```

BAM files work too: hicstuff's, kept with `-n` (`tmp/<prefix>.for.bam` and
`tmp/<prefix>.rev.bam`), or another aligner's with both mates
(`bwa mem -5SP`). metaTOR keeps the primary alignments of both mates with a
mapping quality of at least `--min-quality`, neither duplicate nor QC-failed.
