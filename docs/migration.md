# Migrating to metaTOR 2.0

metaTOR 2.0 no longer aligns reads. It starts from a `.metacool` contact map
of the contigs, which hicstuff makes from the reads, metacooler from pairs
files, or `network` and `pipeline` from pairs or BAM files, and writes the
contact map of the clusters.

This page lists the changes since 1.3.10, the last 1.x release (or since
1.4.0, where a row says so). The [changelog](changelog.md) has the details.

## From 1.x to 2.0

hicstuff (4.0.0 or later) aligns the reads and writes the map; metaTOR
clusters the contigs from it.

```bash
# 1.x
metator pipeline --assembly assembly.fa \
    --forward R1.fq.gz --reverse R2.fq.gz \
    --enzyme DpnII,HinfI --outdir out/

# 2.0
hicstuff pipeline --genome assembly.fa --enzyme DpnII,HinfI \
    --matrix metacool --prefix sample --outdir hic/ R1.fq.gz R2.fq.gz
metator pipeline hic/sample.metacool --assembly assembly.fa \
    --enzyme DpnII,HinfI --outdir out/
```

Several libraries, which 1.x took as comma-separated lists: one hicstuff run
per library, then sum the maps.

```bash
metacooler merge hic_A/A.metacool hic_B/B.metacool -o AB.metacool
metator pipeline AB.metacool --assembly assembly.fa \
    --enzyme DpnII,HinfI --outdir out/
```

MGE clusters: 1.x read the tables and pairs files of a run; 2.0 reads the map
of its clusters. `metator dump` writes the MGE clusters and their hosts in
`clusters.txt` and `cluster_members.txt`.

```bash
# 1.x
metator mge --contigs out/contig_data_final.txt \
    --binning out/binning.txt --mges mges.txt --assembly assembly.fa \
    --outdir mge/ out/alignment_0_sorted.pairs.gz

# 2.0
metator mge out/sample.metacool --mges mges.txt \
    --assembly assembly.fa --outdir mge/
metator host mge/sample.metacool --outdir mge/
metator dump mge/sample.metacool --outdir tables/
```

### Same network, equivalent clusters

On a mock community of 39 genomes (12.8 M read pairs, 174 Mb assembly):

- **Same network.** hicstuff kept the same 10,894,611 read pairs as 1.3.10,
  and 2.0 built the same 263,961 edges (weights equal to rounding) and hits
  on all 20,566 contigs.
- **Equivalent clusters.** Louvain is random in both versions. 2.0 agreed
  with three 1.3.10 runs as closely as they agreed with each other
  (length-weighted adjusted Rand index 0.977–0.992, against 0.971–0.992).
- **Same quality.** 24 HQ and 7 MQ clusters, as in each 1.3.10 run; 29 of the
  39 genomes had at least 90% of their sequence in one cluster of at least 95%
  purity (1.3.10: 28–29).
- **Several libraries.** On a two-library metagenome (355 Mb assembly), the
  map summed by `metacooler merge` gave the network 1.3.10 built from both,
  and clusters of the same quality.
- **`metator mge`.** On the 169 contigs of the mock's 13 plasmids and 5
  phages: the same contacts and 127 MGE clusters, but for the score of one
  whose contigs tied (see [What a run writes](#what-a-run-writes)).

### Comparing with the outputs of 1.x

- Count a 1.x run's quality classes from `plot/bins_size_distribution.png`,
  not its log, which counted HQ clusters as MQ (see the
  [changelog](changelog.md#fixed)).
- Take its clusters from `final_bin_unscaffold/` or `contig_data_final.txt`:
  its `binning.txt` left out the contigs of clusters partitioned again but
  kept whole (27% of the binned bases on the mock).
- To run 1.3.10 again, pin `pandas<3`: it stops in the partition under
  pandas 3, which bioconda installs.

## The command line

<div class="mt-versions" markdown>

| metaTOR 1.x | metaTOR 2.0 | Why |
|---|---|---|
| `metator pipeline -1 R1 -2 R2 -a ASSEMBLY` | `hicstuff pipeline -g ASSEMBLY -e ENZYME -M metacool -o hicstuff -P sample R1 R2`, then `metator pipeline hicstuff/sample.metacool -a ASSEMBLY` | The first argument is the contact map, `MAP`: a `.metacool` file or a `.mmetacool` entry (`file.mmetacool::entries/<name>`). See [Other tools](other-tools.md#from-the-reads-to-the-map). |
| `metator pipeline -1 FOR.bam -2 REV.bam -S bam -a ASSEMBLY`, `metator pipeline -1 L1.bam,L2.bam -S bam -b bwa -a ASSEMBLY`, `metator pipeline -1 L1.pairs,L2.pairs -S pair -a ASSEMBLY` | `metator pipeline FOR.bam REV.bam -a ASSEMBLY`, `metator pipeline L1.bam,L2.bam -a ASSEMBLY`, `metator pipeline L1.pairs,L2.pairs -a ASSEMBLY` | Pairs or BAM files are arguments (lists or glob patterns), their type detected; BAM files of one mate each: the R1 files, then the R2 files. The run makes the map, named after `--prefix` or the assembly (see [The inputs of a run](cli.md#inputs)). |
| `metator network -1 R1 -2 R2 -a ASSEMBLY` | `metator network MAP -a ASSEMBLY`, or `metator network PAIRS_OR_BAM... -a ASSEMBLY` | As `pipeline`. The network is stored in the map the run writes (`-F` replaces one; `--inplace` writes to `MAP`). |
| `metator partition -c contig_data_network.txt -n network.txt -a ASSEMBLY` | `metator partition MAP` | It reads the network from the map, and no assembly. `-c`, `-n` and `-a` are unknown options. |
| `metator validation -c contig_data_partition.txt -n network.txt -f DIR -a ASSEMBLY` | `metator validation MAP -a ASSEMBLY` | It reads the network and the partition from the map. `-n`, `-c` and `-f` are unknown options. |
| `-s, --size` of `partition` | Removed | `partition` keeps every cluster; `validation --size` chooses those miComplete checks. |
| A network of your own (`partition -n`, `validation -n`; HiCzin weights, for example) | Removed | metaTOR partitions only the network that `metator network` stored. metacooler has other normalisations (`metacooler normalize hiczin`, `normcc`) for its own commands. |
| `-t, --threads` of `network`: to sort the network | To read BAM files and bin pairs into the map | The network is summed in memory. |
| `metator partition`, `metator validation` | `partition` stores its partition in the map, which `validation` reads | No tables or FASTA files pass between them (see [What a run writes](#what-a-run-writes)). |
| `-1, --forward` and `-2, --reverse` of `network` and `pipeline` | Removed: alignments and pairs are arguments | hicstuff makes the map from reads. Several libraries: one map each, summed by `metacooler merge`, or one hicstuff run on comma-separated R1 and R2 lists. |
| `-S, --start fastq`, `bam` or `pair` | Removed | `fastq`: use hicstuff. `bam`, `pair`: give the files as arguments. |
| `-b, --aligner bowtie2` or `bwa` | Removed | hicstuff uses bowtie2 by default; `-m 5spm` runs `bwa mem -5SPM`, close to `--aligner bwa` (`bwa mem -5SP`). |
| `-B, --aligner-mode normal`, `iterative` or `cutsite` | Removed | hicstuff's `--mapping` (`-m`). |
| `-q, --min-quality` | Mapping quality of both mates in BAM inputs (30); refused with a map or pairs files | For reads: hicstuff's `--quality-min` (`-q`, 30). |
| `--edge N` (`-E`) of `network` and `pipeline` | Removed | The map holds no read positions: every pair counts (see [What a run counts](#what-a-run-counts)). |
| `-a, --assembly` as a bowtie2 or bwa index basename | Refused: a FASTA file, plain or gzipped | metaTOR aligns no reads. |
| `pipeline --no-metacool` (1.4.0) | Removed | A run always writes the map of the clusters. |
| | `pipeline --inplace` | New: writes the map of the clusters over `MAP`, not to the output directory. |
| | `--seed` of `partition`, `validation`, `pipeline` and `mge` | New: the same map and seed give the same clusters (Louvain, Leiden, `mge --random`). |
| | `metator partition MAP`, `metator validation MAP` | New: the steps one by one, updating the map as `pipeline` does. |
| | `--inplace` of `network`, `partition`, `validation` and `mge` | New: update `MAP` in place. Implied when the output directory already holds `MAP`. |
| `metator mge -c contig_data_final.txt -b binning.txt -m MGES -a ASSEMBLY PAIRSFILE...` | `metator mge MAP -m MGES -a ASSEMBLY` | `mge` reads contacts and contig data from the map, which needs the network of `metator network` and metacooler's contact statistics (0.4.0+). Pairs files, `-c`, `-b`, pypairix and pyfastx are gone. |
| | `mge --inplace`, `mge -F` | New: `mge` stores its MGE clusters in the map (round `mge`): that of its output directory (replaced only with `-F`), or `MAP` with `--inplace`. |
| | `mge` refuses `-n/--network`, `-f/--fasta`, `-c/--contigs-data`, `-M/--method`, `--checkv-db`, `-d/--depth`, `-p/--plot`, `-t/--threads` and `-S/--threshold-asso` of 1.3.10 | 1.4.0 removed the MetaBAT input, host association, CheckV and plots of `mge`, and renamed `--fasta` to `-a/--assembly` and `--contigs-data` to `-c/--contigs`. `metator host` finds the hosts. |
| `mge --random`: written as the real clustering | The round `mge_random` and `fasta/mges_random/` | The real clustering is left untouched. |
| `metator qc` | Removed | `metacooler qc MAP --per-genome qc.tsv`; `hicstuff pipeline -e ENZYME -f --plot` (see [Other tools](other-tools.md)). |
| `metator contactmap` | Removed | `metacooler extract MAP --genomes NAME -o NAME.cool`. |
| `metator scaffold`, and `-v, --scaffold` and `-J, --junctions` of `pipeline` | Removed | `hicstuff pipeline -M graal` on a cluster's FASTA file, then instaGRAAL. |
| `metator pairs` | Removed | hicstuff's pairs are ready for `pairix -p pairs`. |
| `metator metacool` (1.4.0) | Removed | `pipeline` and `validation` write the map of the clusters; by hand, `metacooler reassign` and `annotate` (see [Other tools](other-tools.md)). |
| `-C, --cluster-matrix` of `partition`, `validation` and `pipeline` (`clustering_matrix*.npz`) | Removed | No replacement. `clusters.txt`, `cluster_members.txt` and the map's `/clusters` group give the contigs of every cluster of every round. |
| `-s, --self-contacts` of `network` | Removed | No self-edges. The contacts within each contig are in the map's `qc_intra_contig_short` and `qc_intra_contig_long` columns. |
| `python -m metator.main`, `metator.main.main()` | Removed: use `metator`, or `python -m metator` | An alias since 1.4.0. |
| `metator host -n NETWORK -b BINNING -m MGES -c CONTIG_DATA [-t 0.1]`, and the association of `mge` (`-S/--threshold-asso`, 0.1) | `metator host MAP [-s 10] [-c 5]` | Broken in 1.3.10, removed in 1.4.0, back on the map: it finds the hosts of each MGE cluster (round `mge`) among the MAGs, and stores them in the map. A host has at least 10% of the MGE cluster's signal with the MAGs (`-s`) and at least 5 contigs in contact with it (`-c`, `--min-contigs`). See [the guide](guides/mge.md#the-hosts). |
| `-A spinglass` (`validation`, `pipeline`) | Refused | It never worked. |
| `-p, --prefix` of `partition`, `validation` and `pipeline`: any characters | Letters, digits, `_` and `-` only | miComplete drops dots from cluster names. |
| Abbreviated long options (`--thread`) | Refused | The command line is rich-click's, not docopt's. |
| Command names in any case (`metator Pipeline`) | Refused | Same. |
| `-T, --tmpdir DIR`, removed at the end of the run | A private subdirectory of DIR is removed; DIR is left | 1.x removed DIR, whatever it held. `-N` keeps `DIR/<random>/`. |
| A usage error exits with status 1 | Status 2 | Unknown options or values, non-numbers, or a missing `MAP`. |
| A `text_credentials.json` file in the working directory sent every error message to the web API it named | Removed | |

</div>

## What a run refuses

| What | Why |
|---|---|
| A `--depth` file that lacks a contig of the assembly, or gives one another length | It is read by contig name. 1.x read it line by line, and mixed up contigs listed in another order. |
| An assembly whose FASTA file starts with blank or `#` lines, with Biopython 1.87 or later | Biopython refuses it: start the file with its first `>` line. |

## What a run writes

<div class="mt-versions" markdown>

| metaTOR 1.x | metaTOR 2.0 |
|---|---|
| | The contact map of the clusters: `pipeline MAP`, `partition MAP` and `validation MAP` write it to the output directory under the name of `MAP`, unless `--inplace` (or an output directory that holds `MAP`) updates `MAP` itself. |
| | `clusters.txt` and `cluster_members.txt`, in the temporary directory (`-N` keeps it): the clusters of every round, from the map's `/clusters` group. |
| `contig_data_final.txt`, `bin_summary.txt`, `binning.txt`, miComplete's results (`overlapping_micomplete_results.txt`, `recursive_micomplete_step_<n>.txt`) and `recursive_bin/`, in the output directory | In the temporary directory (`-N`). The output directory holds the map, `fasta/`, `plots/`, the log and, from `pipeline`, `cool/`, `tables/` and the report. `metator dump MAP -o DIR` writes the tables from the map, in map order; its `contig_data_final.txt` has no `ID` column or row number. |
| | `mge`: the round `mge` of the map and the `mge` entry of its metadata (`mge_random` with `--random`). |
| `alignment_<i>.bam`, `alignment_<i>_for.bam`, `alignment_<i>_rev.bam`: the alignments of library `i` | None: hicstuff aligns (and keeps its BAM files with `-n`). As inputs, BAM files give the same network: each mate is taken at its 5' end, where 1.x took its leftmost base. |
| `alignment_<i>_sorted.pairs.gz` and its index: the pairs of library `i` | None: hicstuff writes `<prefix>.valid_idx.pairs.gz`. A run from BAM files writes `library_<i>.pairs` in its temporary directory (`-N`). |
| `network_<i>.txt`, `hit_data_alignment.txt`: per library | None: the map sums the libraries. |
| `network.txt` and `contig_data_network.txt` (`network`, `pipeline`) | None in the output directory. `network` stores the network in the map: the contig columns `GC_content`, `Hit`, `Shotgun_coverage`, `Restriction_site`, and the bin weights `metator_<normalization>`. `metator dump` writes `network.txt`, with contig names where 1.x wrote IDs (see [The network](formats.md#the-network)). |
| Rows of `contig_data_partition.txt` and `contig_data_final.txt`, and contig IDs, in assembly order, every contig of the assembly | In map order (longest first), the map's contigs only; clusters and FASTA records are numbered in that order. The shotgun abundance of a cluster is its share of the map's contigs. `mge` leaves out, with a warning, MGE contigs the map lacks. |
| 1.4.0: `<outdir>/metator.metacool` (`<prefix>.metacool`) | `<outdir>/<MAP>` (`sample.metacool` from `sample.metacool`). `--prefix` names only the clusters, except from pairs or BAM files: the map is then `<prefix>.metacool`, or `<assembly>.metacool`. |
| `bin_summary.txt` | `cluster_summary.txt`, in the temporary directory (`-N`) and from `metator dump`. |
| `overlapping_bin/` (`partition`, `pipeline`) | None: the round `overlapping` of the map. The FASTA files that miComplete checks are in `overlapping_clusters/` of the temporary directory (`-N`). |
| `contig_data_partition.txt` (`partition`) | None: the rounds `core` and `overlapping` of the map. The validation writes it in its temporary directory (`-N`). |
| `recursive_bin/` | `recursive_clusters/`, in the temporary directory (`-N`). |
| `final_bin/` (`validation`), `final_bin_unscaffold/` (`pipeline`) | `fasta/mags_with-mges/`; from `pipeline`, also `fasta/mags/` (contigs and stitched) and `cool/` (see [Outputs](cli.md#conventions)). |
| | The results of a 1.x run (tables, `*_bin/` folders, `plot/`) in the output directory of `validation` or `pipeline` are refused, even with `-F`, and never removed. Move them, or use another `-o`. |
| `scaffold_bin/` (`pipeline --scaffold`) | None (see [Other tools](other-tools.md)). |
| `mges_bin_summary.tsv` (`mge`) | None: the round `mge` of the map holds the MGE clusters (see [MGE clusters](formats.md#mge-clusters)), unrounded, without the empty host columns. `metator dump` writes them in `clusters.txt` (`BinLength`, `ContigsNumber` as `length`, `n_contigs`). |
| `mges_binned.fa` (`mge`) | `fasta/mges/`: per MGE cluster, `MetaTOR_MGE_<ID>.stitched.fa` (contigs longest first, joined by 200 `N`) and `MetaTOR_MGE_<ID>.fa` (one record per contig). |
| The order of the MGE contigs (`mge`): the assembly order | Longest first, ties by name: equal scores are merged, and single-contig clusters numbered, in that order. For a length-sorted assembly (SPAdes), the result matches 1.x except among contigs of equal length. |
| `<tmp>/recursive_clustering/recursive_bins/` (`validation`, `pipeline`; `-N`) | None: subnetworks are partitioned in memory. |
| Contig data columns `Core_bin_ID`, `Core_bin_contigs`, `Core_bin_size`, the same for `Overlapping_bin_*` and `Recursive_bin_*`, and `Final_bin` | `bin` becomes `cluster` (`Core_cluster_ID`, ..., `Final_cluster`). No step reads a 1.x table: rename its columns with `metator.io.rename_1x_columns`. For metacooler, `--contigs-format Name,Final_cluster` reads a 2.0 table. |
| 1.4.0: a map of the clusters that could not be written gave a warning, and exit status 0 | The run stops with an error, keeping its other outputs and the map it would have replaced. |
| The log, `metator_<command>_<date>.log` | `<name>.metator_<date>.log`, `<name>` being the map's name without `.metacool` (`sample.metator_20261123225838.log`). Its header has `## cmd:` and one line per option, as hicstuff's; `pipeline` numbers its steps; the last line gives what the run wrote and how long it took (see the [CLI reference](cli.md#conventions)). |
| `plot/` (`validation`, `pipeline`): `bins_distribution.png`, `bins_size_distribution.png`, `MAGs_GC_distribution.png`, `MAGs_HiC_cov_distribution.png`, `MAGs_SG_cov_distribution.png` and `network_heatmap.png` | `plots/bins_distribution.pdf` and `plots/bins_size_distribution.pdf` (see [Viewing the clusters](guides/binning.md#viewing-the-clusters)): contamination from 0 to 105% (a cluster without markers at 0, not −100%), labels in clusters, not bins, and the total size of the map's contigs, not the assembly. Earlier plots are replaced only with `-F`. The other plots are gone (see [Other tools](other-tools.md)). |
| 1.4.0: the parameters of the run in the map's metadata (`metator` key), `edge`, `scaffold`, `junctions` and `cluster_matrix` among them | One entry per step under `metator` (`network`, `partition`, `validation`): version, command line, parameters, and what it read of the map; `partition` records the network it partitioned (`network_digest`), which `validation` checks, and `validation` how `cluster_summary.txt` is written. Gone: `edge`, `scaffold`, `junctions`, `cluster_matrix`, `forward`, `reverse`, `start`, `aligner`, `aligner_mode`, and 1.4.0's top-level `version`, `command` and `parameters`. |

</div>

## What a run counts

<div class="mt-versions" markdown>

| metaTOR 1.x | metaTOR 2.0 |
|---|---|
| At `--edge 0`, without `--self-contacts`, the network left out a pair with a read at position 0 or on the last base of its contig. | Such a pair counts, as in the map: on `tests_data/alignment.pairs`, 360 of 13,415 edges change. |
| A pair with a read past the end of its contig: left out (but with `--self-contacts`), its hits counted. | Left out, hits included. |
| With `--self-contacts` (`-s`), an edge from each contig to itself | No such edge. |
| `mge`, from pairs files: a read at position 0 and its mate 1001 bp away, in one contig, counted as more than 1 kb apart | The map counts the read on the first base of its contig: 1000 bp apart, not counted. |

</div>

See [The contact map of the clusters](guides/metacool.md) for a run from a
map and the map a run writes.
