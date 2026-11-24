# Migrating to metaTOR 2.0

metaTOR 2.0 starts from a `.metacool` contact map of the contigs, which
hicstuff makes from the reads (`hicstuff pipeline -M metacool`) and metacooler
from pairs files (`metacooler create pairs`), and writes the contact map of
the clusters over it. It no longer reads reads, alignments or pairs files,
`mge` included, which reads the map. This page lists every
change to the command line since metaTOR 1.3.10, the last release of 1.x
(and since 1.4.0, where a row says so), with what
replaced what was removed; the inputs that 1.x took and a run now refuses; and
what a run writes and counts differently. The
[changelog](changelog.md) has the details, the errors that changed, and the
other changes to the clusters and outputs (those of 1.4.0, for example:
Louvain and Leiden run through igraph).

## From 1.x to 2.0

metaTOR 2.0 no longer maps reads: hicstuff (from 4.0.0) maps them and writes
a `.metacool` contact map, from which metaTOR clusters the contigs.

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
per library, then sum the maps:

```bash
metacooler merge hic_A/A.metacool hic_B/B.metacool -o AB.metacool
metator pipeline AB.metacool --assembly assembly.fa \
    --enzyme DpnII,HinfI --outdir out/
```

MGE clusters, from the outputs of a run: in 1.x, its tables and its pairs
files, which `pipeline` removes without `--no-clean-up`; in 2.0, the map of
its clusters, whose MGE clusters, and their hosts, `metator dump` writes in
`clusters.txt` and `cluster_members.txt`:

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

- hicstuff kept the same 10,894,611 read pairs as 1.3.10, each on the same two
  contigs, and metaTOR 2.0 built the same network: 263,961 edges, their
  weights equal to floating-point rounding, and the same hits on all 20,566
  contigs.
- The clusters differ from run to run in both versions (Louvain is random, and
  2.0 runs it through igraph): 2.0 agreed with three 1.3.10 runs as closely as
  these agreed with each other (adjusted Rand index over the contigs that
  either run clusters, weighted by length: 0.977–0.992, against 0.971–0.992).
- Same quality: 24 HQ and 7 MQ clusters in 2.0, as in each 1.3.10 run, and the
  same genomes recovered: 29 of the 39 had at least 90% of their sequence in
  one cluster of at least 95% purity, against 28–29 for 1.3.10.
- Several libraries: on a metagenome sequenced in two libraries (355 Mb
  assembly), each library gave the same pairs in both workflows, the map that
  `metacooler merge` sums gave the network that 1.3.10 partitioned from the
  two, and the clusters were of the same quality, alone or merged.
- `metator mge`, on the 169 contigs of the mock's 13 plasmids and 5 phages:
  the same contacts between them and within each, and the same 127 MGE
  clusters, but for the score of one, whose contigs tied: 2.0 merges two
  equal scores in the order of the contigs, the longest first, where 1.x
  followed the assembly (see [What a run writes](#what-a-run-writes)).

### Comparing with the outputs of 1.x

- Count the quality classes of a 1.x run from its plot,
  `plot/bins_size_distribution.png`, not from its log, which counted the HQ
  clusters as MQ (see the [changelog](changelog.md#fixed)): on the mock, the
  log of a run said 7 HQ and 24 MQ, its plot 24 and 7.
- Take its clusters from `final_bin_unscaffold/` or `contig_data_final.txt`:
  its `binning.txt` left out the contigs of the clusters partitioned again but
  kept whole (27% of the binned bases on the mock).
- To run 1.3.10 again, pin `pandas<3`: bioconda installs pandas 3, under which
  1.3.10 stops in the partition.

## The command line

<div class="mt-versions" markdown>

| metaTOR 1.x | metaTOR 2.0 | Why |
|---|---|---|
| `metator pipeline -1 R1 -2 R2 -a ASSEMBLY` | `hicstuff pipeline -g ASSEMBLY -e ENZYME -M metacool -o hicstuff -P sample R1 R2`, then `metator pipeline hicstuff/sample.metacool -a ASSEMBLY` | `pipeline` and `network` take the contact map of the contigs, `MAP`, as their first argument, which they require: a `.metacool` file, or an entry of a `.mmetacool` collection (`file.mmetacool::entries/<name>`). hicstuff makes it from the reads (`-M metacool` comes with hicstuff 4.0.0), metacooler from pairs files: see [Other tools](other-tools.md#from-the-reads-to-the-map). |
| `metator network -1 R1 -2 R2 -a ASSEMBLY` | `metator network MAP -a ASSEMBLY` | As `pipeline`. The network is stored in the map that the run writes to its output directory (`-F` to replace one there; `--inplace`: in `MAP` itself): `partition` and `validation` read it there. |
| `metator partition -c contig_data_network.txt -n network.txt -a ASSEMBLY` | `metator partition MAP` | The map, which `metator network MAP` stored the network in, is the input of `partition`, which requires it: its `-c` and `-n` are unknown options, and so is `-a`: it reads no assembly, and writes no FASTA file. |
| `metator validation -c contig_data_partition.txt -n network.txt -f DIR -a ASSEMBLY` | `metator validation MAP -a ASSEMBLY` | As `partition`: it requires the map, which holds the network and the partition that `metator partition MAP` stored there; `-n`, `-c` and `-f` are unknown options. |
| `-s, --size` of `partition`: the size of the overlapping clusters it wrote | Removed: an unknown option | `partition` keeps every cluster. `validation --size` chooses the overlapping clusters that miComplete checks, as `pipeline --size` does. |
| A network of your own, given to `partition -n` and `validation -n` (the weights of HiCzin, for example) | Removed | `partition` and `validation` read the network that `metator network` stored in the map, with its normalization (`--normalization`). metacooler computes such normalisations as weights of a map (`metacooler normalize hiczin`, `normcc`), which its own commands read; metaTOR partitions its own only. |
| `-t, --threads` of `network`: the threads that sorted the network | Removed: an unknown option | The network is no file to sort: it is summed in memory. |
| `metator partition`, `metator validation` | `partition` stores its partition in the map, which `validation` reads | `partition` writes no table and no FASTA file, and `validation` reads none: see [What a run writes](#what-a-run-writes). |
| `-1, --forward` and `-2, --reverse` of `network` and `pipeline`: the reads, their alignments or their pairs, in comma-separated lists for several libraries | Removed: unknown options | hicstuff makes the map from the reads: `hicstuff pipeline -M metacool R1 R2`, from hicstuff 4.0.0 (see [Other tools](other-tools.md#from-the-reads-to-the-map)). For several libraries, the map of each, which `metacooler merge` sums, as metaTOR summed their contacts; or one hicstuff run, which reads comma-separated lists of R1 and R2 files as one library. |
| `-S, --start fastq`, `bam` or `pair` | Removed: an unknown option | `fastq`: hicstuff's `pipeline`. `bam`: hicstuff's `pipeline -S bam`, on the BAM files of the two mates. `pair`: `metacooler create pairs --pairs FILE` (`--pairs` once per file, their contacts summed). |
| `-b, --aligner bowtie2` or `bwa` | Removed: an unknown option | hicstuff aligns the reads: with bowtie2 by default; `-m 5spm` (hicstuff 4.0.0) aligns the mates with `bwa mem -5SPM`, close to `--aligner bwa`, which ran `bwa mem -5SP`. |
| `-B, --aligner-mode normal`, `iterative` or `cutsite` | Removed: an unknown option | hicstuff's `--mapping` (`-m`), whose modes these were. |
| `-q, --min-quality` | Removed: an unknown option | hicstuff's `--quality-min` (`-q`), 30 by default, as here. |
| `--edge N` (`-E`) of `network` and `pipeline` | Removed: an unknown option | The network counts every pair, as its contact map does, which holds no read positions. Even at `--edge 0`, the default, 1.x left some pairs out: see [What a run counts](#what-a-run-counts). |
| `-a, --assembly` given as the basename of a bowtie2 or bwa index | Refused: `--assembly` is the FASTA file of the assembly, plain or compressed with gzip | metaTOR aligns no reads: hicstuff does, from its own index of the assembly. |
| `pipeline --no-metacool` (1.4.0): no contact map of the clusters | Removed: an unknown option | The map of the clusters is what a run writes, to its output directory, under the name of `MAP`, which it leaves as it was; `--inplace` writes it over `MAP` instead. |
| | `pipeline --inplace` | New: writes the map of the clusters over `MAP`, in place, rather than to the output directory. |
| | `metator partition MAP`, `metator validation MAP` | New: the steps run one by one on a `.metacool` map, their first argument, which they update as `pipeline MAP` does: `partition MAP` stores its partition (every core and overlapping cluster) in the clustering rounds of the map, and `validation MAP` reads it there and writes the map of the clusters. |
| | `--inplace` of `network`, `partition`, `validation` and `mge` | New, with `MAP`: update `MAP` in place, rather than write the map of the run to the output directory. A run whose output directory holds `MAP` already updates it in place without it: the steps of a run, one by one, in one output directory (`network MAP -o DIR`, then `partition DIR/MAP -o DIR`...). |
| `metator mge -c contig_data_final.txt -b binning.txt -m MGES -a ASSEMBLY PAIRSFILE...` | `metator mge MAP -m MGES -a ASSEMBLY` | `mge` takes the contact map of the run, `MAP`, as its first argument, which it requires: a map that `metator network` stored its network in (after `network`, `pipeline` or `validation`; it reads no cluster of the map). The contacts between the MGE contigs come from it, and those within each, between reads more than 1 kb apart, from its contact statistics (the QC block that metacooler counts from the pairs, at 1001 bp: hicstuff's `pipeline -M metacool` and `metacooler create pairs`, with metacooler 0.4.0 or later). Its pairs files, `-c` and `-b` are gone: the size, GC content and hits of each contig come from the map, which `metator network` stored them in. A map without those contact statistics (or with them measured from its pixels, or at another distance) or without the network of `network` is refused before anything is written, and so is an assembly that is not a FASTA file, plain or compressed with gzip. pairix and pypairix are no longer needed, and metaTOR no longer depends on pypairix, nor on pyfastx: metacooler writes its FASTA files. |
| | `mge --inplace`, `mge -F` | New: `mge` stores its MGE clusters in the map, as the round `mge` of its clustering rounds: in the map of its output directory (replaced only with `-F`), or in `MAP` itself with `--inplace`. The summary of the MGE clusters that an earlier version wrote in the output directory is refused, with or without `-F`, and never removed (see [What a run writes](#what-a-run-writes)). |
| | `mge` refuses `-n/--network`, `-f/--fasta`, `-c/--contigs-data`, `-M/--method`, `--checkv-db`, `-d/--depth`, `-p/--plot`, `-t/--threads` and `-S/--threshold-asso` of 1.3.10: unknown options | 1.4.0 removed the MetaBAT input, the association with the hosts, CheckV and the plots of `mge`, and renamed its options (`-a/--assembly` for `--fasta`, `-c/--contigs` for `--contigs-data`); 2.0 takes the contacts and the contig data from the map, and `metator host` finds the hosts of the MGE clusters (see below). |
| `mge --random`: a random clustering, written as the real one | The round `mge_random` of the map, and `mges_binned_random.fa` | The real clustering, the round `mge` and `mges_binned.fa` of an earlier run, are left as they were. |
| `metator qc` | Removed: an unknown command | metacooler and hicstuff measure a library: `metacooler qc MAP --per-genome qc.tsv` for the contacts within contigs, within clusters and between them, and hicstuff's `pipeline -e ENZYME -f --plot` for the 3C events (see [Other tools](other-tools.md)). |
| `metator contactmap` | Removed: an unknown command | `metacooler extract MAP --genomes NAME -o NAME.cool` writes the contact map of a cluster; hicstuff, on the FASTA file of the cluster, a map of restriction fragments. |
| `metator scaffold`, and `-v, --scaffold` and `-J, --junctions` of `pipeline` | Removed: an unknown command, unknown options | hicstuff's `pipeline -M graal` on the FASTA file of a cluster, then instaGRAAL. |
| `metator pairs` | Removed: an unknown command | hicstuff's `pipeline` writes its pairs sorted and compressed with bgzip, ready for `pairix -p pairs`; pairix's tools sort, compress and index other pairs files. |
| `metator metacool` (1.4.0) | Removed: an unknown command | `pipeline` and `validation` write the contact map of the clusters; by hand, `metacooler reassign` and `annotate` (see [Other tools](other-tools.md)). |
| `-C, --cluster-matrix` of `partition`, `validation` and `pipeline`: the clustering matrices, `clustering_matrix_partition.npz`, `clustering_matrix_recursive.npz` and `clustering_matrix.npz` | Removed: an unknown option | No replacement for the share of the iterations that put each two contigs together. The clustering rounds of the run give the contigs of each cluster of each round: `clusters.txt` and `cluster_members.txt`, and the `/clusters` group of the map. |
| `-s, --self-contacts` of `network` | Removed: an unknown option | The network has no edge from a contig to itself. The contacts within each contig are in the map, in its `qc_intra_contig_short` and `qc_intra_contig_long` columns, when it holds metacooler's QC block (the default; not with `--no-qc`). |
| `python -m metator.main`, `metator.main.main()` | Removed: `metator`, or `python -m metator` | The former entry point, an alias of the command line since 1.4.0. |
| `metator host -n NETWORK -b BINNING -m MGES -c CONTIG_DATA [-t 0.1]`, and the association of `mge` (`-S/--threshold-asso`, 0.1) | `metator host MAP [-s 10] [-c 5]` | `host` failed from 1.3.10, as it called a function that the rewrite of the `host` module removed, and 1.4.0 removed it. It is back, on the map that `metator mge` stored its MGE clusters in (the round `mge`), whose genomes are the MAGs (the final clusters, or the bins of another binner made its genomes with `metacooler reassign`), with the network of `metator network`: it finds the hosts of each MGE cluster, not of each MGE contig, and stores them in the map (see [Clustering mobile genetic elements](guides/mge.md#the-hosts)). Its rule is that of the `host` module since 1.3.10, not that of the association of 1.3.x: a MAG is a host with at least 10% of the signal of the MGE cluster with the MAGs (`-s`), where 1.3.x counted 0.1 of its signal with every other contig, and at least 5 of its contigs in contact with the MGE cluster (`-c`), which 1.3.x did not require. Its options, which read the tables of 1.x, are unknown options, and `-c` is `--min-contigs`: a contig data file given to it is refused as a number of contigs (a usage error). |
| `-A spinglass` (`validation`, `pipeline`) | Refused: not a value of `--algorithm` | It never worked: its code was commented out. |
| `-p, --prefix` of `partition`, `validation` and `pipeline`: any characters | Letters, digits, `_` and `-` only | miComplete drops the dots from the names of the clusters. |
| Long options abbreviated (`--thread` for `--threads`) | Refused: give each option in full | The command line is rich-click's, where it was docopt's. |
| Command names in any case (`metator Pipeline …`) | Refused: an unknown command | The command line is rich-click's, where metaTOR capitalised the name it gave docopt. |
| `-T, --tmpdir DIR`: the temporary files in DIR, which the run removed at its end | In a private subdirectory of DIR, which the run removes, leaving DIR. `validation` and `pipeline` refuse a DIR in `final_clusters/` of their output directory, or in a directory of the results of an earlier run there (`recursive_clusters/`, `recursive_bin/`...: a usage error) | 1.x removed the directory it was given, whatever else it held. With `-N` (`--no-clean-up`), the files kept are in `DIR/<random>/`. The run makes `final_clusters/` anew, and `-F` removed DIR with it. |
| A usage error exits with status 1 | Status 2 | An unknown option, an unknown value of `--algorithm` or `--normalization`, a value that is not a number for a numeric option, or a missing `MAP`, is a usage error, refused before the command starts. |
| A `text_credentials.json` file in the working directory sent every error message to the web API it named | Removed | |

</div>

## What a run refuses

A run now refuses these inputs, which 1.x took, with the reason logged as an
error, before any work unless said otherwise:

| What | Why |
|---|---|
| A `--depth` file that lacks a contig of the assembly, or gives one another length | The file is read by contig name. 1.x read it line by line, in the order of the assembly, and gave the contigs of a file in another order the values of others. |
| An assembly whose FASTA file starts with blank or `#` lines, with Biopython 1.87 or later | Biopython refuses it: start the file with its first `>` header line. |

## What a run writes

<div class="mt-versions" markdown>

| metaTOR 1.x | metaTOR 2.0 |
|---|---|
| | The contact map of the contigs, grouped by cluster: `pipeline MAP`, and the steps `partition MAP` and `validation MAP`, write it to the output directory, under the name of `MAP`, which they leave as it was, through metacooler (in an edit session: a copy of `MAP` staged next to the map they write, edited, validated, then renamed into place), and lock it while they run (metacooler's `<map>.lock`); `--inplace` updates `MAP` itself instead (through a symbolic link, the file it points to), as does a run whose output directory holds `MAP` already. |
| | `clusters.txt` and `cluster_members.txt`, in the temporary directory of the run (`-N` keeps it): the clusters of every clustering round, and their contigs, which the map stores in its `/clusters` group. |
| `contig_data_final.txt`, `bin_summary.txt`, `binning.txt`, miComplete's (`overlapping_micomplete_results.txt`, `recursive_micomplete_step_<n>.txt`) and `recursive_bin/`, in the output directory (`validation`, `pipeline`) | In the temporary directory of the run, in a directory of its own that the log names, which `-N` keeps: the output directory holds the map of the clusters, `final_clusters/`, `plot/` and the log only. The map of the clusters holds what the tables hold, and `metator dump MAP -o DIR` writes `contig_data_final.txt`, `cluster_summary.txt` and `binning.txt` from it, with `clusters.txt`, `cluster_members.txt` and `network.txt`: its `contig_data_final.txt` has no `ID` column and no row number, and its rows, with those of `binning.txt` and `cluster_members.txt`, come in the order of the map of the clusters. An earlier run's tables and `recursive_clusters/` in the output directory are refused, with or without `-F`, and never removed: left there, they would pass for the run's. |
| | `mge`: the round `mge` of the clustering rounds of the map, its MGE clusters, after those of the run, and the `mge` entry of its metadata; with `--random`, the round and the entry `mge_random`. `metator dump MAP -o DIR` writes the rounds in `clusters.txt` and `cluster_members.txt`. |
| `alignment_<i>.bam` (`--aligner bwa`), `alignment_<i>_for.bam` and `alignment_<i>_rev.bam` (bowtie2): the alignments of library `i` | None: hicstuff aligns the reads (and keeps its alignments with `-n`). |
| `alignment_<i>_sorted.pairs.gz` and its index (`network`; `pipeline` with `-N`): the pairs of library `i` | None: hicstuff writes the pairs of a library, `<prefix>.valid_idx.pairs.gz`, sorted and compressed, which it makes the map from. |
| `network_<i>.txt`: the network of library `i`; `hit_data_alignment.txt`: the hits of each contig in each library (with several libraries) | None: the map of a run sums its libraries. |
| `network.txt` and `contig_data_network.txt` (`network`, `pipeline`) | None in the output directory: `contig_data_network.txt`, the contig data that `network` stores in the map, is in its temporary directory, which `-N` keeps. `network` stores the network in the map (the `GC_content`, `Hit`, `Shotgun_coverage` and `Restriction_site` columns of its contigs table, the weights `metator_<normalization>` of its bins, and the `network` entry of its metadata), which `partition` and `validation` read. `metator dump MAP -o DIR` writes `network.txt` from the map of the clusters of a run: the names of the two contigs of each edge, where 1.x wrote their IDs (see [The network](formats.md#the-network)). `metacooler network MAP --level contig --balance metator_<normalization>` writes the network of any map that holds it, by contig names, with the contacts within each contig, which are no edge: drop them. |
| The rows of `contig_data_partition.txt` and `contig_data_final.txt`, and the IDs of the contigs, in the order of the assembly, every contig of the assembly | In the order of the map (longest first, as hicstuff and metacooler make it), its contigs only: the ID of a contig is its position in the map. The clusters are numbered in that order, and so are the records of their FASTA files. The shotgun abundance of a cluster is its share of the sequenced bases of the contigs of the map: the clusters of a map made with a contig filter add up to 100%, where 1.x counted every contig of the assembly. `mge` reads the map: an MGE contig that the map lacks is left out, with a warning that names it. |
| 1.4.0: `<outdir>/metator.metacool` (`<prefix>.metacool`), the contact map of the clusters | `<outdir>/<MAP>`: the map of the clusters goes to the output directory under the name of `MAP` (`sample.metacool` from `sample.metacool`), or over `MAP` with `--inplace`; `--prefix` names the clusters only. |
| `bin_summary.txt` | `cluster_summary.txt`, in the temporary directory of the run (`-N`), which `metator dump` writes from the map |
| `overlapping_bin/` (`partition`, `pipeline`) | None: the overlapping clusters are in the map, the round `overlapping` of its clustering rounds. The validation writes the FASTA files of those it checks in `overlapping_clusters/` of its temporary directory (`-N` keeps it). |
| `contig_data_partition.txt` (`partition`) | None: the partition is in the map, the rounds `core` and `overlapping` of its clustering rounds and the `partition` entry of its metadata. The validation writes the table in its temporary directory (`-N` keeps it); the core and overlapping clusters of each contig are in the `contig_data_final.txt` that `metator dump` writes. |
| `recursive_bin/` | `recursive_clusters/`, in the temporary directory of the run (`-N`) |
| `final_bin/` (`validation`), `final_bin_unscaffold/` (`pipeline`) | `final_clusters/` |
| | The results of a run of 1.x in the output directory of `validation` or `pipeline` (its tables, `contig_data_*.txt`, `bin_summary.txt`, `binning.txt`, `network.txt`, miComplete's and the clustering matrices, and `recursive_bin/`, `final_bin/`, `final_bin_unscaffold/`, `overlapping_bin/` and `scaffold_bin/`, and the plots of `plot/` that 2.0 no longer draws, `MAGs_*_distribution.png` and `network_heatmap.png`) are refused, with or without `-F`, and never removed, as the tables of 2.0's development builds are: left there, they would pass for the results of the run. `-F` replaces only what the run writes itself, `final_clusters/`, the two plots of `plot/` and the map of the output directory. Move them away, or give the run another output directory with `-o`. |
| `scaffold_bin/` (`pipeline --scaffold`) | None: `--scaffold` is gone (see [Other tools](other-tools.md)). |
| `mges_bin_summary.tsv` (`mge`) | None: the round `mge` of the map holds the MGE clusters, with their length, GC content, hits and score (see [MGE clusters](formats.md#mge-clusters)), unrounded where the table rounded `GC` and `BinningScore` to two decimals, and without the three columns of their hosts, which `mge` left empty (`None`, `NA`, `NA`). `metator dump MAP -o DIR` writes the round in `clusters.txt` (`BinLength` and `ContigsNumber` as `length` and `n_contigs`) and `cluster_members.txt`. An earlier run's `mges_bin_summary.tsv` in the output directory is refused, with or without `-F`, and never removed: left there, it would pass for the run's. |
| `mges_binned.fa` (`mge`): the lines of each contig as in the assembly, joined by 200 `N` | The same records, each cluster's contigs the longest first (ties by name), joined by 200 `N`, in lines of 80 bases (metacooler writes it). |
| The order of the MGE contigs (`mge`): that of the rows of `contig_data_final.txt`, the order of the assembly | The longest first, ties by name, whatever the order of the map or of the assembly: two equal scores are merged in that order, and the clusters of one contig numbered in it. For an assembly sorted by length (SPAdes), the MGE clusters and their names are those of 1.x, but for contigs of the same length; for another, two clusters of one contig may swap names, and a tie of scores may merge other contigs. |
| `<tmp>/recursive_clustering/recursive_bins/`, in the temporary directory (`validation`, `pipeline`; kept with `-N`): the subnetwork of each cluster that the validation partitioned again | None: the validation partitions the subnetwork of each cluster in memory, and writes none. |
| The columns of the contig data tables `Core_bin_ID`, `Core_bin_contigs`, `Core_bin_size`, `Overlapping_bin_ID`, `Overlapping_bin_contigs`, `Overlapping_bin_size`, `Recursive_bin_ID`, `Recursive_bin_contigs`, `Recursive_bin_size` and `Final_bin` | `Core_cluster_ID`, `Core_cluster_contigs`, `Core_cluster_size`, `Overlapping_cluster_ID`, `Overlapping_cluster_contigs`, `Overlapping_cluster_size`, `Recursive_cluster_ID`, `Recursive_cluster_contigs`, `Recursive_cluster_size` and `Final_cluster`. No step renames the columns of a table of 1.x, and none reads a table of 1.x: they read the map. Rename the columns of a table of 1.x first (`metator.io.rename_1x_columns` does). metacooler's `--contigs-format metator` names `Final_bin`: give `--contigs-format Name,Final_cluster` for a table of 2.0. |
| 1.4.0: a map of the clusters that could not be written, a warning, and exit status 0 | The run stops with the reason, logged as an error. Its other outputs are kept, and so is the map it would have replaced. |
| The log, `metator_<command>_<date>.log`: `## Command used:` in its header, an empty line on each side of its line of dashes, and a line of dashes and the date at its end | `## cmd:`, then the value of each option (`## <option>: <value>`), and the line of dashes alone, as hicstuff's header; its last line says what the run wrote and how long it took (`Final contact map written after 0h 8m 31s: <map>`). `pipeline` marks its steps (`[1/3] network`...), the partition counts its iterations, the validation its recursive steps (`Recursive step 1/9`), and the bookkeeping of a run is logged at DEBUG, which the log does not show (see the [CLI reference](cli.md#conventions)). |
| `plot/` (`validation`, `pipeline`): `bins_distribution.png`, `bins_size_distribution.png`, `MAGs_GC_distribution.png`, `MAGs_HiC_cov_distribution.png`, `MAGs_SG_cov_distribution.png` and `network_heatmap.png` | `plot/` with `bins_distribution.png` and `bins_size_distribution.png`, the plots of the quality of the clusters (see [Viewing the clusters](guides/binning.md#viewing-the-clusters)): the contamination drawn from 0 to 105%, a cluster without a marker gene at 0, where 1.x drew it at −100%; the pie chart names clusters where 1.x named bins (`Contaminated clusters`, `Other clusters`, `Clusters threshold`), gives the `Total size of the contigs`, those of the map, where 1.x gave the `Total size of the assembly`, and `Other contigs` and `Completeness/Contamination (%)` are spelt so (1.x: `Others contigs`, `Completness`). The plots of an earlier run are refused without `-F`, and replaced with it, where 1.x wrote over them; the other files of `plot/` are left as they are, but for the other plots of 1.x, refused with or without `-F` (see above). The other plots are gone: `metacooler show` draws the contacts between the clusters from the map, and its contigs table holds the GC content, hits and coverage of each contig (see [Other tools](other-tools.md)). |
| 1.4.0: the parameters of the run in the metadata of its map (its `metator` key), `edge`, `scaffold`, `junctions` and `cluster_matrix` among them | One entry per step under the `metator` key: `network`, `partition` and `validation`, each the run that made it (the step, or `pipeline`), with the version of metaTOR, its command line and parameters, and what it read of the map (`read`: the columns of the contigs table, the weights of the bins and the clustering rounds). No `edge`, `scaffold`, `junctions` or `cluster_matrix`, nor `forward`, `reverse`, `start`, `aligner`, `aligner_mode` or `min_quality`. `map`: the map that the run started from. `network`: its normalization, and under `read`, the columns of the map that it read in place of `--enzyme` and `--depth` (`Restriction_site`, `Shotgun_coverage`; `columns_from_map` in a development build). `partition`: the normalization and the identity of the network that it partitioned (`network_digest`), which `validation` checks. `validation`: its normalization, and how the summary of the clusters, `cluster_summary.txt`, is written from the genomes table: `summary_rows`, the clusters in order; `summary_columns`, its columns; `summary_formats`, the printf format of each (`%.2f`, `%.4f`, `%d`, `%s`). The keys that 1.4.0 wrote at the top of `metator` (`version`, `command`, `parameters`; `columns_from_map` in a development build of 2.0) are dropped when a validation writes its entry. |

</div>

## What a run counts

<div class="mt-versions" markdown>

| metaTOR 1.x | metaTOR 2.0 |
|---|---|
| At `--edge 0`, the default, without `--self-contacts` (`-s`), the network left out a pair with a read at position 0 or on the last base of its contig. | Such a pair is a contact, as in the contact map: on `tests_data/alignment.pairs`, 360 of 13,415 edges change. |
| A pair with a read past the end of its contig: left out of the network, but with `--self-contacts`, which tested no positions; its hits counted. | Left out, and so are its hits. |
| With `--self-contacts` (`-s`), an edge from each contig to itself, its contacts within itself | No such edge: `--self-contacts` is gone. |
| `mge`: a pair within a contig with a read at position 0 (which a 1-based pairs file should not hold, but some do), from the pairs files: with its mate 1001 bp away, a contact between reads more than 1 kb apart | Read from the map, whose contact statistics count such a read on the first base of its contig (metacooler moves it there): its mate 1001 bp away is 1000 bp away, and the contact is not counted. |

</div>

See [the guide to the contact map](guides/metacool.md) for a run from a map
and the map a run writes.
