# Migrating to metaTOR 2.0

metaTOR 2.0 no longer aligns reads. It starts from a `.metacool` contact map
of the contigs, made by hicstuff from the reads (`hicstuff pipeline -M
metacool`), by metacooler from pairs files (`metacooler create pairs`), or by
`network` and `pipeline` themselves from pairs or BAM files. It writes the
contact map of the clusters, and `mge` reads that map.

This page lists the changes since 1.3.10, the last 1.x release (or since
1.4.0, where a row says so): the command line, the inputs a run now refuses,
and what a run writes and counts. The [changelog](changelog.md) has the
details.

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

Several libraries, which 1.x took as comma-separated lists: run hicstuff once
per library, then sum the maps.

```bash
metacooler merge hic_A/A.metacool hic_B/B.metacool -o AB.metacool
metator pipeline AB.metacool --assembly assembly.fa \
    --enzyme DpnII,HinfI --outdir out/
```

MGE clusters: 1.x read the tables and pairs files of a run (which `pipeline`
removes without `--no-clean-up`); 2.0 reads the map of its clusters.
`metator dump` writes the MGE clusters and their hosts in `clusters.txt` and
`cluster_members.txt`.

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
  on the same contigs, and 2.0 built the same 263,961 edges (weights equal to
  rounding) and the same hits on all 20,566 contigs.
- **Equivalent clusters.** Louvain is random in both versions (2.0 runs it
  through igraph). 2.0 agreed with three 1.3.10 runs as closely as they agreed
  with each other (length-weighted adjusted Rand index 0.977–0.992, against
  0.971–0.992).
- **Same quality.** 24 HQ and 7 MQ clusters, as in each 1.3.10 run; 29 of the
  39 genomes had at least 90% of their sequence in one cluster of at least 95%
  purity, against 28–29 for 1.3.10.
- **Several libraries.** On a metagenome sequenced in two libraries (355 Mb
  assembly), the map summed by `metacooler merge` gave the network 1.3.10
  built from both, and clusters of the same quality.
- **`metator mge`.** On the 169 contigs of the mock's 13 plasmids and 5 phages:
  the same contacts and the same 127 MGE clusters, but for the score of one
  whose contigs tied (2.0 merges ties longest contig first, 1.x in assembly
  order; see [What a run writes](#what-a-run-writes)).

### Comparing with the outputs of 1.x

- Count the quality classes of a 1.x run from its plot,
  `plot/bins_size_distribution.png`, not its log, which counted HQ clusters as
  MQ (see the [changelog](changelog.md#fixed)). On the mock, the log said 7 HQ
  and 24 MQ, the plot 24 and 7.
- Take its clusters from `final_bin_unscaffold/` or `contig_data_final.txt`:
  its `binning.txt` left out the contigs of clusters partitioned again but
  kept whole (27% of the binned bases on the mock).
- To run 1.3.10 again, pin `pandas<3`: under pandas 3, which bioconda installs,
  1.3.10 stops in the partition.

## The command line

<div class="mt-versions" markdown>

| metaTOR 1.x | metaTOR 2.0 | Why |
|---|---|---|
| `metator pipeline -1 R1 -2 R2 -a ASSEMBLY` | `hicstuff pipeline -g ASSEMBLY -e ENZYME -M metacool -o hicstuff -P sample R1 R2`, then `metator pipeline hicstuff/sample.metacool -a ASSEMBLY` | `pipeline` and `network` take the contact map, `MAP`, as first argument: a `.metacool` file or a `.mmetacool` entry (`file.mmetacool::entries/<name>`). hicstuff 4.0.0+ makes it from reads, metacooler from pairs files (see [Other tools](other-tools.md#from-the-reads-to-the-map)); pairs or BAM files also work (below). |
| `metator pipeline -1 FOR.bam -2 REV.bam -S bam -a ASSEMBLY`, `metator pipeline -1 L1.bam,L2.bam -S bam -b bwa -a ASSEMBLY`, `metator pipeline -1 L1.pairs,L2.pairs -S pair -a ASSEMBLY` | `metator pipeline FOR.bam REV.bam -a ASSEMBLY`, `metator pipeline L1.bam,L2.bam -a ASSEMBLY`, `metator pipeline L1.pairs,L2.pairs -a ASSEMBLY` | Pairs or BAM files are the arguments (comma-separated lists or glob patterns), and their type is detected. BAM files hold both mates (as `-b bwa`) or one mate each (bowtie2, as `-S bam`): give the R1 files, then the R2 files. The run makes the map, one bin per contig, named after `--prefix` or the assembly (see [The inputs of a run](cli.md#inputs)). |
| `metator network -1 R1 -2 R2 -a ASSEMBLY` | `metator network MAP -a ASSEMBLY`, or `metator network PAIRS_OR_BAM... -a ASSEMBLY` | As `pipeline`. The network is stored in the map the run writes to its output directory (`-F` replaces one there; `--inplace` writes to `MAP`), where `partition` and `validation` read it. |
| `metator partition -c contig_data_network.txt -n network.txt -a ASSEMBLY` | `metator partition MAP` | `partition` reads the network from the map. `-c`, `-n` and `-a` are unknown options: it reads no assembly and writes no FASTA file. |
| `metator validation -c contig_data_partition.txt -n network.txt -f DIR -a ASSEMBLY` | `metator validation MAP -a ASSEMBLY` | `validation` reads the network and the partition from the map. `-n`, `-c` and `-f` are unknown options. |
| `-s, --size` of `partition`: the size of the overlapping clusters it wrote | Removed | `partition` keeps every cluster. `validation --size`, like `pipeline --size`, chooses the overlapping clusters that miComplete checks. |
| A network of your own, given to `partition -n` and `validation -n` (HiCzin weights, for example) | Removed | `partition` and `validation` read the network that `metator network` stored in the map. metacooler computes other normalisations (`metacooler normalize hiczin`, `normcc`) for its own commands. |
| `-t, --threads` of `network`: the threads that sorted the network | The threads that read BAM files and bin pairs into the map | The network is summed in memory: there is no file to sort. |
| `metator partition`, `metator validation` | `partition` stores its partition in the map, which `validation` reads | No tables or FASTA files pass between them: see [What a run writes](#what-a-run-writes). |
| `-1, --forward` and `-2, --reverse` of `network` and `pipeline`: reads, alignments or pairs, as comma-separated lists for several libraries | Removed. Alignments and pairs are the arguments of the run | hicstuff 4.0.0+ makes the map from reads (`hicstuff pipeline -M metacool R1 R2`; see [Other tools](other-tools.md#from-the-reads-to-the-map)). Several libraries: one map each, summed by `metacooler merge`, or one hicstuff run on comma-separated R1 and R2 lists. BAM and pairs files: see [The inputs of a run](cli.md#inputs). |
| `-S, --start fastq`, `bam` or `pair` | Removed | `fastq`: use hicstuff. `bam` and `pair`: give the files as arguments; their type is read from their first bytes. |
| `-b, --aligner bowtie2` or `bwa` | Removed | hicstuff aligns with bowtie2 by default; `-m 5spm` (4.0.0) runs `bwa mem -5SPM`, close to `--aligner bwa` (`bwa mem -5SP`). |
| `-B, --aligner-mode normal`, `iterative` or `cutsite` | Removed | hicstuff's `--mapping` (`-m`). |
| `-q, --min-quality` | Minimum mapping quality of both mates in BAM inputs, 30 by default; refused with a map or pairs files | For reads, hicstuff's `--quality-min` (`-q`), also 30 by default. |
| `--edge N` (`-E`) of `network` and `pipeline` | Removed | The network counts every pair: the map holds no read positions. Even at `--edge 0`, 1.x left some pairs out: see [What a run counts](#what-a-run-counts). |
| `-a, --assembly` as the basename of a bowtie2 or bwa index | Refused: `--assembly` is a FASTA file, plain or gzipped | metaTOR aligns no reads. |
| `pipeline --no-metacool` (1.4.0) | Removed | A run always writes the map of the clusters, to its output directory under the name of `MAP` (with `--inplace`, over `MAP`). |
| | `pipeline --inplace` | New: writes the map of the clusters over `MAP` instead of to the output directory. |
| | `--seed` of `partition`, `validation`, `pipeline` and `mge` | New: seeds Louvain or Leiden, and `mge --random`. The same map and seed give the same clusters. |
| | `metator partition MAP`, `metator validation MAP` | New: the steps run one by one on a map, which they update as `pipeline MAP` does. |
| | `--inplace` of `network`, `partition`, `validation` and `mge` | New: update `MAP` in place. A run whose output directory already holds `MAP` does so without it (steps run one by one in one directory). |
| `metator mge -c contig_data_final.txt -b binning.txt -m MGES -a ASSEMBLY PAIRSFILE...` | `metator mge MAP -m MGES -a ASSEMBLY` | `mge` reads the contacts and the contig data from the map. The map needs the network of `metator network`, and the contact statistics within contigs that hicstuff (`-M metacool`) or `metacooler create pairs` count at 1001 bp (metacooler 0.4.0+). Pairs files, `-c` and `-b` are gone; so are the pypairix and pyfastx dependencies. |
| | `mge --inplace`, `mge -F` | New: `mge` stores its MGE clusters in the map, as the round `mge`: in the map of its output directory (replaced only with `-F`), or in `MAP` with `--inplace`. |
| | `mge` refuses `-n/--network`, `-f/--fasta`, `-c/--contigs-data`, `-M/--method`, `--checkv-db`, `-d/--depth`, `-p/--plot`, `-t/--threads` and `-S/--threshold-asso` of 1.3.10: unknown options | 1.4.0 removed the MetaBAT input, host association, CheckV and plots of `mge`, and renamed options (`-a/--assembly` for `--fasta`, `-c/--contigs` for `--contigs-data`). 2.0 reads the map; `metator host` finds the hosts. |
| `mge --random`: a random clustering, written as the real one | The round `mge_random` of the map, and `fasta/mges_random/` | The real clustering (round `mge`, `fasta/mges/`) is left untouched. |
| `metator qc` | Removed | `metacooler qc MAP --per-genome qc.tsv`; hicstuff's `pipeline -e ENZYME -f --plot` for the 3C events (see [Other tools](other-tools.md)). |
| `metator contactmap` | Removed | `metacooler extract MAP --genomes NAME -o NAME.cool`; or hicstuff on the FASTA file of the cluster. |
| `metator scaffold`, and `-v, --scaffold` and `-J, --junctions` of `pipeline` | Removed | hicstuff's `pipeline -M graal` on the FASTA file of a cluster, then instaGRAAL. |
| `metator pairs` | Removed | hicstuff writes its pairs sorted and bgzipped, ready for `pairix -p pairs`; pairix's tools handle other pairs files. |
| `metator metacool` (1.4.0) | Removed | `pipeline` and `validation` write the map of the clusters; by hand, `metacooler reassign` and `annotate` (see [Other tools](other-tools.md)). |
| `-C, --cluster-matrix` of `partition`, `validation` and `pipeline`: `clustering_matrix_partition.npz`, `clustering_matrix_recursive.npz`, `clustering_matrix.npz` | Removed | No replacement for the share of iterations that put two contigs together. `clusters.txt`, `cluster_members.txt` and the `/clusters` group of the map give the contigs of every cluster of every round. |
| `-s, --self-contacts` of `network` | Removed | The network has no self-edges. The contacts within each contig are in the map's `qc_intra_contig_short` and `qc_intra_contig_long` columns (metacooler's QC block, unless `--no-qc`). |
| `python -m metator.main`, `metator.main.main()` | Removed: use `metator`, or `python -m metator` | An alias of the command line since 1.4.0. |
| `metator host -n NETWORK -b BINNING -m MGES -c CONTIG_DATA [-t 0.1]`, and the association of `mge` (`-S/--threshold-asso`, 0.1) | `metator host MAP [-s 10] [-c 5]` | `host` was broken in 1.3.10 and removed in 1.4.0. It is back: it reads the MGE clusters (round `mge`), the MAGs (the genomes) and the network from the map, finds the hosts of each MGE cluster (not each contig), and stores them in the map (see [the guide](guides/mge.md#the-hosts)). A MAG is a host with at least 10% of the MGE cluster's signal with the MAGs (`-s`; 1.3.x: 0.1 of its signal with any contig) and at least 5 contigs in contact with it (`-c`, `--min-contigs`; new). The 1.x options are unknown options. |
| `-A spinglass` (`validation`, `pipeline`) | Refused | It never worked. |
| `-p, --prefix` of `partition`, `validation` and `pipeline`: any characters | Letters, digits, `_` and `-` only | miComplete drops dots from cluster names. |
| Abbreviated long options (`--thread` for `--threads`) | Refused: give each option in full | The command line is rich-click's, not docopt's. |
| Command names in any case (`metator Pipeline …`) | Refused | Same. |
| `-T, --tmpdir DIR`: temporary files in DIR, which the run removed | In a private subdirectory of DIR, which the run removes, leaving DIR | 1.x removed DIR, whatever it held. With `-N` (`--no-clean-up`), the files are kept in `DIR/<random>/`. |
| A usage error exits with status 1 | Status 2 | An unknown option or value, a non-numeric number, or a missing `MAP`, refused before the command starts. |
| A `text_credentials.json` file in the working directory sent every error message to the web API it named | Removed | |

</div>

## What a run refuses

A run now refuses these inputs, which 1.x took:

| What | Why |
|---|---|
| A `--depth` file that lacks a contig of the assembly, or gives one another length | It is read by contig name. 1.x read it line by line, and gave wrong values to contigs listed in another order. |
| An assembly whose FASTA file starts with blank or `#` lines, with Biopython 1.87 or later | Biopython refuses it: start the file with its first `>` line. |

## What a run writes

<div class="mt-versions" markdown>

| metaTOR 1.x | metaTOR 2.0 |
|---|---|
| | The contact map of the contigs, grouped by cluster: `pipeline MAP`, `partition MAP` and `validation MAP` write it to the output directory under the name of `MAP`, which they leave unchanged. `--inplace`, or an output directory that already holds `MAP`, updates `MAP` itself. |
| | `clusters.txt` and `cluster_members.txt`, in the temporary directory (`-N` keeps it): the clusters of every round and their contigs, stored in the map's `/clusters` group. |
| `contig_data_final.txt`, `bin_summary.txt`, `binning.txt`, miComplete's results (`overlapping_micomplete_results.txt`, `recursive_micomplete_step_<n>.txt`) and `recursive_bin/`, in the output directory (`validation`, `pipeline`) | In the temporary directory of the run, which the log names and `-N` keeps. The output directory holds the map, `fasta/mags_with-mges/`, `plots/` and the log; `pipeline` adds `fasta/mags/`, `cool/`, `tables/` and the report. `metator dump MAP -o DIR` writes the tables from the map: `contig_data_final.txt` (no `ID` column or row number), `cluster_summary.txt`, `binning.txt`, `clusters.txt`, `cluster_members.txt` and `network.txt`, rows in the order of the map. |
| | `mge`: the round `mge` of the map, after the rounds of the run, and the `mge` entry of its metadata; with `--random`, `mge_random`. `metator dump` writes the rounds in `clusters.txt` and `cluster_members.txt`. |
| `alignment_<i>.bam` (`--aligner bwa`), `alignment_<i>_for.bam` and `alignment_<i>_rev.bam` (bowtie2): the alignments of library `i` | None: hicstuff aligns the reads (and keeps its BAM files with `-n`). BAM files are inputs of a run, which takes the same contacts from them (primary alignments of both mates, mapping quality at least `--min-quality`, no duplicates or QC failures), each mate at its 5' end where 1.x took its leftmost base. The network is the same. |
| `alignment_<i>_sorted.pairs.gz` and its index (`network`; `pipeline` with `-N`): the pairs of library `i` | None: hicstuff writes `<prefix>.valid_idx.pairs.gz`. A run from BAM files writes `library_<i>.pairs` (from 1) in its temporary directory (`-N`). |
| `network_<i>.txt`: the network of library `i`; `hit_data_alignment.txt`: the hits of each contig in each library | None: the map sums the libraries. |
| `network.txt` and `contig_data_network.txt` (`network`, `pipeline`) | None in the output directory; `contig_data_network.txt` is in the temporary directory (`-N`). `network` stores the network in the map: the contig columns `GC_content`, `Hit`, `Shotgun_coverage` and `Restriction_site`, the bin weights `metator_<normalization>`, and the `network` metadata entry. `metator dump` writes `network.txt` with contig names where 1.x wrote IDs (see [The network](formats.md#the-network)). `metacooler network MAP --level contig --balance metator_<normalization>` writes it too, with the contacts within each contig, which are no edge: drop them. |
| The rows of `contig_data_partition.txt` and `contig_data_final.txt`, and the contig IDs, in assembly order, every contig of the assembly | In map order (longest first), the contigs of the map only: a contig's ID is its position in the map, and clusters and FASTA records are numbered in that order. The shotgun abundance of a cluster is its share of the map's contigs, not of the whole assembly. `mge` leaves out, with a warning, an MGE contig the map lacks. |
| 1.4.0: `<outdir>/metator.metacool` (`<prefix>.metacool`), the contact map of the clusters | `<outdir>/<MAP>`: named after `MAP` (`sample.metacool` from `sample.metacool`), or over `MAP` with `--inplace`. `--prefix` names only the clusters, except for a run from pairs or BAM files, whose map is `<prefix>.metacool`, or `<assembly>.metacool` without `--prefix`. |
| `bin_summary.txt` | `cluster_summary.txt`, in the temporary directory (`-N`); `metator dump` writes it from the map. |
| `overlapping_bin/` (`partition`, `pipeline`) | None: the round `overlapping` of the map. The validation writes the FASTA files it checks in `overlapping_clusters/` of its temporary directory (`-N`). |
| `contig_data_partition.txt` (`partition`) | None: the rounds `core` and `overlapping` and the `partition` metadata entry of the map. The validation writes the table in its temporary directory (`-N`); `metator dump`'s `contig_data_final.txt` holds the core and overlapping clusters of each contig. |
| `recursive_bin/` | `recursive_clusters/`, in the temporary directory (`-N`). |
| `final_bin/` (`validation`), `final_bin_unscaffold/` (`pipeline`) | `fasta/mags_with-mges/`, one FASTA file per final cluster; from `pipeline`, also `fasta/mags/` (each MAG as its contigs and stitched, contigs joined by 200 `N`) and `cool/`, a balanced cool file per MAG (see [Outputs](cli.md#conventions)). |
| | The results of a 1.x run in the output directory of `validation` or `pipeline` (its tables, `*_bin/` folders, and the PNG plots of `plot/`) are refused, even with `-F`, and never removed: they would pass for the results of the run. `-F` replaces only what the run writes. Move them away, or use another `-o`. |
| `scaffold_bin/` (`pipeline --scaffold`) | None: `--scaffold` is gone (see [Other tools](other-tools.md)). |
| `mges_bin_summary.tsv` (`mge`) | None: the round `mge` of the map holds the MGE clusters, with their length, GC content, hits and score (see [MGE clusters](formats.md#mge-clusters)), unrounded, without the three host columns that `mge` left empty. `metator dump` writes it in `clusters.txt` (`BinLength`, `ContigsNumber` as `length`, `n_contigs`) and `cluster_members.txt`. |
| `mges_binned.fa` (`mge`): the lines of each contig as in the assembly, joined by 200 `N` | `fasta/mges/`, two files per MGE cluster: `MetaTOR_MGE_<ID>.stitched.fa` (its contigs longest first, ties by name, joined by 200 `N`, in lines of 80 bases) and `MetaTOR_MGE_<ID>.fa` (one record per contig). |
| The order of the MGE contigs (`mge`): the order of the assembly | Longest first, ties by name: equal scores are merged, and single-contig clusters numbered, in that order. For an assembly sorted by length (SPAdes), the MGE clusters and names match 1.x except among contigs of equal length. Otherwise, two single-contig clusters may swap names, and a tie of scores may merge other contigs. |
| `<tmp>/recursive_clustering/recursive_bins/`: the subnetwork of each cluster partitioned again (`validation`, `pipeline`; kept with `-N`) | None: the validation partitions the subnetworks in memory. |
| The contig data columns `Core_bin_ID`, `Core_bin_contigs`, `Core_bin_size`, the same for `Overlapping_bin_*` and `Recursive_bin_*`, and `Final_bin` | `bin` becomes `cluster`: `Core_cluster_ID`, ... `Recursive_cluster_size`, and `Final_cluster`. No step reads a 1.x table: rename its columns first (`metator.io.rename_1x_columns`). For metacooler, `--contigs-format Name,Final_cluster` reads a 2.0 table. |
| 1.4.0: a map of the clusters that could not be written gave a warning, and exit status 0 | The run stops with an error. Its other outputs, and the map it would have replaced, are kept. |
| The log, `metator_<command>_<date>.log`, with `## Command used:` in its header | `<name>.metator_<date>.log`, `<name>` being the map's name without `.metacool` (`sample.metator_20261123225838.log`). Its header has `## cmd:` and one `## <option>: <value>` line per option, as hicstuff's; its last line says what the run wrote and how long it took. `pipeline` numbers its steps (`[1/4] network` ... `[4/4] export`), and the validation its recursive steps (see the [CLI reference](cli.md#conventions)). |
| `plot/` (`validation`, `pipeline`): `bins_distribution.png`, `bins_size_distribution.png`, `MAGs_GC_distribution.png`, `MAGs_HiC_cov_distribution.png`, `MAGs_SG_cov_distribution.png` and `network_heatmap.png` | `plots/bins_distribution.pdf` and `plots/bins_size_distribution.pdf` (see [Viewing the clusters](guides/binning.md#viewing-the-clusters)). Contamination is drawn from 0 to 105%, a cluster without markers at 0 (1.x: −100%). The labels say clusters, not bins; the total size is that of the map's contigs, not the assembly; typos are fixed (`Others contigs`, `Completness`). Earlier plots in `plots/` are replaced only with `-F`. The other plots are gone: use `metacooler show` and the contigs table of the map (see [Other tools](other-tools.md)). |
| 1.4.0: the parameters of the run in the map's metadata (its `metator` key), `edge`, `scaffold`, `junctions` and `cluster_matrix` among them | One entry per step under `metator` (`network`, `partition`, `validation`), each with the version, command line and parameters of the run that made it, and what it read of the map (`read`); `map` names the map the run started from. `network` adds its normalization and the columns read in place of `--enzyme` and `--depth`; `partition` its normalization and `network_digest`, which `validation` checks; `validation` its normalization and how `cluster_summary.txt` is written from the genomes table (`summary_rows`, `summary_columns`, `summary_formats`). Gone: `edge`, `scaffold`, `junctions`, `cluster_matrix`, `forward`, `reverse`, `start`, `aligner`, `aligner_mode`, and 1.4.0's top-level `version`, `command` and `parameters`. |

</div>

## What a run counts

<div class="mt-versions" markdown>

| metaTOR 1.x | metaTOR 2.0 |
|---|---|
| At `--edge 0`, the default, without `--self-contacts`, the network left out a pair with a read at position 0 or on the last base of its contig. | Such a pair is a contact, as in the map: on `tests_data/alignment.pairs`, 360 of 13,415 edges change. |
| A pair with a read past the end of its contig: left out of the network (but with `--self-contacts`), its hits counted. | Left out, hits included. |
| With `--self-contacts` (`-s`), an edge from each contig to itself | No such edge: `--self-contacts` is gone. |
| `mge`, from pairs files: a pair within a contig with a read at position 0 and its mate 1001 bp away counted as a contact more than 1 kb apart | The map counts such a read on the first base of its contig, so the mates are 1000 bp apart and the contact is not counted. |

</div>

See [The contact map of the clusters](guides/metacool.md) for a run from a
map and the map a run writes.
