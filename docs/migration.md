# Migrating to metaTOR 2.0

metaTOR 2.0 starts from a `.metacool` contact map of the contigs, which
hicstuff makes from the reads (`hicstuff pipeline -M metacool`) and metacooler
from pairs files (`metacooler create pairs`), and writes the contact map of
the clusters over it. It no longer reads reads, alignments or pairs files
(but `mge`, which still counts contacts from the pairs). This page lists every
change to the command line since metaTOR 1.x (commit `4253551`), with what
replaced what was removed; the inputs that 1.x took and a run now refuses; and
what a run writes and counts differently. The
[changelog](changelog.md) has the details, the errors that changed, and the
other changes to the clusters and outputs (those of 1.4.0, for example:
Louvain and Leiden run through igraph).

## The command line

| metaTOR 1.x | metaTOR 2.0 | Why |
|---|---|---|
| `metator pipeline -1 R1 -2 R2 -a ASSEMBLY` | `hicstuff pipeline -g ASSEMBLY -e ENZYME -M metacool -o hicstuff -P sample R1 R2`, then `metator pipeline hicstuff/sample.metacool -a ASSEMBLY` | `pipeline` and `network` take the contact map of the contigs, `MAP`, as their first argument, which they require: a `.metacool` file, or an entry of a `.mmetacool` collection (`file.mmetacool::entries/<name>`). hicstuff makes it from the reads (`-M metacool` comes with hicstuff's next release), metacooler from pairs files: see [Other tools](other-tools.md#from-the-reads-to-the-map). |
| `metator network -1 R1 -2 R2 -a ASSEMBLY` | `metator network MAP -a ASSEMBLY` | As `pipeline`. The map is only read. |
| `metator partition`, `metator validation` | The same, but for the names of their outputs | A script that gives the overlapping clusters of `partition` to `validation -f` names `overlapping_clusters/`: see [What a run writes](#what-a-run-writes). |
| `-1, --forward` and `-2, --reverse` of `network` and `pipeline`: the reads, their alignments or their pairs, in comma-separated lists for several libraries | Removed: unknown options | hicstuff makes the map from the reads: `hicstuff pipeline -M metacool R1 R2`, from hicstuff's next release (see [Other tools](other-tools.md#from-the-reads-to-the-map)). For several libraries, the map of each, which `metacooler merge` sums, as metaTOR summed their contacts; or one hicstuff run, which reads comma-separated lists of R1 and R2 files as one library. |
| `-S, --start fastq`, `bam` or `pair` | Removed: an unknown option | `fastq`: hicstuff's `pipeline`. `bam`: hicstuff's `pipeline -S bam`, on the BAM files of the two mates. `pair`: `metacooler create pairs --pairs FILE` (`--pairs` once per file, their contacts summed). |
| `-b, --aligner bowtie2` or `bwa` | Removed: an unknown option | hicstuff aligns the reads: with bowtie2 by default; `-m 5spm` (hicstuff's next release) aligns the mates with `bwa mem -5SPM`, close to `--aligner bwa`, which ran `bwa mem -5SP`. |
| `-B, --aligner-mode normal`, `iterative` or `cutsite` | Removed: an unknown option | hicstuff's `--mapping` (`-m`), whose modes these were. |
| `-q, --min-quality` | Removed: an unknown option | hicstuff's `--quality-min` (`-q`), 30 by default, as here. |
| `--edge N` (`-E`) of `network` and `pipeline` | Removed: an unknown option | The network counts every pair, as its contact map does, which holds no read positions. Even at `--edge 0`, the default, 1.x left some pairs out: see [What a run counts](#what-a-run-counts). |
| | `pipeline --output-map FILE` | New: writes the map of the clusters of the run to FILE, and leaves the map of the run as it was. |
| | `metator partition MAP`, `metator validation MAP` | New: the steps run one by one on a `.metacool` map, their first argument, which they update as `pipeline MAP` does: `partition MAP` stores its overlapping clusters in the clustering rounds of the map, and `validation MAP` writes the map of the clusters. Without `MAP`, they run as before. |
| | `partition --output-map FILE`, `validation --output-map FILE` | New, with `MAP`: write the map to FILE, and leave `MAP` as it was. |
| | `mge --map MAP`, `mge --output-map FILE`, `mge -F` | New: `mge` also stores its MGE clusters in the map of the run, as the round `mge` of its clustering rounds, in place or in FILE (replaced only with `-F`). It still reads the pairs files, and refuses a map in their place. |
| | `pipeline --no-metacool` | New: writes no contact map of the clusters, and leaves the map of the run as it was. |
| `metator qc` | Removed: an unknown command | metacooler and hicstuff measure a library: `metacooler qc MAP --per-genome qc.tsv` for the contacts within contigs, within clusters and between them, and hicstuff's `pipeline -e ENZYME -f --plot` for the 3C events (see [Other tools](other-tools.md)). |
| `metator contactmap` | Removed: an unknown command | `metacooler extract MAP --genomes NAME -o NAME.cool` writes the contact map of a cluster; hicstuff, on the FASTA file of the cluster, a map of restriction fragments. |
| `metator scaffold`, and `-v, --scaffold` and `-J, --junctions` of `pipeline` | Removed: an unknown command, unknown options | hicstuff's `pipeline -M graal` on the FASTA file of a cluster, then instaGRAAL. |
| `metator pairs` | Removed: an unknown command | hicstuff's `pipeline` writes its pairs sorted and compressed with bgzip, ready for `pairix -p pairs`; pairix's tools sort, compress and index other pairs files. |
| `metator metacool` (1.4.0) | Removed: an unknown command | `pipeline` and `validation` write the contact map of the clusters; by hand, `metacooler reassign` and `annotate` (see [Other tools](other-tools.md)). |
| `-C, --cluster-matrix` of `partition`, `validation` and `pipeline`: the clustering matrices, `clustering_matrix_partition.npz`, `clustering_matrix_recursive.npz` and `clustering_matrix.npz` | Removed: an unknown option | No replacement for the share of the iterations that put each two contigs together. The clustering rounds of the run give the contigs of each cluster of each round: `clusters.txt` and `cluster_members.txt`, and the `/clusters` group of the map. |
| `-s, --self-contacts` of `network` | Removed: an unknown option | The network has no edge from a contig to itself. The contacts within each contig are in the map, in its `qc_intra_contig_short` and `qc_intra_contig_long` columns, when it holds metacooler's QC block (the default; not with `--no-qc`). |
| `python -m metator.main`, `metator.main.main()` | Removed: `metator`, or `python -m metator` | The former entry point, an alias of the command line since 1.4.0. |
| `metator host` | Removed: an unknown command | It always failed since 1.3.10: it called a function that the rewrite of the `host` module removed. |
| `-A spinglass` (`validation`, `pipeline`) | Refused: not a value of `--algorithm` | It never worked: its code was commented out. |
| `-p, --prefix` of `partition`, `validation` and `pipeline`: any characters | Letters, digits, `_` and `-` only | miComplete drops the dots from the names of the clusters. |
| Long options abbreviated (`--thread` for `--threads`) | Refused: give each option in full | The command line is rich-click's, where it was docopt's. |
| Command names in any case (`metator Pipeline …`) | Refused: an unknown command | The command line is rich-click's, where metaTOR capitalised the name it gave docopt. |
| `-T, --tmpdir DIR`: the temporary files in DIR, which the run removed at its end | In a private subdirectory of DIR, which the run removes, leaving DIR | 1.x removed the directory it was given, whatever else it held. With `-N` (`--no-clean-up`), the files kept are in `DIR/<random>/`. |
| A usage error exits with status 1 | Status 2 | An unknown option, an unknown value of `--algorithm` or `--normalization`, a value that is not a number for a numeric option, or a missing `MAP`, is a usage error, refused before the command starts. |
| A `text_credentials.json` file in the working directory sent every error message to the web API it named | Removed | |

## What a run refuses

A run now refuses these inputs, which 1.x took, with the reason logged as an
error, before any work unless said otherwise:

| What | Why |
|---|---|
| A `--depth` file that lacks a contig of the assembly, or gives one another length | The file is read by contig name. 1.x read it line by line, in the order of the assembly, and gave the contigs of a file in another order the values of others. |
| An assembly whose FASTA file starts with blank or `#` lines, with Biopython 1.87 or later | Biopython refuses it: start the file with its first `>` header line. |

## What a run writes

| metaTOR 1.x | metaTOR 2.0 |
|---|---|
| | The contact map of the contigs, grouped by cluster: `pipeline MAP`, and the steps `partition MAP` and `validation MAP`, update the map of the run in place (staged, validated, then renamed over it), and lock it while they run (`<map>.metator.lock`); `--output-map FILE` writes FILE instead. |
| | `clusters.txt` and `cluster_members.txt`: the clusters of every clustering round, and their contigs, also stored in the `/clusters` group of the map. |
| | `mge --map`: the round `mge` of the clustering rounds of the map, its MGE clusters, after those of the run. |
| `alignment_<i>.bam` (`--aligner bwa`), `alignment_<i>_for.bam` and `alignment_<i>_rev.bam` (bowtie2): the alignments of library `i` | None: hicstuff aligns the reads (and keeps its alignments with `-n`). |
| `alignment_<i>_sorted.pairs.gz` and its index (`network`; `pipeline` with `-N`): the pairs of library `i` | None: hicstuff writes the pairs of a library, `<prefix>.valid_idx.pairs.gz`, sorted and compressed, which `mge` reads once indexed with pairix (see [Clustering mobile genetic elements](guides/mge.md#running-it)). |
| `network_<i>.txt`: the network of library `i`; `hit_data_alignment.txt`: the hits of each contig in each library (with several libraries) | None: the map of a run sums its libraries. |
| 1.4.0: `<outdir>/metator.metacool` (`<prefix>.metacool`), the contact map of the clusters | None: the map of the clusters goes over `MAP`, or to `--output-map`; `--prefix` names the clusters only. |
| `bin_summary.txt` | `cluster_summary.txt` |
| `overlapping_bin/` | `overlapping_clusters/`: a script that gives it to `validation -f` names the new directory. |
| `recursive_bin/` | `recursive_clusters/` |
| `final_bin/` (`validation`), `final_bin_unscaffold/` (`pipeline`) | `final_clusters/` |
| | The outputs of a run of 1.x in the output directory (`overlapping_bin/`, `recursive_bin/`, `final_bin/`, `final_bin_unscaffold/`, `bin_summary.txt`) and the `plot/` of an earlier run are left as they are, with or without `-F`: they are no output of 2.0. Remove them, or run in another directory. |
| `scaffold_bin/` (`pipeline --scaffold`) | None: `--scaffold` is gone (see [Other tools](other-tools.md)). |
| `mges_bin_summary.tsv` (`mge`) | `mges_cluster_summary.tsv` |
| `<tmp>/recursive_clustering/recursive_bins/`, in the temporary directory (`validation`, `pipeline`; kept with `-N`) | `<tmp>/recursive_clustering/recursive_clusters/` |
| The columns of the contig data tables `Core_bin_ID`, `Core_bin_contigs`, `Core_bin_size`, `Overlapping_bin_ID`, `Overlapping_bin_contigs`, `Overlapping_bin_size`, `Recursive_bin_ID`, `Recursive_bin_contigs`, `Recursive_bin_size` and `Final_bin` | `Core_cluster_ID`, `Core_cluster_contigs`, `Core_cluster_size`, `Overlapping_cluster_ID`, `Overlapping_cluster_contigs`, `Overlapping_cluster_size`, `Recursive_cluster_ID`, `Recursive_cluster_contigs`, `Recursive_cluster_size` and `Final_cluster`. No step renames the columns of a table of 1.x: `partition` keeps them as extra columns, which `validation` carries into `contig_data_final.txt`; `validation` refuses a table without the columns of a partition (`Overlapping_cluster_ID`...), naming them; `mge` reads only `Name`, `Size`, `Hit` and `GC_content`, the same in both versions; `metator.host` still reads the tables of 1.x. Rename the columns of a table of 1.x first (`metator.io.rename_1x_columns` does). metacooler's `--contigs-format metator` names `Final_bin`: give `--contigs-format Name,Final_cluster` for a table of 2.0. |
| 1.4.0: a map of the clusters that could not be written, a warning, and exit status 0 | The run stops with the reason, logged as an error. Its other outputs are kept, and so is the map it would have replaced. |
| `plot/` (`validation`, `pipeline`): `bins_distribution.png`, `bins_size_distribution.png`, `MAGs_GC_distribution.png`, `MAGs_HiC_cov_distribution.png`, `MAGs_SG_cov_distribution.png` and `network_heatmap.png` | None: metaTOR draws no plot. The log of the validation counts the clusters in each quality class, and `metacooler show` draws the contacts between them from the map (see [Other tools](other-tools.md)). |
| 1.4.0: the parameters of the run in the metadata of its map (its `metator` key), `edge`, `scaffold`, `junctions` and `cluster_matrix` among them | No `edge`, `scaffold`, `junctions` or `cluster_matrix`, nor `forward`, `reverse`, `start`, `aligner`, `aligner_mode` or `min_quality`. `map`: the map that the run started from. `columns_from_map`: the columns of the map that the run read (`Restriction_site`, `Shotgun_coverage`). |

## What a run counts

| metaTOR 1.x | metaTOR 2.0 |
|---|---|
| At `--edge 0`, the default, without `--self-contacts` (`-s`), the network left out a pair with a read at position 0 or on the last base of its contig. | Such a pair is a contact, as in the contact map: on `tests_data/alignment.pairs`, 360 of 13,415 edges change. |
| A pair with a read past the end of its contig: left out of the network, but with `--self-contacts`, which tested no positions; its hits counted. | Left out, and so are its hits. |
| With `--self-contacts` (`-s`), an edge from each contig to itself, its contacts within itself | No such edge: `--self-contacts` is gone. |

See [the guide to the contact map](guides/metacool.md) for a run from a map
and the map a run writes.
