# Migrating to metaTOR 2.0

metaTOR 2.0 makes the `.metacool` contact map of a run its primary input and
output. The command of the original workflow is spelled as before: a run
from reads, alignments or pairs takes the same command line. This page lists
every other change to the command line since metaTOR 1.x (commit `4253551`),
with what replaced what was removed; the inputs that 1.x took and a run now
refuses; and what a run writes and counts differently. The
[changelog](changelog.md) has the details, the errors that changed, and the
other changes to the clusters and outputs (those of 1.4.0, for example:
Louvain and Leiden run through igraph).

## The command line

| metaTOR 1.x | metaTOR 2.0 | Why |
|---|---|---|
| `metator pipeline -1 R1 -2 R2 -a ASSEMBLY`, with `-S fastq`, `bam` or `pair` | The same | The original workflow is unchanged. |
| `metator network -1 R1 -2 R2 -a ASSEMBLY`, `metator partition`, `metator validation` | The same, but for the names of their outputs | A script that gives the overlapping clusters of `partition` to `validation -f` names `overlapping_clusters/`: see [What a run writes](#what-a-run-writes). |
| | `metator pipeline MAP -a ASSEMBLY`, `metator network MAP -a ASSEMBLY` | New: a run from a `.metacool` contact map of the contigs, the first argument (a file, or an entry `file.mmetacool::entries/<name>`), in place of `-1`, `-2` and `--start`, which a run from a map refuses. |
| `metator network -1 …`: `-1` required | `-1` needed only without `MAP` | A run takes a map, or reads, alignments or pairs. |
| `--edge N` (`-E`) of `network` and `pipeline` | Removed: refused, with the reason | The network counts every pair, as its contact map does, which holds no read positions. Even at `--edge 0`, the default, 1.x left some pairs out: see [What a run counts](#what-a-run-counts). |
| | `pipeline --output-map FILE` | New: writes the map of the clusters of the run to FILE, and leaves the input map of a run from a map as it was. |
| | `metator partition MAP`, `metator validation MAP` | New: the steps run one by one on a `.metacool` map, their first argument, which they update as `pipeline MAP` does: `partition MAP` stores its overlapping clusters in the clustering rounds of the map, and `validation MAP` writes the map of the clusters. Without `MAP`, they run as before. |
| | `partition --output-map FILE`, `validation --output-map FILE` | New, with `MAP`: write the map to FILE, and leave `MAP` as it was. |
| | `metator qc MAP`, `metator contactmap MAP` | New: the `.metacool` contact map of the run in place of the pairs files, for what a map gives: the noise of the library (`qc`, without the events of the pairs, which need their reads), and the contact map of the contigs of a cluster, as `metacooler extract` writes it (`contactmap`, a `.cool` file). Their options that need the read positions are refused then. |
| | `mge --map MAP`, `mge --output-map FILE`, `mge -F` | New: `mge` also stores its MGE clusters in the map of the run, as the round `mge` of its clustering rounds, in place or in FILE (replaced only with `-F`). It still reads the pairs files, and refuses a map in their place. |
| | `pipeline --no-metacool` | New: keeps no contact map (a run from reads still builds the map of its contigs, and removes it at the end). |
| | `metator metacool` | New command: writes the contact map of a run, from its pairs files or from the map it started from, and its final tables. |
| `metator host` | Removed: an unknown command | It always failed since 1.3.10: it called a function that the rewrite of the `host` module removed. |
| `-A spinglass` (`validation`, `pipeline`) | Refused: not a value of `--algorithm` | It never worked: its code was commented out. |
| `-p, --prefix` of `partition`, `validation` and `pipeline`: any characters | Letters, digits, `_` and `-` only (`qc -p` still takes any) | miComplete drops the dots from the names of the clusters. |
| Long options abbreviated (`--thread` for `--threads`) | Refused: give each option in full | The command line is rich-click's, where it was docopt's. |
| Command names in any case (`metator Pipeline …`) | Refused: an unknown command | The command line is rich-click's, where metaTOR capitalised the name it gave docopt. |
| `-T, --tmpdir DIR`: the temporary files in DIR, which the run removed at its end | In a private subdirectory of DIR, which the run removes, leaving DIR | 1.x removed the directory it was given, whatever else it held. With `-N` (`--no-clean-up`), the files kept are in `DIR/<random>/`. |
| A usage error exits with status 1 | Status 2 | An unknown value of `--algorithm`, `--aligner`, `--aligner-mode`, `--normalization`, `--start` or `--mat-fmt`, or a value that is not a number for a numeric option, is a usage error, refused before the command starts. |
| A `text_credentials.json` file in the working directory sent every error message to the web API it named | Removed | |

## What a run refuses

A run now refuses these inputs, which 1.x took, with the reason logged as an
error, before any work unless said otherwise:

| What | Why |
|---|---|
| From reads, alignments or pairs: an assembly whose contig names are not ASCII | The contact map of the run needs ASCII names. Rename the contigs. |
| A `--depth` file that lacks a contig of the assembly, or gives one another length | The file is read by contig name. 1.x read it line by line, in the order of the assembly, and gave the contigs of a file in another order the values of others. |
| An assembly whose FASTA file starts with blank or `#` lines, with Biopython 1.87 or later | Biopython refuses it: start the file with its first `>` header line. |
| From pairs: an empty pairs file, or a zip archive | The run builds its contact map from the pairs files, with metacooler, which reads them plain or compressed with gzip, bzip2, zstd or LZ4, whatever their name. |
| In the output directory of a run of 1.x, without `-F`: its outputs named after bins. `partition` refuses `overlapping_bin/`; `validation` refuses `recursive_bin/`, `final_bin/`, `final_bin_unscaffold/` and `bin_summary.txt`; `pipeline` refuses all of these and, with `--scaffold`, `scaffold_bin/` | The run would leave them next to its own outputs, which are named otherwise. `-F` removes them. |
| From pairs: a pairs file whose rows have different numbers of fields | The run stops with metacooler's reason, before its network. |

## What a run writes

| metaTOR 1.x | metaTOR 2.0 |
|---|---|
| | `<outdir>/metator.metacool` (`<prefix>.metacool`): the contact map of the contigs, grouped by cluster, written by every `pipeline` run from reads, alignments or pairs. |
| | A run from a map updates that map in place (staged, validated, then renamed over it), and locks it while it runs (`<map>.lock`): `pipeline MAP`, and the steps `partition MAP` and `validation MAP`. |
| | `clusters.txt` and `cluster_members.txt`: the clusters of every clustering round, and their contigs, also stored in the `/clusters` group of the map. |
| | `mge --map`: the round `mge` of the clustering rounds of the map, its MGE clusters, after those of the run. |
| `pipeline` from reads or alignments: `alignment_<i>_sorted.pairs.gz` and its index, sorted in the output directory, removed at the end of the run (kept with `-N`) | Written to the temporary directory of the run as they are made, `alignment_<i>.pairs`, and removed once the network is built; sorted, indexed and kept in the output directory with `-N` only (`--scaffold` sorts them too, and removes them at the end). Without `-N` or `--scaffold`, a run that stops before its network is built leaves them in its temporary directory; one that stops later (in the partition or the validation) leaves none, where 1.x left them, sorted, in the output directory: restart it from the contact map of its contigs, which it leaves (`metator pipeline <outdir>/metator.metacool … -F`). |
| `bin_summary.txt` | `cluster_summary.txt`. `qc -O` still finds the `bin_summary.txt` of a 1.x run. |
| `overlapping_bin/` | `overlapping_clusters/`: a script that gives it to `validation -f` names the new directory. |
| `recursive_bin/` | `recursive_clusters/` |
| `final_bin/` (`validation`), `final_bin_unscaffold/` (`pipeline`) | `final_clusters/` |
| `scaffold_bin/` (`pipeline --scaffold`) | `scaffolded_clusters/` |
| `mges_bin_summary.tsv` (`mge`) | `mges_cluster_summary.tsv` |
| `<tmp>/recursive_clustering/recursive_bins/`, in the temporary directory (`validation`, `pipeline`; kept with `-N`) | `<tmp>/recursive_clustering/recursive_clusters/` |
| The columns of the contig data tables `Core_bin_ID`, `Core_bin_contigs`, `Core_bin_size`, `Overlapping_bin_ID`, `Overlapping_bin_contigs`, `Overlapping_bin_size`, `Recursive_bin_ID`, `Recursive_bin_contigs`, `Recursive_bin_size` and `Final_bin` | `Core_cluster_ID`, `Core_cluster_contigs`, `Core_cluster_size`, `Overlapping_cluster_ID`, `Overlapping_cluster_contigs`, `Overlapping_cluster_size`, `Recursive_cluster_ID`, `Recursive_cluster_contigs`, `Recursive_cluster_size` and `Final_cluster`. Every command still reads the tables of 1.x, their columns under the new names. The Python functions that take a table of contig data (`metator.figures.plot_figures`, `metator.validation.correct_final_bin`...) need the new names: read a table of 1.x with `metator.io.read_contig_data`, or give it to `metator.io.rename_1x_columns`. metacooler's `--contigs-format metator` names `Final_bin`: give `--contigs-format Name,Final_cluster` for a table of 2.0. |
| 1.4.0: a map of the clusters that could not be written, a warning, and exit status 0 | The run stops with the reason, logged as an error. Its other outputs are kept, and so is the map it would have replaced. |
| `plot/network_heatmap.png`: drawn from the pairs files, at 2000 dpi whatever its size; a heatmap that could not be drawn stopped the run | Drawn from the contact map of the run, whatever it started from: the contacts of a bin of the map are spread over the 50 kb windows that it covers, as a map holds no read positions. At about one pixel per window. A heatmap that cannot be drawn is a warning. Without a map, `validation` draws it from the pairs files of its output directory, as before. |
| 1.4.0: the parameters of the run in the metadata of its map (its `metator` key), `edge` among them | No `edge`. `map`: the map that the run started from (none for a run from reads, alignments or pairs), with `start` none for `pipeline MAP` and `network MAP` (`partition` and `validation` have no `--start`). `columns_from_map`: the columns of a map that a run from it read (`Restriction_site`, `Shotgun_coverage`). |

## What a run counts

| metaTOR 1.x | metaTOR 2.0 |
|---|---|
| At `--edge 0`, the default, without `--self-contacts` (`-s`), the network left out a pair with a read at position 0 or on the last base of its contig. | Such a pair is a contact, as in the contact map: on `tests_data/alignment.pairs`, 360 of 13,415 edges change. |
| A pair with a read past the end of its contig: left out of the network, but with `--self-contacts`, which tested no positions; its hits counted. | Left out, with or without `--self-contacts`, hits included. |

See [the guide to the contact map](guides/metacool.md) for a run from a map
and the map a run writes.
