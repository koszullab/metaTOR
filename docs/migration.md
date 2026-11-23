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
| `metator network -1 R1 -2 R2 -a ASSEMBLY`, `metator partition`, `metator validation` | The same | |
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
| From pairs: a pairs file whose rows have different numbers of fields | The run stops with metacooler's reason, before its network. |

## What a run writes

| metaTOR 1.x | metaTOR 2.0 |
|---|---|
| | `<outdir>/metator.metacool` (`<prefix>.metacool`): the contact map of the contigs, grouped by cluster, written by every `pipeline` run from reads, alignments or pairs. |
| | A run from a map updates that map in place (staged, validated, then renamed over it), and locks it while it runs (`<map>.metator.lock`): `pipeline MAP`, and the steps `partition MAP` and `validation MAP`. |
| | `clusters.txt` and `cluster_members.txt`: the clusters of every clustering round, and their contigs, also stored in the `/clusters` group of the map. |
| | `mge --map`: the round `mge` of the clustering rounds of the map, its MGE clusters, after those of the run. |
| 1.4.0: a map of the clusters that could not be written, a warning, and exit status 0 | The run stops with the reason, logged as an error. Its other outputs are kept, and so is the map it would have replaced. |

## What a run counts

| metaTOR 1.x | metaTOR 2.0 |
|---|---|
| At `--edge 0`, the default, without `--self-contacts` (`-s`), the network left out a pair with a read at position 0 or on the last base of its contig. | Such a pair is a contact, as in the contact map: on `tests_data/alignment.pairs`, 360 of 13,415 edges change. |
| A pair with a read past the end of its contig: left out of the network, but with `--self-contacts`, which tested no positions; its hits counted. | Left out, with or without `--self-contacts`, hits included. |

See [the guide to the contact map](guides/metacool.md) for a run from a map
and the map a run writes.
