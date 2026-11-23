# Change Log

All notable changes to this project will be documented in this file.

## Unreleased

metaTOR 2.0 clusters contigs, and nothing else. It starts from a `.metacool`
contact map of the contigs, which hicstuff makes from the reads and metacooler
from pairs files, and every step reads that map and stores what it finds in it:
the map is what the steps share, and it becomes the contact map of the
clusters. The commands that did what hicstuff, metacooler, instaGRAAL or pairix
do are gone. [Migrating to metaTOR 2.0](migration.md) lists every change to the
command line and the outputs; [Other tools](other-tools.md) says what replaces
each removed command.

### Breaking changes
- **A contact map is the only input.** `network`, `partition`, `validation`,
  `pipeline` and `mge` take a `.metacool` contact map of the contigs, `MAP`,
  as their first argument, which they require: a file, or an entry of a
  `.mmetacool` collection (`file.mmetacool::entries/<name>`), at any bin size.
  hicstuff makes it from the reads (`hicstuff pipeline -M metacool`, from
  hicstuff 4.0.0), metacooler from pairs files (`metacooler create
  pairs`) or a `.cool` file (`create cool`), and `metacooler merge` sums the
  maps of several libraries. Gone with the reads, alignments and pairs
  inputs: `-1/--forward`, `-2/--reverse`, `-S/--start`, `-b/--aligner`,
  `-B/--aligner-mode`, `-q/--min-quality` and `-E/--edge`; the alignments,
  pairs files and per-library networks a run wrote (`alignment_<i>*.bam`,
  `alignment_<i>_sorted.pairs.gz`, `network_<i>.txt`,
  `hit_data_alignment.txt`); `--assembly` as the basename of an aligner's
  index (it is the FASTA file of the assembly, plain or compressed with
  gzip); and the dependencies on hicstuff and pysam, and on bowtie2, bwa and
  samtools. A removed option is an unknown option, a usage error (exit
  status 2).
- **The map is the interface between the steps.**
  - `network MAP` stores the network in the map: what it measured on each
    contig in the contigs table (`GC_content`, `Hit`, and `Shotgun_coverage`
    and `Restriction_site` when it has them), the normalization of the
    network as the weights `metator_<normalization>` of the bins (none for
    `None`), and its run under the `network` entry of the `metator` key of
    the metadata. It writes no `network.txt`, and writes
    `contig_data_network.txt` in its temporary directory, which `-N` keeps,
    not in the output directory; it loses `-t/--threads`: the network is in
    memory, and the external `sort` is no longer run.
  - `partition MAP` reads the network from the map, and stores its
    partition there: every core cluster (`core_00001`...) and every
    overlapping cluster, whatever its size, with their contigs, as the
    clustering rounds `core` and `overlapping`, and its run under the
    `partition` entry, with the identity of the network it partitioned
    (`network_digest`). It writes no `contig_data_partition.txt` and no
    `overlapping_bin/`, and loses `-a/--assembly`, `-c/--contigs`,
    `-n/--network` and `-s/--size`.
  - `validation MAP` reads the network and the partition from the map, and
    refuses a partition made on another network than the one stored there.
    Its `--size` chooses the overlapping clusters that miComplete checks. It
    writes the FASTA files of the final clusters in the output directory,
    `final_clusters/`, and stores the rest in the map: the final cluster of
    each contig, as its genome, the summary of the clusters in the genomes
    table, every clustering round, and its run under the `validation` entry,
    with how the summary is written from the map (`summary_rows`,
    `summary_columns`, `summary_formats`). It loses `-c/--contigs`,
    `-f/--fasta` and `-n/--network`.
  - The tables of the validation (`contig_data_final.txt`,
    `cluster_summary.txt`, `binning.txt`, `clusters.txt`,
    `cluster_members.txt` and miComplete's results) and `recursive_clusters/`
    are written in the temporary directory of the run, which
    `-N/--no-clean-up` keeps and the log names: `metator dump` writes them
    from the map (see Added).
  - **The results of an earlier run in the output directory are refused,
    with or without `-F`, and never removed**: the tables and recursive
    clusters of a development build of 2.0, and the outputs of 1.x (its
    tables and networks, `bin_summary.txt`, `network.txt`, the clustering
    matrices, `recursive_bin/`, `final_bin/`, `final_bin_unscaffold/`,
    `overlapping_bin/`, `scaffold_bin/`, `plot/`). Left there, they would pass
    for the results of the run; a file of one of these names, a user's own
    too, is refused and left as it is. `-F` replaces only what a run writes
    itself: `final_clusters/` and the map of the output directory.
  - `pipeline MAP` runs the three steps on the map, and writes the map of
    the clusters once, at its end: the steps run one by one in one output
    directory give the map of `pipeline`. It goes to the output directory,
    under the name of `MAP` (`<outdir>/sample.metacool` from
    `sample.metacool`), where 1.4.0 wrote `<outdir>/metator.metacool`.
- **The map of a run goes to its output directory, and `MAP` is left as it
  was.** A run writes its map to its output directory, under the name of
  `MAP` (`NAME.metacool` for an entry of a collection): metacooler copies
  `MAP` there, the run edits the copy, and metacooler renames it into place.
  `--inplace`, on every command that writes the map, updates `MAP` itself
  instead, and so does a run whose output directory holds `MAP` already: the
  steps of a run, one by one, in one output directory (`network MAP -o DIR`,
  then `partition DIR/MAP -o DIR`...). A run writes the map through
  metacooler, in an edit session: every edit of a step goes into one copy of
  the input map, staged next to the map that the run writes
  (`.<map>.<uuid>.<n>.d/`), validated in full, then renamed onto it, so a
  run that stops for any reason leaves the map as it was, and a map updated
  run after run does not grow. A run holds metacooler's lock of the map it
  writes (`<map>.lock`) from its checks to its end. A map of an earlier run
  in the output directory is replaced only with `-F`. Refused before any
  work: a map that another process holds, a map updated in place without
  write permission, a directory that is not writable, a hard link of `MAP`
  in the output directory, an entry of a collection with `--inplace`. A map
  updated in place through a symbolic link is updated where the link points;
  a symbolic link in the output directory is replaced, not its target. A map that changed since the
  run read it is refused when the run writes its map, and a read of a map
  that another run replaces is refused. Neither metaTOR nor metacooler
  handles signals: a run killed by SIGTERM (a job scheduler's time limit) or
  SIGKILL leaves its lock and its staged copy, which the next writer on the
  same host takes over and removes; on another host, delete the files that
  the error names. miComplete runs in the process group of the run, so a
  Ctrl-C, or a scheduler's signal to the job, reaches it. A map of the
  clusters that cannot be written stops the run with the reason, where 1.4.0
  warned and exited 0.
- **Clusters, not bins.** The outputs and columns of a run are named after
  clusters: `bin_summary.txt` is `cluster_summary.txt`, `final_bin/`
  (`validation`) and `final_bin_unscaffold/` (`pipeline`) are
  `final_clusters/`, `recursive_bin/` is `recursive_clusters/`, and in the
  contig data tables `Core_bin_ID`, `Core_bin_contigs`, `Core_bin_size`,
  their overlapping and recursive equivalents and `Final_bin` are
  `Core_cluster_ID`... and `Final_cluster`. `binning.txt` and the names of
  the clusters are kept. `metator.host` still reads the tables of 1.x
  (`metator.io.rename_1x_columns`).
- **The contigs of a run are those of its map, in its order.** The ID of a
  contig is its position in the map, and the rows of the tables, the numbers
  of the clusters and the records of their FASTA files follow that order
  (the longest first, as hicstuff and metacooler make a map), where they
  followed the assembly. The contigs of the assembly that the map lacks are
  in no table, and the shotgun abundance of a cluster is its share of the
  contigs of the map.
- **The network.** It counts every pair, as the map does: a pair with a read
  on the last base of a contig, or at position 0, is a contact, which
  `--edge 0` left out (on `tests_data/alignment.pairs`, 360 of 13,415 edges
  change), and a pair with a read past the end of its contig is none. Each
  edge weighs the contacts of two contigs times the weights of each,
  `n × w_i × w_j`, which differ from 1.x's `n / √(h_i × h_j)` in their last
  bits; the edges of weight 0 are left out.
- **`mge` reads the map.** `metator mge MAP -m MGES -a ASSEMBLY` takes the
  contacts between the MGE contigs from the map, and those within each
  contig, between reads more than 1 kb apart, from its contact statistics
  (metacooler's QC block, counted from the pairs at 1001 bp, with metacooler
  0.4.0 or later); the size, GC content and hits of each contig are those
  that `metator network` stored. A map without those statistics, or without
  the network, is refused before anything is written, and an MGE contig that
  the map lacks is left out, with a warning. Gone: its pairs files, `-c` and
  `-b`, and the dependencies on pypairix and pyfastx. It takes
  the MGE contigs the longest first, ties by name, whatever the order of the
  map, which numbers some clusters otherwise and may break a tie of scores
  otherwise. It stores the MGE clusters in the map, as the round `mge` (and
  `--random`, the round `mge_random`), where it wrote
  `mges_bin_summary.tsv`, which is refused in the output directory, with or
  without `-F`, and never removed; `mges_binned.fa` is written by metacooler, its contigs joined by
  200 `N`, in lines of 80 bases.
- **Removed commands and features**, with what replaces them in
  [Other tools](other-tools.md): `metator qc`, `contactmap`, `metacool`,
  `scaffold` and `pairs`; `pipeline -v/--scaffold` and `-J/--junctions`; the
  plots of `validation` and `pipeline` (`plot/`); `-C/--cluster-matrix` of
  `partition`, `validation` and `pipeline`, and the clustering matrices it
  wrote; `-s/--self-contacts` of `network`; `pipeline --no-metacool` of
  1.4.0, as the map of the clusters is what a run writes (to its output
  directory, or over `MAP` with `--inplace`); and the module `metator.main` (`python -m metator.main`):
  run `metator`, or `python -m metator`.
- **Seeded runs.** The temporary directory of a command is named from the
  random numbers of the system, not from Python's `random` module, which the
  partition draws from: a run seeded as before gives other clusters.
- **Python API.** Removed: the modules `metator.align`, `metator.contact_map`,
  `metator.figures` (`bin_quality` is now `metator.validation.bin_quality`),
  `metator.main`, `metator.quality_check`, `metator.regions` and
  `metator.scaffold`; the readers and writers of the files 2.0 no longer reads
  or writes (`metator.io.sort_pairs`, `write_contigs`, `import_anvio_binning`,
  `read_contig_data`...; `metator.network.normalize_pair`, `write_edge`,
  `alignment_to_contacts`...; `metator.partition.read_network`,
  `generate_fasta`, `update_contigs_data`...;
  `metator.validation.update_contigs_data_recursive`, `write_bins_contigs`...;
  `metator.mge.build_matrix`, `generate_mges_fasta`...); the map writers of
  the development builds (`reassign_metacool`, `store_network`,
  `write_rounds`, `write_round`), now `metator.metacool.write_map`; and the
  unused CheckM functions and `metator.partition.spinglass_partition`.
  `metator.partition.partition` returns the partition, a
  `metator.clusters.Partition`, and writes nothing;
  `metator.validation.recursive_decontamination` takes the partition and the
  map; `metator.mge.mge_binning` takes the map.

### Added
- `metator dump MAP -o DIR` writes the tables of a run from the map of its
  clusters, read only: `contig_data_final.txt` (without its `ID` column),
  `cluster_summary.txt`, `binning.txt`, `clusters.txt` and
  `cluster_members.txt` (every clustering round, `mge` included) and
  `network.txt` (by contig names, where 1.x wrote IDs). The values are those
  of the tables that the run kept in its temporary directory; the rows of
  `contig_data_final.txt`, `binning.txt` and `cluster_members.txt` follow
  the map of the clusters. A map without the clusters of a validation, whose
  clusters were edited after the run (by `metacooler reassign`), whose network
  is not the one its clusters were made from, or that changes while it is
  read, is refused, and the tables are replaced all at once.
- The clustering rounds of a run: every cluster of every round (`core`,
  `overlapping`, `recursive_<k>`, `final`, `mge`, `mge_random`), with its
  parent and miComplete's columns, and the contigs of each, in the
  `/clusters` group of the map, which `clusters.txt` and
  `cluster_members.txt` write out.
- The normalization of the network is stored in the map, as the weights
  `metator_<normalization>` of its bins, the way cooler stores balancing
  weights: `cooler` (`matrix(balance="metator_empirical_hit")`) and
  `metacooler show`, `dump` and `network` (`--balance`) read them. They record
  the columns they are computed from, so that metacooler drops them when an
  edit changes one.
- A run reads the restriction sites and the shotgun coverage of the contigs
  from the `Restriction_site` and `Shotgun_coverage` columns of the map when
  `--enzyme` and `--depth` are not given; given, they replace them, with a
  warning that counts the contigs whose value differs.
- The metadata of the map record, under the `metator` key, one entry per
  step (`network`, `partition`, `validation`, `mge`): the version of
  metaTOR, the command line and parameters of the run, and what it read of
  the map (`read`).
- `partition` and `validation` refuse a map whose network an edit changed
  since `network` stored it, with the step to run again: a column, a value or
  the weights dropped, or contacts whose hits are no longer its `Hit` column
  (as after `metacooler merge`, which keeps the contigs table of its first
  map).
- `--inplace` on every command that writes the map (`network`, `partition`,
  `validation`, `pipeline`, `mge`): the run updates `MAP` in place, rather
  than write its map to the output directory; and `-F` on `network` and
  `mge`, to replace the map of an earlier run in the output directory.
  `partition -F`, which replaced the outputs of an earlier partition (it
  writes none now), replaces the map of the output directory.
- Python API: `metator.dump`; `metator.metacool.write_map`, the one writer of
  the maps of the steps, and `map_rounds`, `read_partition`; `metator.network`
  `map_network`, `stored_network`, `contig_contacts`, `contig_hits`,
  `normalization_weights` and `relabel_network`; `metator.clusters`
  (`Partition`, `Rounds`, `write_round_tables`); `metator.mge.check_mge_map`,
  `map_matrix` and `write_mges_fasta`; `metator.validation.contaminated_bins`;
  `metator.partition.iterate`, and the `progress` argument of
  `algo_partition`, `louvain_iterations` and `leiden_iterations`, which logs
  the iterations as they go; `metator.log.duration`.

### Changed
- Validation runs Prodigal and HMMER through pyrodigal and pyhmmer, their
  Python bindings, which pip installs with metaTOR, in place of the
  prodigal and hmmsearch programs that miComplete runs
  (`metator._micomplete`): metaTOR runs no program that pip does not
  install, with pip, with conda (`metator.yaml`) or in the Docker image.
  miComplete gives the clusters of the tests the rows that prodigal 2.6.3
  and HMMER 3.4 gave them, and a cluster of less than 20 kb still gets no
  gene, as prodigal finds none in it.
- metaTOR requires metacooler 0.4.0: its contact statistics of each contig
  (`qc_intra_contig_short`, `qc_intra_contig_long`, split at 1001 bp), its
  clustering rounds, its locked edit sessions and its FASTA writer.
- The FASTA files of the clusters are written by metacooler
  (`metacooler.write_fasta`), the records as they are in the assembly; a
  gzip-compressed assembly is decompressed once per run, in its temporary
  directory, which needs room for it. No pyfastx index is made.
- metaTOR reads back the numbers of the tables it wrote exactly, under
  pandas 2 as under pandas 3: the map stores the run's numbers.
- The partition computes the distances between the core clusters in
  `--threads` threads, where it forked as many processes: Python warned that a
  fork of the run, which runs threads, may deadlock. SciPy computes the
  distances without holding the GIL, so the threads run in parallel.
- Python's warnings (cooler gives the problems of a map with them) are shown
  and written to the log of a run, as metaTOR's own messages.
- **The log of a run.** Its header gives, as hicstuff's does, the version,
  the date, the command line (`## cmd:`, where it said `## Command used:`)
  and the value of each option (`## <option>: <value>`), then a line of
  dashes, without the empty lines around it; the log ends with what the run
  wrote and how long it took (`Final contact map written after 0h 8m 31s:
  <map>`), where it ended with a line of dashes and the date. `pipeline`
  marks its three steps (`[1/3] network`, `[2/3] partition`,
  `[3/3] validation`), and no longer logs its parameters, which the header
  gives. The partition counts its iterations as they go, at each tenth of
  them, 10 seconds apart at least, and how long they took; the validation
  says how long miComplete took and how many clusters are contaminated, then
  logs each recursive step (`Recursive step 2/9: ...`), how many clusters it
  partitions again and replaces, and when no further step is needed, where
  it printed each step in a banner of dashes. The counts of the contacts, of
  the core and overlapping clusters, and of the clusters checked are one
  line each. The bookkeeping of a run (the tables written in its temporary
  directory, the temporary directory removed, the staged copy of the map,
  and metacooler's accounts of what it reads and writes) is logged at the
  DEBUG level, which the log does not show.
- The documentation moves from Read the Docs to GitHub Pages, at
  <https://koszullab.github.io/metator/>: Markdown pages built with Zensical,
  a guide per step and command, the CLI reference (which `tests/test_docs.py`
  checks against the CLI), the file formats and the Python API reference.
  The Sphinx sources and the scripts of the curation tutorial
  (`docs/example/scripts/`, which ran CheckM) are gone: the curation guide
  checks the curated clusters with miComplete.
- The docstrings of the Python API are shorter: a summary, then numpy
  sections (Parameters, Returns, Raises), which the API reference renders as
  tables.

### Fixed
- The depth file (`--depth`) is read by contig name: read line by line, a
  file in another order than the assembly gave the contigs the values of
  others. A file that lacks a contig, or gives it another length, is refused.
- `network` and `pipeline` read an assembly compressed with gzip, which they
  refused as no FASTA file.
- `validation` reads the contig names of its tables as text: contigs named
  by numbers (`000017`, `1e5`), or `NA`, kept their clusters, where they lost
  them or the map was not written. A clustered contig that the map lacks, and
  final clusters that are not those of their FASTA files, are refused.
- `binning.txt` holds every contig of the final clusters: those of a cluster
  partitioned again but kept whole were left out.
- The log of the validation counts the clusters in the quality classes of
  `metator.validation.bin_quality`, the contamination being the weighted
  redundancy − 1; it divided it by the completeness.
- The `recursive` column of `cluster_summary.txt` says whether the cluster
  comes from the recursive partition.
- `metator validation --help` gives the rule of the recursive partition: a
  weighted completeness of 33%, not 50%.
- `pipeline` and `validation` check their arguments and their assembly before
  `-F` replaces the clusters of an earlier run; `-F` replaces a directory of
  clusters that is a symbolic link, not its target; a `--tmpdir` in
  `final_clusters/`, which `-F` replaces, or in a directory of the results of
  an earlier run, which the run refuses, is refused as a usage error.
- `pipeline` and `validation` stop with a clear error when no cluster of the
  partition reaches `--size`, where miComplete failed.
- An output path that holds `::`, which metacooler reads as `file::group`, is
  refused; a map given as `file::/` is the file.
- The validation leaves no pyfastx index (`.fxi`) next to the assembly or the
  FASTA files of the clusters.
- The FASTA files of the recursive and final clusters (`validation`,
  `pipeline`) and of the MGE clusters (`mge`) were written by the `pyfastx`
  command through the shell (`pyfastx extract … >`): without `pyfastx` on the
  `PATH` (a virtual environment used without activating it), every file came
  out empty, without an error. metacooler writes them now.

## [1.4.0] - 2026-11-17

### Added
- `pipeline` writes the contact map of all the contigs, grouped by bin, as
  `metator.metacool` in its output directory (`<prefix>.metacool` with
  `--prefix`), once the bins are final; `--no-metacool` skips it. A
  `.metacool` file is the contact map of
  [metacooler](https://github.com/js2264/metacooler): a cooler whose
  chromosomes are the contigs, which `cooler` reads as any `.cool` file and
  `metacooler` by genome. It has one bin per contig and one genome per final
  bin, plus `__unbinned__` for the contigs in no bin, which keep all their
  contacts (all the contigs of a run with no final bin). Each pair whose two
  ends lie on contigs of the assembly is one contact, and the libraries are
  summed. The genomes table holds the columns of `bin_summary.txt`; the
  contigs table the `GC_content`, `Hit`, `Shotgun_coverage` and
  `Restriction_site` columns of `contig_data_final.txt`, except those that
  are `-` on every contig (`Shotgun_coverage` without `--depth`,
  `Restriction_site` without `--enzyme`), and never the bins; and the file
  metadata the version, command line and parameters of the run, under the
  `metator` key. The map is written last: when it cannot be written (for
  example, from a misaligned bin summary, see below), a warning says why
  and the run keeps its other outputs and its exit status. A map left in
  the output directory by an earlier run is removed first, so that it
  cannot pass for the map of this run. See the README.
- `metator metacool` writes the same file from the pairs files,
  `contig_data_final.txt` and `bin_summary.txt` of a run. Python API:
  `metator.metacool.write_metacool`. It refuses bins named as metacooler
  marks a contig in no bin (`NA`, `NaN`, `N/A`, `N/D`, `-1` or
  `__unbinned__` in any case, and `nd` in any case other than `ND`, which
  marks a contig in no bin), whose contigs metacooler would put in
  `__unbinned__` without a word; metaTOR's own bin names never are. When
  it cannot read the contig data (for example, a `Size` that is not a
  whole number) or metacooler cannot write the map, the command logs the
  reason as an error, then raises it.
- `bin_summary.txt` goes to metacooler as it is, its first column naming
  the bins. Its rows must have as many fields as its header, as metaTOR
  writes them: metacooler refuses a misaligned table, such as
  `tests_data/outdir/bin_summary.txt` (18 header fields, 15 on each row:
  its header splits three column names at their space), whose values used
  to be stored under the wrong column names. The pipeline then writes no
  map, with a warning saying why, and `metator metacool` stops with that
  reason.
- metaTOR depends on `metacooler` >= 0.2.0.

### Changed
- metaTOR supports Python 3.10 to 3.14: 3.13 and 3.14 are new, and
  3.9 is no longer supported (see Removed). The Docker image runs
  Python 3.13 instead of 3.11, and takes miComplete from PyPI instead of
  bioconda, whose package requires Python < 3.12.
- `Biopython` is no longer pinned <= 1.80: miComplete is now run through `metator._micomplete`, which restores `Bio.SeqUtils.GC` (removed in Biopython 1.82) and `distutils.spawn` (removed in Python 3.12) before running it.
- With Biopython >= 1.87, FASTA inputs must start with a `>` header line:
  leading blank or `#` lines are rejected. `network` and `pipeline` then
  refuse the assembly, and `qc` stops with an error naming it.
- Louvain and Leiden now run through `igraph`: the bundled Louvain binaries,
  the Leiden jar and the Java requirement are gone, and the wheel is pure
  Python.
- `metator mge`: flags renamed for consistency with the other commands;
  MetaBAT input and the host, CheckV and plotting steps removed (#33).
- `--prefix` (partition, validation, pipeline) accepts only letters, digits,
  `_` and `-`: miComplete drops dots from bin names.
- `validation` and `pipeline` no longer accept `-A spinglass`.
- Temporary files go to a private subdirectory of `--tmpdir`.
- Declare `matplotlib`, and drop `checkv` and `cooler`, in the package
  dependencies.
- metaTOR no longer depends on `scikit-learn`, `looseversion`,
  `packaging` or `requests`.
- Figures are drawn with matplotlib alone: metaTOR no longer depends on
  `seaborn`, whose boxplots pass an argument that matplotlib 3.13
  removes, and requires `matplotlib` >= 3.9. The MAG GC and coverage
  boxplots and the background-noise violins show the same data; their
  colours are no longer muted, and each half violin now stops at the
  smallest and largest values instead of extending past them.
- `networkx` is no longer installed, and `igraph` must be 0.10.3 or later:
  `scaffold` and the MGE host detection build their graphs with `igraph`,
  and the recursive step of `validation` takes the edges of each bin
  straight from `network.txt`, every row as `partition` does (networkx
  kept the last row of an edge listed several times, as in a
  `network.txt` written by 1.3.10 from several libraries). Otherwise the
  outputs do not change, with one exception: `scaffold` (and
  `pipeline --scaffold`) now reads a linear scaffold from the end whose
  contig comes first in the bin fasta, and a circular one from its first
  contig. networkx chose the end from the order of a Python set for
  scaffolds holding less than half of the bin's contigs: such a scaffold
  can come out reverse-complemented or starting from another contig, with
  the same contigs, junctions and relative orientations. Python API:
  `metator.host.build_contig_graph` returns an `igraph.Graph`, and the
  `network` of `metator.validation.recursive_clustering` and
  `recursive_clustering_worker` is the network table (a
  `pandas.DataFrame` with the columns contig1, contig2 and weight).
- The command line uses `rich-click` instead of `docopt`: the help pages
  look different, with the options grouped in panels, but the commands,
  their options, short aliases and defaults are unchanged. An unknown
  value of `--algorithm`, `--aligner`, `--aligner-mode`,
  `--normalization`, `--start` or `--mat-fmt`, or a value that is not a
  number for a numeric option, is refused before the command starts,
  with exit status 2 (`network --start network` included: it was
  accepted, then failed while aligning). The other usage errors (a
  missing option or argument, an unknown command or option) also exit
  with status 2 instead of 1. Long options can no longer be
  abbreviated (`--thread` for `--threads`). The log files record the
  command line as typed. `python -m metator` runs the command line.
- Pairs files are sorted by hicstuff (`hicstuff.io.sort_pairs`, so
  hicstuff >= 3.3.0 is required): in `metator pairs`, after the
  alignment in `network` and `pipeline`, and when `scaffold`, `qc` or
  `contactmap` read a pairs file that has no pypairix index. The sorted
  files keep their names, their bgzip compression and their pypairix
  index, and `metator pairs` its options. The pairs keep the order of
  pairtools sort, contig names byte by byte (`chr1`, `chr10`, `chr2`)
  then positions, so their blocks come in the same order as with 1.3.10
  and the pairs of contigs whose names differ only by leading zeros
  (`contig_01` and `contig_1`) are still indexed like any others; only
  pairs with the same two contigs and the same two positions can come
  in another order among themselves. The header says
  `#sorted: chr1-chr2-pos1-pos2` instead of keeping the `#sorted: readID`
  of the input. Python API: `metator.io.sort_pairs_pairtools` is now
  `metator.io.sort_and_index_pairs`, whose `tmp_dir` defaults to the
  temporary directory of `sort` (`$TMPDIR`, else `/tmp`) instead of the
  working directory.
- `pairtools` is no longer installed: neither metaTOR nor hicstuff >= 3.3.0
  uses it.

### Removed
- Python 3.9 support, a breaking change for Python 3.9 users: metaTOR
  now requires Python 3.10 or later, and pip or uv on Python 3.9 keep
  installing 1.3.10.
- `metator host`, which always failed since 1.3.10: it called
  `metator.host.host_detection`, which the rewrite of the `host` module
  removed. Python API: `metator.io.import_network`, used only by it.
- Text notifications: when the working directory held a
  `text_credentials.json` file, every error message was sent to the web
  API named in it. Python API: `metator.log.setup_text_logging`,
  `TEXT_CREDENTIALS_DEFAULT_PATH` and `text_formatter` are gone.
- Python API: the CheckV plot helpers of `metator.figures`
  (`barplot_bins_number`, `barplot_bins_size`, `build_vmags_summary` and
  `pie_bins_size_distribution`), unused since `mge` lost its CheckV step.
- The `docopt` dependency.
- Python API: `metator.commands`, the docopt command classes. Each
  command is a function of `metator.cli`, whose `cli` is the `metator`
  command line; `metator.main.main` is an alias of it. Like any click
  command, it ends by raising `SystemExit` (status 0 after a successful
  run) instead of returning.

### Fixed
- A fresh install of 1.3.10 did not start: it never declared `docopt`,
  which `hicstuff` stopped providing in 3.2.5. The command line no longer
  needs it.
- `network`, `partition`, `validation` and `pipeline` deleted the whole
  `--tmpdir` directory, including files that were already in it.
- `pipeline` (fastq or bam start) crashed while plotting at the end, after
  writing its bins: the network heatmap read the compressed pairs files as
  text. Runs with `--depth` crashed at the same step on a missing column.
- `pipeline` (fastq or bam start, without `--no-clean-up`) deleted the
  pairs files it had made from the reads but left their pypairix index
  (`alignment_N_sorted.pairs.gz.px2`) in the output directory, as 1.3.10
  did. The index now goes with its pairs file; pairs files given with
  `--start pair` are still kept, with their index.
- `pipeline --scaffold` left the messages of the processes scaffolding the
  bins (`Scaffolding <bin>...`) out of its log file wherever these
  processes are not forked: on macOS, and on Linux from Python 3.14, whose
  default start method is forkserver. They only reached the terminal.
- With pandas 3, when a contaminated bin was split into several bins, only
  one of them was kept in the final bins.
- With `--aligner bwa`, every pair carried the first mate's position for
  both mates.
- With `--aligner bowtie2` (the default), `network` and `pipeline` paired
  the forward and reverse reads assuming an order of the read names other
  than the one `samtools sort -n` sorts them in: the pair `r5` was lost
  when the reverse reads were `r.5` (unpaired) then `r5`, and names such
  as `r_1` and `r5` stopped the run with `TypeError`. The reads are now
  paired in the samtools order.
- `--aligner-mode iterative` (`network`, `pipeline`, and `cutsite` mode
  without `--enzyme`, which falls back to it) only ever aligned the first
  40 bp of each read: it told hicstuff that the reads were 20 bp long, so
  after a round at 20 bp (which aligns nothing with bowtie2's local mode)
  the reads were cut to 40 bp and never extended. The reads that do not
  align are now extended by 20 bp at a time up to the length of the
  longest read of the file (not of the first one, hicstuff's default,
  which would still cut the reads of a trimmed file whose first read is
  short). A read is kept at the first length at which it aligns with a
  mapping quality of at least `--min-quality` (30 by default; the Python
  API's `metator.align.align` takes it as `min_qual`): hicstuff was given
  40, whatever `--min-quality`. Networks built in iterative mode change: on
  the test data, 56,810 pairs instead of 33,358 (41,563 in normal mode).
- `--aligner bwa` in iterative mode crashed once the reads were aligned,
  with `UnboundLocalError`. The forward and reverse reads, which this mode
  aligns separately, are now paired as with bowtie2, keeping one alignment
  per read: hicstuff writes a read once per round when bwa also gives it a
  supplementary alignment.
- With several Hi-C libraries, `network.txt` listed each edge once per
  library; steps reading it with networkx kept only the last library's
  contacts.
- Two runs started from the same directory deleted each other's miComplete
  files.
- `contactmap` selected no contig for a numeric bin id (`--name 8`), and
  failed without `--outdir`.
- `contactmap --mat-fmt cool` wrote no `.cool` file: hicstuff 3.1.4 to
  3.2.5 swap the chr2 and pos2 columns in `pairs2cool`. hicstuff >= 3.3.0,
  now required, does not.
- `mge --random` always failed.
- A bowtie2 or bwa index given as `../path` was looked up at `/path`, and
  so, in iterative mode, was an index built from the fasta in a `--tmpdir`
  given as `../path`. In iterative mode, hicstuff looks the index up again
  and strips the leading dots of its path: it is now given the absolute
  path of the index.
- pandas 3 compatibility.
- Pairs sorting and `bgzip` failures are reported where they happen.
- `pairs` failed on a pairs file without a `#columns` header line, where
  pairtools >= 1.1 looked up the chr1, chr2, pos1 and pos2 columns, or
  without any header line, which pairtools refused; so did `scaffold`,
  `qc` and `contactmap` when they had to sort such a file.
- Core bin detection is no longer quadratic in the number of core bins.
- `scaffold` and `pipeline --scaffold` never scored a junction between the
  starts of two contigs that follow each other in the bin fasta, so such
  contigs could come out in the wrong order or orientation.
- `scaffold` and `pipeline --scaffold` could reverse the first or last
  contig of a linear scaffold (or close the scaffold into a circle): the
  outer end of that contig, which has no partner, took the contig next to
  its other end.
- Python API: `metator.regions.Bin.insert`, `invert` and `delete` raised
  `AttributeError` after they had changed the contigs; the length of a
  `Bin` now follows these edits.
- Python API: when the contigs of a `metator.regions.Bin` came back to the
  layout they had when its boundaries were first read (a contig inserted,
  then deleted), the bin kept the boundaries of the layout in between: a
  wrong length, and wrong or failing contig lookups.
- Python API: `metator.regions.Bin.get_breakpoints` raised `ValueError` for
  any bin of two contigs or more.
- Python API: `metator.regions.Bin.get_contigs_between` named every piece
  after the bin instead of its contig, so the pieces could not be read from
  the fasta.
- `validation` failed at the very end, after writing its bins, when no bin
  was contaminated: it removed `<assembly>.fxi`, which pyfastx only writes
  when a bin is re-partitioned.
- `partition` failed the same way at the end when no bin reached `--size`.
- `HiC_abundance` in `bin_summary.txt` (`validation`, `pipeline`) was half
  the share of the Hi-C reads in the bin: the values are now twice as
  large, and those of a fully binned metagenome sum to 100%, like
  `SG_abundance`.
- `mge` with a pairs file that has no pypairix index: the contacts of a
  pair whose first read lies on the later of the two MGE contigs (in the
  contig table order) were left out of the MGE contact matrix. The matrix,
  and so the MGE bins, now count them, as with an indexed pairs file.

## [1.3.10] - 2025-05-21
- `host` module rewritten to annotate the bacterial hosts of each MGE (#27).
  Known issue: the `metator host` command was not updated and fails.
- `pairtools sort` writes its temporary files to the temporary directory.
- Docker image: install `procps`.

## [1.3.9] - 2025-03-18
- Pin `biopython<=1.80`, required by `micomplete` 1.1.1; better management
  of the dependencies and of the CI (#24).

## [1.3.8] - 2025-03-01
- Release after `metator` is added to `bioconda`.

## [1.3.7] - 2025-02-22
- Auto-release upon tag creation.
- Generate a Docker container hosted on `ghcr.io/koszullab/metator` upon release.

## [1.3.6] - 2025-02-21
- Prepare for bioconda release.
- Only rely on python `pypairix` package (installable from pip). 
- Binaries are now put into `./bin` instead of `./external/artifacts/`. 
- micomplete is now available from `bioconda`. 

## [1.3.4] - 2025-02-19
- Package now relies on `pyproject.toml` for build configuration with `hatch`. 
- Binaries for `louvain` and `leiden` clustering algorithms are now embedded in the package.
- Uses pre-commit hooks for code formatting and linting.
- Fix deprecated Bio.SeqUtils.GC to Bio.SeqUtils.gc_fraction.
- `Biopython` is pinned <= 1.80 to work with `micomplete 1.1.1`.

## [1.3.3] - 2023-11-27
- Improve ci.
- Add pairix as requirements.
- Correct pipeline bug.

## [1.3.2] - 2023-11-27
- Add log files for all modules.

## [1.3.1] - 2023-11-20
- Switch pairs as an argument for `metator mge`.
- Correct bug in QC.
- Add checkV in the requirements.

## [1.3.0] - 2023-11-09
- Merge MetaVIR (https://github.com/ABignaud/MetaVir) module to MetaTOR.
- Automatic generation of the compressed and sorted pairs files in the main 
pipeline.
- Correct QC issues with low-quality libraries.
- Add a conda environment file.

## [1.2.8] - 2023-04-12
- Correct issue with anvio binning file.

## [1.2.7] - 2023-04-07
- Correct bug during final plot generation.

## [1.2.6] - 2023-04-05
- Correct bug of contigs data final bin value from multiple recursive steps.

## [1.2.5] - 2023-04-05
- Remove recursive parallel which seems to failed in some case.
- Correct bug of contigs data final bin value from multiple recursive steps.
- Add exception to avoid error while mutliple metator are started from the same 
working directory.

## [1.2.4] - 2023-04-04
- Add a prefix option for fasta files.

## [1.2.3] - 2023-04-03
- Modify pairs option as an argument for QC module.
- Correct bug if only one MAG corrected in the recursive step.
- Correct figures to work with miComplete output.
- Improve quality check metrics based on mock community results.
- Update QC tests.

## [1.2.2] - 2023-03-29
- Repair the pipeline module.

## [1.2.1] - 2023-03-29
- Make scaffolding optional in pipeline.
- Remove the option to start with network file.
- Allow parallel scaffolding at the end of the pipeline.
- Force pairix index used in the contact map module.
- Change pairs from option to arguments to handle easier multiple files as 
arguments.
- Debug scaffold module default parameters.
- Update docs.

## [1.2.0] - 2023-03-17
- Remove the skip validation option.
- Use miComplete instead of checkM for the recursive step.
- Make multiple recursive step if necessary. 
- Sort pairs and do the scaffolding in the pipeline module.
- Add QC, pairs and scaffold modules.
- Correct pairs indexed loading.
- Add option to remove reads mapping at the edges of the contigs.
- Add possibility to used compressed pairs as input.
- Parallelization of the recursive step.
- Change bin id to make them sortable.
- Correct reads extraction bugs for contact map module.
- Correct bugs from cutsite mapping mode.
- Remove Looseversion sorting Warning.
- Update requirements file.
- Change output plot colors.
- Add tests for QC module.
- Update tests of network, partition and validation modules.
- Update documentation.

## [1.1.6] - 2022-10-06
- Add iterative and cutsite alignment mode.

## [1.1.5] - 2022-07-05
- Add possibility to start with bam files.
- Add SCN normalization to balance the heatmap.
- Add conda environnement and docs.
- Add prefix name for metator contact map.

## [1.1.4] - 2021-07-23

- Correct broken release 1.1.3

## [1.1.3] - 2021-07-23

- Avoid raising error if no pairs extracted in contact map module.
- Add network heatmap plot in the ouput.

## [1.1.2] - 2021-07-13

- Add the possibility to read pairix index files in the contact map module to retrieve faster pairs from a 2D region.

## [1.1.1] - 2021-06-30

Correct some minor issues:

- Correct pairs files:
  - Correct position index (previously 0-based but should be 1-based).
  - Change pairs entry to have an upper triangle shape.
- Correct recursive bin issue in the contig data file were unbinned contigs or binned contigs in bins smaller than the size threshold given were given an id of an existing bin.
- Correct issue in the contact map were contigs needed to be size sorted if a size threshold were given.
- Add some multiple network files in the output if multiple input samples were given.

## [1.1.0] - 2021-05-26

- Add some Figures in the output:
  - Completion/contamination distribution of the MAGs.
  - Distribution of the quality of the MAGs in the assembly.
  - Distribution of the GC and coverage inside each MAG.
- Add test for the new figure module.
- Add a module contactmap to build HiC matrix from some metaTOR objects such as contigs, bins, or other objects with contigs labels from the user. It could be use with metator contactmap [arguments]. The output could be use then with hicstuff view or cooler show to visualize matrix or with instagraal to try to scaffold MAGs.
- Modify the pairs files:
  - Change pairs file extension to .pairs.
  - Reorder the columns of .pairs file: readID-chr1-pos1-chr2-pos2-strand1-strand2 instead of readID-chr1-pos1-strand1-chr2-pos2-strand2.
  - Use tabulation to separate column instead of space to match the official specification.
- Correction of some minor issues:
  - Stop the duplication of the logging message.
  -Silence the HTSlib index error (from pysam.AlignementFile) which was meaningless.
  - Remove pyfastx index.
  - Put as optional the output of the clustering matrix which needs high memory usage for large assembly and optimize memory usage to build it.
  - Modify Assignment method to dataframe to avoid assignment to copy of the dataframe (SettingWithCopyWarning).

## [1.0.4] - 2021-05-14

- Change the format of the pairs file to match with the official specification.
- Add the HiC coverage (mean number of contacts in 1kb region) in the final bin summary.
- Add possibility to start with the pairs file.
- Add info about general alignment if multiples files are given.
- Organize temporary folder.

## [1.0.3] - 2021-05-04
  
- Add clustering matrix in the output of the partition, validation and pipeline modules.
- Add output files documentation.
- Add codeco.

## [1.0.2] - 2021-05-01

- Modification of the ouput files:
  - Add bin coverage based on the depth file given in the bin_summary file.
  - Add a final bin column in the contig data file after the validation to specify the name of the final bin of the contig. A "ND" value is given if the contig isn't binned.
  - Add the pairs file of the alignment to be able to produce HiC matrices from it using hicstuff for example.
  - Add a log file with all information of the run.
  - With the pipeline workflow remove intermediate contig data files.
- Modify call to input data to use only header names and not order of the column.
- Change the parameter of the bowtie2 alignement back to --very-sensitive-local.
- Add the total size of the different category of bins in the log.
- Correct minor issues in the doc.

## [1.0.1] - 2021-04-27

Correct bug in validation module:

- louvain_recursif: correct type error to identify contigs of a subnetwork which could triggers empty subnetwork error.

## [1.0.0] - 2021-04-26

New version of metaTOR in Python3
