# Change Log

All notable changes to this project will be documented in this file.

## Unreleased

metaTOR 2.0 clusters contigs, and nothing else: every step reads a
`.metacool` contact map of the contigs and stores its results in it.
[Migrating to metaTOR 2.0](migration.md) details every change; [Other
tools](other-tools.md) says what replaces each removed command. Below,
`<name>` is the name of the map without `.metacool`.

### Breaking changes
- **A contact map is the input.** `network`, `partition`, `validation`,
  `pipeline` and `mge` require a `.metacool` map of the contigs, `MAP`, as
  first argument: a file or a collection entry
  (`file.mmetacool::entries/<name>`), at any bin size. hicstuff
  (`pipeline -M metacool`, 4.0.0+) or metacooler (`create pairs`, `create
  cool`, `merge`) makes it; so do `network` and `pipeline` from pairs or BAM
  files.
- **Removed options:** `-1/--forward`, `-2/--reverse`, `-S/--start`,
  `-b/--aligner`, `-B/--aligner-mode` and `-E/--edge`, now usage errors (exit
  status 2). `--assembly` is the FASTA file (plain or gzip), not an index.
- **Removed outputs:** `alignment_<i>*.bam`, `alignment_<i>_sorted.pairs.gz`,
  `network_<i>.txt` and `hit_data_alignment.txt`. Dropped dependencies:
  hicstuff, bowtie2, bwa, samtools, pypairix and pyfastx.
- **The map is the interface between the steps.**
  - `network MAP` stores the contig values (`GC_content`, `Hit`,
    `Shotgun_coverage`, `Restriction_site`), the normalization (bin weights
    `metator_<normalization>`) and its run. No `network.txt`;
    `contig_data_network.txt` goes to the temporary directory. The network is
    summed in memory: `-t/--threads` no longer sorts it.
  - `partition MAP` stores every core and overlapping cluster (rounds `core`,
    `overlapping`); no `contig_data_partition.txt` or `overlapping_bin/`.
    Removed: `-a/--assembly`, `-c/--contigs`, `-n/--network`, `-s/--size`.
  - `validation MAP` reads the partition from the map (refused if made on
    another network), stores the final clusters as its genomes, and writes
    `fasta/mags_with-mges/` and `plots/`. `--size` picks the overlapping
    clusters miComplete checks. Removed: `-c/--contigs`, `-f/--fasta`,
    `-n/--network`.
  - The tables (`contig_data_final.txt`, `cluster_summary.txt`,
    `binning.txt`, `clusters.txt`, `cluster_members.txt`, miComplete's) and
    `recursive_clusters/` go to the temporary directory (`-N` keeps it);
    `metator dump` writes them from the map.
  - `pipeline MAP` writes the map of the clusters once, at its end, as
    `<outdir>/<name>.metacool` (1.4.0: `<outdir>/metator.metacool`): the map
    that the steps give when run one by one in one output directory.
- **Where the map goes.** A run writes its map to its output directory under
  the name of `MAP`, leaving `MAP` as it was; `--inplace`, or an output
  directory that holds `MAP`, updates `MAP` itself. A run that stops leaves
  the map as it was (see [Conventions](cli.md#conventions)); one that cannot
  write it fails (1.4.0 warned and exited 0).
- **Earlier outputs.** The outputs of 1.x left in the output directory
  (tables, networks, `*_bin/` folders, clustering matrices, old plots) are
  refused, with or without `-F`, and never removed. `-F` replaces only what
  the run writes: its folders, cluster plots, tables, report and map. Earlier
  plots are refused without `-F` (1.4.0 overwrote them).
- **Renamed outputs: clusters, not bins.** `bin_summary.txt` →
  `cluster_summary.txt`; `final_bin/`, `final_bin_unscaffold/` →
  `fasta/mags_with-mges/`; `recursive_bin/` → `recursive_clusters/`;
  `plot/*.png` → `plots/*.pdf`; columns `Core_bin_ID`... `Final_bin` →
  `Core_cluster_ID`... `Final_cluster` (`metator.io.rename_1x_columns`
  converts a 1.x table). Plot labels say clusters.
- **Renamed log:** `<name>.metator_<date>.log` (1.x:
  `metator_<command>_<date>.log`).
- **Contigs follow the map.** A contig's ID is its position in the map
  (longest first), which orders the table rows, cluster numbers and FASTA
  records. Contigs missing from the map are in no table; shotgun abundance is
  a share of the map's contigs.
- **The network counts every pair**, as the map does: a read at position 0
  or on a contig's last base counts (360 of 13,415 edges change on
  `tests_data/alignment.pairs`); one past the contig end does not. Edges
  weigh `n × w_i × w_j` (1.x: `n / √(h_i × h_j)`, different in the last
  bits); edges of weight 0 are dropped.
- **`mge` reads the map:** `metator mge MAP -m MGES -a ASSEMBLY` needs the
  network of `metator network` and metacooler's contact statistics (0.4.0+)
  in the map. Removed: the pairs inputs, `-c` and `-b`. MGE contigs missing
  from the map are skipped with a warning; the others are taken longest
  first, which can renumber clusters. The MGE clusters go to the map (rounds
  `mge`, `mge_random`), not `mges_bin_summary.tsv`, and `fasta/mges/`
  replaces `mges_binned.fa`; both old files are refused in the output
  directory.
- **Removed commands and features** (replacements in [Other
  tools](other-tools.md)): `metator qc`, `contactmap`, `metacool`, `scaffold`
  and `pairs`; `pipeline -v/--scaffold`, `-J/--junctions` and `--no-metacool`;
  `-C/--cluster-matrix` and its matrices; `network -s/--self-contacts`; the
  plots `network_heatmap.png`, `MAGs_GC_distribution.png`,
  `MAGs_HiC_cov_distribution.png` and `MAGs_SG_cov_distribution.png`;
  `python -m metator.main` (use `metator` or `python -m metator`).
- **Seeded runs** give other clusters than before: the temporary directory
  name no longer draws from Python's `random` module.
- **Python API.**
  - Removed: the modules `metator.align`, `contact_map`, `main`,
    `quality_check`, `regions` and `scaffold`; the readers and writers of
    files 2.0 no longer handles (`metator.io.sort_pairs`...);
    the unused CheckM functions and `metator.partition.spinglass_partition`.
    The development map writers (`reassign_metacool`, `store_network`,
    `write_rounds`, `write_round`) are now `metator.metacool.write_map`.
  - `metator.figures` keeps the quality plots only:
    `figures_bins_distribution`, `figures_bins_size_distribution` and
    `plot_figures` (PDF files in `plots/`); `bin_quality` moved to
    `metator.validation`.
  - `metator.partition.partition` returns a `metator.clusters.Partition` and
    writes nothing; `metator.validation.recursive_decontamination` and
    `metator.mge.mge_binning` take the map.
  - `metator.host` reads the map (`find_hosts`). Removed: `Bin`, `Mag`,
    `MgeMag`, `create_bins`, `map_contigs_to_mags`, `build_contig_graph`,
    `compute_mge_mag_interactions`, `annotate_hosts`, and the unused noise
    estimate (`estimate_noise`, `plot_mag_noise`, `classify_mags`,
    `background_noise.png`).

### Added
- **`network` and `pipeline` from pairs or BAM files:** `INPUT...` is the
  map, or the pairs or BAM/SAM files of the libraries (paths, comma-separated
  lists or globs). The run first makes the map in its output directory
  (`<prefix>.metacool` with `pipeline --prefix`, else `<assembly>.metacool`;
  `[1/6] map` in the log). BAM pairs: primary alignments of both mates, MAPQ
  ≥ `-q/--min-quality` (30), no duplicates or QC failures, at their 5' ends.
  pysam is a dependency again; `network -t` reads the files. See [The inputs
  of a run](cli.md#inputs).
- **`pipeline --mges FILE`** runs `mge --reassign` and `host` after the
  step `markers` (`[5/7] mge`, `[6/7] host`, `[7/7] export`), with
  `--mge-threshold-bin`, `--mge-min-share` and `--mge-min-contigs`.
- **`metator markers MAP -a ASSEMBLY`** finds miComplete's marker genes on
  each contig, as miComplete finds them in its genome (pyrodigal trained on
  the contigs of the genome in map order, pyhmmer with miComplete's
  `--evalue`, `--bias` and `--domain-cutoff`; metagenome mode for the contigs
  of `__unbinned__` and of genomes under 20 kb), and stores them in the map:
  the contigs columns `markers_Bact105` and `markers_Arch131` (text,
  `Ribosomal_L22,SecY*2`; `*0` for a marker miComplete sees but drops, which
  counts in the weighted completeness only) and the `markers` metadata entry,
  with how the genes were found and the markers and weights of both sets.
  Summed over any group of contigs they give miComplete's completeness,
  redundancy, weighted completeness and weighted redundancy
  (`metator.markers.marker_quality`), those of a genome exactly.
  `pipeline` runs it after the validation (`[4/5] markers`); every run keeps
  them. See [Marker genes of the contigs](formats.md#marker-genes-of-the-contigs).
- **`metator host MAP`** finds the MAGs that host the MGE clusters stored by
  `mge`: at least `--min-share` (10%) of the MGE cluster's signal with the
  MAGs, and `--min-contigs` (5) contigs in contact. It writes
  `tables/mge_hosts.txt` and `plots/mge_hosts.pdf`, and stores the hosts in
  the map. Bins of another binner work once `metacooler reassign` makes them
  genomes.
- **`metator dump MAP -o DIR`** writes the tables of a run from its map,
  `network.txt` included (by contig name). See [`dump`](cli.md#dump).
- **`metator export MAP [-a ASSEMBLY] -o DIR`** writes every output of a run
  again from its map: `cool/`, `tables/`, `plots/`, the report, and `fasta/`
  with the assembly. `pipeline` runs it last. See [`export`](cli.md#export).
- **`metator report MAP -o DIR`** writes the report of a run:
  `<name>.report.html` (a page with its MAGs, MGE hosts, parameters and
  plots), `<name>.report.txt` (a versioned text format) and
  `<name>.report_mqc.yaml` (a MultiQC General Statistics row). `export` and
  `pipeline` write it too. New dependencies: `jinja2`, `pyyaml`. See [The
  report of a run](formats.md#the-report-of-a-run).
- **Output folders** (see [Conventions](cli.md#conventions)):
  `fasta/mags_with-mges/`, `fasta/mags/` (each MAG as contigs and stitched),
  `fasta/mges/`, `fasta/mges_random/`, `cool/` (a balanced cool file per
  MAG, with its MGE clusters), `tables/` and `plots/`.
- **`mge --reassign`** makes each MGE cluster a genome of the map;
  `binning.txt` then lists the MAGs and the MGE genomes.
- **`--seed`** on `partition`, `validation`, `pipeline` and `mge`: the same
  map and seed give the same clusters.
- **`--inplace`** on `network`, `partition`, `validation`, `pipeline` and
  `mge`; **`-F`** on `network` and `mge` (and `partition -F`) replaces the
  map in the output directory.
- **The map stores** every clustering round (`core`, `overlapping`,
  `recursive_<k>`, `final`, `mge`, `mge_random`) with its contigs, the
  network normalization as bin weights that `cooler` and `metacooler
  --balance` read, the cluster summary at full precision, and one metadata
  entry per step (version, command line, parameters, what it read).
- Without `--enzyme` and `--depth`, a run reads the map's `Restriction_site`
  and `Shotgun_coverage` columns.
- `partition` and `validation` refuse a map whose network changed since
  `network` stored it (e.g. after `metacooler merge`).
- Python API: `metator.dump`; `metator.clusters`; `metator.host`
  (`find_hosts`...); `metator.markers` (`contig_markers`, `group_markers`,
  `marker_quality`, `quality_of`, `encode_markers`, `decode_markers`...);
  `metator.metacool.write_map`, `map_rounds`, `read_partition`,
  `holds_hosts`, `hosts_network`, `holds_markers`;
  `metator.network.map_network`, `stored_network`, `contig_contacts`,
  `contig_hits`, `normalization_weights`, `relabel_network`;
  `metator.mge.check_mge_map`, `map_matrix`, `write_mges_fasta`;
  `metator.validation.contaminated_bins`; `metator.partition.iterate`, and
  `progress` on its iteration functions; `metator.log.duration`.
- The documentation page [The end-to-end pipeline](end-to-end.md): from the
  reads to MAGs and MGE hosts, with MEGAHIT, hicstuff, geNomad and metaTOR.

### Changed
- The validation runs Prodigal and HMMER through pyrodigal and pyhmmer:
  metaTOR needs no program that pip does not install.
- metaTOR requires metacooler 0.4.0, which writes the FASTA files of the
  clusters. A gzip assembly is decompressed once per run.
- Numbers read back from the tables are exact, under pandas 2 and 3.
- The partition computes core-cluster distances in threads, not processes.
- Python warnings are shown and logged.
- The quality plots are drawn without pyplot; one that fails is a warning
  (1.4.0 stopped).
- **The log:** a hicstuff-style header (`## cmd:`, `## <option>: <value>`)
  and a last line with what the run wrote and how long it took. `pipeline`
  numbers its steps (`[1/5] network`...); the partition logs its iterations
  as they go, the validation each recursive step; bookkeeping is at DEBUG
  level.
- `metator --help` lists the commands in three panels: Pipeline, Steps and
  Outputs.
- The documentation moves from Read the Docs to GitHub Pages
  (<https://koszullab.github.io/metator/>), in Markdown built with Zensical,
  with a CLI reference that `tests/test_docs.py` checks. The Sphinx sources
  and the CheckM curation scripts (`docs/example/scripts/`) are gone.
- The documentation was trimmed: shorter pages, each fact said once.
- Shorter numpy-style docstrings.

### Fixed
- `mge` refuses an unreadable `--mges` list before any work (it stopped with
  a traceback), and skips blank lines (`IndexError`).
- `--depth` is read by contig name: a file in another order than the
  assembly gave contigs wrong values. Missing contigs or wrong lengths are
  refused.
- `--depth` lengths in scientific notation (contigs of 1 Mb or more) stopped
  `network` and `pipeline` (#26).
- `network` and `pipeline` read a gzip-compressed assembly.
- `--enzyme` counts sites in lower-case (soft-masked) bases too:
  `Restriction_site`, and the `RS` and `theoritical_hit` weights, were wrong.
- An assembly naming a contig in several records is refused (the last one was
  kept).
- `validation` keeps the clusters of contigs named like numbers (`000017`,
  `1e5`) or `NA`.
- `binning.txt` lists the contigs of a cluster partitioned again but kept
  whole.
- The validation log counts the quality classes as
  `metator.validation.bin_quality` does: from 1.2.0 to 1.4.0 it swapped HQ
  and MQ, and divided the contamination by the completeness.
- `plots/bins_distribution.pdf` draws the contamination from 0 to 105%
  (a cluster without marker genes was at −100%); typos fixed in the labels
  of the quality plots.
- The `recursive` column of `cluster_summary.txt` says whether a cluster
  comes from the recursive partition.
- `metator validation --help` gives the recursive partition threshold as
  33% weighted completeness, not 50%.
- `pipeline` and `validation` check their arguments before `-F` removes
  anything, and refuse a `--tmpdir` inside a folder they replace.
- `pipeline` and `validation` stop with a clear error when no cluster
  reaches `--size` (miComplete failed).
- An output path containing `::` is refused; a map given as `file::/` is the
  file.
- No pyfastx index (`.fxi`) is left next to the assembly or the FASTA files.
- The FASTA files of the clusters were empty, without error, when `pyfastx`
  was not on the `PATH`.

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
