# Change Log

All notable changes to this project will be documented in this file.

## Unreleased

### Added
- `network` and `pipeline` start from a `.metacool` contact map with
  `--start metacool`, the map given to `-1` (a file, or an entry
  `file.mmetacool::entries/<name>`), at any bin size. The network is that
  of the pairs the map was made of, byte for byte at `--edge 0` (see
  Changed). A map may lack contigs of the assembly, which get no
  contacts; a map of other contigs or lengths, or that stores both
  triangles, is refused, and so are `--edge` above 0, several maps and
  `pipeline --scaffold`. A map sums its libraries: no `network_{i}.txt`
  or `hit_data_alignment.txt`. Python API:
  `metator.network.metacool_to_contacts`, `read_metacool_contacts`.
- A `pipeline` run from a map writes the map of its bins from it, with
  `metacooler reassign`: its bins and contacts, grouped by final bin as
  from pairs, with this run's genome and contig columns. Its contact
  statistics are recomputed from the pixels (at one bin per contig, every
  intra-contig contact then counts as short-range). The output may be
  the input map, replaced only once the new map is validated. Python API:
  `metator.metacool.reassign_metacool`. The README shows the route from
  the pairs of `hicstuff` through `metacooler create pairs`.
- `metator metacool` writes this map from a `.metacool` map given in place
  of the pairs files, so that a run step by step (`network -S metacool`,
  `partition`, `validation`) can write the map that `pipeline -S
  metacool` writes.

### Changed
- Python API: `metator.partition.partition` no longer takes
  `temp_directory`, nor `generate_fasta` and
  `metator.validation.update_contigs_data_recursive` `tmpdir`, nor
  `metator.mge.mge_binning`, `generate_mges_fasta` and
  `MetatorObject.write_fasta` `tmp_dir`: they only held the contig lists
  of the `pyfastx` command. `partition` and `mge` write no temporary file
  any more.
- At `--edge 0` (the default), the network of pairs files counts the pairs
  as metacooler counts them in a map: a pair with a read on the last base
  of a contig, or at position 0 (which a 1-based pairs file should not
  hold, but some do), is a contact; a pair with a read past the end of its
  contig is none, at any `--edge`, and its hits no longer count. The
  network of a `.metacool` map is then that of its pairs, byte for byte.
  On `tests_data/alignment.pairs`, 360 of 13,415 edges change (298 new,
  62 of another weight), from its 472 pairs between two contigs with a
  read at position 0.
- The documentation moves from Read the Docs to GitHub Pages, at
  <https://koszullab.github.io/metator/>: Markdown pages in `docs/`, built
  by Zensical with chromokit, the look of metacooler's documentation
  (`uv sync --extra docs`, then `uv run zensical serve`). It holds a
  quickstart, the concepts of the binning, a guide per step and command,
  the CLI reference, the file formats, and the Python API reference,
  written from the docstrings. The tutorials of `docs/example/` are guides
  in `docs/guides/`, brought up to date with the commands and outputs of
  this version, the changelog is `docs/changelog.md`, and `CONTRIBUTING.md`
  points to the contributing guide of the site. `tests/test_docs.py` checks
  the CLI reference against the CLI. The docstrings that disagreed with
  their code (`metator.align.get_contact_pairs`, `merge_alignment`,
  `process_bamfile` and `process_bwa_bamfile`,
  `metator.scaffold.parallel_scaffold` and
  `metator.validation.correct_final_bin`) are corrected.
- Python API: `metator.validation.recursive_clustering` takes the
  clustering matrix of the earlier recursive steps (`clustering_matrix`)
  and returns it, updated, instead of writing
  `clustering_matrix_recursive.npz`, which `recursive_decontamination`
  writes once, at the end.

### Removed
- Python API: the CheckM functions (`metator.validation.checkm` and
  `checkm_compare_bins`, `metator.io.check_checkm` and
  `read_results_checkm`), unused since miComplete replaced CheckM, and
  `metator.partition.spinglass_partition`, whose body was commented out.
- The Sphinx documentation of Read the Docs (`docs/conf.py`,
  `docs/index.rst`, `docs/api/*.rst`, `docs/requirements.txt`...), and the
  scripts of the curation tutorial, `docs/example/scripts/`:
  `anvio_figures.py` read CheckM's results with `read_results_checkm`,
  removed above, and `checkM.sh` ran CheckM. The curation guide checks the
  curated bins with miComplete, through the Python API, instead.

### Fixed
- `network` and `pipeline --start bam --aligner bwa` no longer fail without
  `--reverse` (`AttributeError: 'NoneType' object has no attribute
  'split'`): a BAM file of bwa holds both mates. `pipeline` requires
  `--forward`, as `network` does: a run without it failed with the same
  error.
- `validation` and `pipeline --cluster-matrix` no longer fail as soon as a
  bin is partitioned again (`ValueError: axis 1 index … exceeds matrix
  dimension 2`): the clustering matrix of the recursive partition was built
  2 columns wide. It now holds the pairs of contigs of every bin partitioned
  again, with the values of its latest partition: each recursive step
  overwrote the matrix of the steps before, usually with an empty one, as
  the last step partitions no bin.
- `binning.txt` holds every contig of the final bins, as
  `contig_data_final.txt` does: the contigs of a bin partitioned again but
  kept whole, which had gone to a recursive part, were left out. Python
  API: `metator.validation.write_bins_contigs` writes the `Final_bin` of the
  contigs, and no longer takes `prefix`.
- The log of the validation counts the bins in the quality classes of
  `plot/bins_size_distribution.png` (`metator.figures.bin_quality`, the
  contamination being the weighted redundancy − 1). It counted as HQ the
  bins of 90% completeness or more with a contamination above 5%, and as MQ
  those of 5% or less, and divided the contamination by the completeness,
  where the plot does not.
- `metator validation --help` and the docstrings of `metator.validation`
  give the rule of the recursive partition: a weighted completeness of 33%,
  not 50%, and the completeness its parts must keep.
- The `recursive` column of `bin_summary.txt` says whether the bin comes
  from the recursive partition: every bin read `True` when no bin was
  partitioned again, and a bin made at an earlier recursive step `False`.
- The validation, and `pipeline --scaffold`, leave no pyfastx index
  (`<bin>.fa.fxi`) next to the fasta files of the final bins.
- The depth file (`--depth`) is read by contig name: it was read line by
  line, in the order of the assembly, and a depth file in another order
  gave the contigs the length and coverage of others, without a word. A
  depth file that lacks a contig of the assembly, or gives it another
  length, stops the run with an error, and the lengths of the contigs come
  from the assembly.
- `network` and `pipeline` read an assembly compressed with gzip, which
  they refused as not being a FASTA file.
- The index of the assembly is that of its own name only
  (`assembly.fa.amb`... for bwa, `assembly.fa.1.bt2`... for bowtie2): with
  bwa, any bwa index of the directory of the assembly, that of another
  assembly for instance, was taken for its index, and with bowtie2 any
  index whose name starts with that of the assembly.
- `network --self-contacts` leaves out the contacts between contigs within
  `--edge` of a contig end, as without `--self-contacts`: `--edge` was
  ignored with it.
- `contactmap --object recursive_bin` selects the recursive bin within its
  overlapping bin: recursive bins are numbered within their overlapping
  bin, and the recursive bins of the same number of every overlapping bin
  were mapped together. It needs the name of the bin
  (`metator_00012_00003`), and the names of bins with a `--prefix` holding
  `_` are read right. The fasta file of the map is named after `--name`.
- `qc` says why it checks nothing: it logged `No MAGs have been found.`
  whatever the reason (no high-quality MAG, no contig of more than 100 kb,
  or no pair between them), and logs the thresholds given with
  `--threshold` too.
- `metator metacool` and the map of `pipeline` refuse a
  `contig_data_final.txt` whose rows have more or fewer fields than its
  header, which they read shifted or padded, without a word, and one
  compressed otherwise than with gzip as `.gz`, which metacooler failed to
  read once the map was written.
- `pipeline` checks its arguments and its assembly before replacing the
  bins of an earlier run with `--force`: a refused run no longer empties
  `overlapping_bin/`, `recursive_bin/` and `final_bin_unscaffold/` first.
- The fasta files of the bins (`partition`, `validation`, `pipeline`), of
  the MGE bins (`mge`) and of `contactmap` are written with pyfastx's
  Python API (`metator.io.write_contigs`), byte for byte as before. The
  `pyfastx` command wrote them through the shell: without it on the `PATH`
  (a virtual environment used without activating it), every file was
  empty, without an error. A contig missing from the assembly now stops
  the run with `KeyError`, where the command left it out of its bin.
- `pipeline`, `contactmap` and `mge` remove the pyfastx index of the
  assembly only when it is there: they failed at their very end, with
  `FileNotFoundError`, when another run sharing the assembly had removed
  it.
- The network heatmap of `pipeline` and `validation` reads only the files
  of the output directory whose name ends in `.pairs.gz`: an unescaped dot
  also took, for example, `alignment_sorted_pairs.gz` for a pairs file.
- `qc` no longer fails with `ZeroDivisionError` when none of the pairs
  between its high-quality MAGs falls within one of them, or none is
  informative: the ratios that have nothing to divide are 0, and a pie
  chart of no pair is not drawn. The ratio of informative contacts
  between contigs is logged as a percentage, as its `%` says.
- `pipeline` and `validation` stop with a clear error when no bin of the
  partition reached `--size` (a small assembly with the default 500 kb,
  say). They failed in miComplete with `EmptyDataError: No columns to
  parse from file`, leaving its working directory in the output
  directory.

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
  example, contig names that are not ASCII, which cooler cannot read, or a
  misaligned bin summary, see below), a warning says why and the run keeps
  its other outputs and its exit status. A map left in the output directory
  by an earlier run is removed first, so that it cannot pass for the map of
  this run. See the README.
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
