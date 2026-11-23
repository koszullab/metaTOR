# Change Log

All notable changes to this project will be documented in this file.

## Unreleased

### Added
- A `pipeline` run with several libraries builds the map of each library,
  the entries `library_0`, `library_1`... of `<prefix>.mmetacool`, and the
  map of their sum, `<prefix>.metacool`, merged from them with metacooler
  merge (contacts and contact statistics summed; its counts are stored as
  64-bit integers, where a map written from all the pairs files at once
  stores 32-bit ones). It reads `network.txt` and `contig_data_network.txt`
  from the sum, and the network of each library (`network_N.txt`) and
  `hit_data_alignment.txt` from its map, instead of parsing the pairs
  files again: the same files, byte for byte, as from the pairs files. So
  zstd and LZ4 pairs files work with several libraries too (not with
  `--scaffold`, which reads the pairs files), and a pair on a contig that
  the assembly lacks is left out, as metacooler warns, where the run
  stopped on a `KeyError`. `--no-metacool` removes both files at the end.
  Python API: `metator.metacool.write_library_maps`, and `libraries` of
  `metator.network.metacool_to_contacts`.
- The clustering rounds of a run: the validation writes `clusters.txt`,
  every cluster of every round (`overlapping`, `recursive_<k>`, `final`)
  with the round and name of its parent and miComplete's columns, and
  `cluster_members.txt`, the contigs of each. The map of the clusters
  stores both in its `/clusters` group (metacooler 0.4.0), so that the
  intermediate clusters are kept without being written into its genomes
  and contigs tables, which hold the final ones. `metator metacool` stores
  none. Python API: `metator.clusters`, and `rounds` of
  `metator.metacool.reassign_metacool`, which stores them (without them,
  the map keeps none).
- `pipeline --output-map FILE` writes the map of the clusters to FILE
  (replaced only with `-F`). A run from reads, alignments or pairs writes
  the maps of its libraries next to FILE, `<name>.mmetacool`, with several
  (and removes an earlier run's there, with one): one that is there is
  refused without `-F` too, before any work. A FILE named `.mmetacool`,
  where the maps of the libraries would be written over the map, is
  refused.
- `network` and `pipeline` start from a `.metacool` contact map, their
  first argument (`metator pipeline MAP -a ASSEMBLY`: a file, or an entry
  `file.mmetacool::entries/<name>`), at any bin size, in place of `-1`,
  `-2` and `--start`, which a run from a map refuses. The network is that
  of the pairs the map was made of, byte for byte, as both count every
  pair (see Changed). A map may lack contigs of the assembly, which get
  no contacts; a map of other contigs or lengths, or that stores both
  triangles, is refused, and so are a first argument that is a file but
  no map, and `pipeline --scaffold`. A map sums its libraries: no
  `network_{i}.txt` or `hit_data_alignment.txt`. `-1` is no longer
  required: a run takes a map or reads, alignments or pairs. Python API:
  `metator.network.metacool_to_contacts`, `read_metacool_contacts`.
- A `pipeline` run from a map writes the map of its bins from it, with
  `metacooler reassign`: its bins and contacts, grouped by final bin as
  from pairs, with this run's genome and contig columns (the contig
  columns last, as from pairs, whatever columns the input map carried),
  and the contact statistics of the input map, carried over to the bins
  (a map written by metacooler 0.2.0 has them computed again from the
  pixels: at one bin per contig, every intra-contig contact then counts
  as short-range). The output may be the input map, replaced only once
  the new map is validated. Python API:
  `metator.metacool.reassign_metacool`. The README shows the route from
  the pairs of `hicstuff` through `metacooler create pairs`.
- `metator metacool` writes this map from a `.metacool` map given in place
  of the pairs files, so that a run step by step from a map (`network
  MAP`, then `partition` and `validation`) can write the map that
  `pipeline MAP` writes, but for its clustering rounds.
- A run from a map (`network MAP` and `pipeline MAP`)
  reads the restriction sites and the shotgun coverage of the contigs from
  the `Restriction_site` and `Shotgun_coverage` columns of the map, which
  metaTOR stores in the maps it writes with `--enzyme` and `--depth`: the
  same network and tables as with these options, without them (a run
  from such a map wrote `-` in their place, and the map of its bins
  dropped the columns). The options fill in a column that the map lacks
  and replace one that it has, with a warning that counts the contigs
  whose value differs; the map of the bins stores what the run used. A
  column that does not give every contig of the assembly a value is not
  read, with a warning that names the contigs it lacks, and a
  normalization that needs a column that neither the map nor the options
  give is refused before the run replaces anything. The metadata of the
  map of the bins that `pipeline` writes says which columns the run read
  (`columns_from_map`). With pandas 2, a coverage of more than 15
  significant digits, or of 16 or more digits after the decimal point,
  can be stored one unit off in its last place, which changes the last
  digits of the contig data files and of a network normalized by it. A
  run from pairs reads its contig data before it replaces the bins of an
  earlier run (`-F`): a depth file of another assembly, which it refuses
  then, leaves them (they were removed first). Python API:
  `metator.network.contig_data_from_map`, `map_columns`,
  `check_map_contigs`.
- `partition` and `validation` take the map of a run step by step as
  their first argument too (`metator partition MAP ...`, a file or an
  entry of a collection), and update it as `pipeline MAP` does: the steps
  run one by one on a map (`network MAP`, `partition MAP`, `validation
  MAP`) write the map that `pipeline MAP` writes from the same clusters,
  but for the command line and parameters it stores (Louvain draws random
  numbers: steps run as separate commands differ from `pipeline` as two
  `pipeline` runs do). `partition MAP` stores its overlapping clusters
  that reach `--size`, with their contigs, as the round `overlapping` of
  the clustering rounds of the map, in place of those it held, and leaves
  its genomes and contigs tables and the `metator` key of its metadata as
  they were (an earlier run's, until `validation MAP` writes the new
  ones); the round records the step that made it. `validation MAP` writes
  the map of the clusters, with every round, and draws the heatmap of the
  clusters from the map.
  Each writes the map in place (staged, validated in full, then renamed
  over it, a symbolic link replaced, not its target, the mode kept) or to
  `--output-map FILE` (replaced only with `-F`), and holds its lock while
  it runs. Refused before any work, the reason logged as an error, and
  before `-F` replaces the directories of clusters of an earlier run:
  what `pipeline MAP` refuses (a map that is not there, is no map, or that
  metacooler cannot open, such as a `.mmetacool` collection, an entry of a
  collection without `--output-map`, an `--output-map` that exists without
  `-F`, a map that another run holds), a map that is no map of the
  contigs of their contig data (`-c`: other contigs or lengths), an
  `--output-map` that is one of their inputs, `--output-map` without a
  map, and an invalid `--prefix`. Without a map, they run as before, but
  for an invalid `--prefix`, now refused before `-F` replaces their
  directories of clusters; `pipeline` refuses an `--output-map` that is
  its assembly or depth file, which the map replaced. Python API:
  `metator.metacool.write_rounds`, `metator.clusters.overlapping_round`,
  `Rounds.tables`, and the `names` of `metator.metacool.store_rounds` (a
  round with no cluster).
- `qc` and `contactmap` take the `.metacool` contact map of the run (a
  file, or an entry of a collection) in place of their pairs files, for
  what a map gives: the contacts between contigs, those of the pairs it
  was made of. `metator qc MAP` measures the noise of the library, the
  contacts within contigs, within MAGs and between MAGs, the noise
  contact ratio and the noise score, each equal to what `qc` measures
  from the pairs files, when their pairs are in the upper triangle, as
  metaTOR's alignment writes them; the events of the pairs within a
  contig, which need the positions and strands of their reads, are not
  measured, which the log says. `metator contactmap MAP` writes the
  contact map of its contigs as `metacooler extract` writes it,
  `<name>.cool`, at the bins of the map (one per contig in a map of
  metaTOR): hicstuff's map of its restriction fragments from the pairs
  files, summed by pair of contigs. Refused before any work, the reason
  logged as an error: a map with pairs files or several maps, a map that
  metacooler cannot open, and the options that need the read positions
  (`qc --threshold` and `--plot`, `contactmap --filter`, `--pcr-dup` and
  a `--mat-fmt` other than `cool`); then, before `contactmap` writes or
  replaces an output, a map that gives a contig another length than the
  contig data table, or that stores both triangles. `--enzyme` is not
  used, and a contig that the map lacks has no contacts, with a warning.
  Python API: `metator.quality_check.extract_map_pairs`, `map_quality`,
  `noise_metrics`; `metator.contact_map.extract_contact_map`;
  `metator.metacool.contig_contacts`, `held_contigs`, `extract_contigs`.
- `mge --map MAP` stores the MGE clusters of the run in the contact map
  of the run, as the round `mge` of its clustering rounds: one cluster
  per row of `mges_cluster_summary.tsv`, with its contigs and its columns
  (its length and number of contigs computed by metacooler), and the
  version, command line and parameters of the run, after the rounds of
  the binning run, which are kept, in place of the round `mge` of an
  earlier run. A column of whole numbers, or of booleans, of the rounds
  of the binning run, which the round `mge` lacks, is stored as
  decimals, or as text, from then on, as any column with a missing value
  is. The genomes and contigs tables, the final clusters, and the rest
  of the map are left as they were. The map is written as
  `partition MAP` writes it, in place or to `--output-map FILE`
  (replaced only with `-F`, new to `mge`), and locked while it is
  written; what `partition MAP` refuses is refused before any work (a
  map of other contigs or lengths than the contig data, an
  `--output-map` that is one of the inputs of the run), and so are a map
  given in place of the pairs files, which `mge` reads for the contacts
  within each MGE contig, and `--random` with `--map`. A run that
  replaces the clusters of the map (`partition MAP`, `validation MAP`,
  `pipeline MAP`, `metator metacool MAP`) drops the round `mge`, made
  with the clusters it replaces, with a warning. Python API:
  `metator.metacool.write_round`, `metator.clusters.mge_round`.

### Changed
- **Breaking:** metaTOR calls its groups of contigs clusters, as its contact
  map does: the bins are those of the map. The outputs of a run are renamed:
  `bin_summary.txt` is `cluster_summary.txt`; `overlapping_bin/` and
  `recursive_bin/` are `overlapping_clusters/` and `recursive_clusters/`;
  `final_bin/` (`validation`) and `final_bin_unscaffold/` (`pipeline`) are
  `final_clusters/`; `scaffold_bin/` is `scaffolded_clusters/`;
  `mges_bin_summary.tsv` (`mge`) is `mges_cluster_summary.tsv`. In the
  contig data tables, `Core_bin_ID`, `Core_bin_contigs` and `Core_bin_size`,
  and the columns of the overlapping and recursive bins, are
  `Core_cluster_ID`, `Core_cluster_contigs`, `Core_cluster_size` and so on,
  and `Final_bin` is `Final_cluster`. Every command still reads the contig
  data tables of metaTOR 1.x, their columns under the new names
  (`metator.io.rename_1x_columns`), and `qc -O` finds the `bin_summary.txt`
  of a 1.x run. `binning.txt`, which anvi'o imports, and the names of the
  clusters are kept. A script that gives the overlapping clusters of
  `partition` to `validation -f` names `overlapping_clusters/` now: see
  [Migrating to metaTOR 2.0](migration.md). A run in the output directory of
  a run of 1.x refuses its outputs named after bins without `-F`, as those
  of 2.0, and removes them with it. The Python functions that take a table
  of contig data need the 2.0 names (read a table of 1.x with
  `metator.io.read_contig_data`). `contactmap --object Final_bin` reads
  `Final_cluster`, unless the table has a `Final_bin` column. metacooler's
  `--contigs-format metator` names `Final_bin`: give
  `--contigs-format Name,Final_cluster` for a table of 2.0. The help, the
  logs, the texts of the plots and the docs call them clusters too, and the
  help groups the three steps and `pipeline` under Clustering.
- `pipeline` no longer sorts, compresses and indexes the pairs files of
  a run from reads or alignments, which it removed at its end: the maps
  of the contigs read them once, in any order, and the networks come
  from the maps. Without `--scaffold` or `-N`, which keep them
  sorted in the output directory as before, it writes them to its
  temporary directory as they are made (`alignment_<i>.pairs`), and
  removes them once its network is built: a run that stops later leaves
  none, but leaves the contact map of its contigs, which a run restarts
  from. The tables of the alignments of the mates of each library, as
  large as its pairs, are removed once its pairs are written, where they
  stayed in the temporary directory until the end of the run (`-N` keeps
  them, in `pipeline` and `network`). The map and the network are those
  of the sorted pairs. Python API: `metator.align.get_contact_pairs` takes
  `sort_pairs` and `keep_temp` (by default, it removes those tables).
- `qc` logs the pairs within contigs, within MAGs and between MAGs that
  its noise contact ratio and noise score are computed from, from the
  pairs files as from a map.
- The heatmap of the network of the final clusters
  (`plot/network_heatmap.png`) is drawn from the contact map of a
  `pipeline` run, whatever it started from: a run from a map drew none.
  The contacts of a bin of the map, with another bin or within itself,
  are spread over the 50 kb windows that it covers, as a map holds no
  read positions; those within a window are left out, as from the pairs.
  A heatmap that cannot be drawn is a warning, and the run goes on: it
  stopped the run. The image is saved at about one pixel per window,
  where it was saved at 2000 dpi whatever its size, with lines at least
  a pixel wide between the clusters: on two libraries of MK01 (3,263
  windows), the matrix takes about 5 s where the pairs files took 51 s,
  and the heatmap step peaks at 1.8 GB of memory where the image at 2000
  dpi took it to 5.8 GB. The SCN normalization makes two copies of the
  matrix, where hicstuff's made several at each iteration. `metator
  validation MAP` draws it from the map; without a map, `metator
  validation` draws it from the pairs files of its output directory.
  Python API: `metator.figures.build_matrix_map`, `heatmap_windows`,
  `normalize_scn`, `plot_heatmap`; `plot_figures` and
  `metator.validation.recursive_decontamination` take the map
  (`metacool`).
- A `pipeline` run from a `.metacool` map updates that map
  in place, where it wrote a new map at `<outdir>/<prefix>.metacool`: the
  map of the clusters is written next to it, validated in full, then
  renamed over it (a symbolic link is replaced, not its target), so a run
  that stops, for any reason, leaves it as it was. `--output-map FILE`
  writes it elsewhere and leaves the input map alone; a run from an entry
  of a `.mmetacool` collection needs it. While a run writes a map,
  `<map>.metator.lock` holds it: a second run that would write the same map stops
  before any work, naming the run that holds it, and a lock left by a run
  that died on the same host is removed with a warning; `metator metacool`
  locks the map it writes too. So the directory of the map must be
  writable, and a map whose directory is not is refused before any work:
  write it elsewhere with `--output-map`. A map updated in place keeps its
  mode. A run that updates a map in place (`pipeline MAP`, `partition
  MAP`, `validation MAP`, `mge --map`) warns, before any work, of the
  clusters and the clustering rounds of an earlier run (or of another
  binner) that it replaces there, of which the map keeps no record.
  SIGTERM (a job scheduler's time limit) stops a run as Ctrl-C does: it
  removes its lock and the files of the map it was writing next to the
  map, and ends miComplete, which runs in a process group of its own,
  with the processes it started; its temporary directory stays, as after
  any failure. A SIGTERM or a Ctrl-C that Python would drop, as a pool of
  processes starts or ends, kills the run at once, which leaves the map
  as it was. `--output-map` refuses, before any work: `--no-metacool`,
  which writes no map; an entry of a collection, a `.mmetacool` file, a
  directory, and the collection of an input entry. Only the runs that
  write a map take its lock: a read of a map (`network MAP`, `qc MAP`,
  `contactmap MAP`, the heatmap, the input map of a run with
  `--output-map`) stops with an error when another run replaces it, or
  writes it, as it reads it: metacooler reads a map one chunk of pixels at
  a time, and such a read mixed the bins of one map with the pixels of
  the other, which gave wrong contacts. Python API:
  `metator.metacool.unchanged_while_read`.
- The map of the clusters of a run, and the map of `metator metacool`, are
  written through a fresh file: the edits that make it (metacooler's
  relabelling, `annotate`, the clustering rounds) left the space of what
  they replaced in it, which HDF5 does not reclaim, so a map updated in
  place grew with each run (from 129 kB to 272 kB over four updates of the
  map of the synthetic test data; within 8 kB now, the history of the
  runs). The map of the contigs that a stopped run leaves is not
  rewritten.
- **Breaking:** a map of the clusters that cannot be written stops
  `pipeline` with the reason, logged as an error (it was a warning, and
  the run went on to exit 0). The other outputs of the run are kept, and
  so is the map it would have replaced.
- Validation runs Prodigal and HMMER through pyrodigal and pyhmmer, their
  Python bindings, which pip installs with metaTOR, in place of the
  prodigal and hmmsearch programs that miComplete runs
  (`metator._micomplete`): prodigal and hmmer no longer need installing,
  with pip, with conda (`metator.yaml`) or in the Docker image. miComplete
  gives the bins of the tests the rows that prodigal 2.6.3 and HMMER 3.4
  gave them, and a bin of less than 20 kb still gets no gene, as prodigal
  finds none in it.
- `pipeline` from reads, alignments or pairs reads its pairs files into
  its `.metacool` map once, as soon as it has them: the map of the
  contigs, every contig in `__unbinned__`, with their contacts, the
  contact statistics measured on the reads and the contig columns measured
  on the assembly. At the end, it writes the map of its bins over it with
  `metacooler reassign`, where `write_metacool` read every pairs file
  again: the map is the one it wrote before, and that `metator metacool`
  writes from the pairs, but for its history (the `annotate` and
  `reassign` keys of its metadata) and the type of its counts with several
  libraries (64-bit integers, which `metator metacool` stores as 32-bit
  ones).
  `--no-metacool` builds it all the same and removes it at the end. An
  earlier run's map at its path is removed before the build,
  `--no-metacool` or not (it was removed at the end, and left there with
  `--no-metacool`): an input of the run that is that map, and from pairs a
  contact map, are refused before any work. A map of the contigs that
  cannot be built stops the run before its network, `--no-metacool` or
  not, whatever the number of libraries, with the reason logged as an
  error, and leaves none; a map of the bins that cannot be written over
  it, and a run that stops with an error once it is built, leave the map
  of the contigs, in no bin, which a new run can start from
  (`metator pipeline MAP`, with `-F` in the same output directory). Python
  API: `metator.metacool.write_contig_map`.
- `pipeline` checks what the map of the contigs needs before any work. A
  run that is not from a map refuses an assembly whose contig names are
  not ASCII, which it gives, before it builds its index or aligns a read:
  the map, which needs ASCII names, could only be built once the reads
  were aligned. From pairs, it refuses an empty pairs file, and one that
  metacooler, which builds the map, or metaTOR, which reads the pairs files
  of `--scaffold`, could not read (both go by the first bytes of the file,
  whatever its name): a zip archive and, with `--scaffold`, a zstd or LZ4
  file, which the scaffolding cannot read: without it, they work, as the
  networks come from the maps. A map that metacooler cannot build from the
  pairs files of a run from pairs stops the run with what it needs of them:
  the same number of fields on every row (their header lines may come
  anywhere, so files of one assembly concatenated are one library), and one
  file per library, given comma-separated.
- `pipeline` from reads, alignments or pairs with one library reads
  `network.txt` and `contig_data_network.txt` from its map of the
  contigs, as a run from a map does, instead of parsing its pairs files
  again: the same files, byte for byte. It still writes `network_0.txt`,
  the network of its library, as a copy of `network.txt`, the bytes its
  pairs files give. A run with several libraries reads them from one map
  per library, in a `.mmetacool` collection (see Added). A pair on a
  contig that the assembly lacks is left out, as metacooler warns, where
  the run stopped on a `KeyError`, and pairs that make no contact stop
  the run with an error that says so, where it stopped on a
  `ZeroDivisionError`. Its log names the map ("Information of
  <outdir>/metator.metacool:") instead of the pairs file, with the same
  counts, and no longer has the line that counted the pairs with a read
  past the end of their contig ("N pairs with a read past the end of its
  contig are left out."): metacooler's warning, which stays, counts them
  with the other pairs it leaves out.
- The temporary directory of a command is named from the random numbers
  of the system (`metator.io.generate_temp_dir`), not from Python's
  `random` module, which the partition draws from (igraph): the clusters
  of a run no longer depend on how many of its temporary directories it
  made, so that the steps run one by one in one process, after one seed,
  draw the stream that `pipeline` draws (the alignment of hicstuff still
  draws from it). A run seeded as before gives other clusters.
- metaTOR requires metacooler 0.4.0, whose `/clusters` group holds the
  clustering rounds of the map. The `.metacool` maps it writes also
  store the contact statistics of each contig (the `qc_intra_contig_short`
  and `qc_intra_contig_long` columns of the contigs table), which
  `metacooler reassign` carries over to the bins.
- Python API: `metator.partition.partition` no longer takes
  `temp_directory`, nor `generate_fasta` and
  `metator.validation.update_contigs_data_recursive` `tmpdir`, nor
  `metator.mge.mge_binning`, `generate_mges_fasta` and
  `MetatorObject.write_fasta` `tmp_dir`: they only held the contig lists
  of the `pyfastx` command. `partition` and `mge` write no temporary file
  any more.
- The network of pairs files counts the pairs as metacooler counts them
  in a map: a pair with a read on the last base of a contig, or at
  position 0 (which a 1-based pairs file should not hold, but some do), is
  a contact, which `--edge 0`, the default, left out; a pair with a read
  past the end of its contig is none, and its hits no longer count. The
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
- **Breaking:** `--edge` (`-E`) of `network` and `pipeline`: the network
  counts every pair, as its contact map does, which holds the contacts of
  bins, not the positions of the reads (see Changed for the pairs that
  `--edge 0`, the default, left out). A command line that still gives it
  is refused with this reason. See the
  [migration guide](migration.md).
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
- The commands that write a map refuse an output path that holds `::`,
  which metacooler reads as `file::group`: `metator metacool -o
  'X.metacool::sub'`, without `-F`, wrote the map into a group of
  `X.metacool`, failed, and destroyed the map that `X.metacool` held (a
  partial group was all that was left of it), as `-o 'X.metacool::/' -F`
  did. A map given as
  `file::/`, the root of the file, is the file itself, and a map whose
  file name holds `::` is refused, as metacooler would read a group of
  it; so is a `::` in the `-o` or `--prefix` of a `pipeline` run, which
  names its map.
- `validation -F` refuses clusters to validate (`--fasta`) in a directory
  that `-F` removes (`final_clusters/`, or a directory of 1.x such as
  `final_bin/`), before it removes anything. It removed them, then failed.
- `network` and `pipeline --start bam --aligner bwa` no longer fail without
  `--reverse` (`AttributeError: 'NoneType' object has no attribute
  'split'`): a BAM file of bwa holds both mates. A `pipeline` run given
  no input, neither `--forward` nor a map, failed with the same error: it
  is refused, with what to give.
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
  API: `metator.validation.write_bins_contigs` writes the `Final_cluster` of the
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
- The `recursive` column of `cluster_summary.txt` says whether the bin comes
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
- Python API: `metator.network.alignment_to_contacts` with
  `self_contacts` leaves out the contacts between contigs within `edge`
  of a contig end, as without `self_contacts`: `edge` was ignored with
  it.
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
- `metator metacool`, and the map that `pipeline` writes, no longer leave
  metacooler's partial file (`<outfile>.partial`) next to the outfile
  when metacooler cannot write the map: metacooler removes it and says so
  at INFO, where it logged, at ERROR, that it left it for inspection.
- `metator metacool` and the map of `pipeline` refuse a
  `contig_data_final.txt` whose rows have more or fewer fields than its
  header, which they read shifted or padded, without a word, and one
  compressed otherwise than with gzip as `.gz`, which metacooler failed to
  read once the map was written.
- `pipeline` checks its arguments and its assembly before replacing the
  bins of an earlier run with `--force`: a refused run no longer empties
  `overlapping_clusters/`, `recursive_clusters/` and `final_clusters/` first.
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
