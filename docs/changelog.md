# Change Log

All notable changes to this project will be documented in this file.

## Unreleased

### Added
- `metator dump MAP -o DIR` writes the tables of a run from the map of its
  clusters, read only: `contig_data_final.txt` (one row per contig of the
  map, in its order, without the `ID` column), `cluster_summary.txt` (the
  run's, from the genomes table of the map, as its `validation` entry
  records it: `summary_rows`, `summary_columns` and `summary_formats`),
  `binning.txt`, `clusters.txt` and `cluster_members.txt` (every clustering
  round of the map, the rounds `mge` and `mge_random` of `metator mge`
  included), and `network.txt` (the network that `metator network` stored,
  by contig names: 1.x wrote their IDs). The values are those of the tables that the
  run keeps in its temporary directory with `-N`; the rows of
  `contig_data_final.txt`, `binning.txt` and `cluster_members.txt` come in
  the order of the map of the clusters. A map without the clusters of a
  validation (the `validation` entry of its `metator` metadata, and the
  round `final`), validated by a development build (without
  `summary_rows` and `summary_formats`), whose clusters were edited after
  the run (by `metacooler reassign`), whose network is not the one its
  clusters were made from (`metator network` ran on it after the run: the
  `network_digest` of its `partition` entry), or that changes while it is
  read, is refused, and so is an output directory that holds the tables of
  an earlier dump, without `-F`. The tables are replaced all at once: a
  dump that fails leaves them as they were. Python API: `metator.dump`; `metator.clusters.write_round_tables`,
  which `Rounds.write` calls too, and `metator.metacool.map_rounds`, the
  rounds of a map as `Rounds` gives them; `stored_network`,
  `map_contig_table` and `stored_weights` of `metator.network` read through
  a map already opened (`mc`); `metator.metacool.MapReplacedError`.
- metaTOR reads back the numbers of the tables it wrote exactly (pandas'
  `round_trip` parser, under pandas 2; pandas 3 reads them exactly by
  default): the clustering rounds, miComplete's results and the contig
  data of the partition. pandas 2's
  default parser could read the last digit of one off, and the map would
  store another number than the run's tables.
- **Breaking:** what `validation` finds lives in the map, and `pipeline`
  runs the three steps on one map. `validation` and `pipeline` write the
  FASTA files of the final clusters in the output directory
  (`final_clusters/`) and store the rest in the map: the final cluster of
  each contig, the summary of the clusters in the genomes table, every
  clustering round, and the `validation` entry of the `metator` key of the
  metadata (version, command line, parameters, normalization, what it read
  of the map, `read`, and how `cluster_summary.txt` is written from the
  genomes table: `summary_rows`, the cluster names in the order of its
  lines; `summary_columns`, its columns in order; `summary_formats`, the
  printf format of each column, `%.2f`, `%.4f` for the abundances, `%d` or
  `%s`, a missing value an empty field). The map alone writes the summary
  again, byte for byte. The tables of the validation
  (`contig_data_final.txt`, `cluster_summary.txt`, `binning.txt`,
  `clusters.txt`, `cluster_members.txt`, miComplete's results) and
  `recursive_clusters/` are written in the temporary directory of the run,
  which `-N` keeps and the log names, not in the output directory;
  `metator dump` writes the tables from the map. An earlier run's tables
  and `recursive_clusters/` in the output directory are refused without
  `-F`, before any work, and `-F` removes them: left there, they passed for
  the run's. The map of the clusters is written anew
  from the input map by `metacooler reassign` (what its regroup does not
  carry, an extra group or attribute, is left out with its warning), then
  the network is stored last in the contigs and bins tables: `pipeline`
  stores what the three steps store one by one, in one edit session, and
  the steps one by one still give its map, `validation` storing again the
  network that `network` stored. The FASTA files of the recursive clusters
  are written by `metacooler.write_fasta`, as those of the overlapping
  ones, from an assembly decompressed once per run when it is compressed:
  the validation no longer indexes the assembly with pyfastx, and no longer
  removes an `<assembly>.fxi`. `pipeline --no-metacool` is gone. Every
  step's entry of the `metator` key records what it read of the map under
  `read` (`columns_from_map` of the `network` entry is now its
  `read.contigs`). Python API: `metator.metacool.write_map`, the one writer
  of the maps of the steps, which refuses clusters given with rounds or
  `round_names`, and a clustered contig that the map lacks; `FinalClusters`;
  `Network.entry` None keeps the map's entry; `CLUSTER_COLUMNS` and
  `METATOR_TABLE_COLUMNS`, metaTOR's own columns, which a run drops from
  the contigs table of its map, 1.x names included;
  `metator.io.SummaryLayout`, which `write_bin_summary` returns;
  `metator.validation.write_recursive_fasta`, `name_recursive_clusters` and
  `plain_assembly`; `correct_final_bin(contigs_data, written, bin_summary)`
  takes the contigs written in each FASTA file (no directory), and refuses
  one that the contig data lacks; `recursive_clustering` returns the
  contigs of the new clusters too; `recursive_decontamination` takes
  `metacool=` (the map) and `written=` (the contigs of each FASTA file
  written), and returns the contig data and the `SummaryLayout` of the
  summary; `write_overlapping_fasta` returns the contigs of each file;
  `metator.cli.check_final_clusters`, `refuse_earlier_outputs`,
  `remove_earlier_outputs` and `EARLIER_OUTPUTS`.
- **Breaking:** the map is the network that the steps share: `network MAP`
  stores the network of the contigs in the map, which `partition` and
  `validation` read. It stores what it measured on each contig in the
  contigs table of the map (`GC_content`, `Hit`, and `Shotgun_coverage` and
  `Restriction_site` when it has them), the normalization of the network as
  the weights `metator_<normalization>` of its bins (none for `None`), and
  its command line, parameters, normalization and the columns of the map
  that it read, under the `network` entry of the `metator` key of the
  metadata of the map; they replace those of an earlier network. The map is
  updated in place, under its lock, or written to the new `--output-map`
  (`-F` to replace it), which an entry of a collection needs. `pipeline`
  stores the network of its contigs in the map of its clusters as `network`
  does, so that the steps run one by one on a map give the map of
  `pipeline`, normalization included. Python API:
  `metator.metacool.store_network` and `Network`, and `network` of
  `reassign_metacool`; `metator.network.contig_contacts`, `contig_hits`,
  `log_contacts`, `contig_data_table`, `contig_table`, `map_contig_table`,
  `map_network`, `network_weights`, `network_normalization`,
  `stored_network` and `StoredNetwork`, `stored_weights` and
  `relabel_network`; `metator.partition.network_graph`.
- The map of the clusters of a `pipeline` run stores the normalization of its
  network (`--normalization`, `empirical_hit` by default), as the weights
  `metator_<normalization>` of its bins, the way cooler stores balancing
  weights: `cooler` reads them (`matrix(balance="metator_empirical_hit")`),
  and so do `metacooler show`, `dump` and `network` with `--balance`. The
  contacts of the map normalized by them and summed by pair of contigs are
  the edges of the network, but those of weight 0. They are computed from the `Hit`,
  `Shotgun_coverage` and `Restriction_site` columns of its contigs table,
  which the column records (`depends_on`), so that metacooler drops it when
  an edit changes one of them. A map without a column that the
  normalization reads (`Restriction_site` without `--enzyme`) stores none,
  with a warning; a run from the map of an earlier run replaces the
  normalization of that run. A normalization that metacooler fails to store
  fails the map, which is then not written, with the reason. Python API:
  `metator.metacool.store_normalization`,
  `metator.network.normalization_weights`, and `normalization` of
  `metator.metacool.reassign_metacool`.
- The clustering rounds of a run: the validation writes `clusters.txt`,
  every cluster of every round (`core`, `overlapping`, `recursive_<k>`, `final`)
  with the round and name of its parent and miComplete's columns, and
  `cluster_members.txt`, the contigs of each. The map of the clusters
  stores both in its `/clusters` group (metacooler 0.4.0), so that the
  intermediate clusters are kept without being written into its genomes
  and contigs tables, which hold the final ones. Python API:
  `metator.clusters`, and `rounds` of `metator.metacool.reassign_metacool`,
  which stores them (without them, the map keeps none).
- `pipeline --output-map FILE` writes the map of the clusters to FILE
  (replaced only with `-F`), and leaves the map of the run as it was.
- `network` and `pipeline` start from a `.metacool` contact map, their
  first argument (`metator pipeline MAP -a ASSEMBLY`: a file, or an entry
  `file.mmetacool::entries/<name>`), at any bin size: a map that hicstuff
  makes from the reads, or metacooler from pairs files (see Removed). The
  contacts of the network are those of the pairs the map was made of, as
  both count every pair (see Changed). A map may lack contigs of the
  assembly, which get no contacts; a map of other contigs or lengths, or
  that stores both triangles, is refused, and so is a first argument that
  is a file but no map. Python API: `metator.network.contig_contacts`.
- A `pipeline` run from a map writes the map of its bins from it, with
  `metacooler reassign`: its bins and contacts, grouped by final bin,
  with this run's genome and contig columns (the contig columns last,
  whatever columns the input map carried),
  and the contact statistics of the input map, carried over to the bins
  (a map written by metacooler 0.2.0 has them computed again from the
  pixels: at one bin per contig, every intra-contig contact then counts
  as short-range). The output may be the input map, replaced only once
  the new map is validated. Python API:
  `metator.metacool.reassign_metacool`.
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
  digits of the contig data files and of a network normalized by it.
  Python API:
  `metator.network.contig_data_from_map`, `map_columns`,
  `check_map_contigs`.
- `partition` and `validation` take the map of a run step by step as
  their first argument too (`metator partition MAP ...`, a file or an
  entry of a collection), and update it as `pipeline MAP` does: the steps
  run one by one on a map (`network MAP`, `partition MAP`, `validation
  MAP`) write the map that `pipeline MAP` writes from the same clusters,
  but for the command line and parameters it stores (Louvain draws random
  numbers: steps run as separate commands differ from `pipeline` as two
  `pipeline` runs do). `partition MAP` stores its partition, every core
  cluster (`core_00001`...) and every overlapping cluster, with their
  contigs, as the rounds `core` and `overlapping` of the clustering rounds
  of the map, in place of those it held, and what the run was (its command
  line, parameters, and the normalization and the identity of the network
  it partitioned, `network_digest`: the SHA-256 of what the network is
  computed from, by contig name, which a regroup of the map leaves as it
  was) as the `partition` entry of the `metator` key of its metadata, and
  leaves its genomes and contigs tables and the rest of its metadata as
  they were (an earlier run's, until `validation MAP` writes the new ones);
  the rounds record the step that made them. `validation MAP` reads the
  partition there, and writes the map of the clusters, with every round,
  and keeps the `partition` entry, as `pipeline MAP` stores it.
  `network MAP` keeps the clusters of the map, its partition among them,
  and warns when the network that it stores is not the one the partition
  was made on.
  Each writes the map in place (through metacooler: staged, validated in
  full, then renamed over it, the mode kept; through a symbolic link, the
  file it points to) or to `--output-map FILE` (replaced only with `-F`), and
  holds its lock while it runs. Refused before any work, the reason logged as an error, and
  before `-F` replaces the directories of clusters of an earlier run:
  what `pipeline MAP` refuses (a map that is not there, is no map, or that
  metacooler cannot open, such as a `.mmetacool` collection, an entry of a
  collection without `--output-map`, an `--output-map` that exists without
  `-F`, a map that another run holds), a map without the partition that
  `partition MAP` stores, or whose partition was made on another network
  than the one stored in it (`validation`: run `partition MAP` again), an
  `--output-map` that is one of their inputs, `--output-map` without a
  map, and an invalid `--prefix`;
  `pipeline` refuses an `--output-map` that is its assembly or depth file,
  which the map replaced. Python API: `metator.metacool.write_rounds` (and
  its `entries`), `read_partition` (and its `network`), `holds_partition`,
  `partition_network` and `NETWORK_DIGEST`; `metator.network.network_digest`
  and `StoredNetwork.digest`; `metator.clusters.Partition`,
  `partition_rounds`, `record_partition`, `core_name`, `overlapping_name`
  and `cluster_number`, and `Rounds.tables`; the `names` of
  `metator.metacool.store_rounds` (a round with no cluster), and the
  `partition` of `reassign_metacool`; `metator.partition.partition_table`;
  `metator.validation.write_overlapping_fasta`.
- **Breaking:** `mge` reads the contact map of the run, `MAP`, its first
  argument, which it requires (`metator mge MAP -m MGES -a ASSEMBLY`: a map
  that `metator network` stored its network in, after `network`,
  `pipeline` or `validation`, or an entry of a collection with
  `--output-map`): the contacts between the MGE contigs are those of the
  map, summed by pair of contigs, whatever its bin size, and those within
  each contig, between reads more than 1 kb apart, are the contact
  statistics of the map (metacooler's QC block, counted from the pairs at
  1001 bp: `qc_intra_contig_long`). The size, GC content and hits of each
  contig are those that `metator network` stored in the map; `mge` reads no
  cluster of it. A map without those contact statistics, with them
  measured from its pixels or at another distance than 1001 bp, or without
  the columns of `network`, is refused before anything is written, with
  how to make one (after a run with `--output-map`: the map that the run
  wrote); so is a map that holds none of the MGE contigs, and an assembly
  that is not a FASTA file, plain or compressed with gzip, as `network`,
  `validation` and `pipeline` refuse it. An MGE contig that the map lacks
  is left out, with a warning that names it. On the mock
  community, at one bin per contig, 1 kb and 5 kb, the contacts, MGE
  clusters and scores are those of the pairs. A pair within a contig with
  a read at position 0 (which a 1-based pairs file should not hold) and its
  mate 1001 bp away is no longer counted: metacooler counts that read on
  the first base of its contig, 1000 bp from its mate.
- **Breaking:** `mge` takes the MGE contigs the longest first, ties by
  name, whatever the order of the map, where it took them in the order of
  the rows of the `contig_data_final.txt` that it was given: that of the
  map that the validation read, in the table of the run (which `-N`
  keeps), or that of the map of the clusters, in the table that
  `metator dump` writes, which numbered some clusters otherwise. Two equal
  scores are merged in that order, and the clusters of one contig numbered
  in it. The tables of a run from a map of the contigs
  that hicstuff or metacooler made, longest first, gave the same order but
  for contigs of the same length; those of 1.x, in the order of the
  assembly, gave it for an assembly sorted by length (SPAdes). From another
  order, two clusters of one contig may swap names, and a tie of scores
  may merge other contigs. The contigs of each record of `mges_binned.fa`
  come in that order too.
- `mge` stores its MGE clusters in the map, as the round `mge` of its
  clustering rounds: one cluster per MGE cluster, named as its record of
  `mges_binned.fa`, with its contigs and its columns (its length and number
  of contigs computed by metacooler, its GC content, hits and clustering
  score), after the other rounds, which are kept, in place of the round
  `mge` of an earlier run; `metator dump` writes it with them, in
  `clusters.txt` and `cluster_members.txt`. The `mge` entry of the
  `metator` key of the metadata of the map records the version, command
  line and parameters of the run, and what it read of the map under `read`
  (the `GC_content`, `Hit` and `qc_intra_contig_long` columns of its
  contigs table, and the `qc` key of its metadata). A column of whole
  numbers, or of booleans, of the rounds of the binning run, which the
  round `mge` lacks, is stored as decimals, or as text, from then on, as
  any column with a missing value is. The genomes and contigs tables and
  the rest of the map are left as they were. The map is written as
  `partition MAP` writes it, in place or to `--output-map FILE` (replaced
  only with `-F`, new to `mge`), and locked while it is written; what
  `partition MAP` refuses is refused before any work (an `--output-map`
  that is one of the inputs of the run), and so is the summary of the MGE
  clusters that an earlier version wrote in the output directory
  (`mges_bin_summary.tsv`, or `mges_cluster_summary.tsv` of a development
  build), without `-F`, which removes it: left there, it passed for the
  run's. `--random` stores its random clustering as the round
  `mge_random`, with its own entry, and writes `mges_binned_random.fa`:
  the round `mge` and `mges_binned.fa` are left as they were. A run that
  replaces the clusters of the map (`partition MAP`, `validation MAP`,
  `pipeline MAP`) drops the rounds `mge` and `mge_random`, and their
  entries, with a warning, as it drops every round that it does not make:
  run `mge` again on its map. Python API: `metator.mge.check_mge_map`,
  `mge_contigs`, `map_matrix` and `write_mges_fasta`;
  `metator.mge.mge_binning` takes the map and the MGE contigs, and returns
  the summary and the clusters; the `mge` of `metator.metacool.write_map`
  stores the entry of its round in the metadata too;
  `metator.clusters.mge_round` takes a summary, and the name of the round;
  `metator.clusters.MGE_RANDOM`; `metator.cli.MGE_EARLIER_OUTPUTS`, and the
  `outputs` of `refuse_earlier_outputs` and `remove_earlier_outputs`.

### Changed
- **Breaking:** `partition` and `validation` require the map, `MAP`, and read
  the network that `metator network` stored in it, and `partition` the data of
  its contigs too. They read it in full before they replace the outputs of an
  earlier run (`-F`), and refuse a map without it, or whose network an edit
  changed since (a column, a value or the weights dropped, or contacts whose
  hits are no longer its `Hit` column, as after `metacooler merge`, which
  keeps the contigs table of the first map), with the step to run again.
  `partition` loses `-c, --contigs` and `-n, --network`,
  `validation` loses `-n, --network`; `network` loses `-t, --threads`, as
  there is no network file to sort. The network is in memory: each edge
  weighs the contacts of two contigs times the weights of the normalization
  of each, `n × w_i × w_j`, which differ from `n / √(h_i × h_j)` in their last
  bits, and the edges of weight 0 are left out.
- **Breaking:** the map is the partition that `validation` checks:
  `partition` stores it in the map (see Added), and writes neither
  `contig_data_partition.txt` nor `overlapping_clusters/`; `validation`
  reads it there, and writes the FASTA files of the overlapping clusters
  that reach its `--size`, through metacooler (`metacooler.write_fasta`,
  the records as they are in the assembly, in the order of the partition),
  and the contig data table of the partition, in its temporary directory
  (`-N` keeps them), which needs room for them (at most the size of the
  assembly) and, from a gzip-compressed assembly, for its decompressed
  copy, which metacooler writes there while it writes them (pyfastx read
  it compressed); and so does `pipeline`, which wrote
  `overlapping_clusters/` and `contig_data_partition.txt` in its output
  directory. `partition` loses `-a, --assembly` and `-s, --size`: it reads no
  assembly, and keeps every cluster; `validation` loses `-c, --contigs` and
  `-f, --fasta`, and refuses a map without the partition (the rounds `core`
  and `overlapping`, their clusters named as `partition` names them),
  naming the step to run. Its `--size` chooses the overlapping clusters that
  miComplete checks, as that of `pipeline` does. The round `overlapping` of
  the map of the clusters holds every overlapping cluster, those that
  miComplete checked with its columns, after the round `core`; so do
  `clusters.txt` and `cluster_members.txt`. Python API:
  `metator.partition.partition(algorithm, network, iterations,
  overlapping_parameter, resolution_parameter, threads, prefix)` returns
  the partition, a `metator.clusters.Partition`, and writes nothing: it
  took the assembly, the contig data, the output and FASTA directories and
  the size, and returned the path of the table that it wrote;
  `metator.validation.recursive_decontamination` takes the partition
  (`partition`, after `prefix`, required), whose core and overlapping
  clusters it records as the first rounds of the run.
- **Breaking:** the contigs of a run are those of its map, in its order
  (longest first, as hicstuff and metacooler make it), and the ID of a contig
  is its position in the map: the rows of `contig_data_final.txt`, the
  numbers of the clusters and the order of the
  records of their FASTA files follow it, where they followed the
  assembly. The contigs of the assembly that the map lacks are in no table
  of the run, and the shotgun abundance of a cluster is its share of the
  sequenced bases of the contigs of the map: the clusters of a map made with
  a contig filter add up to 100%, where 1.x counted every contig of the
  assembly. `metator.network.create_contig_data`
  gives no `id` and `hit` of a contig: its ID is its position in the map,
  and its hits are counted on the map. `validation` numbers them in
  the order of the map as it reads it: a validation run again on the map of
  the clusters, whose contigs are grouped by cluster, numbers them in that
  order.
- `validation MAP` stores the normalization of the network that `network`
  stored in the map in the map of its clusters, where it dropped it.
- The columns of the map that a run read (`columns_from_map`) are recorded
  under the `network` entry of the `metator` key of its metadata, where
  `pipeline` recorded them under the `metator` key.
- **Breaking:** metaTOR calls its groups of contigs clusters, as its contact
  map does: the bins are those of the map. The outputs of a run are renamed:
  `bin_summary.txt` is `cluster_summary.txt`; `recursive_bin/` is
  `recursive_clusters/`, and `overlapping_bin/` is gone (the partition is
  in the map);
  `final_bin/` (`validation`) and `final_bin_unscaffold/` (`pipeline`) are
  `final_clusters/`; `mges_bin_summary.tsv` (`mge`) is gone (see
  Removed). In the contig data tables, `Core_bin_ID`,
  `Core_bin_contigs` and `Core_bin_size`, and the columns of the overlapping
  and recursive bins, are `Core_cluster_ID`, `Core_cluster_contigs`,
  `Core_cluster_size` and so on, and `Final_bin` is `Final_cluster`.
  `binning.txt`, which anvi'o imports, and the names of the clusters are
  kept: see [Migrating to metaTOR 2.0](migration.md). No step renames the
  columns of a contig data table of 1.x, and `partition` and `validation`
  read none; and `metator.host`
  still reads the tables of 1.x (`metator.io.rename_1x_columns`). A run in the
  output directory of a run of 1.x leaves the outputs that 1.x named after
  bins as they are, with or without `-F`. metacooler's
  `--contigs-format metator` names `Final_bin`: give `--contigs-format
  Name,Final_cluster` for a table of 2.0. The help, the logs and the docs
  call them clusters too, and the help groups the three steps and `pipeline`
  under Clustering.
- A `pipeline` run from a `.metacool` map updates that map in place, where
  it wrote a new map at `<outdir>/<prefix>.metacool`. Every map is written
  through metacooler (0.4.0), in an edit session (`metacooler.staged_edit`):
  the edits of a step go into one copy of the map, staged in a private
  directory next to it (`.<map>.<uuid>.<n>.d/`), which metacooler validates
  in full, then renames over the map, so a run that stops, for any reason,
  leaves it as it was; the copy starts fresh, so the map holds at most the
  space that the edits of the last run left behind, and a map updated run
  after run does not grow. A map given through a symbolic link
  is updated where the link points, and the link keeps pointing at it. The
  map keeps its mode, owner, group and extended attributes where the run may
  set them, and a warning names those it cannot keep. `--output-map FILE`
  writes it elsewhere and leaves the input map alone, unless FILE names the
  file of the input map, which the run then updates in place: a symbolic
  link at FILE is replaced (with `-F`, as any FILE that exists), not its
  target, and a hard link of the input map there (the same file under
  another name) is refused before any work. A run from an entry of a
  `.mmetacool` collection needs `--output-map`. A run holds metacooler's
  lock of the map it writes, `<map>.lock` (FILE's own, with
  `--output-map`), from its checks to its end, and reads the map under it:
  a second run, or a metacooler command that writes the map in place
  (`annotate --inplace`, `normalize`...), stops before any work, naming the
  process that holds it; runs from one map, each with its own
  `--output-map`, run side by side. So the directory of the map must be
  writable, with room, beside the map, for its copy and for the new map that
  `metacooler reassign` writes from it, plus cooler's temporary files, and a
  map whose directory is not, or a map without write permission, is refused
  before any work: write it elsewhere with
  `--output-map`. A map that changed since the run read it (another run
  replaced it, or a program that takes no metacooler lock, `cooler balance`
  or a script, wrote it) is refused when the run writes its map, rather
  than mixed with it: the run stops with the reason. Python API:
  `reassign_metacool`, `write_rounds` and `write_round` of
  `metator.metacool` take `planned=`, the map as the run read it, which
  they refuse to write from if it changed since; `is_input_map` says
  whether a write updates the map in place, `is_hard_link_of` refuses the
  other names of its file, and `make_map_directory` makes the directory of
  a map. A run that updates a map in
  place (`pipeline MAP`, `partition MAP`, `validation MAP`, `mge MAP`)
  warns, before any work, of the clusters and the clustering rounds of an
  earlier run (or of another binner) that it replaces there, of which the
  map keeps no record. `--output-map` refuses, before any work:
  `--no-metacool`, which writes no map; an entry of a collection, a
  directory, and the collection of an input entry. Only the runs that write
  a map take its lock: a read of a map (`network MAP`, the input map of a
  run with `--output-map`, and the checks of a run that read its contigs
  and its columns) stops with an error, logged as one, when another run
  replaces it, or writes it, as it reads it: metacooler reads a map one
  chunk of pixels at a time, and checks its file before and after each, as
  such a read would mix the bins of one map with the pixels of the other.
  Python API: `metator.metacool.read_errors` raises metacooler's refusal of
  such a read (`MapChangedError`) as a `MetacoolError`.
- Neither metaTOR nor metacooler handles signals. A run stopped by Ctrl-C
  removes its staged copy and its lock, and leaves the map as it was; its
  temporary directory stays, as after any failure. A run killed by SIGTERM
  (a job scheduler's time limit, or a cancelled job) or SIGKILL leaves its
  lock and its staged copy next to the map: the next run (or metacooler
  command) that writes the map on the same host takes the lock over and
  removes what other writers left beside the map, with a warning each; a
  lock of another host, or of another machine or PID namespace, is never
  taken over: the error names the files to delete. miComplete runs in the
  process group of the run, so a Ctrl-C in the terminal, or the signal that
  a job scheduler sends to the job, reaches it too.
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
- The temporary directory of a command is named from the random numbers
  of the system (`metator.io.generate_temp_dir`), not from Python's
  `random` module, which the partition draws from (igraph): the clusters
  of a run no longer depend on how many of its temporary directories it
  made, so that the steps run one by one in one process, after one seed,
  draw the stream that `pipeline` draws. A run seeded as before gives
  other clusters.
- metaTOR requires metacooler 0.4.0, whose `/clusters` group holds the
  clustering rounds of the map. The `.metacool` maps it writes also
  store the contact statistics of each contig (the `qc_intra_contig_short`
  and `qc_intra_contig_long` columns of the contigs table), which
  `metacooler reassign` carries over to the bins.
- Python API: `metator.partition.partition` no longer takes
  `temp_directory`, nor
  `metator.validation.update_contigs_data_recursive` `tmpdir`: they only
  held the contig lists of the `pyfastx` command. `partition` writes no
  temporary file any more.
- `mges_binned.fa`: metacooler writes it, one record per MGE cluster as
  before, its contigs the longest first, ties by name (see Added), joined
  by 200 `N`, in lines of 80 bases, where the lines of each contig were
  those of the assembly. Nothing is written next to the assembly: no
  pyfastx index.
- The network counts the pairs as metacooler counts them in a map: a pair with a read on the last base of a contig, or at
  position 0 (which a 1-based pairs file should not hold, but some do), is
  a contact, which `--edge 0`, the default, left out; a pair with a read
  past the end of its contig is none, and its hits no longer count. The
  contacts of a `.metacool` map are then those of its pairs, exactly.
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
  the CLI reference against the CLI. The docstring of
  `metator.validation.correct_final_bin`, which disagreed with its code, is
  corrected.

### Removed
- `metator.metacool.reassign_metacool`, `store_network`, `write_rounds` and
  `write_round` (all now `write_map`), `read_contig_data` and
  `PANDAS_COMPRESSIONS`; `metator.validation.update_contigs_data_recursive`
  and `write_bins_contigs`; `metator.io.write_contigs`;
  `metator.cli.run_metadata`'s `columns_from_map` (see `network_entry`);
  `pipeline --no-metacool`.
- **Breaking:** the pairs route of `mge`: its pairs files (`PAIRSFILE`),
  `-c, --contigs` (`contig_data_final.txt`), `-b, --binning`
  (`binning.txt`, or an anvi'o collection, which only set columns that no
  output read) and `--map`, which the map, its first argument, replaces;
  the summary of the MGE clusters, `mges_bin_summary.tsv`
  (`mges_cluster_summary.tsv` in a development build), which the round
  `mge` of the map holds (unrounded, where the table rounded `GC` and
  `BinningScore` to two decimals), but for the three columns of the hosts
  of the MGE clusters (`MetagenomicBin`, `AssociationScore`,
  `AssociatedBins`), which `mge` always left empty; the dependency on
  pypairix, whose index of the pairs it read, and the one on pyfastx, which
  no step uses any more. Python API: `metator.mge.build_matrix` (by
  `map_matrix`), `generate_mge_bins_pairs`, `generate_mges_fasta` (by
  `write_mges_fasta`), and the `outfile` of `generate_bin_summary`, which
  returns the summary; `metator.io.import_anvio_binning`,
  `import_contig_data_mges`, `write_mge_data` and `remove_pyfastx_index`;
  `metator.metacool.check_map_of_contigs`.
- **Breaking:** `network.txt` and `contig_data_network.txt`: the network is
  in the map (see Added). `metacooler network MAP --level contig --balance
  metator_<normalization>` writes it, by contig names, with the contacts
  within each contig, which are no edge: drop them (see
  [The network](formats.md#the-network)). A network of your
  own, given to `partition -n` and `validation -n`, is gone with them. The
  network is no longer sorted with the external `sort`, nor written to
  temporary files, and the validation no longer writes the subnetwork of
  each cluster it partitions again to its temporary directory. Python API:
  `metator.network.normalize_pair`, `write_edge` and `write_contig_data`;
  `metator.partition.read_network`; `metator.io.sort_pairs`.
- **Breaking:** `contig_data_partition.txt` and `overlapping_clusters/` of
  `partition` and `pipeline`: the partition is in the map (see Changed).
  `metacooler dump MAP -t members` lists the contigs of each core and
  overlapping cluster. `validation` writes the FASTA files of the
  overlapping clusters that reach its `--size`, those that miComplete
  checks, in its temporary directory (`-N` keeps them), as
  `overlapping_clusters/` held them; `metacooler fasta MAP --fasta ASSEMBLY
  --round overlapping -o DIR` writes every overlapping cluster, whatever its
  size, its contigs in the order of the map: not the same files.
  Python API: `metator.partition.generate_fasta` and
  `update_contigs_data`, whose table `metator.partition.partition_table`
  makes; `metator.clusters.overlapping_round`;
  `metator.validation.check_partition_table` and `PARTITION_COLUMNS`.
- **Breaking:** the commands that did what hicstuff, metacooler, instaGRAAL
  or pairix do, and the plots: metaTOR clusters contigs, and nothing else.
  The [Other tools](other-tools.md) page says what replaces each one:
  - `metator qc`: `metacooler qc MAP --per-genome qc.tsv`, and hicstuff's
    `pipeline -e ENZYME -f --plot` for the 3C events of a library;
  - `metator contactmap`: `metacooler extract MAP --genomes NAME -o
    NAME.cool`;
  - `metator scaffold`, and the `-v/--scaffold` and `-J/--junctions`
    options of `pipeline`, with its `scaffold_bin/`: hicstuff's
    `pipeline -M graal` on the FASTA file of a cluster, then instaGRAAL;
  - `metator pairs`: hicstuff's `pipeline` writes its pairs sorted and
    compressed with bgzip, ready for `pairix -p pairs`;
  - `metator metacool`: `pipeline` and `validation` write the map of the
    clusters; `metacooler reassign` and `annotate` by hand;
  - the plots of `validation` and `pipeline`, `plot/` (the quality, size,
    GC content and coverage of the clusters, and the heatmap of their
    network): the log counts the clusters in each quality class, and
    `metacooler show` draws the map.

  A command line that names one of them is refused as an unknown command or
  option, and a run leaves the `plot/` of an earlier run as it is. Python
  API:
  the modules `metator.figures` (but `bin_quality`, now
  `metator.validation.bin_quality`), `metator.quality_check`,
  `metator.contact_map`, `metator.scaffold` and `metator.regions`;
  `metator.metacool.write_metacool`; and `metator.io.get_pairs_data` and
  `read_bin_summary`.
- **Breaking:** the reads, alignments and pairs inputs of `network` and
  `pipeline`: a run starts from a `.metacool` contact map of the contigs,
  `MAP`, their first argument, which they now require. hicstuff makes it
  from the reads (`hicstuff pipeline -M metacool`, from hicstuff's next
  release), and metacooler from
  pairs files (`metacooler create pairs`) or a `.cool` file (`metacooler
  create cool`); `metacooler merge` sums the maps of several libraries (see
  [Other tools](other-tools.md#from-the-reads-to-the-map) and the
  [migration guide](migration.md)). On the test reads of 1.x, the map of
  hicstuff and the map that metacooler makes from `tests_data/alignment.pairs`
  give the same network and contig data, byte for byte. Removed with them:
  - the options `-1/--forward`, `-2/--reverse`, `-S/--start`,
    `-b/--aligner`, `-B/--aligner-mode` and `-q/--min-quality`: a command
    line that gives one is refused by click as an unknown option, a usage
    error (exit status 2), as is `--edge` (below); click may suggest an
    option of a near name (`--force` for `--forward`);
  - the outputs of the alignment: `alignment_<i>.bam`,
    `alignment_<i>_for.bam` and `alignment_<i>_rev.bam`; the pairs files
    `alignment_<i>_sorted.pairs.gz` and their index (`network`, and
    `pipeline -N`); the network of each library, `network_<i>.txt`
    (`network_0.txt` included), and `hit_data_alignment.txt`. `pipeline` no
    longer writes `<outdir>/metator.metacool`: it writes the map of the
    clusters over `MAP`, or to `--output-map`, and `--prefix` names the
    clusters only;
  - `--assembly` as the basename of a bowtie2 or bwa index: it is the FASTA
    file of the assembly, plain or compressed with gzip, and one compressed
    with bzip2 or zip, which pyfastx cannot read, is refused before any work;
  - the dependencies on hicstuff and pysam, and the programs bowtie2, bwa
    and samtools: prodigal and hmmer stay, for miComplete.

  Python API: the module `metator.align`; `metator.network.alignment_to_contacts`,
  `precompute_network`, `compute_network` and `write_hit_data`, and
  `nb_alignment` and `hit_data`, of `create_contig_data`, which returns the
  contig data alone; and `metator.io.check_fasta_index`,
  `generate_fasta_index`, `retrieve_fasta`, `sort_and_index_pairs` and
  `check_pypairix`.
- Python's warnings (`warnings.warn`, which cooler gives the problems of a
  map with) are shown, and written to the log file of a run, as metaTOR's
  own messages (`WARNING :: ...`): logging dropped them once metaTOR no
  longer imported hicstuff, which printed them.
- **Breaking:** `--edge` (`-E`) of `network` and `pipeline`: the network
  counts every pair, as its contact map does, which holds the contacts of
  bins, not the positions of the reads (see Changed for the pairs that
  `--edge 0`, the default, left out). A command line that still gives it
  is refused as an unknown option. See the
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
- **Breaking:** the small features that nothing else used:
  - `--cluster-matrix` (`-C`) of `partition`, `validation` and `pipeline`,
    and the clustering matrices that it wrote, `clustering_matrix_partition.npz`,
    `clustering_matrix_recursive.npz` and `clustering_matrix.npz`: the share
    of the iterations that put each two contigs together has no replacement.
    The clustering rounds of a run (`clusters.txt`, `cluster_members.txt`
    and the `/clusters` group of the map) give the contigs of each cluster
    of each round;
  - `--self-contacts` (`-s`) of `network`: the network has no edge from a
    contig to itself. The contacts within each contig are in the map, in
    its `qc_intra_contig_short` and `qc_intra_contig_long` columns, when it
    holds metacooler's QC block (the default, but with `--no-qc`);
  - the module `metator.main` (`python -m metator.main`), the former entry
    point: run `metator`, or `python -m metator`.

  A command line that gives one of these options is refused as an unknown
  option. Python API: `metator.partition.build_clustering_matrix`; the
  `cluster_matrix` argument of `metator.partition.partition`, which returned
  the path of the clustering matrix too (it returns the partition now: see
  Changed), of `metator.validation.recursive_clustering`, which returns
  `(contamination, contigs_data, parent_dict)` (without the path of the
  matrix), and of `recursive_decontamination`, which returns nothing (it
  returned the path of the matrix); `metator.io.read_contig_data` and
  `save_sparse_matrix`; and `metator.cli.check_output_path`.

### Fixed
- `validation` reads the contig names of its tables as text: contigs named
  by numbers (`000017`, `1e5`) kept their names, where they were read as
  numbers and lost their clusters. A contig of a final FASTA file that the
  contig data lacks, final clusters that are not those of the FASTA files
  of `final_clusters/`, and a clustered contig that the map lacks are
  refused, named, rather than written into a wrong map.
- `-F` replaces a directory of clusters of an earlier run that is a
  symbolic link (`validation`, `pipeline`): the link is removed, not its
  target, where the run failed with `OSError`, some of the directories
  already replaced. A link to nothing is refused without `-F`, as a
  directory is.
- The commands that write a map refuse an output path that holds `::`,
  which metacooler reads as `file::group`. A map given as `file::/`, the
  root of the file, is the file itself, and a map whose file name holds `::`
  is refused, as metacooler would read a group of it; so is a `::` in an
  `--output-map`.
- `validation` and `pipeline` refuse a `--tmpdir` in a directory that they
  make anew, which `-F` removes (`recursive_clusters/`, `final_clusters/`),
  as a usage error, before they create or remove anything. `-F` removed it
  with them, and the run then failed on its temporary files.
- A contig named `NA` (or `N/A`, `null`...) is a contig in the validation:
  read as a missing value, it made the `Size` column of
  `contig_data_final.txt` a float, and the run wrote no map of its clusters.
- `binning.txt` holds every contig of the final bins, as
  `contig_data_final.txt` does: the contigs of a bin partitioned again but
  kept whole, which had gone to a recursive part, were left out. Python
  API: `metator.validation.write_bins_contigs` writes the `Final_cluster` of the
  contigs, and no longer takes `prefix`.
- The log of the validation counts the bins in the quality classes of
  `metator.validation.bin_quality`, the contamination being the weighted
  redundancy − 1. It counted as HQ the bins of 90% completeness or more with
  a contamination above 5%, and as MQ those of 5% or less, and divided the
  contamination by the completeness.
- `metator validation --help` and the docstrings of `metator.validation`
  give the rule of the recursive partition: a weighted completeness of 33%,
  not 50%, and the completeness its parts must keep.
- The `recursive` column of `cluster_summary.txt` says whether the bin comes
  from the recursive partition: every bin read `True` when no bin was
  partitioned again, and a bin made at an earlier recursive step `False`.
- The validation leaves no pyfastx index (`<bin>.fa.fxi`) next to the fasta
  files of the final bins.
- The depth file (`--depth`) is read by contig name: it was read line by
  line, in the order of the assembly, and a depth file in another order
  gave the contigs the length and coverage of others, without a word. A
  depth file that lacks a contig of the assembly, or gives it another
  length, stops the run with an error, and the lengths of the contigs come
  from the assembly.
- `network` and `pipeline` read an assembly compressed with gzip, which
  they refused as not being a FASTA file.
- The map that `pipeline` writes no longer leaves metacooler's partial file
  (`<outfile>.partial`) next to the outfile when metacooler cannot write the
  map: metacooler removes it and says so at INFO, where it logged, at ERROR,
  that it left it for inspection.
- The map of `pipeline` refuses a `contig_data_final.txt` whose rows have
  more or fewer fields than its header, which it read shifted or padded,
  without a word, and one compressed otherwise than with gzip as `.gz`,
  which metacooler failed to read once the map was written.
- `pipeline` checks its arguments and its assembly before replacing the
  bins of an earlier run with `--force`: a refused run no longer empties
  `recursive_clusters/` and `final_clusters/` first.
- The fasta files of the recursive and final bins (`validation`,
  `pipeline`) and of the MGE bins (`mge`) are written with pyfastx's Python API
  (`metator.io.write_contigs`), byte for byte as before. The `pyfastx`
  command wrote them through the shell: without it on the `PATH` (a virtual
  environment used without activating it), every file was empty, without an
  error. A contig missing from the assembly now stops the run with
  `KeyError`, where the command left it out of its bin.
- `pipeline` and `mge` remove the pyfastx index of the assembly only when it
  is there: they failed at their very end, with `FileNotFoundError`, when
  another run sharing the assembly had removed it.
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
