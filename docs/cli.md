# CLI reference

The `metator` command has one subcommand per task. This page documents every
command and option; `-h/--help` prints the same in a terminal.

```console
metator --help
metator <command> --help
metator --version
```

## Conventions

| Convention | Meaning |
|---|---|
| Several libraries | Comma-separated lists, in the same order in `--forward` and `--reverse`: `-1 L1_R1.fq.gz,L2_R1.fq.gz -2 L1_R2.fq.gz,L2_R2.fq.gz`. |
| Several enzymes | Separated by commas: `DpnII,HinfI`. |
| Pairs files | `qc`, `contactmap`, `metacool`, `scaffold` and `mge` take the pairs files of a run as arguments: `alignment_<i>_sorted.pairs.gz`, which `metator network` writes and `metator pipeline` keeps with `--no-clean-up`. |
| `-F, --force` | Replace the outputs of an earlier run. Without it, a command refuses to overwrite them. |
| `-T, --tmpdir` | Each run works in a directory of its own inside it, removed at the end, unless `--no-clean-up` or the run fails. In `scaffold`, `-T` is `--threshold`. |
| `-N, --no-clean-up` | Keep the temporary files. In `qc`, it is `-n`. |
| Log files | Each command writes its log, `metator_<command>_<date>.log`, into its output directory; `metacool`, `scaffold`, `pairs` and `mge`, into the current directory. |

New to metaTOR? Start with the [Quickstart](quickstart.md); for hits, core,
overlapping and recursive bins, see [Concepts](concepts.md).

## Binning

### `pipeline` { #pipeline }

Run the three steps of the binning, `network`, `partition` and `validation`,
in one command, and write the contact map of the bins. **Guide:**
[Running the pipeline](guides/pipeline.md).

```console
metator pipeline [OPTIONS]
```

**Recommended:** `--aligner bwa`. See
[Aligning the reads](guides/network.md#aligning-the-reads).

| Option | Default | Description |
|---|---|---|
| `-1, --forward` | *required* | The R1 reads (FASTQ), or their alignments (BAM), or the pairs files, of each library, separated by commas. With `--start metacool`, one `.metacool` contact map, or an entry of a `.mmetacool` collection (`file.mmetacool::entries/<name>`). |
| `-2, --reverse` | none | The R2 reads, or their alignments, in the same order. Not needed from pairs files, a contact map, or the BAM files of bwa, which hold both mates. |
| `-a, --assembly` | *required* | The assembly: a FASTA file, or the basename of its bowtie2 or bwa index. |
| `-S, --start` | `fastq` | What `--forward` and `--reverse` are: `fastq`, `bam`, `pair` or `metacool`. |
| `-d, --depth` | none | The coverage of the contigs, from `jgi_summarize_bam_contig_depths` of MetaBAT2. |
| `-b, --aligner` | `bowtie2` | `bwa`, the recommended one, or `bowtie2`. |
| `-B, --aligner-mode` | `normal` | `normal`, `iterative` or `cutsite`. `cutsite` needs `--enzyme`: without it, the reads are aligned iteratively. |
| `-e, --enzyme` | none | The restriction enzymes of the library. |
| `-q, --min-quality` | `30` | The minimum mapping quality of both mates of a pair. |
| `-m, --normalization` | `empirical_hit` | `None`, `abundance`, `length`, `RS`, `empirical_hit` or `theoritical_hit`. `RS` and `theoritical_hit` need `--enzyme`; `abundance` and `theoritical_hit` need `--depth`. See [Normalising the network](guides/network.md#normalising-the-network). |
| `-E, --edge` | `0` | Only 0 in `pipeline`, whatever the start: it counts every pair, as its contact map does (`metator network` applies `--edge` to pairs files). |
| `-A, --algorithm` | `louvain` | The partition algorithm: `louvain` or `leiden`. |
| `-i, --iterations` | `100` | The iterations of the partition. |
| `-O, --overlap` | `80` | The percentage of iterations in which two core bins must fall together to be merged into one overlapping bin. |
| `-r, --res-param` | `1.0` | The resolution of the Leiden algorithm. |
| `-s, --size` | `500000` | The size, in bp, a bin must reach to be kept. |
| `-C, --cluster-matrix` | off | Write the clustering matrices (`clustering_matrix*.npz`). |
| `-j, --rec-iter` | `10` | The iterations of the recursive partition of the contaminated bins. |
| `-P, --rec-overlap` | `90` | The overlap threshold of the recursive partition, in percent. |
| `-v, --scaffold` | off | Scaffold the final bins, into `scaffold_bin/`. Refused with `--start metacool`. |
| `-J, --junctions` | `NNNNN` | The sequence between two contigs of a scaffold. |
| `-o, --outdir` | the current directory | The output directory. |
| `-p, --prefix` | none | Name the bins `<prefix>_metator_…` and the map `<prefix>.metacool`. Letters, digits, `_` and `-` only. |
| `-F, --force` | off | Replace the directories of bins of an earlier run. |
| `-T, --tmpdir` | `./tmp` | The temporary directory, where the index of the assembly is built. |
| `-t, --threads` | `1` | Threads for the alignment, the partition and miComplete. |
| `-N, --no-clean-up` | off | Keep the temporary files, and the pairs files of the run. |
| `--no-metacool` | off | Do not keep the contact map of the contigs, grouped by bin, `metator.metacool` (`<prefix>.metacool` with `--prefix`): a run from reads, alignments or pairs still builds the map of its contigs from the pairs files, and removes it at the end. |

The validation is the step that needs the most memory: miComplete runs
prodigal and hmmer on every bin.

```bash
# Reads to bins, with the recommended aligner
metator pipeline -b bwa -e DpnII -d depth.txt -t 16 -1 R1.fq.gz -2 R2.fq.gz \
    -a assembly.fa -o binning

# Two libraries, keeping the pairs files
metator pipeline -b bwa -e DpnII -t 16 -1 L1_R1.fq.gz,L2_R1.fq.gz \
    -2 L1_R2.fq.gz,L2_R2.fq.gz -a assembly.fa -N -o binning

# From a contact map of the contigs
metator pipeline -S metacool -1 sample.metacool -a assembly.fa -e DpnII -o binning
```

### `network` { #network }

Align the reads, pair the mates, and build the network of contigs: the first
step of the pipeline. Writes `network.txt`, `contig_data_network.txt` and the
pairs files. **Guide:**
[Aligning reads and building the network](guides/network.md).

```console
metator network [OPTIONS]
```

| Option | Default | Description |
|---|---|---|
| `-1, --forward` | *required* | As in `pipeline`. |
| `-2, --reverse` | none | As in `pipeline`. |
| `-a, --assembly` | *required* | As in `pipeline`. |
| `-S, --start` | `fastq` | `fastq`, `bam`, `pair` or `metacool`. |
| `-d, --depth` | none | As in `pipeline`. |
| `-b, --aligner` | `bowtie2` | `bwa`, the recommended one, or `bowtie2`. |
| `-B, --aligner-mode` | `normal` | `normal`, `iterative` or `cutsite`. |
| `-e, --enzyme` | none | The restriction enzymes of the library. |
| `-q, --min-quality` | `30` | The minimum mapping quality of both mates of a pair. |
| `-n, --normalization` | `empirical_hit` | As `--normalization` (`-m`) in `pipeline`. |
| `-E, --edge` | `0` | Leave out the pairs with a mate within this many bp of the end of its contig. Not available with `--start metacool`. |
| `-s, --self-contacts` | off | Count the contacts of each contig with itself. |
| `-o, --outdir` | the current directory | The output directory. |
| `-T, --tmpdir` | `./tmp` | The temporary directory. |
| `-t, --threads` | `1` | Threads for the alignment. |
| `-N, --no-clean-up` | off | Keep the temporary files. |

```bash
metator network -1 R1.fq.gz -2 R2.fq.gz -a assembly.fa -b bwa -e DpnII -t 16 -o binning
```

### `partition` { #partition }

Partition the network, many times, into core and overlapping bins: the second
step. Writes `overlapping_bin/` and `contig_data_partition.txt`. **Guide:**
[Partition and validation](guides/binning.md).

```console
metator partition [OPTIONS]
```

| Option | Default | Description |
|---|---|---|
| `-a, --assembly` | *required* | The assembly, as a FASTA file. |
| `-c, --contigs` | *required* | The contig data table of the network step, `contig_data_network.txt`. |
| `-n, --network` | *required* | The network, `network.txt`. |
| `-A, --algorithm` | `louvain` | `louvain` or `leiden`. |
| `-i, --iterations` | `100` | The iterations of the partition. |
| `-O, --overlap` | `80` | The overlap threshold, in percent. |
| `-r, --res-param` | `1.0` | The resolution of the Leiden algorithm. |
| `-s, --size` | `500000` | The size, in bp, an overlapping bin must reach to be written. |
| `-C, --cluster-matrix` | off | Write `clustering_matrix_partition.npz`. |
| `-o, --outdir` | `./` | The output directory. |
| `-p, --prefix` | none | Name the bins `<prefix>_metator_…`. |
| `-F, --force` | off | Replace `overlapping_bin/`. |
| `-T, --tmpdir` | `./tmp` | The temporary directory. |
| `-t, --threads` | `1` | Threads for the partition. |
| `-N, --no-clean-up` | off | Keep the temporary files. |

```bash
metator partition -a assembly.fa -c binning/contig_data_network.txt -n binning/network.txt \
    -t 16 -o binning
```

### `validation` { #validation }

Check the bins with miComplete, and partition the contaminated ones again:
the third step. Writes `final_bin/`, `recursive_bin/`, `bin_summary.txt`,
`contig_data_final.txt` and `binning.txt`. **Guide:**
[Partition and validation](guides/binning.md#validation).

```console
metator validation [OPTIONS]
```

A bin is partitioned again when its weighted completeness reaches 33% and its
contamination, (weighted redundancy − 1) / weighted completeness, reaches 5%.
Its parts replace it when the most complete of them keeps enough of its
completeness.

| Option | Default | Description |
|---|---|---|
| `-a, --assembly` | *required* | The assembly, as a FASTA file. |
| `-c, --contigs` | *required* | The contig data table of the partition, `contig_data_partition.txt`. |
| `-f, --fasta` | *required* | The directory of the bins to check, `overlapping_bin/`. |
| `-n, --network` | *required* | The network, `network.txt`. |
| `-A, --algorithm` | `louvain` | `louvain` or `leiden`. |
| `-i, --iterations` | `10` | The iterations of the recursive partition. |
| `-O, --overlap` | `90` | The overlap threshold of the recursive partition, in percent. |
| `-r, --res-param` | `1.0` | The resolution of the Leiden algorithm. |
| `-s, --size` | `500000` | The size, in bp, a bin must reach to be kept. |
| `-C, --cluster-matrix` | off | Write `clustering_matrix_recursive.npz`. |
| `-o, --outdir` | `./` | The output directory. |
| `-p, --prefix` | none | Name the bins `<prefix>_metator_…`. |
| `-F, --force` | off | Replace `recursive_bin/` and `final_bin/`. |
| `-T, --tmpdir` | `./tmp` | The temporary directory. |
| `-t, --threads` | `1` | Threads for the partition and miComplete. |
| `-N, --no-clean-up` | off | Keep the temporary files. |

```bash
metator validation -a assembly.fa -c binning/contig_data_partition.txt -n binning/network.txt \
    -f binning/overlapping_bin -t 16 -o binning
```

## After the binning

### `metacool` { #metacool }

Write the contact map of all the contigs, grouped by bin, as a `.metacool`
file, as `pipeline` does at the end of a run. **API:**
`metator.metacool.write_metacool`. **Guide:**
[The contact map of the bins](guides/metacool.md).

```console
metator metacool [OPTIONS] PAIRSFILE...
```

`PAIRSFILE` is the pairs files of the run, whose contacts are summed, or, for
a run that started from a contact map, that `.metacool` map.

| Option | Default | Description |
|---|---|---|
| `-c, --contig-data` | *required* | `contig_data_final.txt`: the contigs, their final bins, and what metaTOR measured on them. |
| `-b, --bin-summary` | none | `bin_summary.txt`, whose columns are added to the genomes table. |
| `-o, --outfile` | *required* | The `.metacool` file to write. |
| `-n, --assembly-name` | none | The name of the assembly, stored in the file. Refused for a map written from a map, which keeps its own. |
| `-t, --threads` | `1` | Threads binning the pairs. |
| `-F, --force` | off | Replace the output file. |

```bash
metator metacool -c binning/contig_data_final.txt -b binning/bin_summary.txt \
    -o sample.metacool binning/alignment_0_sorted.pairs.gz
```

### `contactmap` { #contactmap }

Build the contact map of one contig or bin, with hicstuff. **API:**
`metator.contact_map.generate_contact_map`. **Guide:**
[Viewing and scaffolding MAGs](guides/visualisation.md#the-contact-map-of-a-bin).

```console
metator contactmap [OPTIONS] PAIRSFILE...
```

| Option | Default | Description |
|---|---|---|
| `-a, --assembly` | *required* | A FASTA file holding the contigs of the map: the assembly, or the bin. |
| `-c, --contig-data` | *required* | `contig_data_final.txt`, or a contig data table with a column of other bins. |
| `-e, --enzyme` | *required* | The restriction enzymes of the library: the map is at the resolution of their fragments. |
| `-n, --name` | *required* | The name of the contig or bin, or the numeric ID of a core or overlapping bin. A recursive bin needs its name (`metator_00012_00003`): recursive bins are numbered within their overlapping bin. |
| `-O, --object` | `final_bin` | What `--name` is: `contig`, `core_bin`, `overlapping_bin`, `recursive_bin`, `final_bin`, or the header of a column of the contig data table. |
| `-m, --mat-fmt` | `graal` | The format of the map: `graal`, `bg2` or `cool`. |
| `-f, --filter` | off | Leave out the uncuts and loops. |
| `-D, --pcr-dup` | off | Leave out the PCR duplicates. |
| `-s, --min-size` | `5000` | Leave out the contigs shorter than this, in bp. |
| `-o, --outdir` | `./contact_map_<name>` | The output directory. |
| `-F, --force` | off | Replace the outputs. |
| `-T, --tmpdir` | `./tmp` | The temporary directory. |
| `-t, --threads` | `1` | Threads. |
| `-N, --no-clean-up` | off | Keep the temporary files. |

```bash
metator contactmap -a assembly.fa -c binning/contig_data_final.txt -e DpnII \
    -n metator_00007_00000 binning/alignment_0_sorted.pairs.gz
```

### `scaffold` { #scaffold }

Order and orient the contigs of a bin into scaffolds, from the contacts
between their ends. **API:** `metator.scaffold.get_scaffolds`. **Guide:**
[Viewing and scaffolding MAGs](guides/visualisation.md#scaffolding-a-bin).

```console
metator scaffold [OPTIONS] PAIRSFILE...
```

| Option | Default | Description |
|---|---|---|
| `-b, --bin-name` | *required* | The name of the bin. |
| `-i, --input-fasta` | *required* | The FASTA file of the bin. |
| `-T, --threshold` | `0.05` | The score above which the ends of two contigs are joined. |
| `-w, --window-size` | `5000` | The size, in bp, of the ends of the contigs. |
| `-j, --junctions` | `NNNNN` | The sequence between two contigs of a scaffold. |
| `-o, --out-fasta` | `./<bin>_scaffolded.fa` | The FASTA file of the scaffolds. |
| `-O, --out-frags` | `./<bin>_info_frags.txt` | The table of the contigs of each scaffold. |
| `-t, --threads` | `1` | Threads to sort pairs files that are not sorted yet. |

```bash
metator scaffold -b metator_00007_00000 -i binning/final_bin_unscaffold/metator_00007_00000.fa \
    binning/alignment_0_sorted.pairs.gz
```

### `pairs` { #pairs }

Sort pairs files by contig and position, compress them with bgzip and index
them with pairix, for the commands that read them. **API:**
`metator.io.sort_and_index_pairs`.

```console
metator pairs [OPTIONS] PAIRSFILE...
```

Each output goes next to its input, named after it without its last
extension: `sample.pairs` gives `sample_sorted.pairs.gz` and its index
`sample_sorted.pairs.gz.px2`; `sample.pairs.gz` gives
`sample.pairs_sorted.pairs.gz`.

| Option | Default | Description |
|---|---|---|
| `-F, --force` | off | Replace the outputs. |
| `-r, --remove` | off | Delete each input once it is sorted and indexed. |
| `-t, --threads` | `1` | Threads for the sort. |

### `qc` { #qc }

Measure the quality of the library on the large contigs of the good MAGs:
uncuts, loops, informative contacts, noise. **API:**
`metator.quality_check.quality_check`. **Guide:**
[Checking a library](guides/qc.md).

```console
metator qc [OPTIONS] PAIRSFILE...
```

| Option | Default | Description |
|---|---|---|
| `-a, --assembly` | *required* | The assembly, as a FASTA file. |
| `-e, --enzyme` | *required* | The restriction enzymes of the library. |
| `-O, --metator-dir` | the current directory | The output directory of the run, where `bin_summary.txt` and `contig_data_final.txt` are. |
| `-b, --bin-summary` | `<metator-dir>/bin_summary.txt` | The bin summary. |
| `-c, --contig-data` | `<metator-dir>/contig_data_final.txt` | The contig data table. |
| `-t, --threshold` | set from the library | The uncut and loop thresholds, in restriction sites: `uncuts,loops`. |
| `-P, --plot` | off | Write the plots of the events. |
| `-p, --prefix` | none | The name of the sample, on the plots and their files. |
| `-o, --outdir` | `./metator_qc` or `./metator_qc_<prefix>` | The output directory. |
| `-T, --tmpdir` | `./tmp` | The temporary directory. |
| `-n, --no-clean-up` | off | Keep the temporary files. |

```bash
metator qc -a assembly.fa -e DpnII -O binning -p sample -P binning/alignment_0_sorted.pairs.gz
```

### `mge` { #mge }

Bin the contigs of mobile genetic elements from the contacts between them.
**API:** `metator.mge.mge_binning`. **Guide:**
[Binning mobile genetic elements](guides/mge.md).

```console
metator mge [OPTIONS] PAIRSFILE...
```

| Option | Default | Description |
|---|---|---|
| `-a, --assembly` | *required* | The assembly, as a FASTA file. |
| `-c, --contigs` | *required* | The contig data table of the run. |
| `-b, --binning` | *required* | The bins of the contigs: `binning.txt`, or a collection exported by anvi'o. |
| `-m, --mges` | *required* | The MGE contigs, one name per line. |
| `-s, --threshold-bin` | `0.8` | The score above which two contigs, or bins, are merged. |
| `-r, --random` | off | Shuffle the contigs between the bins: a random binning of the same shape. |
| `-o, --outdir` | the current directory | The output directory. |
| `-T, --tmpdir` | `./tmp` | The temporary directory. |
| `-N, --no-clean-up` | off | Keep the temporary files, and the index of the assembly. |

```bash
metator mge -a assembly.fa -c binning/contig_data_final.txt -b binning/binning.txt \
    -m mges.txt -o mge binning/alignment_0_sorted.pairs.gz
```
