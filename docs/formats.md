# File formats

The files metaTOR reads and writes. Every table is tab-separated.

## Inputs

### Assembly

An uncompressed FASTA file of the contigs (`network` and `pipeline` refuse a
compressed one), or the basename of its bowtie2
or bwa index (`--assembly`). Contig names must be ASCII for the contact map
of the bins (see [Names and failures](guides/metacool.md#names-and-failures)).

### Depth

The output of `jgi_summarize_bam_contig_depths` of MetaBAT2 (`--depth`), a
header line, then one line per contig:

```text
contigName	contigLen	totalAvgDepth	sample.bam	sample.bam-var
NODE_522	22786	4.76595	4.76595	15.5294
```

metaTOR reads its lines in the order of the contigs of the assembly, not by
name: they must list the same contigs, in the same order, as
`jgi_summarize_bam_contig_depths` does from reads aligned on that assembly.
It takes the size of each contig from the second column and its coverage
from the third.

### Pairs

The pairs files that metaTOR writes, and that `--start pair` reads, follow the
[pairs format](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md)
of the 4D Nucleome, with 7 columns: read name, contig and position of the
first mate, contig and position of the second, and the strand of each.

```text
## pairs format v1.0
#columns: readID chr1 pos1 chr2 pos2 strand1 strand2
#sorted: chr1-chr2-pos1-pos2
#shape: upper triangle
#chromsize: NODE_522 22786
...
H9:1:HGJYYBBXY:8:1102:3792:13693	NODE_100063	161	NODE_100063	402	-	-
```

metaTOR writes them sorted, compressed with bgzip and indexed with pairix
(`alignment_<i>_sorted.pairs.gz` and `.px2`), as `metator pairs` does for
other pairs files.

### Lists of contigs

The MGE contigs of `metator mge --mges`: one contig name per line.

The bins of `metator mge --binning`: a contig and its bin per line, without
header, as `binning.txt` or the collection exported by anvi'o.

## The contig data tables

`contig_data_network.txt`, `contig_data_partition.txt` and
`contig_data_final.txt` hold one row per contig of the assembly, in its order,
with a header. Each step adds its columns to those of the step before:

| Column | Written by | What it holds |
|---|---|---|
| `ID` | network | The ID of the contig, from 1, in the order of the assembly: the ID of the network |
| `Name` | network | Its name |
| `Size` | network | Its length, in bp |
| `GC_content` | network | Its GC content, in percent |
| `Hit` | network | Its hits: the ends of the pairs that lie on it, in all libraries (a pair within the contig counts twice) |
| `Shotgun_coverage` | network | Its coverage in the depth file, `-` without `--depth` |
| `Restriction_site` | network | Its number of restriction sites plus one (its restriction fragments), `-` without `--enzyme` |
| `Core_bin_ID`, `Core_bin_contigs`, `Core_bin_size` | partition | Its core bin: the ID, number of contigs and size in bp of the bin |
| `Overlapping_bin_ID`, `Overlapping_bin_contigs`, `Overlapping_bin_size` | partition | Its overlapping bin, likewise |
| `Recursive_bin_ID`, `Recursive_bin_contigs`, `Recursive_bin_size` | validation | Its recursive bin, likewise: `00000`, `-`, `-` when its bin was not partitioned again |
| `Final_bin` | validation | Its final bin, or `ND` (not determined) |

```text
	ID	Name	Size	GC_content	Hit	Shotgun_coverage	Restriction_site	Core_bin_ID	Core_bin_contigs	Core_bin_size	Overlapping_bin_ID	Overlapping_bin_contigs	Overlapping_bin_size	Recursive_bin_ID	Recursive_bin_contigs	Recursive_bin_size	Final_bin
0	1	NODE_522	22786	62.0775915035548	1683	4.76595	162	00001	2	23301	00001	15	49021	00000	-	-	metator_00001_00000
1	2	NODE_1404	15015	64.2823842823842	952	4.00585	118	00002	1	15015	00002	5	19711	00000	-	-	ND
```

`contig_data_final.txt` starts with a column without header, the row number
from 0. A contig with no edge in the network (no hits, or no contact with
another contig) has `-` in the columns of its core and overlapping bins. As
inputs, the tables are read by the headers of their columns: other columns
may be added, in any order, such as the bins of another binning for
`metator contactmap --object`.

With several libraries, `hit_data_alignment.txt` holds the hits of each
contig in each library: its ID, its name, then one column per library,
without header.

## The network

`network.txt` holds the network of contigs, one line per pair of contigs in
contact, without header: the IDs of the two contigs, and the normalised
weight of their edge. `network_<i>.txt` holds the network of library `i`
alone.

```text
1	2	0.014220405247506712
1	3	0.008170769116255661
```

## The bins

### Bin summary

`bin_summary.txt` holds one row per final bin, named in its first column. Its
first columns are those of miComplete:

| Column | What it holds |
|---|---|
| `Length` | The size of the bin, in bp |
| `GC-content` | Its GC content, in percent |
| `Present Markers` | The single-copy marker genes found in it |
| `Completeness`, `Redundancy` | The share of the markers found, and their mean number of copies (1 for no duplicate) |
| `Weighted completeness`, `Weighted redundancy` | The same, with miComplete's weights of the markers of the set |
| `Contigs`, `N50`, `L50`, `N90`, `L90` | Its number of contigs, and the N50, L50, N90 and L90 of their sizes |
| `CDs` | Its number of coding sequences, found by prodigal |
| `recursive` | Whether the bin comes from the last recursive step that split bins: in metaTOR 1.4.0, every bin reads `True` when no bin is partitioned again, and a bin made at an earlier step reads `False` |
| `step` | The recursive step that made the bin, 0 for an overlapping bin |
| `parent` | The bin it was split from (an overlapping bin at step 1, a recursive bin after), empty otherwise |
| `HiC_abundance` | The share of the hits of all the contigs that are on the bin, in percent |
| `SG_abundance` | The share of the shotgun reads (size × coverage) of all the contigs that are on the bin, in percent; with `--depth` only |

```text
	Length	GC-content	Present Markers	Completeness	Redundancy	Weighted completeness	Weighted redundancy	Contigs	N50	L50	N90	L90	CDs	recursive	step	parent	HiC_abundance	SG_abundance
metator_00001_00000	49021.00	62.03	1.00	0.01	1.00	0.02	1.00	15.00	9019.00	2.00	1390.00	9.00	70.00	True	0		2.4736	2.6314
```

Each marker set, bacterial (105 markers) and archaeal (131), is tried on
every bin, which keeps the one that gives it the higher weighted
completeness. The plots take the weighted redundancy minus one as the
contamination of a bin; the log of the validation, and the recursive
partition, divide it by the weighted completeness (see
[Bin quality](concepts.md#bin-quality)).

### binning.txt

Each contig of a final bin and its bin, without header: the collection that
anvi'o imports (see [Curating MAGs with anvi'o](guides/curation.md)).

!!! warning "`binning.txt` can miss contigs"

    When a contaminated bin is partitioned again but kept whole, its contigs
    that went to a recursive part are left out of `binning.txt` in
    metaTOR 1.4.0, though `contig_data_final.txt` gives their bin. The full
    binning is in `contig_data_final.txt`:

    ```bash
    awk -F '\t' 'NR > 1 && $NF != "ND" {print $3 "\t" $NF}' contig_data_final.txt > bins.txt
    ```

```text
NODE_522	metator_00001_00000
NODE_2487	metator_00007_00000
```

### miComplete results

`overlapping_micomplete_results.txt` and `recursive_micomplete_step_<n>.txt`
hold miComplete's results on the overlapping bins and on the recursive bins
of each step: three comment lines (`##`), then the columns of the bin summary
up to `CDs`, and the marker set kept, `Bacteria` or `Archaea`.

### Clustering matrices

With `--cluster-matrix`, `clustering_matrix_partition.npz`,
`clustering_matrix_recursive.npz` and, in `pipeline`, `clustering_matrix.npz`
(the mean of both) hold how often each two contigs were binned together, from
0 to 1, as a sparse matrix of scipy (`scipy.sparse.load_npz`). In metaTOR
1.4.0, only the matrix of the partition is right: see the
[warning](guides/binning.md#validation) on the recursive ones. Row and
column *i* are the contig of ID *i* (row and column 0 are empty), and only the
upper triangle is stored. See
[The clustering matrices](guides/python.md#the-clustering-matrices).

### FASTA files of the bins

`overlapping_bin/`, `recursive_bin/step_<n>/` (the recursive bins of step
`n`) and `final_bin/` (`final_bin_unscaffold/` in `pipeline`) hold one FASTA
file per bin, `<bin>.fa`, with the records of its contigs as they are in the
assembly. The final bins also hold the pyfastx
index of each, `<bin>.fa.fxi`, which the validation leaves there.

## The contact map of the bins

`metator.metacool` is a `.metacool` file of
[metacooler](https://github.com/js2264/metacooler): a cooler, one bin per
contig, whose contigs are grouped into one genome per final bin. See
[The contact map of the bins](guides/metacool.md).

## Scaffolds

`<bin>_scaffolded.fa` holds one record per scaffold, named
`metator_scaffold_0001` and on, with its shape and the length of its contigs:

```text
>metator_scaffold_0001 linear length: 12167
```

`<bin>_info_frags.txt` lists the contigs of each scaffold, in order, after
the same header line: the name, start (0), end and orientation (`1`, or `-1`
for a reversed contig) of each.

```text
>metator_scaffold_0001 linear length: 12167
NODE_2487	0	11273	1
NODE_202465	0	894	1
```

## MGE bin summary

`mges_bin_summary.tsv` holds one row per MGE bin, with a header:

| Column | What it holds |
|---|---|
| `BinName` | `MetaTOR_MGE_00001` and on |
| `BinLength`, `GC`, `Hit` | The size, GC content and hits of the bin |
| `BinningScore` | The lowest score at which its contigs were merged, `NA` for a single contig |
| `ContigsNumber`, `Contigs` | Its number of contigs, and their names, separated by commas |
| `MetagenomicBin`, `AssociationScore`, `AssociatedBins` | For its host MAGs: `None`, `NA` and `NA`, as `metator mge` does not look for them |

```text
BinName	BinLength	GC	Hit	BinningScore	ContigsNumber	Contigs	MetagenomicBin	AssociationScore	AssociatedBins
MetaTOR_MGE_00064	30610.00	68.93	18383.00	1.20	2.00	NODE_8199_len_17842,NODE_10116_len_12768	None	NA	NA
MetaTOR_MGE_00105	23587.00	37.36	2662.00	NA	1.00	NODE_10949_len_23587	None	NA	NA
```

## Logs

Each command writes its log, `metator_<command>_<date>.log`: the version of
metaTOR, the date and the command line, then what the run did. hicstuff's
own log goes next to it when a command runs hicstuff (`contactmap`).
