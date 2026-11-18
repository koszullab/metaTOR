# The contact map of the bins

`metator pipeline` writes the contact map of all the contigs, grouped by bin,
as `metator.metacool` in the output directory (`<prefix>.metacool` with
`--prefix`). It reads the pairs files into the map once, as soon as it has
them: the map of the contigs, in no bin yet (every contig in `__unbinned__`,
the longest first), with their contacts, the contact statistics that metacooler
measures on the reads and what metaTOR measures on the assembly (`GC_content`,
and `Shotgun_coverage` and `Restriction_site` with `--depth` and `--enzyme`). A
run with one library reads its network from this map, as a run from a map does
(see [Starting from a contact map](#starting-from-a-contact-map)):
`network.txt` and `contig_data_network.txt` are those of its pairs files, byte
for byte, as `pipeline` counts every pair (`--edge` must be 0, whatever the
start). It still writes `network_0.txt`, the network of its library, as a copy
of `network.txt`: its pairs files give the same bytes. A run with several
libraries builds the map of each library, the entries `library_0`,
`library_1`... of `metator.mmetacool` (`<prefix>.mmetacool`), and the map of
their sum, `metator.metacool`, merged from them (their contacts and contact
statistics summed; its counts are stored as 64-bit integers): it reads
`network.txt` from the sum, and the network of each library (`network_N.txt`)
and `hit_data_alignment.txt` from its map, the files that its pairs files give,
byte for byte. So a pair on a contig that the assembly lacks is left out, as
metacooler warns (it counts them, with the pairs that have a read past the end
of their contig), and the run goes on, whatever the number of libraries. At the
end of the run, `pipeline` writes the map of the bins over the map of the
contigs, from the final tables, with `metacooler reassign`; `--no-metacool`
removes the maps instead (the `.mmetacool` collection too). The maps of the
libraries stay maps of the contigs, in no bin, without the columns measured on
the assembly, which the sum holds. Its file is larger than that of a map
written from the pairs in one go (by 12 to 32% on our test data), for the same
data: the columns measured on the assembly, dropped and written again at the
end, leave free space in it. A `.metacool` file is the contact map format of
[metacooler](https://github.com/js2264/metacooler): a
[cooler](https://github.com/open2c/cooler) whose chromosomes are the contigs,
grouped by genome. `cooler` and the tools built on it read it as any `.cool`
file, and metacooler adds the genomes.

## What it holds

- **One bin per contig, and one genome per final bin.** The contigs that are
  in no final bin (`ND` in `contig_data_final.txt`) are kept, with all their
  contacts, in the `__unbinned__` genome: all the contigs, for a run with no
  final bin. The contigs come genome by genome, the longest first.
- **The contacts.** Each pair whose two ends lie on contigs of the assembly is
  one contact between these two contigs (an intra-contig contact when both
  ends are on the same contig). These are the raw contacts of the pairs
  files, without the filters and the normalisation of the network, and those
  of several libraries are summed.
- **The genomes table** holds the columns of `bin_summary.txt` (completeness,
  redundancy, Hi-C and shotgun abundances…), next to the length, contig and bin
  counts of each genome and the contact statistics that metacooler computes
  while reading the pairs.
- **The contigs table** (the chromosomes table of `cooler`) holds what
  metaTOR measured on each contig: the `GC_content`, `Hit`,
  `Shotgun_coverage` and `Restriction_site` columns of
  `contig_data_final.txt`, except those that are `-` on every contig
  (`Shotgun_coverage` of a run without `--depth`, `Restriction_site` of a run
  without `--enzyme`, when the map that a run starts from does not give them
  either), and the contact statistics that metacooler computes
  within each contig (`qc_intra_contig_short`, `qc_intra_contig_long`). The
  genome of a contig is its final bin.
- **The metadata** hold the version of metaTOR, the command line and its
  parameters, under the `metator` key.

```python
import cooler
import metacooler

clr = cooler.Cooler("binning/metator.metacool")
mc = metacooler.MetaCooler("binning/metator.metacool")
mc.genomes  # one row per bin, and __unbinned__, with the columns of bin_summary.txt
mc.contigs  # one row per contig, with its GC content and Hi-C hits
```

```bash
metacooler info binning/metator.metacool
metacooler dump -t genomes -H binning/metator.metacool
```

## `metator metacool`

`metator metacool` writes the same file from the pairs files,
`contig_data_final.txt` and `bin_summary.txt` of a run, for a run made step by
step, or one that skipped the map, but for the history that it stores (the
`annotate` and `reassign` keys of its metadata) and the type of its counts with
several libraries (64-bit integers, which `metator metacool` stores as 32-bit
ones):

```bash
metator metacool -c binning/contig_data_final.txt -b binning/bin_summary.txt \
    -o sample.metacool binning/alignment_0_sorted.pairs.gz binning/alignment_1_sorted.pairs.gz
```

For a run that started from a contact map step by step (`network -S
metacool`, `partition`, `validation`), it writes the map from that map, as
`pipeline -S metacool` does; the map keeps its assembly name, so
`--assembly-name` is refused:

```bash
metator metacool -c contig_data_final.txt -b bin_summary.txt -o binning.metacool sample.metacool
```

## Names and failures

The contig and bin names must be ASCII, as `cooler` reads them as ASCII text.
No bin may be named as metacooler marks a contig in no bin: `NA`, `NaN`,
`N/A`, `N/D`, `-1` or `__unbinned__` in any case, and `nd` in any case other
than `ND`, which marks a contig in no bin. metaTOR's own bin names never are,
but those of a hand-edited `contig_data_final.txt` can be. Each row of
`bin_summary.txt` must have as many fields as its header, as metaTOR writes
it.

When the map of the clusters cannot be written, `metator pipeline` and
`metator metacool` log the reason as an error and stop with it; `pipeline`
keeps its other outputs. `pipeline` removes a map that
an earlier run left in the output directory before it builds the map of the
contigs in its place (`--no-metacool` or not), as it no longer matches the
bins: from then on, the map there is this run's, or there is none. So a run
from reads, alignments or pairs refuses, before any work, an input that is
that map (or a symbolic link to it or from it), and, from pairs, a contact
map: start from a map with `--start metacool`. A map of the contigs that
cannot be built (for example, from pairs files that metacooler cannot read,
see below) stops the run before its network, `--no-metacool` or not, whatever
the number of libraries: the reason is logged as an error, in the log file of
the run, and no map is left. When the map of the bins cannot be written over
it (for example, from a misaligned `bin_summary.txt`), the map of the contigs
stays, in no bin. So it does when the run stops with an error once it is
built (in the partition or the validation, for example), `--no-metacool` or
not: a new run can start from it, without aligning the reads again, with
`--start metacool -1 <outdir>/metator.metacool` (`<outdir>/<prefix>.metacool`
with `--prefix`) and `-F`, which replaces the bin directories that the
stopped run left in the output directory; `--scaffold` is not available then
(see [Starting from a contact map](#starting-from-a-contact-map)).

`pipeline` checks what the map of the contigs needs before any work. Its
contig names must be ASCII: a run that is not from a map refuses an assembly
with other contig names, which it gives, before it builds the index of the
assembly or aligns a read; rename these contigs to ASCII names, in the
assembly and in the alignments or pairs files made from it (and build its
index again, for an assembly given as an index). metacooler, which builds the
map, and metaTOR, which reads the pairs files of `--scaffold`, both find the
compression of a pairs file from its first bytes, whatever its name:
metacooler reads plain, gzip, bzip2, zstd and LZ4 files, metaTOR plain, gzip,
bzip2 and zip ones. So a run from pairs refuses, before any work, an empty
pairs file; a zip archive, which metacooler cannot read; and, with
`--scaffold`, a zstd or LZ4 file, which metaTOR cannot read: decompress the
file, or compress it with gzip. zstd and LZ4 files work without `--scaffold`:
the networks come from the maps. metacooler also needs the same number of
fields on every row of a pairs file; its header lines (those that start with
`#`) may come anywhere, so pairs files of one assembly concatenated with `cat`
are one library. Give the pairs files of several libraries comma-separated,
one file per library. A run from pairs that metacooler cannot read stops with
the reason and these needs.

## Starting from a contact map

`metator network` and `metator pipeline` can build the network from a
`.metacool` contact map of the contigs of the assembly, instead of reads,
alignments or pairs: `--start metacool`, with the map given to `--forward`
(one file, or an entry of a `.mmetacool` collection,
`file.mmetacool::entries/<name>`):

```bash
metator pipeline -S metacool -1 sample.metacool -a assembly.fa -d depth.txt -e DpnII -o binning
```

- **Any bin size**, from one bin per contig (as the maps written by metaTOR)
  to fixed bins. The contacts between two contigs are those between their
  bins, and the hits of a contig those of its bins: the network and
  `contig_data_network.txt` are those of the pairs the map holds, byte for
  byte, at `--edge 0`.
- **Its contigs must be contigs of the assembly**, with their lengths: the
  contig IDs and GC contents still come from the assembly. It may lack some
  of them, as a map made with a contig filter does (`metacooler create pairs
  --filter 'contigs.length >= 1kb'`, which also leaves out the pairs with a
  read on the contigs that it drops): those get no contacts and no hits, and
  the log says how many they are. A map of other contigs or lengths, or that
  stores both triangles of the matrix, stops the run with an error.
- **Its restriction sites and coverage**: the `Restriction_site` and
  `Shotgun_coverage` columns of its contigs table, which metaTOR stores in
  the maps it writes with `--enzyme` and `--depth`, give the restriction
  sites and the shotgun coverage of the contigs. So a run from such a map
  needs neither option: it gives the same network and tables, byte for byte
  (with pandas 2, for coverages of at most 15 significant digits and 15
  digits after the decimal point, as jgi's depth files have: metacooler can
  store a longer one one unit off in its last place, which changes the last
  digits of the contig data files and of a network normalized by it, and
  has `--depth` warn that it replaces the map's value).
  - `--enzyme` and `--depth` fill in a column that the map lacks, and replace
    one that it has, with a warning that counts the contigs whose value
    differs. The map of the bins stores what the run used, and its metadata
    that `pipeline` writes which columns it read (`columns_from_map`).
  - A column is read only if it gives every contig of the assembly a value: a
    number of restriction sites (plus one) of at least 1, a coverage of at
    least 0. A contig of no length, which no map holds, gets 1 and 0. A
    column with no valid value, such as one of `-`, as metaTOR writes what it
    did not measure, is as no column. One that lacks some contigs is not
    read, with a warning that names them, and the map of the bins does not
    carry it.
  - A normalization that needs a column that neither the map nor the options
    give (`RS`, `abundance`, `theoritical_hit`) is refused before the run
    replaces anything: before any work when the map has no such column, once
    the assembly is read when the column lacks contigs.
- **Its libraries are already summed**: no `network_<i>.txt` is written.
  `--edge` must be 0, as a map holds no read positions (`pipeline` refuses
  any other value, whatever the start: it counts every pair, as its map does;
  `metator network` applies it to pairs files), and `pipeline --scaffold`,
  which reads the pairs files, is refused.
- **The heatmap** of the bins' network is drawn from the `.pairs.gz` files of
  the output directory: a run from a map draws none, and warns that no pairs
  file was found, unless an earlier run left its pairs files there.

At the end of such a run, `pipeline` updates the input map in place with the
map of its clusters, with `metacooler reassign` (`--output-map FILE` writes it
to FILE instead):

- the bins and contacts of the input map, at its bin size; one genome per
  final bin and `__unbinned__` for the other contigs, grouped genome by genome,
  the longest first, as in the map of a run from pairs;
- the columns of `bin_summary.txt` and `contig_data_final.txt`, last, as in
  the map of a run from pairs: they replace those the input map had (the map
  of an earlier binning, for example);
- the contact statistics (the `qc` block of the metadata and the `qc_*`
  columns of the genomes and contigs tables) of the input map, when it has
  them, carried over to the bins: metacooler keeps the counts within each
  contig (measured on the reads, for a map made from pairs files) and counts
  the contacts between contigs again. A map written by metacooler 0.2.0 holds
  no count per contig: its statistics are computed again from the pixels,
  which hold no read positions (at one bin per contig, every contact within a
  contig then counts as short-range, and the numbers of pairs read and left
  out are lost, which metacooler warns about for a map made from pairs files);
- the assembly name and the history of the commands of the input map.

The input map is replaced only once the map of the clusters is written next to
it and validated in full. When it cannot be, the run stops with the reason, and
the input map is left as it was, as it is when the run stops for any other
reason (Ctrl-C included). A symbolic link is replaced, not its target.
`--output-map FILE` leaves the input map as it was: a run from an entry of a
`.mmetacool` collection needs it, as a run does not update a collection in
place, and a FILE that exists is replaced only with `-F`.

While a run writes a map, a lock file next to it (`<map>.lock`) holds it: a
second run that would write the same map stops before any work, naming the run
that holds it. A lock that a run left behind on the same host when it died is
removed, with a warning; one from another host is not: remove it once that
run is over.

## From hicstuff to metaTOR through metacooler

A map to start from can be made from the pairs files of
[hicstuff](https://koszullab.github.io/hicstuff/), aligned in mapping mode
`5spm`, with no binning yet:

```bash
# 1. Align each library, without -f or -D: they leave out pairs (3C events that
#    hicstuff filters, PCR duplicates) that metaTOR's own alignment keeps.
hicstuff pipeline --mapping 5spm -g assembly.fa -e DpnII -o hicstuff_L1 -P L1 \
    L1_R1.fq.gz L1_R2.fq.gz
hicstuff pipeline --mapping 5spm -g assembly.fa -e DpnII -o hicstuff_L2 -P L2 \
    L2_R1.fq.gz L2_R2.fq.gz

# 2. One map of both libraries, one bin per contig, every contig in __unbinned__.
samtools faidx assembly.fa
awk -v OFS='\t' '{print $1, 0, $2}' assembly.fa.fai > contigs.bed
metacooler create pairs --pairs hicstuff_L1/L1.valid_idx.pairs.gz \
    --pairs hicstuff_L2/L2.valid_idx.pairs.gz --bins contigs.bed --no-ps \
    --assembly-name assembly -o sample.metacool

# 3. Cluster the contigs: sample.metacool becomes the map of the clusters, in place.
metator pipeline -S metacool -1 sample.metacool -a assembly.fa -e DpnII -o binning
```

!!! note "Mapping mode `5spm` needs a hicstuff release after 3.3.0"

    hicstuff 3.3.0, the latest release so far, has no mapping mode `5spm`:
    with it, leave `--mapping 5spm` out.

- `--bins 1000` makes bins of 1 kb instead, taking the lengths of the contigs
  from the `#chromsize` lines of hicstuff's pairs files, or from `--fasta`.
- `--no-ps`: otherwise metacooler also writes `sample.ps.tsv`, a P(s) table of
  the single `__unbinned__` genome.
- `--assembly-name`: otherwise the assembly of the map, and of the map of the
  bins, is `unknown`.
- `--contig-order fasta` keeps the contigs in the order of the assembly (by
  default, the longest first). It needs `--fasta` or the `#chromsize` lines of
  the pairs files: metacooler does not keep the order of a BED file.
  `--fasta` writes the index of the assembly (`.fai`) next to it, in a
  directory that must be writable.

On the test library of the [Quickstart](../quickstart.md), this route finds
the same 28,510 contacts between contigs as metaTOR's own alignment with
`--aligner bwa`.
