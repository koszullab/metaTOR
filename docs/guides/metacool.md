# The contact map of the bins

At the end of a run, `metator pipeline` writes the contact map of all the
contigs, grouped by bin, as `metator.metacool` in the output directory
(`<prefix>.metacool` with `--prefix`); `--no-metacool` skips it. A
`.metacool` file is the contact map format of
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
  without `--enzyme`). The genome of a contig is its final bin.
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
step, or one that skipped the map:

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

When the map cannot be written (contig names that are not ASCII, or a
misaligned `bin_summary.txt`), `metator pipeline` logs a warning saying why,
and keeps its other outputs and its exit status; `metator metacool` logs the
reason as an error and stops with it. A map that an earlier run left in the
output directory is removed before the new one is written, as it no longer
matches the bins: a run whose map cannot be written leaves none.

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
  contig IDs, GC contents, coverages and restriction sites still come from the
  assembly, `--depth` and `--enzyme`. It may lack some of them, as a map made
  with a contig filter does (`metacooler create pairs --filter 'contigs.length
  >= 1kb'`, which also leaves out the pairs with a read on the contigs that it
  drops): those get no contacts and no hits, and the log says how many they
  are. A map of other contigs or lengths, or that stores both triangles of the
  matrix, stops the run with an error.
- **Its libraries are already summed**: no `network_<i>.txt` is written.
  `--edge` must be 0, as a map holds no read positions, and
  `pipeline --scaffold`, which reads the pairs files, is refused.
- **The heatmap** of the bins' network is drawn from the `.pairs.gz` files of
  the output directory: a run from a map draws none, and warns that no pairs
  file was found, unless an earlier run left its pairs files there.

At the end of such a run, `pipeline` writes the map of its bins from the input
map, with `metacooler reassign`, as a new file:

- the bins and contacts of the input map, at its bin size; one genome per
  final bin and `__unbinned__` for the other contigs, grouped genome by genome,
  the longest first, as in the map of a run from pairs;
- the columns of `bin_summary.txt` and `contig_data_final.txt`, which replace
  those the input map had (the map of an earlier binning, for example);
- the contact statistics (the `qc` block of the metadata and the `qc_*`
  columns of the genomes table), computed again from the pixels when the input
  map has them: at one bin per contig, every contact within a contig counts
  as short-range, and the numbers of pairs read and left out are not known;
- the assembly name and the history of the commands of the input map.

The map of the bins may go where the input map is, to bin its contigs again in
the same output directory (`-S metacool -1 binning/metator.metacool -o binning
-F`): the new map is written next to it and validated before it replaces it,
and the old map is left as it was, with a warning, when the new one cannot be
written. A symbolic link there is replaced, not its target.

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

# 3. Bin the contigs; binning/metator.metacool is the map of the bins.
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
