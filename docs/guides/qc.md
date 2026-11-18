# Checking a library

`metator qc` measures the quality of a metaHiC library on the genomes that
the binning recovered best. Within the contigs of good MAGs, each pair can be
told apart as a real contact or an artefact of the protocol, and a pair that
joins two of these MAGs can only be noise.

## What it measures

1. **The MAGs.** It takes the bins of `bin_summary.txt` with a weighted
   completeness above 70% and a weighted redundancy under 1.15, and their
   contigs longer than 100 kb.
2. **The pairs.** It reads, from the pairs files, the pairs between these
   contigs, and finds the restriction fragment of each mate, from the
   restriction sites of `--enzyme`.
3. **The events.** Each pair within a contig is one of hicstuff's
   [3C events](https://koszullab.github.io/hicstuff/guides/filter/):
   - an **uncut** (*religated* in the log): the two mates face each other
     (`+-`) across a few restriction sites, a fragment that was not cut, or
     was cut and ligated again;
   - a **loop**: the mates face away from each other (`-+`) across a few
     restriction sites, the two ends of a fragment ligated together;
   - a **weird**: the mates are on the same fragment and the same strand;
   - **informative** otherwise.

   A pair between two contigs of the same MAG is informative, and one between
   two MAGs is noise.

"A few" is a threshold for the uncuts and one for the loops, in restriction
sites: hicstuff sets them from the library (as `hicstuff filter` does), or
`--threshold` (`-t`) gives them, as `uncuts,loops`: `-t 3,3`.

## Running it

```bash
metator qc -a assembly.fa -e DpnII -O binning -p sample -P binning/alignment_0_sorted.pairs.gz
```

It reads `bin_summary.txt` and `contig_data_final.txt` in `--metator-dir`
(`-O`, the current directory by default), or where `--bin-summary` (`-b`) and
`--contig-data` (`-c`) say, and the pairs files of the run (kept with
`metator pipeline --no-clean-up`).

| Option | Default | Effect |
|---|---|---|
| `--enzyme` (`-e`) | *required* | The enzymes of the library, separated by commas |
| `--threshold` (`-t`) | set from the library | The uncut and loop thresholds |
| `--prefix` (`-p`) | none | The name of the sample, on the plots and their files |
| `--plot` (`-P`) | off | Write the plots |
| `--outdir` (`-o`) | `metator_qc`, or `metator_qc_<prefix>` | Where the plots and the log go |
| `--no-clean-up` (`-n`) | off | Keep the temporary files: the pairs extracted, with their fragments |

With `--plot`, it writes `<prefix>_event.pdf`, the distribution of the events
from which the thresholds are set (not drawn when `--threshold` gives them),
and `<prefix>_camembert_plot.pdf`, the share of each kind of pair.

## The results

The results are in the log, `metator_qc_<date>.log`. On a simulated
metagenome of two genomes, in contigs larger than 100 kb:

```console
INFO :: 30000 pairs have been extracted.
INFO :: Filtering with thresholds: uncuts=2 loops=2
INFO :: Number of MAGs: 2
INFO :: Religated ratio: 4.77%.
INFO :: Loop ratio: 4.73%.
INFO :: Weirds ratio: 3.77%.
INFO :: Informative contacts estimation: 85.97%.
INFO :: Ratio inter/intra contigs: 9.88%.
INFO :: Noise contact ratio: 0.88%
INFO :: Noise score: 1.74E-02
```

| Line | What it is |
|---|---|
| `Number of MAGs` | The MAGs selected |
| `Religated ratio`, `Loop ratio`, `Weirds ratio` | The uncuts, loops and weirds, in percent of the pairs within the MAGs (within a contig, or between two contigs of the same MAG) |
| `Informative contacts estimation` | The informative pairs, in percent of all the pairs extracted |
| `Ratio inter/intra contigs` | The share of the informative pairs that join two contigs |
| `Noise contact ratio` | The pairs between two MAGs, in percent of all the pairs extracted |
| `Noise score` | The contacts per pair of MAGs, against those per MAG: *x* / (*y* + *x*), with *x* the pairs between MAGs divided by the number of pairs of MAGs, and *y* the pairs within MAGs divided by the number of MAGs. 0 for a library with no noise |

When no MAG passes, or none has a contig longer than 100 kb, as on the small
test library of the [Quickstart](../quickstart.md), no pair is extracted, and
the log says why, such as `The 2 high-quality MAGs have no contig longer than
100 kb: no quality check.`
