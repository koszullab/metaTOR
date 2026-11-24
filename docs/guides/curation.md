# Curating MAGs with anvi'o

metaTOR's clusters can keep some contamination, more so when the assembly
keeps small contigs: they make the clusters more complete, and noisier.
[anvi'o](https://anvio.org) removes it by hand, from the coverage of the
contigs across libraries and their sequence composition. This page goes from
metaTOR to anvi'o and back; anvi'o's
[tutorial on refining bins](https://merenlab.org/2015/05/11/anvi-refine/)
explains the refinement itself.

The commands are for anvi'o 7. They use the output directory `binning` of a
`pipeline` run on the map of two libraries, `L1` and `L2`, each made with
`hicstuff pipeline -M metacool -n` (hicstuff 4.0.0; see
[Other tools](../other-tools.md#from-the-reads-to-the-map)). `-n` keeps
hicstuff's alignments, which give the coverage of the contigs. The tables of
the run are in `binning/tables` (for a run step by step,
`metator dump binning/sample.metacool -o binning/tables` writes them).

## 1. The contigs database

The database holds only the contigs of the final clusters: those of
`binning.txt`, whose records are in `fasta/mags_with-mges/`. With `--mges`,
`binning.txt` also lists the MGE genomes, whose contigs are in
`fasta/mges/`.

```bash
mkdir -p anvio
cut -f 1 binning/tables/binning.txt > anvio/contigs.txt
cat binning/fasta/mags_with-mges/*.fa > anvio/contigs.fa
```

anvi'o only accepts contig names of ASCII letters, digits and underscores:
rename the contigs of the assembly before the run if needed (see anvi'o's
[note on FASTA files](https://merenlab.org/2016/06/22/anvio-tutorial-v2/#take-a-look-at-your-fasta-file)).
Build the database without splitting the contigs (`-L 0`), then find the
single-copy genes, which give completeness and redundancy:

```bash
anvi-gen-contigs-database -f anvio/contigs.fa -o anvio/contigs.db -n metator -L 0 -T 16
anvi-run-hmms -c anvio/contigs.db -T 16
```

Optional annotations that help curation:

```bash
anvi-setup-scg-taxonomy                    # once
anvi-run-scg-taxonomy -c anvio/contigs.db -T 16
anvi-setup-ncbi-cogs                       # once
anvi-run-ncbi-cogs -c anvio/contigs.db -T 16 --sensitive
anvi-display-contigs-stats anvio/contigs.db
```

## 2. The coverage profiles

anvi'o compares the coverage of the contigs across samples: profile each
library on its own (shotgun libraries can be added the same way). With `-n`,
hicstuff keeps the alignments of each mate in its temporary directory
(`<outdir>/tmp` by default), `<prefix>.for.bam` and `<prefix>.rev.bam`.
`anvi-init-bam` sorts and indexes them. The first mates are enough:

```bash
for lib in L1 L2; do
    anvi-init-bam hicstuff_${lib}/tmp/${lib}.for.bam -o anvio/${lib}.bam
    anvi-profile -i anvio/${lib}.bam -c anvio/contigs.db -o anvio/HiC_${lib} -S HiC_${lib} \
        -M 500 -T 16 --skip-hierarchical-clustering --skip-SNV-profiling \
        --contigs-of-interest anvio/contigs.txt
done
anvi-merge -c anvio/contigs.db -o anvio/merged -S HiC -W anvio/HiC_*/PROFILE.db
```

- `--skip-hierarchical-clustering`: anvi'o clusters each cluster when it is
  refined.
- `--skip-SNV-profiling` saves much time.
- `-M 500` keeps the contigs of 500 bp or more, as in the assembly.

## 3. The clusters of metaTOR

`binning.txt` is an anvi'o collection, a contig and its cluster per line:

```bash
anvi-import-collection binning/tables/binning.txt -c anvio/contigs.db -p anvio/merged/PROFILE.db \
    -C metator --contigs-mode
```

## 4. Refining the clusters

`anvi-refine` opens a cluster in anvi'o's interactive interface, in a web
browser:

```bash
anvi-refine -c anvio/contigs.db -p anvio/merged/PROFILE.db -C metator -b metator_00001_00000
```

Two clusters of a human gut sample: the first cleaned of its contamination,
the second split in two.

|A cluster cleaned|A cluster split in two|
|:-:|:-:|
|![A cluster cleaned in anvi'o](../images/anvio_example_12.png)|![A cluster split in two in anvi'o](../images/anvio_example_186.png)|

| Cluster | Before: completeness – contamination (%) | After |
|:-:|:-:|:-:|
| Cleaned | 97.2 – 19.7 | 95.8 – 9.9 |
| Split | 95.8 – 98.6 | 94.4 – 0.0, and 93.0 – 0.0 |

## 5. The curated clusters

Export the collection, a contig (split) and its cluster per line:

```bash
anvi-export-collection -p anvio/merged/PROFILE.db -C metator -O anvio/curated
```

From `anvio/curated.txt`, metaTOR's Python API writes one FASTA file per
cluster, checks them with miComplete as `metator validation` does, counts
them by quality class, and draws the quality plots that the run draws in
`binning/plots/`:

```python
import metacooler
import metator.figures as mtf
import metator.validation as mtv
import pandas as pd

# The collection exported by anvi'o: each split and its cluster, a split read
# as its contig (NODE_1_split_00001 is NODE_1).
collection = pd.read_csv("anvio/curated.txt", sep="\t", header=None, names=["split", "cluster"], dtype=str)
binning = collection.set_index(collection["split"].str.split("_split").str[0])["cluster"]
binning = binning[~binning.index.duplicated(keep="last")]

# One FASTA file per curated cluster, written by metacooler from the map of the run.
groups = {bin_name: list(contigs.index) for bin_name, contigs in binning.groupby(binning)}
metacooler.write_fasta("binning/sample.metacool", "assembly.fa", "anvio/fasta", groups=groups)

# miComplete on every cluster, with the bacterial and the archaeal markers.
mtv.micomplete_quality("anvio/fasta", "anvio/micomplete.txt", threads=16)

# The quality class of each curated cluster (HQ, MQ, LQ, Contaminated, Other).
summary = pd.read_csv("anvio/micomplete.txt", sep="\t", comment="#", index_col=0)
quality = summary.apply(lambda row: mtv.bin_quality(row["Weighted completeness"], row["Weighted redundancy"]), axis=1)
print(quality.value_counts())

# The quality plots of the validation, for the curated clusters.
total = int(metacooler.MetaCooler("binning/sample.metacool").contigs["length"].sum())
mtf.figures_bins_distribution(summary, "anvio/bins_distribution.pdf")
mtf.figures_bins_size_distribution(summary, total, 500000, "anvio/bins_size_distribution.pdf")
```

Compare them with `binning/plots/bins_distribution.pdf` and
`binning/plots/bins_size_distribution.pdf`, from before the curation.

To look at the contacts of a curated cluster, add the curated clusters to the
contig data table as a column, regroup the map of the run with
`metacooler reassign`, then extract a cluster with `metacooler extract` (see
[Other tools](../other-tools.md)):

```python
# binning: the curated cluster of each contig, read from the collection as above.
contigs = pd.read_csv("binning/tables/contig_data_final.txt", sep="\t", dtype={"Name": str})
contigs["Anvio"] = contigs["Name"].map(binning).fillna("ND")
contigs.to_csv("binning/contig_data_anvio.txt", sep="\t", index=False)
```

```bash
metacooler reassign binning/sample.metacool --contigs binning/contig_data_anvio.txt \
    --contigs-format Name,Anvio -o anvio/curated.metacool
metacooler extract anvio/curated.metacool --genomes metator_00001_00000 -o metator_00001_00000.cool
```
