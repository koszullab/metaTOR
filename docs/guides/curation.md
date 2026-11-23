# Curating MAGs with anvi'o

The clusters of metaTOR can keep some contamination, the more so as the assembly
keeps small contigs: they make the clusters more complete, and noisier.
[anvi'o](https://anvio.org) can remove that contamination by hand, from the
coverage of the contigs across libraries and their sequence composition. This
page goes from the outputs of metaTOR to anvi'o and back; anvi'o's own
[tutorial on refining bins](https://merenlab.org/2015/05/11/anvi-refine/)
explains the refinement itself.

The commands below were written for anvi'o 7. They run on the output
directory `binning` of a run on the map of two libraries, `L1` and `L2`, each
made with `hicstuff pipeline -M metacool -n`, from hicstuff 4.0.0 (see
[Other tools](../other-tools.md#from-the-reads-to-the-map)): `-n` keeps the
alignments of hicstuff, which give the coverage of the contigs. The tables of
the run, `binning.txt` and `contig_data_final.txt`, are written from the map
of its clusters, `sample.metacool`, by `metator dump`, in `tables` below.

## 1. The contigs database

The database holds only the contigs of the final clusters, those of
`binning.txt`, whose records the run wrote in `final_clusters/`, one FASTA file
per cluster:

```bash
# The tables of the run, from the map of its clusters
metator dump sample.metacool -o tables
mkdir -p anvio
cut -f 1 tables/binning.txt > anvio/contigs.txt
cat binning/final_clusters/*.fa > anvio/contigs.fa
```

anvi'o reads only contig names made of ASCII letters, digits and underscores:
rename the contigs of the assembly before the run if needed (see anvi'o's
[note on FASTA files](https://merenlab.org/2016/06/22/anvio-tutorial-v2/#take-a-look-at-your-fasta-file)).
Build the database without splitting the contigs (`-L 0`), and find the
single-copy genes, which give the completeness and redundancy of the clusters:

```bash
anvi-gen-contigs-database -f anvio/contigs.fa -o anvio/contigs.db -n metator -L 0 -T 16
anvi-run-hmms -c anvio/contigs.db -T 16
```

More annotations help to curate, and are optional:

```bash
anvi-setup-scg-taxonomy                    # once
anvi-run-scg-taxonomy -c anvio/contigs.db -T 16
anvi-setup-ncbi-cogs                       # once
anvi-run-ncbi-cogs -c anvio/contigs.db -T 16 --sensitive
anvi-display-contigs-stats anvio/contigs.db
```

## 2. The coverage profiles

anvi'o compares the coverage of the contigs across samples: give it each
library on its own (shotgun libraries can be added the same way). Its profiles
are made from BAM files sorted by position and indexed, which `anvi-init-bam`
writes from those of hicstuff: with `-n`, hicstuff keeps the alignments of
each mate in its temporary directory (`<outdir>/tmp` by default), as
`<prefix>.for.bam` and `<prefix>.rev.bam`. Profile those of the first mates,
as the second mates tell the same:

```bash
for lib in L1 L2; do
    anvi-init-bam hicstuff_${lib}/tmp/${lib}.for.bam -o anvio/${lib}.bam
    anvi-profile -i anvio/${lib}.bam -c anvio/contigs.db -o anvio/HiC_${lib} -S HiC_${lib} \
        -M 500 -T 16 --skip-hierarchical-clustering --skip-SNV-profiling \
        --contigs-of-interest anvio/contigs.txt
done
anvi-merge -c anvio/contigs.db -o anvio/merged -S HiC -W anvio/HiC_*/PROFILE.db
```

- `--skip-hierarchical-clustering` leaves anvi'o's hierarchical clustering to
  each cluster, when it is refined.
- `--skip-SNV-profiling` saves much time, without the single-nucleotide
  variants.
- `-M 500` keeps the contigs of 500 bp or more, as in the assembly.

## 3. The clusters of metaTOR

`binning.txt` is an anvi'o collection, a contig and its cluster per line:

```bash
anvi-import-collection tables/binning.txt -c anvio/contigs.db -p anvio/merged/PROFILE.db \
    -C metator --contigs-mode
```

## 4. Refining the clusters

`anvi-refine` opens a cluster in anvi'o's interactive interface, in a web
browser:

```bash
anvi-refine -c anvio/contigs.db -p anvio/merged/PROFILE.db -C metator -b metator_00001_00000
```

Two clusters of a human gut sample: the first cleaned of its contamination, the
second split into two clusters.

|A cluster cleaned|A cluster split in two|
|:-:|:-:|
|![A cluster cleaned in anvi'o](../images/anvio_example_12.png)|![A cluster split in two in anvi'o](../images/anvio_example_186.png)|

| Cluster | Before: completeness – contamination (%) | After |
|:-:|:-:|:-:|
| Cleaned | 97.2 – 19.7 | 95.8 – 9.9 |
| Split | 95.8 – 98.6 | 94.4 – 0.0, and 93.0 – 0.0 |

## 5. The curated clusters

Export the collection, as a contig (split) and its cluster per line:

```bash
anvi-export-collection -p anvio/merged/PROFILE.db -C metator -O anvio/curated
```

From `anvio/curated.txt`, metaTOR's Python API writes one FASTA file per
cluster, checks the clusters with miComplete as `metator validation` does, and
counts them in the quality classes of its log:

```python
import metacooler
import metator.validation as mtv
import pandas as pd

# The collection exported by anvi'o: each split and its cluster, a split read
# as its contig (NODE_1_split_00001 is NODE_1).
collection = pd.read_csv("anvio/curated.txt", sep="\t", header=None, names=["split", "cluster"], dtype=str)
binning = collection.set_index(collection["split"].str.split("_split").str[0])["cluster"]
binning = binning[~binning.index.duplicated(keep="last")]

# One FASTA file per curated cluster, the records of its contigs as they are
# in the assembly, written by metacooler from the contact map of the run.
groups = {bin_name: list(contigs.index) for bin_name, contigs in binning.groupby(binning)}
metacooler.write_fasta("sample.metacool", "assembly.fa", "anvio/fasta", groups=groups)

# miComplete on every cluster, with the bacterial and the archaeal markers.
mtv.micomplete_quality("anvio/fasta", "anvio/micomplete.txt", threads=16)

# The quality class of each curated cluster (HQ, MQ, LQ, Contaminated, Other).
summary = pd.read_csv("anvio/micomplete.txt", sep="\t", comment="#", index_col=0)
quality = summary.apply(lambda row: mtv.bin_quality(row["Weighted completeness"], row["Weighted redundancy"]), axis=1)
print(quality.value_counts())
```

To look at the contacts of a curated cluster, add the curated clusters to the
contig data table as a column, give them to the contact map of the run with
`metacooler reassign`, then extract a cluster with `metacooler extract` (see
[Other tools](../other-tools.md)):

```python
# binning: the curated cluster of each contig, read from the collection as above.
contigs = pd.read_csv("tables/contig_data_final.txt", sep="\t", dtype={"Name": str})
contigs["Anvio"] = contigs["Name"].map(binning).fillna("ND")
contigs.to_csv("binning/contig_data_anvio.txt", sep="\t", index=False)
```

```bash
metacooler reassign sample.metacool --contigs binning/contig_data_anvio.txt \
    --contigs-format Name,Anvio -o anvio/curated.metacool
metacooler extract anvio/curated.metacool --genomes metator_00001_00000 -o metator_00001_00000.cool
```
