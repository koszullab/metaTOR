# The end-to-end pipeline

```bash
## 1. Getting the assembly
# Assemble the shotgun reads of the sample with MEGAHIT (or metaSPAdes:
# spades.py --meta), keeping the contigs of 500 bp or more.
megahit -1 shotgun_R1.fq.gz -2 shotgun_R2.fq.gz --min-contig-len 500 -t 16 -o megahit
ln -s megahit/final.contigs.fa sample.fa

## 2. Mapping the Hi-C pairs
# Align the Hi-C reads on the assembly and make the contact map of its
# contigs: hicstuff/sample.metacool.
hicstuff pipeline -g sample.fa -e DpnII,HinfI -M metacool -t 16 \
    -o hicstuff -P sample hic_R1.fq.gz hic_R2.fq.gz

## 3. Annotating MGEs
# Find the viral and plasmid contigs with geNomad, then list them, one name
# per line (proviruses, which are part of a host contig, are left out).
genomad download-database .
genomad end-to-end --cleanup --splits 8 sample.fa genomad genomad_db
tail -q -n +2 genomad/sample_summary/sample_{virus,plasmid}_summary.tsv \
    | cut -f 1 | grep -v '|provirus' > mges.txt

## 4. Clustering and assigning hosts with metaTOR
# Bin the contigs into MAGs, cluster the MGE contigs, and find their hosts.
metator pipeline hicstuff/sample.metacool -a sample.fa -e DpnII,HinfI \
    --mges mges.txt -t 16 -o metator
```

This takes a metagenome from its reads to MAGs, MGE clusters and their hosts.
It needs a shotgun library (for the assembly) and a Hi-C library of the same
sample.

1. **Assembly.** metaTOR does not assemble. Contigs shorter than 500 bp carry
   too few contacts to be clustered. MEGAHIT refuses an existing `-o`
   directory.
2. **Contact map.** `hicstuff pipeline -M metacool` (hicstuff 4.0.0 or later)
   writes `hicstuff/sample.metacool`. Give it the enzymes of the library.
   metaTOR can also make the map itself from pairs or BAM files: see
   [The inputs of a run](cli.md#inputs).
3. **MGEs.** geNomad writes its tables in `genomad/sample_summary/`, named
   after the assembly. `--mges` takes any list of contig names, one per line,
   so another annotation tool works too.
4. **metaTOR.** `metator pipeline` builds the network, partitions it,
   validates the clusters with miComplete, clusters the MGE contigs and finds
   the MAGs that host them. Its outputs, in `metator/`:

    | Output | What it holds |
    |---|---|
    | `sample.metacool` | The contact map of the clusters |
    | `fasta/` | The MAGs and the MGE clusters, as FASTA files |
    | `cool/` | A cool file per MAG, with the MGE clusters it hosts |
    | `tables/`, `plots/` | The tables and plots of the run, among them `tables/mge_hosts.txt` |
    | `sample.report.html` | The report of the run (with `sample.report.txt` and `sample.report_mqc.yaml`, for MultiQC) |
    | `sample.metator_<date>.log` | The log |

See [Running the pipeline](guides/pipeline.md) for the options, and
[Clustering mobile genetic elements](guides/mge.md) for the MGE steps.
