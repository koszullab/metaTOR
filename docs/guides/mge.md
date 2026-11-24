# Clustering mobile genetic elements

The partition clusters plasmids and phages, the mobile genetic elements
(MGEs) of a community, poorly: their contigs are small, an MGE may contact
several hosts, and MGEs carry none of the marker genes that the validation
counts. `metator mge` clusters them on their own. Given the MGE contigs of
the assembly, found beforehand with an annotation tool, it groups those of
the same element from the Hi-C contacts between them, and stores the MGE
clusters in the contact map of the run. `metator host` then finds the MAGs
that host them.

![The clustering of the MGE contigs](../images/metator_mge.png)

## How it works

1. **Contacts.** Between two MGE contigs: the contacts of the map, summed by
   pair of contigs, whatever its bin size. Within a contig: only pairs whose
   mates are more than 1 kb apart, which the map's contact statistics count
   (metacooler's QC column `qc_intra_contig_long`).
2. **Scores.** The score of contigs *i* and *j* is their contacts divided by
   the geometric mean of the contacts within each:
   *A<sub>i,j</sub>* = *x<sub>i,j</sub>* / √(*x<sub>i,i</sub>* × *x<sub>j,j</sub>*).
3. **Merges.** The two contigs with the highest score are merged if it is
   above `--threshold-bin` (0.8). The merged cluster takes the contacts of
   both, and its scores are computed again. This repeats until no score is
   above the threshold. Ties are broken in the order of the MGE contigs:
   longest first, then by name.

Each contig never merged is an MGE cluster of its own. Clusters are numbered
in the order they are made, merged ones first, so a run gives the same MGE
clusters, with the same names, on the map of the contigs and on the map of
the clusters.

## Running it

```bash
# A map that metator network, pipeline or validation wrote in binning: updated in place
metator mge binning/sample.metacool -a assembly.fa -m mges.txt -o binning
```

`metator pipeline MAP --mges mges.txt` runs it with `--reassign` as its step
`mge`, then [`metator host`](#the-hosts) (see
[Mobile genetic elements and their hosts](pipeline.md#mobile-genetic-elements-and-their-hosts)).

| Input | What it is |
|---|---|
| `MAP` | A map that `metator network` stored its network in (the map that `network`, `pipeline` or `validation` wrote). It must hold the contact statistics of its contigs at 1001 bp, which `hicstuff pipeline -M metacool` and `metacooler create pairs` write with metacooler 0.4.0 or later (`--qc-threshold 1001` with an earlier one) |
| `--assembly` (`-a`) | The FASTA file the map was made from, plain or gzipped |
| `--mges` (`-m`) | The MGE contigs, one name per line. Those missing from the map are left out, with a warning |

| Option | Default | Effect |
|---|---|---|
| `--threshold-bin` (`-s`) | 0.8 | The score above which two contigs or clusters merge. Lower merges more |
| `--random` (`-r`) | off | Shuffle the contigs between the clusters, keeping the number of contigs of each: a random clustering to compare with. Stored as the round `mge_random`, written to `fasta/mges_random/` |
| `--reassign` | off | Make each MGE cluster a genome of the map (see [MGE genomes](#mge-genomes)). Not with `--random` |
| `--outdir` (`-o`) | current directory | Where `fasta/mges/` and the map go |
| `--inplace` | off | Update `MAP` itself rather than the map of the output directory |

A map that `mge` cannot use is refused before anything is written, with the
reason and how to fix it (see [`mge`](../cli.md#mge)).

## Outputs

| Output | What it holds |
|---|---|
| `fasta/mges/` | Two FASTA files per MGE cluster (`MetaTOR_MGE_00001`...): `<name>.fa`, its contigs, and `<name>.stitched.fa`, one record with the contigs joined by 200 `N`. Longest contig first |
| The round `mge` of the map | One cluster per MGE cluster: its contigs, length, GC content, hits and clustering score (`BinningScore`) |
| `sample.metator_<date>.log` | The log, in the current directory (`sample` for `sample.metacool`) |

The clustering score of an MGE cluster is the lowest score at which its
contigs merged (empty for a single contig). The stitched file holds one
sequence, for the tools that check the quality of an MGE genome and expect
one sequence per genome.

`mge` writes no summary table: the round `mge` holds it, and
[`metator dump`](../cli.md#dump) writes it out. `-F` replaces the
`fasta/mges/` of an earlier run.

## In the contact map of the run

`mge` adds the MGE clusters to the
[clustering rounds](metacool.md#clustering-rounds) of the map, as the round
`mge`, and keeps the other rounds (`core`, `overlapping`, `recursive_<k>`,
`final`). The entry `mge` of the `metator` key of the map's metadata records
the run: version, command line, parameters, and the columns it read. The
genomes and contigs tables are left as they were, unless `--reassign`.

```python
import metacooler

mc = metacooler.MetaCooler("binning/sample.metacool")
mc.clusters[mc.clusters["round"] == "mge"]  # one row per MGE cluster
mc.cluster_members("mge")                   # their contigs, in map order
mc.metadata["metator"]["mge"]               # the run
```

- `metator dump` writes the round `mge` (and `mge_random`) in
  `clusters.txt` and `cluster_members.txt`.
- A second run replaces the round `mge` only; `--random` replaces only
  `mge_random`.
- The map goes to the output directory, or `MAP` is updated in place when it
  is already there or with `--inplace` (see
  [Step by step](metacool.md#step-by-step)).
- `mge` reads no cluster, so it runs on the map of the contigs (after
  `network`) as well as on the map of the clusters. But `partition`,
  `validation` and `pipeline` drop the rounds `mge` and `mge_random`, with a
  warning: run `mge` again on their map (`pipeline --mges` does it itself).

## MGE genomes

With `--reassign`, each MGE cluster also becomes a genome of the map, named
after it, with its `GC`, `Hit` and `BinningScore` in the genomes table. Its
contigs leave the genomes that held them (MAGs or `__unbinned__`).

- The MAGs keep their other contigs and their validation values, which still
  describe them with the contigs they lost. The log names each MAG that loses
  contigs, and how many. The marker genes that `metator markers` stored on
  their contigs, found in the MAGs with those contigs, give an estimate of
  their numbers without them (see
  [Marker genes of the contigs](../formats.md#marker-genes-of-the-contigs)).
- A MAG made only of MGE contigs is no longer a genome, with a warning; the
  round `final` keeps it.
- The column `mge` of the genomes table is `True` for MGE genomes and `False`
  for MAGs (as text).

```bash
metator mge binning/sample.metacool -a assembly.fa -m mges.txt -o binning --reassign
metator dump binning/sample.metacool -o tables
```

The tables of `metator dump` follow the genomes: `binning.txt`,
`contig_data_final.txt` (`Final_cluster`) and `cluster_summary.txt` list the
MGE genomes next to the MAGs.

```python
import metacooler

mc = metacooler.MetaCooler("binning/sample.metacool")
mc.genomes[mc.genomes["mge"] == "True"]  # the MGE genomes, with GC, Hit and BinningScore
```

A run that replaces the round `mge` (any run but `--random`) first undoes
the MGE genomes of the map: their contigs go back to their MAG or to
`__unbinned__`. So a second run with `--reassign` replaces the MGE genomes,
and a run without it restores the MAGs (a MAG that had lost all its contigs
comes back without values, with a warning: run `validation` again to measure
them).

## The hosts

`metator host` finds the MAGs that host the MGE clusters: the MAGs their
contigs are in contact with. It reads the map that `metator mge` stored its
MGE clusters in, and stores the hosts there.

```bash
# The MGE clusters of the map of a run in binning, then their hosts: updated in place
metator mge binning/sample.metacool -a assembly.fa -m mges.txt -o binning
metator host binning/sample.metacool -o binning
```

1. **The MAGs** are the genomes of the map: the final clusters of
   `validation` or `pipeline`, or the bins of another binner, made genomes
   with `metacooler reassign`. `__unbinned__` and the MGE genomes are not
   MAGs.
2. **The MGE contigs** count for their MGE cluster only. A MAG that holds
   some (a plasmid binned with its host) is taken without them.
3. **The signal** of an MGE cluster with a MAG is the sum of the weights of
   the network edges between their contigs, with the normalisation that
   `metator network` stored (`empirical_hit` by default, which corrects for
   coverage). The number of contigs of the MAG in contact is counted too.
4. **A MAG is a host** when its signal is at least `--min-share` percent
   (10) of the signal of the MGE cluster with all MAGs, and at least
   `--min-contigs` (5) of its contigs are in contact with the MGE cluster. An
   MGE cluster can have one host, several or none.

| Option | Default | Effect |
|---|---|---|
| `--min-share` (`-s`) | 10 | The share of the MGE cluster's signal with the MAGs, in percent, a host needs |
| `--min-contigs` (`-c`) | 5 | The number of contigs of a MAG in contact with the MGE cluster a host needs |
| `--outdir` (`-o`) | current directory | Where `tables/mge_hosts.txt`, `plots/mge_hosts.pdf`, the log and the map go |
| `--inplace` | off | Update `MAP` itself rather than the map of the output directory |

### Reading the hosts

- **A small share of all contacts.** The share counts only the contacts with
  MAGs. An MGE cluster with few such contacts, such as a phage caught in free
  particles, can get a host from a small part of its contacts: on a mock
  community, a free 42 kb phage got a host from 1% of its contacts.
  `share_of_all` in `tables/mge_hosts.txt` is the share of all its contacts;
  when it is small, look at the contact map.
- **Several hosts are not always several species.** A genome split into
  several clusters gives each a share of the signal: on the mock community,
  12 of the 14 MGE clusters with several hosts had clusters of one genome as
  hosts. Related plasmids, assembled into shared contigs, also give a host for
  each genome that carries one. Check the taxonomy of the hosts, or their
  contacts with each other, before reading a broad host range.
- **No noise filter.** `metator host` does not estimate the noise of the
  library, nor keep a MAG as the only host when it has most of the signal.
  The rule of [Bignaud et al. 2025](../citation.md), which adds both, is left
  to the analysis of `tables/mge_hosts.txt`.

### What host writes

| Output | What it holds |
|---|---|
| `tables/mge_hosts.txt` | One row per MGE cluster and MAG in contact: their signal, the contigs of the MAG in contact, the share of the signal (`share`) and of all contacts (`share_of_all`), and whether the MAG is a host (see [Hosts of the MGE clusters](../formats.md#hosts-of-the-mge-clusters)) |
| `plots/mge_hosts.pdf` | The distribution of the signal of the hosts (log scale) and of their share |
| The round `mge` of the map | The number of hosts of each MGE cluster, their names by decreasing share, and their shares (`n_hosts`, `hosts`, `host_shares`) |
| The genomes table of the map | The number of MGE clusters each MAG hosts, and their names (`n_mges`, `mges`) |
| `sample.metator_<date>.log` | The log, in the output directory |

The entry `host` of the `metator` key of the map's metadata records the run,
with the normalisation of the network it read and its digest
(`network_digest`).

```python
import metacooler

mc = metacooler.MetaCooler("binning/sample.metacool")
mge = mc.clusters[mc.clusters["round"] == "mge"]
mge[["name", "n_hosts", "hosts", "host_shares"]]  # the hosts of each MGE cluster
mc.genomes[["name", "n_mges", "mges"]]              # the MGE clusters of each MAG
mc.metadata["metator"]["host"]                      # the run
```

`metator dump` writes the hosts of each MGE cluster in `clusters.txt`, and
the MGE clusters of each MAG in `cluster_summary.txt`.

A second `host` run replaces the hosts. A run that changes what they were
found from (the round `mge`, the genomes or the network) drops them, with a
warning: `mge` but `--random`, `partition`, `validation`, `pipeline`, and
`network` with another network. Run `host` again then; `pipeline --mges`
finds them again itself.

Marbouty et al., 2021, and Bignaud et al., 2025, study the phages of a
community and their hosts (see [Citing metaTOR](../citation.md)).
