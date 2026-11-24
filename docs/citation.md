# Citing metaTOR

Please cite:

> Lyam Baudry, Théo Foutel-Rodier, Agnès Thierry, Romain Koszul and Martial
> Marbouty. 2019. MetaTOR: A Computational Pipeline to Recover High-Quality
> Metagenomic Bins From Mammalian Gut Proximity-Ligation (meta3C) Libraries.
> *Frontiers in Genetics* 10: 753.
> <https://doi.org/10.3389/fgene.2019.00753>

```bibtex
@article{baudry2019metator,
  author  = {Baudry, Lyam and Foutel-Rodier, Th{\'e}o and Thierry, Agn{\`e}s and
             Koszul, Romain and Marbouty, Martial},
  title   = {{MetaTOR}: A Computational Pipeline to Recover High-Quality
             Metagenomic Bins From Mammalian Gut Proximity-Ligation (meta3C)
             Libraries},
  journal = {Frontiers in Genetics},
  volume  = {10},
  pages   = {753},
  year    = {2019},
  doi     = {10.3389/fgene.2019.00753}
}
```

## Related work

The meta3C method, and studies that use it:

- [Metagenomic chromosome conformation capture (meta3C) unveils the diversity of chromosome organization in microorganisms](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4381813/),
  Martial Marbouty, Axel Cournac, Jean-François Flot, Hervé Marie-Nelly,
  Julien Mozziconacci and Romain Koszul. *eLife*, 2014.
- [Scaffolding bacterial genomes and probing host-virus interactions in gut microbiome by proximity ligation (chromosome capture) assay](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5315449/),
  Martial Marbouty, Lyam Baudry, Axel Cournac and Romain Koszul.
  *Science Advances*, 2017.
- [MetaHiC phage-bacteria infection network reveals active cycling phages of the healthy human gut](https://elifesciences.org/articles/60608),
  Martial Marbouty, Agnès Thierry, Gaël A. Millot and Romain Koszul.
  *eLife*, 2021.
- [Phages with a broad host range are common across ecosystems](https://hal.science/pasteur-05269627),
  Amaury Bignaud, Devon Conti, Agnès Thierry, Jacques Serizay, Karine
  Labadie, Julie Poulain, Olivia Cheny, Maritrini Colón-González, Laurent
  Debarbieux, Marianna Guerrero-Osornio, Sophie Helaine, Peter Hill,
  Gwenaelle Le Tinier, Gael A. Millot, Lucia Morales, Andrés Parada, Nadia
  Riera, Gregorio Iraola, Romain Koszul and Martial Marbouty.
  *Nature Microbiology*, 2025.

## The tools metaTOR runs

Please also cite the tools a run relies on:

- Validation: [miComplete](https://pypi.org/project/micomplete/) (Hugoson,
  Lam and Guy, *Bioinformatics* 2020), with Prodigal (Hyatt et al., 2010) and
  HMMER (Eddy, 2011) through [pyrodigal](https://doi.org/10.21105/joss.04296)
  (Larralde, *JOSS* 2022) and
  [pyhmmer](https://doi.org/10.1093/bioinformatics/btad214) (Larralde et al.,
  *Bioinformatics* 2023).
- Partition: [Louvain](https://doi.org/10.1088/1742-5468/2008/10/P10008)
  (Blondel et al., 2008) or [Leiden](https://doi.org/10.1038/s41598-019-41695-z)
  (Traag, Waltman and van Eck, 2019).
- Contact map from the reads: [hicstuff](https://koszullab.github.io/hicstuff/citation/),
  with [bowtie2](https://github.com/BenLangmead/bowtie2) or
  [bwa](https://github.com/lh3/bwa).
- Contact map of the clusters: [metacooler](https://github.com/js2264/metacooler)
  and [cooler](https://github.com/open2c/cooler) (Abdennur and Mirny,
  *Bioinformatics* 2020).
