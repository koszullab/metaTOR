---
title: metaTOR
hide:
  - navigation
  - toc
  - footer
---

<div class="cx-landing" markdown>

<h1 class="cx-masthead"><img src="assets/mark.svg" alt="" width="64" height="64"><span>metaTOR</span></h1>

<p class="cx-tagline">Metagenome binning from Hi-C contacts.</p>

[metaTOR](https://github.com/koszullab/metator) (Metagenomic Tridimensional
Organisation-based Reassembly) clusters the contigs of a metagenome assembly
into genomes from the contacts of a metaHiC library. Contigs of one genome
touch each other far more often than they touch other genomes, and metaTOR
finds those groups in the contig network, checks each cluster with miComplete,
partitions contaminated clusters again, and writes the contact map of the
clusters.
{ .cx-lede }

<div class="cx-install" markdown>

```bash
conda install bioconda::metator
```

</div>

<div class="cx-actions" markdown>

[Why metaTOR?](why.md){ .cx-button }
[Quickstart :lucide-arrow-right:](quickstart.md){ .cx-button .cx-button--filled }

</div>

<div class="cx-badges" markdown>

[![Bioconda](https://img.shields.io/conda/vn/bioconda/metator?color=blue)](https://anaconda.org/bioconda/metator)
[![PyPI version](https://img.shields.io/pypi/v/metator.svg)](https://pypi.org/project/metator/)
[![Python versions](https://img.shields.io/pypi/pyversions/metator.svg)](https://pypi.org/project/metator/)
[![CI](https://github.com/koszullab/metator/actions/workflows/ci.yml/badge.svg)](https://github.com/koszullab/metator/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/koszullab/metator/branch/master/graph/badge.svg)](https://codecov.io/gh/koszullab/metator)

</div>

<p class="cx-kicker">Start here</p>

<div class="cx-tiles" markdown>

[<span class="cx-glyph">:lucide-rocket:</span> <strong>Quickstart</strong> <span class="cx-tile-text">From pairs and an assembly to a contact map and clusters, on the test data.</span>](quickstart.md){ .cx-tile }
[<span class="cx-glyph">:lucide-lightbulb:</span> <strong>Concepts</strong> <span class="cx-tile-text">The network, core and overlapping clusters, validation.</span>](concepts.md){ .cx-tile }
[<span class="cx-glyph">:lucide-square-terminal:</span> <strong>CLI reference</strong> <span class="cx-tile-text">Every command and option.</span>](cli.md){ .cx-tile }
[<span class="cx-glyph">:lucide-braces:</span> <strong>Python API</strong> <span class="cx-tile-text">The functions behind the commands.</span>](api/network.md){ .cx-tile }

</div>

</div>
