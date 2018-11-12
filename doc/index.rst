.. metator documentation master file, created by
   sphinx-quickstart on Sun Nov 11 23:23:52 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to metator's documentation!
===================================

MetaTOR is a metagenomic binning pipeline based on 3C data. The pipeline is
decomposed into four steps:

- Alignment
- Partitioning
- Annotation
- Binning

See the `user manual
<https://github.com/koszullab/metaTOR/blob/master/meta3c_manual.pdf>`_ for
detailed instruction on how to run the pipelines and the parameters to use.

MetaTOR is also a Python library with various modules and functions devoted to
the handling of meta3C data. The library is used by the above pipeline, and the
present documentation exposes its API for advanced usage.

.. toctree::
   :maxdepth: 3
   :caption: Contents:


Reference API
=============

.. toctree::
    :maxdepth: 3

    api/modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
