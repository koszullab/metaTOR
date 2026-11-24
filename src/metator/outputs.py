#!/usr/bin/env python3
# coding: utf-8

"""The files of a run beside its map: where they go in its output directory,
the FASTA files of its genomes and MGE clusters, and the cool files of its
MAGs.

Every command writes into its output directory (-o), under these names:

- ``fasta/mags_with-mges/``: the FASTA file of each final cluster, as the
  validation made it, with the MGE contigs that it holds (validation,
  pipeline);
- ``fasta/mags/``: each MAG of the map of the run, unstitched and stitched
  (see write_fasta_files): its genomes but ``__unbinned__`` and the MGE
  genomes, which hold no MGE contig once ``metator mge --reassign`` made the
  MGE clusters genomes (pipeline);
- ``fasta/mges/``: each MGE cluster, unstitched and stitched (mge, pipeline
  --mges); ``fasta/mges_random/`` for the random clusters of mge --random;
- ``cool/``: one cool file per MAG, with the MGE clusters that it hosts,
  balanced (see write_cool_files; pipeline);
- ``tables/``: the tables of metator dump (pipeline), and ``mge_hosts.txt``
  (host, pipeline --mges);
- ``plots/``: the plots, as PDF files (see metator.figures and
  metator.host).
"""

import os
import shutil
import tempfile
from os.path import join

import numpy as np

import metator.clusters as mtcl
import metator.metacool as mtmc
from metator.log import logger

FASTA_DIR = "fasta"
MAGS_WITH_MGES_DIR = join(FASTA_DIR, "mags_with-mges")
MAGS_DIR = join(FASTA_DIR, "mags")
MGES_DIR = join(FASTA_DIR, "mges")
MGES_RANDOM_DIR = join(FASTA_DIR, "mges_random")
COOL_DIR = "cool"
TABLES_DIR = "tables"

#: The FASTA file of a genome or a cluster as one record, beside the one of
#: its contigs (<name>.fa): <name>.stitched.fa.
STITCHED = ".stitched.fa"

#: N between two contigs of a stitched record, as metaTOR 1.x joined the
#: contigs of an MGE cluster in mges_binned.fa.
SPACER = 200


def write_fasta_files(metacool, assembly, outdir, genomes=None, groups=None, temp_dir=None):
    """Write the FASTA files of genomes of a map, or of groups of its contigs, in
    outdir (made if missing), two for each, side by side:

    - ``<name>.fa``: its contigs, the records of the assembly as they are
      (headers and lines), in the order of the map for a genome, in the order
      given for a group;
    - ``<name>.stitched.fa``: one record named after it, its contigs (the
      longest first for a genome, ties in the order of the map; in the order
      given for a group) joined by SPACER N, in lines of 80 bases.

    metacooler writes both (metacooler.write_fasta, once each: each reads the
    assembly); the stitched files are written in a directory of their own in
    outdir, then renamed. A file of that name in outdir is replaced.

    Parameters
    ----------
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection: the map whose contigs they are.
    assembly : str
        Path to the assembly the map was made from: FASTA, plain or compressed
        with gzip.
    outdir : str
        The directory to write them in.
    genomes : list of str, optional
        The genomes of the map to write, by name.
    groups : dict, optional
        Name -> the names of its contigs, in their order: groups of contigs
        to write in place of genomes (the clusters of a run that the map does
        not hold yet).
    temp_dir : str, optional
        Where a compressed assembly is copied, plain, as it is read: the
        temporary directory of the run.

    Returns
    -------
    list of str
        The genomes or groups written, in order: none for none given.

    Raises
    ------
    metator.metacool.MetacoolError
        If the files cannot be written: the assembly lacks a contig, gives one
        another length than the map, or is not FASTA (a ValueError, KeyError or
        OSError of metacooler, message kept)."""
    import metacooler

    names = list(genomes) if genomes is not None else None
    groups = dict(groups or {}) if genomes is None else None
    # metacooler writes every genome for none named.
    if not (names if names is not None else groups):
        return []
    os.makedirs(outdir, exist_ok=True)
    with mtmc._metacooler_logging(), mtmc._as_metacool_errors(), mtmc.read_errors(metacool):
        written = metacooler.write_fasta(metacool, assembly, outdir, genomes=names, groups=groups, force=True, temp_dir=temp_dir)
        stitched = tempfile.mkdtemp(prefix=".stitched.", dir=outdir)
        try:
            files = metacooler.write_fasta(
                metacool,
                assembly,
                stitched,
                genomes=names,
                groups=groups,
                stitch=True,
                spacer=SPACER,
                force=True,
                temp_dir=temp_dir,
            )
            for name, path in files.items():
                os.replace(path, join(outdir, f"{name}{STITCHED}"))
        finally:
            shutil.rmtree(stitched, ignore_errors=True)
    return list(written)


def map_mags(mc):
    """The MAGs of a map (mc, opened), in the order of its genomes table: its
    genomes but __unbinned__ and the MGE genomes, those that its column mge
    says are and those named after a cluster of its round mge (see
    metator.host.map_bins)."""
    names = [str(name) for name in mc.genomes["name"]]
    mges = set(mtmc.mge_genomes(mc))
    if mtcl.MGE in [entry.get("name") for entry in mc.cluster_rounds]:
        mges |= set(mc.cluster_members(mtcl.MGE)["cluster"].astype(str))
    return [name for name in names if name != mtmc.UNBINNED_GENOME and name not in mges]


def hosted_mges(mc):
    """The MGE clusters that each MAG of a map hosts (mc, opened), by MAG, as
    metator host stores them in its genomes table (its column mges, names
    separated by commas): none for a map without hosts."""
    genomes = mc.genomes
    if "mges" not in genomes.columns:
        return {}
    return {
        str(name): [mge for mge in str(mges).split(",") if mge and mge != "nan"]
        for name, mges in zip(genomes["name"], genomes["mges"], strict=True)
        if isinstance(mges, str) and mges
    }


def write_cool_files(metacool, outdir, mags, hosted=None):
    """Write the cool file of each MAG of a map in outdir (made if missing):
    ``<MAG>.cool``, the contacts of its contigs and of those of the MGE
    clusters that it hosts (hosted: their names, by MAG, the MGE genomes of
    the map), as metacooler extract writes them, at the bins of the map, then
    balanced as cooler balance balances a cool file, the whole of it at once,
    with its defaults (the first two diagonals left out, bins with fewer than
    10 non-zero pixels masked...): the weight column of its bins. A cool file
    that cannot be balanced is a warning, and is kept without weights; one
    whose every bin is masked (too few bins, or contacts, for cooler balance)
    is a warning too, its weights all NaN.

    Returns
    -------
    list of str
        The cool files written, in the order of mags.

    Raises
    ------
    metator.metacool.MetacoolError
        If metacooler cannot extract a MAG (a ValueError, KeyError or OSError
        of metacooler, message kept)."""
    import cooler
    import metacooler

    os.makedirs(outdir, exist_ok=True)
    written = []
    for mag in mags:
        path = join(outdir, f"{mag}.cool")
        with mtmc._metacooler_logging(), mtmc._as_metacool_errors(), mtmc.read_errors(metacool):
            metacooler.extract(metacool, genomes=[mag, *(hosted or {}).get(mag, [])], output=path)
        try:
            weights, _ = cooler.balance_cooler(cooler.Cooler(path), store=True)
        except Exception as error:
            logger.warning(f"{path} is not balanced, and keeps its raw contacts only: {type(error).__name__}: {error}")
        else:
            if not np.isfinite(weights).any():
                logger.warning(
                    f"Every bin of {path} is masked by its balancing (too few bins, or contacts, for cooler balance): its "
                    "weights are all NaN, and its raw contacts are kept."
                )
        written.append(path)
    return written
