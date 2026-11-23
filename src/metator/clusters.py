#!/usr/bin/env python3
# coding: utf-8

"""The clustering rounds of a run: what the map stores in its /clusters
group (see metacooler's format), and the two tables a run writes them to.

A run of metaTOR reaches its final clusters in rounds:

- ``core``: the core clusters of the partition, the contigs that every
  iteration of Louvain or Leiden put together, named ``core_<number>``
  (:func:`core_name`), with no parent: an overlapping cluster merges several
  of them;
- ``overlapping``: every overlapping cluster of the partition, its core
  clusters merged, named ``<prefix>_<number>_00000``
  (:func:`overlapping_name`); the validation checks with miComplete those
  that reach ``--size``;
- ``recursive_<k>``: the clusters that step k of the validation makes by
  partitioning a contaminated cluster again, its parent;
- ``final``: the final clusters, each the cluster of the same name of the
  round it comes from.

``metator partition MAP`` stores the first two in the map
(:func:`partition_rounds`, :func:`metator.metacool.write_rounds`), and the
validation reads them back (:func:`metator.metacool.read_partition`): the
map is the partition that the validation checks. The validation records
every round as it makes them (:class:`Rounds`) and writes them next to its
other tables: ``clusters.txt``, one row per cluster of every round, with its
parent and the columns of miComplete, and ``cluster_members.txt``, one row
per contig of a cluster. The map of the clusters stores both
(:func:`metator.metacool.reassign_metacool`).

``metator mge`` makes one more round after them, ``mge``: its MGE clusters,
which it stores in the map next to those of the run (:func:`mge_round`,
:func:`metator.metacool.write_round`).
"""

import csv
from dataclasses import dataclass
from typing import Any
from os.path import exists, join

import pandas as pd

from metator.log import logger


#: The two tables of the rounds, in the output directory of a run.
CLUSTERS_FILE = "clusters.txt"
CLUSTER_MEMBERS_FILE = "cluster_members.txt"

#: The keywords of pandas.read_csv with which metaTOR reads back the decimal
#: numbers of the tables it wrote: exactly, each the number that was written.
#: Under pandas 2, its round_trip parser: the default one may read the last
#: digit of a number off (24.621666666666666 as 24.621666666666663), and the
#: map would then store, and metator dump write, another number than the
#: run's tables. Under pandas 3, none: its default parser reads them
#: exactly, and it deprecates the keyword.
EXACT_FLOATS: dict[str, Any] = {"float_precision": "round_trip"} if int(pd.__version__.split(".")[0]) < 3 else {}

#: The names of the rounds.
CORE = "core"
OVERLAPPING = "overlapping"
FINAL = "final"
MGE = "mge"

#: The columns of the summary of metator mge (mges_cluster_summary.tsv) that the
#: round mge does not store as columns of its own: the name of each MGE
#: cluster, and what metacooler computes from its contigs (their total length
#: and their number) or stores as its members.
MGE_SUMMARY_KEYS = {"BinName": "name", "BinLength": "length", "ContigsNumber": "n_contigs", "Contigs": "members"}


def recursive_round(step):
    """The name of the round of recursive step *step* (1, 2...)."""
    return f"recursive_{step}"


def core_name(number):
    """The name of the core cluster *number* (1, 2...) of a partition:
    ``core_00001``..., as the ``Core_cluster_ID`` of the contig data tables
    numbers it (``00001``...)."""
    return f"{CORE}_{number:05d}"


def overlapping_name(prefix, number):
    """The name of the overlapping cluster *number* (1, 2...) of a partition:
    ``<prefix>_<number>_00000``, its FASTA file's name and its name in
    miComplete's results, as the validation names its recursive clusters
    ``<prefix>_<number>_<recursive number>`` (*prefix*: ``metator``, or
    ``<--prefix>_metator``)."""
    return f"{prefix}_{number:05d}_{0:05d}"


def cluster_number(round_name, name):
    """The number of a cluster of the round ``core`` or ``overlapping`` from
    its name (see :func:`core_name` and :func:`overlapping_name`): the number
    that the ``Core_cluster_ID`` or the ``Overlapping_cluster_ID`` of the
    contig data tables gives it, as an int. None for a name of another form:
    a cluster that metaTOR's partition did not name."""
    if round_name == CORE:
        number = name[len(CORE) + 1 :] if name.startswith(f"{CORE}_") else ""
    else:
        parts = name.rsplit("_", 2)
        number = parts[1] if len(parts) == 3 and parts[0] and parts[2] == f"{0:05d}" else ""
    return int(number) if len(number) >= 5 and number.isascii() and number.isdigit() else None


@dataclass
class Partition:
    """The clusters of a partition of the network of a map (see
    :func:`metator.partition.partition`): its core and its overlapping
    clusters, by number (from 1), each the IDs of its contigs (their
    positions in the map, from 1, as in the network: see
    :func:`metator.network.map_network`), in the order that the partition
    lists them, and the name of each overlapping cluster (see
    :func:`overlapping_name`).

    The order of the contigs of an overlapping cluster is that of its core
    clusters, by number, then that of their contigs, by ID: the order of the
    records of the FASTA file of the cluster that miComplete checks. A map
    stores the partition as the rounds ``core`` and ``overlapping`` (see
    :func:`partition_rounds`), which give it back
    (:func:`metator.metacool.read_partition`).
    """

    core: dict
    overlapping: dict
    names: dict

    def sizes(self, lengths):
        """The size in bp of each overlapping cluster, by number: the sum of
        the lengths of its contigs (*lengths*, by contig ID - 1)."""
        return {number: int(sum(lengths[contig - 1] for contig in ids)) for number, ids in self.overlapping.items()}


def record_partition(rounds, partition, names, quality=None):
    """Record the rounds ``core`` and ``overlapping`` of a partition in
    *rounds* (a :class:`Rounds`): every core cluster, named
    :func:`core_name`, then every overlapping cluster, named as the
    partition names it, each with its contigs and no parent, in the order of
    their numbers. *names*: the name of each contig, by ID - 1. *quality*:
    the columns of miComplete of each overlapping cluster that it checked,
    by name (none for the others)."""
    quality = quality or {}
    for number, ids in sorted(partition.core.items()):
        rounds.add(CORE, core_name(number), [names[contig - 1] for contig in ids])
    rounds.close_round(CORE)
    for number, ids in sorted(partition.overlapping.items()):
        name = partition.names[number]
        rounds.add(OVERLAPPING, name, [names[contig - 1] for contig in ids], columns=quality.get(name))
    rounds.close_round(OVERLAPPING)
    return rounds


def partition_rounds(partition, names):
    """The rounds ``core`` and ``overlapping`` of a partition, as
    ``metator partition MAP`` stores them in the map: every cluster, with its
    contigs (by name: *names*, the name of each contig by ID - 1), and no
    column of miComplete (see :func:`record_partition`).

    Returns
    -------
    tuple of pandas.DataFrame
        The clusters and their contigs (see :meth:`Rounds.tables`).
    """
    return record_partition(Rounds(), partition, names).tables()


class Rounds:
    """The clustering rounds of a run, in the order they are made.

    Each cluster is recorded once per round, with its contigs (by name), its
    parent (a cluster of an earlier round, by name) and its columns (those of
    miComplete). The parent of a cluster is looked up in the latest round
    that holds a cluster of that name, as a recursive cluster replaces the
    cluster it was split from.
    """

    def __init__(self):
        self.clusters = []
        self.members = []
        self._latest = {}

    def add(self, round_name, name, contigs, parent=None, columns=None):
        """Record the cluster *name* of *round_name*.

        Parameters
        ----------
        round_name : str
            The round.
        name : str
            The name of the cluster.
        contigs : iterable of str
            The names of its contigs.
        parent : str, optional
            The name of the cluster it was split from, in an earlier round.
        columns : dict, optional
            Its other columns (those of miComplete).
        """
        parent_round = self._latest.get(parent) if parent is not None else None
        if parent is not None and parent_round is None:
            raise ValueError(f"The parent of {name}, {parent}, is not a cluster of an earlier round.")
        self.clusters.append(
            {
                "round": round_name,
                "name": name,
                "parent_round": parent_round,
                "parent": parent if parent_round is not None else None,
                **(columns or {}),
            }
        )
        self.members.extend((round_name, name, contig) for contig in contigs)

    def close_round(self, round_name):
        """Make the clusters of *round_name* the latest ones of their names,
        once the round is recorded: its clusters' parents are of earlier
        rounds."""
        for cluster in self.clusters:
            if cluster["round"] == round_name:
                self._latest[cluster["name"]] = round_name

    def tables(self):
        """The two tables of the rounds: one row per cluster of every round,
        with its parent (missing for none) and its columns, and one row per
        contig of a cluster.

        Returns
        -------
        tuple of pandas.DataFrame
            The clusters and their contigs, with the columns of the tables
            that :meth:`write` writes.
        """
        clusters = pd.DataFrame(self.clusters, columns=_columns(self.clusters))
        members = pd.DataFrame(self.members, columns=["round", "cluster", "contig"])
        return clusters, members

    def write(self, outdir):
        """Write the two tables of the rounds in *outdir*: one row per cluster
        of every round, and one row per contig of a cluster.

        Returns
        -------
        tuple of str
            The paths of the two tables.
        """
        clusters_file = join(outdir, CLUSTERS_FILE)
        members_file = join(outdir, CLUSTER_MEMBERS_FILE)
        write_round_tables(*self.tables(), clusters_file, members_file)
        return clusters_file, members_file


def write_round_tables(clusters, members, clusters_file, members_file):
    """Write the two tables of clustering rounds (see :meth:`Rounds.tables`):
    the clusters to clusters_file, a missing value written NA, and their
    contigs to members_file. :meth:`Rounds.write` writes those of a run so,
    and ``metator dump`` those of a map (see :func:`metator.dump.dump`)."""
    clusters.to_csv(clusters_file, sep="\t", index=False, na_rep="NA", quoting=csv.QUOTE_NONE)
    members.to_csv(members_file, sep="\t", index=False, quoting=csv.QUOTE_NONE)
    logger.info(
        f"{len(clusters)} clusters in {clusters['round'].nunique()} rounds are written to {clusters_file} and {members_file}."
    )


def _columns(clusters):
    """The columns of the clusters table: the round, the name and the parent
    first, then the others, in the order they come."""
    columns = ["round", "name", "parent_round", "parent"]
    for cluster in clusters:
        columns.extend(column for column in cluster if column not in columns)
    return columns


def read_rounds(outdir):
    """The two tables of the rounds that a run wrote in *outdir*, as
    :meth:`Rounds.write` writes them.

    Returns
    -------
    tuple of pandas.DataFrame, or None
        The clusters and their contigs, or None when the run wrote neither.

    Raises
    ------
    ValueError
        If only one of the two tables is there.
    """
    paths = [join(outdir, CLUSTERS_FILE), join(outdir, CLUSTER_MEMBERS_FILE)]
    found = [exists(path) for path in paths]
    if not any(found):
        return None
    if not all(found):
        missing = paths[found.index(False)]
        raise ValueError(f"{missing} is missing: the clustering rounds of the run need both of {', '.join(paths)}.")
    # Names as written: a contig called NA, None or nan is a contig. Only
    # the missing values of the clusters table are NA (Rounds.write). Its
    # numbers are read back exactly (EXACT_FLOATS), as the map stores them.
    clusters = pd.read_csv(
        paths[0],
        sep="\t",
        dtype={"name": str, "parent": str, "round": str, "parent_round": str},
        keep_default_na=False,
        na_values=["NA"],
        **EXACT_FLOATS,
    )
    members = pd.read_csv(paths[1], sep="\t", dtype=str, keep_default_na=False)
    return clusters, members


def mge_round(summary_file):
    """The round ``mge``: the MGE clusters of a run of ``metator mge``, from
    the summary that it wrote (``mges_cluster_summary.tsv``), in its order: each
    named as its record (``MetaTOR_MGE_<ID>``), with its contigs (the
    ``Contigs`` column) and no parent. Its other columns are those of the
    summary, as written (GC, Hit, BinningScore, MetagenomicBin,
    AssociationScore, AssociatedBins), but those that the map computes from
    the contigs (see :data:`MGE_SUMMARY_KEYS`): ``BinLength`` is the length
    of the cluster, and ``ContigsNumber`` its number of contigs. ``NA`` is a
    missing value: no number, or an empty text in the map.

    Parameters
    ----------
    summary_file : str
        Path to ``mges_cluster_summary.tsv``.

    Returns
    -------
    tuple of pandas.DataFrame
        The clusters and their contigs (see :meth:`Rounds.tables`).
    """
    # Names as written: a contig called NA, None or nan is one. None, which
    # mge writes for a host it does not know, is a text; NA, in the other
    # columns, a missing value.
    text = ["BinName", "Contigs", "MetagenomicBin", "AssociatedBins"]
    missing = ["BinLength", "GC", "Hit", "BinningScore", "ContigsNumber", "MetagenomicBin", "AssociationScore", "AssociatedBins"]
    summary = pd.read_csv(
        summary_file,
        sep="\t",
        dtype=dict.fromkeys(text, str),
        keep_default_na=False,
        na_values={column: ["NA"] for column in missing},
        **EXACT_FLOATS,
    )
    clusters = pd.DataFrame({"round": MGE, "name": summary["BinName"], "parent_round": None, "parent": None})
    for column in summary.columns:
        if column not in MGE_SUMMARY_KEYS:
            clusters[column] = summary[column].astype(object) if column in text else summary[column]
    members = pd.DataFrame(
        [
            (MGE, name, contig)
            for name, contigs in zip(summary["BinName"], summary["Contigs"], strict=True)
            for contig in contigs.split(",")
        ],
        columns=["round", "cluster", "contig"],
    )
    return clusters, members
