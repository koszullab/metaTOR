#!/usr/bin/env python3
# coding: utf-8

"""The clustering rounds of a run: what the map stores in its /clusters
group (see metacooler's format), and the two tables a run writes them to.

A run of metaTOR reaches its final clusters in rounds:

- ``overlapping``: the overlapping clusters of the partition that reach
  ``--size``, which miComplete checks;
- ``recursive_<k>``: the clusters that step k of the validation makes by
  partitioning a contaminated cluster again, its parent;
- ``final``: the final clusters, each the cluster of the same name of the
  round it comes from.

The validation records them as it makes them (:class:`Rounds`) and writes
them next to its other tables: ``clusters.txt``, one row per cluster of
every round, with its parent and the columns of miComplete, and
``cluster_members.txt``, one row per contig of a cluster. The map of the
clusters stores both (:func:`metator.metacool.reassign_metacool`).
"""

import csv
from os.path import exists, join

import pandas as pd

from metator.log import logger


#: The two tables of the rounds, in the output directory of a run.
CLUSTERS_FILE = "clusters.txt"
CLUSTER_MEMBERS_FILE = "cluster_members.txt"

#: The names of the rounds.
OVERLAPPING = "overlapping"
FINAL = "final"


def recursive_round(step):
    """The name of the round of recursive step *step* (1, 2...)."""
    return f"recursive_{step}"


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
        clusters = pd.DataFrame(self.clusters, columns=_columns(self.clusters))
        clusters.to_csv(clusters_file, sep="\t", index=False, na_rep="NA", quoting=csv.QUOTE_NONE)
        members = pd.DataFrame(self.members, columns=["round", "cluster", "contig"])
        members.to_csv(members_file, sep="\t", index=False, quoting=csv.QUOTE_NONE)
        logger.info(
            f"{len(clusters)} clusters in {clusters['round'].nunique()} rounds are written to {clusters_file} and {members_file}."
        )
        return clusters_file, members_file


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
    # the missing values of the clusters table are NA (Rounds.write).
    clusters = pd.read_csv(
        paths[0],
        sep="\t",
        dtype={"name": str, "parent": str, "round": str, "parent_round": str},
        keep_default_na=False,
        na_values=["NA"],
    )
    members = pd.read_csv(paths[1], sep="\t", dtype=str, keep_default_na=False)
    return clusters, members
