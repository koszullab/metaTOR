#!/usr/bin/env python3
# coding: utf-8

"""The tables of a run, written from the map of its clusters: metator dump.

A run of metaTOR (``metator pipeline``, or ``network``, ``partition`` and
``validation`` one by one) keeps everything it finds in its map (see
metator.metacool): the network of the contigs, every clustering round and the
summary of the final clusters. It writes the tables of 1.x in its temporary
directory only, which ``-N`` keeps. ``metator dump MAP -o DIR`` writes them
from the map, read only:

- ``contig_data_final.txt``: one row per contig of the map, in its order:
  ``Name``, ``Size``, ``GC_content``, ``Hit``, ``Shotgun_coverage``,
  ``Restriction_site``, then the ID, the number of contigs and the size of
  its core, overlapping and recursive clusters, then its ``Final_cluster``,
  with the values that the validation writes (``-`` for no cluster,
  ``00000`` for no recursive cluster, ``ND`` for no final cluster), but no
  ``ID`` column: the ID of a contig was its position in the map that the
  run read, which a run that writes the map of its clusters changes (see
  :func:`contig_data_final`).
- ``cluster_summary.txt``: one row per final cluster, with the columns of
  the summary that the validation recorded, in its order (see
  :func:`cluster_summary`).
- ``binning.txt``: each contig of a final cluster, and its cluster.
- ``clusters.txt`` and ``cluster_members.txt``: every clustering round of the
  map, as the validation writes them (see metator.clusters).
- ``network.txt``: the network of the contigs, as metator network stored it
  (see :func:`network_table`).
"""

import os
from os.path import join

import numpy as np
import pandas as pd

import metator.clusters as mtcl
import metator.metacool as mtmc
import metator.network as mtn
import metator.partition as mtp
from metator.log import logger

#: The tables that metator dump writes, in this order.
CONTIG_DATA_FILE = "contig_data_final.txt"
SUMMARY_FILE = "cluster_summary.txt"
BINNING_FILE = "binning.txt"
NETWORK_FILE = "network.txt"
TABLES = (CONTIG_DATA_FILE, SUMMARY_FILE, BINNING_FILE, mtcl.CLUSTERS_FILE, mtcl.CLUSTER_MEMBERS_FILE, NETWORK_FILE)

#: metacooler's genome of the contigs in no genome: those in no final cluster.
UNBINNED_GENOME = "__unbinned__"

#: The columns of the summary of the final clusters that the validation
#: writes with four decimals (see metator.io.write_bin_summary); the other
#: decimals have two.
ABUNDANCE_COLUMNS = ("HiC_abundance", "SG_abundance")

#: The column of the summary of the final clusters that the validation writes
#: as a whole number: the recursive step that made the cluster. The genomes
#: table of the map may store it as decimals, next to the missing value of
#: __unbinned__.
STEP_COLUMN = "step"

#: The columns of the clustering rounds that metacooler computes from their
#: members (see metacooler.MetaCooler.clusters), which no run writes.
COMPUTED_CLUSTER_COLUMNS = ("n_contigs", "length", "genome")

#: The columns that a recursive step of the validation sets.
RECURSIVE_COLUMNS = ("Recursive_cluster_ID", "Recursive_cluster_contigs", "Recursive_cluster_size")


def check_dumpable(metacool, mc):
    """Refuse a map that holds no run to dump: one without the validation
    entry of its metator metadata, or without the clustering round final,
    which metator validation and pipeline write. The error names what is
    missing, and the step that writes it.

    Raises
    ------
    metator.metacool.MetacoolError
        If the map lacks either.
    """
    missing = []
    metator = mc.metadata.get("metator")
    validation = metator.get("validation") if isinstance(metator, dict) else None
    if not isinstance(validation, dict):
        missing.append("the validation entry of its metator metadata")
    if mtcl.FINAL not in [entry.get("name") for entry in mc.cluster_rounds]:
        missing.append(f"the clustering round {mtcl.FINAL}")
    if missing:
        raise mtmc.MetacoolError(
            f"{metacool} holds no clusters of a validation to dump: it lacks {' and '.join(missing)}. Run metator "
            "validation (or metator pipeline) on it first, which writes the map of the clusters."
        )
    return validation


def contig_data_final(metacool, mc, table):
    """The contig data of the final clusters of a run (contig_data_final.txt),
    from its map: one row per contig of the map, in its order.

    The columns are those that the validation writes, but its ID: the
    contig data that metator network stored (table, see
    metator.network.map_contig_table), the core and the overlapping cluster
    of each contig, from the partition that the map stores (see
    metator.metacool.read_partition and metator.partition.partition_table),
    its recursive cluster, from the rounds recursive_<k> of the map, and its
    final cluster, its genome in the map (ND for none).

    The values are written as the validation writes them: ``-`` for a
    contig in no cluster, ``00000`` for one in no recursive cluster, the
    IDs of the overlapping and recursive clusters with five digits. The ID
    of a core cluster is a number when every contig is in one, as the
    validation reads its table back (``00001`` reads as 1), and has five
    digits when some contig is in none (``-``).

    Parameters
    ----------
    metacool : str
        The map, for the messages.
    mc : metacooler.MetaCooler
        The map, opened.
    table : pandas.DataFrame
        The contig data of the contigs of the map (see
        metator.network.map_contig_table).

    Returns
    -------
    pandas.DataFrame
        The table, without its ID column.
    """
    data = mtp.partition_table(table, mtmc.read_partition(metacool, mc))
    # As the validation reads contig_data_partition.txt (see
    # metator.validation.recursive_decontamination): the IDs of the core
    # clusters are numbers there, but next to the - of a contig in none.
    if (data["Core_cluster_ID"] != "-").all():
        data["Core_cluster_ID"] = data["Core_cluster_ID"].astype(np.int64)
    data["Recursive_cluster_ID"] = f"{0:05d}"
    for column in RECURSIVE_COLUMNS[1:]:
        data[column] = pd.Series("-", index=data.index, dtype=object)
    rows = pd.Index(data["Name"])
    sizes = data["Size"].to_numpy()
    # The rounds recursive_<k> in the order they were made: a contig ends in
    # the cluster of the last step that split its cluster.
    for name in [entry["name"] for entry in mc.cluster_rounds]:
        if not _is_recursive_round(name):
            continue
        members = mc.cluster_members(name)
        for cluster, contigs in members.groupby("cluster", sort=False)["contig"]:
            positions = rows.get_indexer(contigs.astype(str))
            positions = positions[positions >= 0]
            data.loc[positions, "Recursive_cluster_ID"] = str(cluster).rsplit("_", 1)[-1]
            data.loc[positions, "Recursive_cluster_contigs"] = len(positions)
            data.loc[positions, "Recursive_cluster_size"] = int(sizes[positions].sum())
    data["Final_cluster"] = final_clusters(mc)
    return data.drop(columns=["ID"])


def _is_recursive_round(name):
    """Whether a clustering round is a step of the recursive partition of the
    validation (recursive_<k>)."""
    prefix = mtcl.recursive_round("")
    return name.startswith(prefix) and name[len(prefix) :].isdigit()


def final_clusters(mc):
    """The final cluster of each contig of the map, in its order: the name of
    its genome, ND for metacooler's __unbinned__ (a contig in no final
    cluster)."""
    names = mc.genomes["name"].astype(str).to_numpy()
    final = names[mc.contigs["genome_id"].to_numpy()]
    return np.where(final == UNBINNED_GENOME, "ND", final)


def binning(data):
    """Each contig of a final cluster, and its cluster (binning.txt, the table
    that anvi'o imports), from the contig data of the final clusters (see
    contig_data_final), in its order."""
    binned = data[data["Final_cluster"] != "ND"]
    return binned[["Name", "Final_cluster"]]


def cluster_summary(metacool, mc, validation):
    """The summary of the final clusters of a run (cluster_summary.txt), from
    the genomes table of its map: one row per final cluster, named by the
    first column, in the order of the clustering round final, with the
    columns that the validation recorded (the summary_columns of its entry,
    see metator.cli.validation_entry), in their order.

    The values are written as the validation writes them (see
    metator.io.write_bin_summary): the abundances with four decimals, the
    step a whole number, the other decimal numbers with two, and nothing for
    a missing value. The map stores them as the summary gave them.

    Raises
    ------
    metator.metacool.MetacoolError
        If the validation entry records no columns, or the genomes table
        lacks one of them (an edit of the map dropped it).
    """
    columns = validation.get("summary_columns")
    if not isinstance(columns, list) or not all(isinstance(column, str) for column in columns):
        raise mtmc.MetacoolError(
            f"The validation entry of the metator metadata of {metacool} records no columns of the summary of its "
            "clusters (summary_columns): run metator validation on the map again."
        )
    genomes = mc.genomes.set_index(mc.genomes["name"].astype(str))
    missing = [column for column in columns if column not in genomes.columns]
    if missing:
        raise mtmc.MetacoolError(
            f"The genomes table of {metacool} lacks the column(s) {mtmc._shown(missing)} of the summary of its "
            "clusters, which metator validation stored there (an edit of the map dropped them since): run metator "
            "validation on the map again."
        )
    final = mc.clusters
    ordered = [str(name) for name in final.loc[(final["round"] == mtcl.FINAL).to_numpy(), "name"]]
    clusters = [name for name in ordered if name in genomes.index]
    clusters += [name for name in genomes.index if name != UNBINNED_GENOME and name not in set(clusters)]
    summary = genomes.loc[clusters, columns].copy()
    if STEP_COLUMN in summary.columns and pd.api.types.is_numeric_dtype(summary[STEP_COLUMN]):
        summary[STEP_COLUMN] = summary[STEP_COLUMN].round().astype("Int64")
    for column in ABUNDANCE_COLUMNS:
        if column in summary.columns and pd.api.types.is_numeric_dtype(summary[column]):
            summary[column] = summary[column].map(lambda value: "" if pd.isna(value) else "%.4f" % value)
    summary.index.name = None
    return summary


def round_tables(mc):
    """The two tables of the clustering rounds of the map (clusters.txt and
    cluster_members.txt), as metator.clusters.Rounds writes them: every
    cluster of every round, in the order of the rounds and of their clusters,
    with its parent and the columns its round stored (not those that
    metacooler computes), and every contig of a cluster, in the order of the
    map. A missing value is NA, as the validation writes it: the empty text
    that the map stores where a text column has none, too.

    Returns
    -------
    tuple of pandas.DataFrame
        The clusters and their contigs.
    """
    clusters = mc.clusters.drop(columns=list(COMPUTED_CLUSTER_COLUMNS), errors="ignore").reset_index(drop=True)
    for column in clusters.columns:
        if not pd.api.types.is_numeric_dtype(clusters[column]):
            values = clusters[column].astype(object)
            clusters[column] = values.where(values.notna() & (values.astype(str) != ""), None)
    members = mc.cluster_members().reset_index(drop=True)
    return clusters, members[["round", "cluster", "contig"]]


def network_table(table, network):
    """The network of the contigs of a map (network.txt): one row per edge,
    the names of its two contigs and its weight, without the contacts within
    a contig, in the order of the map: from the first contig of the map to
    the last, each with the contigs after it (see
    metator.network.map_network). The weight is the count of the contacts
    of the two contigs times the weight of each, for the normalization that
    metator network stored, or the count for None.

    Parameters
    ----------
    table : pandas.DataFrame
        The contig data of the contigs of the map, in its order.
    network : pandas.DataFrame
        The network of the map, by the positions of its contigs (see
        metator.network.map_network).
    """
    names = table["Name"].astype(str).to_numpy()
    return pd.DataFrame(
        {
            "contig1": names[network["contig1"].to_numpy() - 1],
            "contig2": names[network["contig2"].to_numpy() - 1],
            "weight": network["weight"].to_numpy(),
        }
    )


def existing_tables(outdir):
    """The tables of metator dump that outdir holds already (see TABLES), by
    name."""
    return [name for name in TABLES if os.path.lexists(join(outdir, name))]


def dump(metacool, outdir, force=False):
    """Write the tables of the run that wrote the map metacool (see the module's
    docstring) in outdir, made if it is not there. The map is read only:
    every table comes from one state of the map, which is refused if it was
    replaced or written while it was read.

    Parameters
    ----------
    metacool : str
        Path to the .metacool file of the clusters of a run, or URI of an
        entry of a .mmetacool collection: ``file.mmetacool::entries/<name>``.
    outdir : str
        The directory of the tables.
    force : bool
        Whether to replace the tables that outdir holds already.

    Returns
    -------
    list of str
        The paths of the tables, in the order of TABLES.

    Raises
    ------
    metator.metacool.MetacoolError
        If outdir holds one of the tables already (without force), or the map
        cannot be read, holds no clusters of a validation (see
        check_dumpable), no network, partition or summary that metaTOR
        stored, or changed while it was read.
    """
    existing = existing_tables(outdir)
    if existing and not force:
        raise mtmc.MetacoolError(
            f"{outdir} holds the table(s) of an earlier dump: {', '.join(existing)}. Give -F to replace them, or "
            "another output directory."
        )
    with mtmc._metacooler_logging(), mtmc.read_errors(metacool):
        mc = mtmc.open_map(metacool)
        validation = check_dumpable(metacool, mc)
    # The network, read once: the contig data of the contigs of the map, and
    # the network weighed with the stored normalization (whose Hit column must
    # be that of its contacts).
    stored = mtn.stored_network(metacool)
    with mtmc._metacooler_logging(), mtmc.read_errors(metacool):
        data = contig_data_final(metacool, mc, stored.table)
        summary = cluster_summary(metacool, mc, validation)
        clusters, members = round_tables(mc)
        # Every table is of the map as it was opened: a read through mc checks
        # its file.
        with mc.open():
            pass
    os.makedirs(outdir, exist_ok=True)
    paths = [join(outdir, name) for name in TABLES]
    contig_file, summary_file, binning_file, clusters_file, members_file, network_file = paths
    for path in paths:
        if os.path.islink(path):
            os.remove(path)
    data.to_csv(contig_file, sep="\t", index=False)
    summary.to_csv(summary_file, sep="\t", float_format="%.2f")
    binning(data).to_csv(binning_file, sep="\t", header=False, index=False)
    mtcl.write_round_tables(clusters, members, clusters_file, members_file)
    network_table(stored.table, stored.network).to_csv(network_file, sep="\t", header=False, index=False)
    logger.info(
        f"The tables of {metacool} are written in {outdir}: {len(data)} contigs, {len(summary)} final clusters, "
        f"{clusters['round'].nunique()} clustering rounds and {len(stored.network)} edges of the network."
    )
    return paths
