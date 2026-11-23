#!/usr/bin/env python3
# coding: utf-8

"""The tables of a run, written from the map of its clusters: metator dump.

A run keeps everything it finds in its map and writes its tables in its
temporary directory only. ``metator dump MAP -o DIR`` writes them from the
map, read only: ``contig_data_final.txt`` (in the map's order, without the
``ID`` column, whose values were positions in the map the run read),
``cluster_summary.txt`` (rows, columns and formats as the validation
recorded them), ``binning.txt``, ``clusters.txt`` (with the number of
contigs and the length of each cluster, which the map computes),
``cluster_members.txt`` and ``network.txt``.

The map must be the map of the clusters as its run wrote it (see
check_dumpable); every table is read through one opened map, and they are
renamed into place all at once.
"""

import os
import uuid
from os.path import join

import numpy as np
import pandas as pd

import metator.clusters as mtcl
import metator.io as mio
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

#: The keys of the validation entry of the metator metadata of a map that
#: record how the run wrote the summary of its final clusters
#: (cluster_summary.txt): its rows, its columns and their formats (see
#: metator.cli.validation_entry and metator.io.SummaryLayout).
SUMMARY_KEYS = ("summary_rows", "summary_columns", "summary_formats")

#: The formats of the columns of the summary that the validation writes (see
#: metator.io.SummaryLayout).
SUMMARY_FORMATS = ("%.2f", "%.4f", "%d", "%s")

#: The columns that a recursive step of the validation sets.
RECURSIVE_COLUMNS = ("Recursive_cluster_ID", "Recursive_cluster_contigs", "Recursive_cluster_size")


def check_dumpable(metacool, mc):
    """Refuse a map that holds no run to dump, or whose clusters were edited since
    the run wrote them.

    A map whose genomes are no longer the final clusters of its run (a
    metacooler reassign, a curation by hand) is refused too: its summary would
    contradict its clusters.

    Returns
    -------
    dict
        The validation entry of the metator metadata of the map.

    Raises
    ------
    metator.metacool.MetacoolError
        If the map lacks either, or its clusters were edited."""
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
    from metacooler.schema import CLUSTER_GENOME_COLUMN

    clusters = mc.clusters
    final = clusters.loc[(clusters["round"] == mtcl.FINAL).to_numpy(), ["name", CLUSTER_GENOME_COLUMN]]
    edited = [
        str(name)
        for name, genome in zip(final["name"], final[CLUSTER_GENOME_COLUMN], strict=True)
        if pd.isna(genome) or str(genome) != str(name)
    ]
    final_names = set(final["name"].astype(str))
    extra = [name for name in mc.genomes["name"].astype(str) if name != UNBINNED_GENOME and name not in final_names]
    if edited or extra:
        changed = []
        if edited:
            changed.append(f"final cluster(s) whose contigs are no longer a genome of the map: {mtmc._shown(edited)}")
        if extra:
            changed.append(f"genome(s) that are no final cluster: {mtmc._shown(extra)}")
        raise mtmc.MetacoolError(
            f"The clusters of {metacool} were edited after the run that wrote it (by metacooler reassign, for example): "
            f"{'; and '.join(changed)}. Its summary of the clusters and its clustering rounds no longer describe "
            "them, and its tables would contradict each other. Run metator validation on the map again, or write "
            "the edited assignment with metacooler dump (metacooler dump -t contigs MAP)."
        )
    return validation


def contig_data_final(metacool, mc, table):
    """The contig data of the final clusters of a run (contig_data_final.txt),
    from its map: one row per contig of the map, in its order.

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
        The table, without its ID column."""
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
    validation (``recursive_<k>``)."""
    prefix = mtcl.recursive_round("")
    return name.startswith(prefix) and name[len(prefix) :].isdigit()


def final_clusters(mc):
    """The final cluster of each contig of the map, in its order: the name of its
    genome, ND for metacooler's __unbinned__ (a contig in no final cluster)."""
    names = mc.genomes["name"].astype(str).to_numpy()
    final = names[mc.contigs["genome_id"].to_numpy()]
    return np.where(final == UNBINNED_GENOME, "ND", final)


def binning(data):
    """Each contig of a final cluster, and its cluster (binning.txt, the table
    that anvi'o imports), from the contig data of the final clusters (see
    contig_data_final), in its order."""
    binned = data[data["Final_cluster"] != "ND"]
    return binned[["Name", "Final_cluster"]]


def summary_layout(metacool, validation):
    """How the run wrote the summary of its final clusters (cluster_summary.txt),
    as the validation entry of its map records it (validation, see
    metator.cli.validation_entry): its rows, its columns and their formats,
    with which the values in the genomes table of the map write it again (see
    cluster_summary).

    Returns
    -------
    metator.io.SummaryLayout
        The rows, the columns and the format of each column.

    Raises
    ------
    metator.metacool.MetacoolError
        If the entry lacks one of them, as that of a map that a development
        build of metaTOR 2.0 validated does, or does not record them as metator
        validation records them: each time, run metator validation on the map
        again."""
    missing = [key for key in SUMMARY_KEYS if key not in validation]
    if missing:
        raise mtmc.MetacoolError(
            f"The validation entry of the metator metadata of {metacool} lacks {', '.join(missing)}, which record how "
            "the summary of its clusters (cluster_summary.txt) is written: a development build of metaTOR 2.0 "
            "validated the map. Run metator validation on the map again."
        )
    rows, columns, formats = (validation[key] for key in SUMMARY_KEYS)
    malformed = []
    if not _texts(rows):
        malformed.append("summary_rows, which must be the names of the clusters")
    if not _texts(columns):
        malformed.append("summary_columns, which must be the names of the columns")
    elif not isinstance(formats, dict) or any(formats.get(column) not in SUMMARY_FORMATS for column in columns):
        malformed.append(f"summary_formats, which must give each column one of the formats {', '.join(SUMMARY_FORMATS)}")
    if malformed:
        raise mtmc.MetacoolError(
            f"The validation entry of the metator metadata of {metacool} does not record the summary of its clusters "
            f"as metator validation records it: {'; '.join(malformed)}. Run metator validation on the map again."
        )
    return mio.SummaryLayout(rows, columns, {column: formats[column] for column in columns})


def _texts(values):
    """Whether values is a list of texts."""
    return isinstance(values, list) and all(isinstance(value, str) for value in values)


def cluster_summary(metacool, mc, layout):
    """The summary of the final clusters of a run (cluster_summary.txt), from the
    genomes table of its map, as the validation wrote it (layout, see
    summary_layout): a line per row of the summary, in its order, the name of
    the cluster then the value of each column in the genomes table, written
    with the format of the column (an empty field for a missing value: NaN,
    None or an empty text), after a header, an empty field then the columns.

    Returns
    -------
    pandas.DataFrame
        The fields of the summary, as text, indexed by the names of the
        clusters (an index without a name, which gives the header its empty
        field).

    Raises
    ------
    metator.metacool.MetacoolError
        If the genomes table lacks a column or a cluster of the summary, or
        holds a value that the format of its column does not write (an edit of
        the map changed them since the run)."""
    genomes = mc.genomes.set_index(mc.genomes["name"].astype(str))
    missing = [column for column in layout.columns if column not in genomes.columns]
    if missing:
        raise mtmc.MetacoolError(
            f"The genomes table of {metacool} lacks the column(s) {mtmc._shown(missing)} of the summary of its "
            "clusters, which metator validation stored there (an edit of the map dropped them since): run metator "
            "validation on the map again."
        )
    absent = [row for row in layout.rows if row not in genomes.index]
    if absent:
        raise mtmc.MetacoolError(
            f"The genomes table of {metacool} lacks the cluster(s) {mtmc._shown(absent)} of the summary of its "
            "clusters (summary_rows): run metator validation on the map again."
        )
    values = genomes.loc[layout.rows, layout.columns]
    fields = {}
    for column in layout.columns:
        written = layout.formats[column]
        try:
            fields[column] = ["" if _missing(value) else written % value for value in values[column]]
        except (TypeError, ValueError) as error:
            raise mtmc.MetacoolError(
                f"The column {column} of the genomes table of {metacool} holds a value that metator validation did "
                f"not write there ({written}: {error}): an edit of the map changed it since. Run metator validation "
                "on the map again."
            ) from error
    return pd.DataFrame(fields, index=pd.Index(layout.rows, dtype=object), columns=layout.columns, dtype=object)


def _missing(value):
    """Whether a value of the genomes table of a map is missing: NaN, None, or an
    empty text, which metacooler stores for a missing text."""
    if isinstance(value, str):
        return value == ""
    return bool(pd.isna(value))


def check_network(metacool, mc, stored):
    """Refuse a map whose network is not the one that its clusters were made from:
    metator partition records the digest of the network that it partitioned
    (see metator.metacool.partition_network), which that of the network stored
    in the map (stored, see metator.network.stored_network) must be. metator
    network, run on the map after the run that wrote its clusters (with another
    normalization, or on the contacts of another library), stores another
    network there and keeps the clusters: the network and the clusters of the
    tables would not go together.

    Raises
    ------
    metator.metacool.MetacoolError
        With the remedies: metator network with the normalization of the run,
        which stores its network again on the contacts that it read, or metator
        partition and validation, which cluster the network stored now."""
    made_on = mtmc.partition_network(mc)
    if made_on == stored.digest:
        return
    if made_on is None:
        raise mtmc.MetacoolError(
            f"{metacool} does not record the network that its clusters were made from (the "
            f"{mtmc.NETWORK_DIGEST} of the partition entry of its metator metadata, which metator partition "
            "records): run metator partition, then metator validation, on the map again (or metator pipeline)."
        )
    partition = mtmc._metator_metadata(mc).get("partition")
    normalization = partition.get("normalization") if isinstance(partition, dict) else None
    run = f" (normalization {normalization})" if isinstance(normalization, str) else ""
    option = f" -n {normalization}" if isinstance(normalization, str) else ""
    raise mtmc.MetacoolError(
        f"The network stored in {metacool} is not the one that its clusters were made from: metator network stored "
        f"another network in the map (normalization {stored.normalization}) after the run that wrote them "
        f"partitioned its own{run}, and kept them. Its tables would give the network of the one and the clusters of "
        f"the other. Run metator network{option} on the map, which stores the network of the run again if its "
        "contacts are those that the run read, or metator partition, then metator validation, on the map again "
        "(or metator pipeline), which cluster the network stored now."
    )


def network_table(table, network):
    """The network of the contigs of a map (network.txt): one row per edge, the
    names of its two contigs and its weight, without the contacts within a
    contig, in the order of the map: from the first contig of the map to the
    last, each with the contigs after it (see metator.network.map_network). The
    weight is the count of the contacts of the two contigs times the weight of
    each, for the normalization that metator network stored, or the count for
    None.

    Parameters
    ----------
    table : pandas.DataFrame
        The contig data of the contigs of the map, in its order.
    network : pandas.DataFrame
        The network of the map, by the positions of its contigs (see
        metator.network.map_network)."""
    names = table["Name"].astype(str).to_numpy()
    return pd.DataFrame(
        {
            "contig1": names[network["contig1"].to_numpy() - 1],
            "contig2": names[network["contig2"].to_numpy() - 1],
            "weight": network["weight"].to_numpy(),
        }
    )


def existing_tables(outdir):
    """The tables of metator dump that outdir holds already (see TABLES), by name."""
    return [name for name in TABLES if os.path.lexists(join(outdir, name))]


#: What to do about a map replaced while metator dump reads it.
DUMP_AGAIN = "No table was written: run metator dump again once the run that writes the map has ended."


def dump(metacool, outdir, force=False):
    """Write the tables of the run that wrote the map metacool (see the module's
    docstring) in outdir, made if it is not there. The map is read only: every
    table comes from one state of the map, read through one opened map, which
    is refused if it was replaced or written while it was read.

    The tables are written under temporary names, then renamed: an earlier
    dump's tables are replaced all at once, or not at all.

    Parameters
    ----------
    metacool : str
        Path to the .metacool file of the clusters of a run, or URI of an entry
        of a .mmetacool collection: ``file.mmetacool::entries/<name>``.
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
        If outdir holds one of the tables already (without force), or a
        directory in the place of one, the map cannot be read, holds no
        clusters of a validation or clusters edited since (see check_dumpable),
        no summary of them that metaTOR 2.0 records (see summary_layout), no
        network, partition or summary that metaTOR stored, or a network other
        than the one its clusters were made from (see check_network), or was
        replaced while it was read, or a table cannot be written."""
    blocked = [name for name in TABLES if os.path.isdir(join(outdir, name)) and not os.path.islink(join(outdir, name))]
    if blocked:
        raise mtmc.MetacoolError(
            f"{outdir} holds a directory where metator dump writes {', '.join(blocked)}: remove it, or give another "
            "output directory."
        )
    existing = existing_tables(outdir)
    if existing and not force:
        raise mtmc.MetacoolError(
            f"{outdir} holds the table(s) of an earlier dump: {', '.join(existing)}. Give -F to replace them, or "
            "another output directory."
        )
    try:
        with mtmc._metacooler_logging(), mtmc.read_errors(metacool, DUMP_AGAIN):
            mc = mtmc.open_map(metacool)
            layout = summary_layout(metacool, check_dumpable(metacool, mc))
        # The network, read once, through mc: the contig data of the contigs of
        # the map, and the network weighed with the stored normalization
        # (whose Hit column must be that of its contacts), which must be the
        # one that the clusters were made from.
        stored = mtn.stored_network(metacool, mc)
        with mtmc._metacooler_logging(), mtmc.read_errors(metacool, DUMP_AGAIN):
            check_network(metacool, mc, stored)
            data = contig_data_final(metacool, mc, stored.table)
            summary = cluster_summary(metacool, mc, layout)
            clusters, members = mtmc.map_rounds(mc)
            # Each cluster with its number of contigs and their length, which
            # the map computes, as the run wrote them in its clusters.txt: the
            # MGE clusters too, which mges_bin_summary.tsv gave as
            # ContigsNumber and BinLength.
            contigs = mc.contigs
            lengths = pd.Series(contigs["length"].to_numpy(), index=contigs["name"].astype(str).to_numpy())
            clusters = mtcl.with_sizes(clusters, members, lengths)
    except mtmc.MapReplacedError as error:
        # A replacement found by a read of stored_network, which names
        # another remedy.
        raise mtmc.MapReplacedError(
            f"{metacool} was replaced or written while metator dump read it, by a run that updated it: what was "
            f"read of it would mix two maps. {DUMP_AGAIN}"
        ) from error
    try:
        os.makedirs(outdir, exist_ok=True)
    except OSError as error:
        raise mtmc.MetacoolError(f"The output directory {outdir} cannot be made: {error}.") from error
    writers = {
        CONTIG_DATA_FILE: lambda path: data.to_csv(path, sep="\t", index=False),
        SUMMARY_FILE: lambda path: summary.to_csv(path, sep="\t"),
        BINNING_FILE: lambda path: binning(data).to_csv(path, sep="\t", header=False, index=False),
        NETWORK_FILE: lambda path: network_table(stored.table, stored.network).to_csv(path, sep="\t", header=False, index=False),
    }
    _write_all(metacool, outdir, writers, (clusters, members))
    logger.info(
        f"The tables of {metacool} are written in {outdir}: {len(data)} contigs, {len(summary)} final clusters, "
        f"{clusters['round'].nunique()} clustering rounds and {len(stored.network)} edges of the network."
    )
    return [join(outdir, name) for name in TABLES]


def _write_all(metacool, outdir, writers, rounds):
    """Write the tables of a dump in outdir (see dump): each under a temporary
    name (``.<table>.<random>.tmp``), then each renamed onto its own, once all are
    written, which replaces a table of an earlier dump, or a symbolic link
    there, not its target."""
    token = uuid.uuid4().hex[:12]
    temporary = {name: join(outdir, f".{name}.{token}.tmp") for name in TABLES}
    clusters_file, members_file = temporary[mtcl.CLUSTERS_FILE], temporary[mtcl.CLUSTER_MEMBERS_FILE]
    replaced = []
    try:
        for name, write in writers.items():
            write(temporary[name])
        mtcl.write_round_tables(*rounds, clusters_file, members_file)
        for name in TABLES:
            os.replace(temporary[name], join(outdir, name))
            replaced.append(name)
    except OSError as error:
        for path in temporary.values():
            if os.path.lexists(path):
                os.remove(path)
        written = f"; already replaced: {', '.join(replaced)}" if replaced else "; no table was replaced"
        raise mtmc.MetacoolError(f"The tables of {metacool} cannot be written in {outdir}: {error}{written}.") from error
