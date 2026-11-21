#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""The metaHiC network of the contigs of a .metacool contact map.

The network sums the contacts of the map by pair of contigs (see
contig_contacts: metacooler's network of the contigs), and weighs each edge
with the normalization of the network (see map_network): the count of the
pair times one weight per contig, which metator network stores in the map as
weights of its bins (see metator.metacool.store_network). Its nodes are the
contigs of the map, by their position in it.

What metaTOR measures on each contig (see create_contig_data) is stored in
the contigs table of the map by metator network too, and read from it by
the steps after it (see map_contig_table): the map is the network that they
partition.

Core functions to build the network are:
    - check_map_contigs
    - contig_contacts
    - contig_data_from_map
    - contig_data_table
    - contig_hits
    - create_contig_data
    - map_columns
    - map_contig_table
    - map_network
    - network_digest
    - network_normalization
    - network_weights
    - normalization_weights
    - relabel_network
    - stored_network
    - stored_weights
"""

import logging
import re
from typing import NamedTuple, cast

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import SeqUtils

import metator.io as mio
import metator.metacool as mtmc
from metator.log import logger

# Pixels of a .metacool file read at a time (about 20 bytes each): the
# contacts are summed chunk by chunk, never from the whole pixel table at once.
METACOOL_CHUNKSIZE = 1_000_000

# The columns of the contig data of a run (see contig_table): the contigs of
# its map, in its order, with what metaTOR measured on them. partition writes
# them in contig_data_partition.txt, with those of its clusters. ID is the
# position of the contig in the map, from 1: the id of its node in the
# network (see map_network).
CONTIG_TABLE_COLUMNS = ["ID", "Name", "Size", "GC_content", "Hit", "Shotgun_coverage", "Restriction_site"]


def contig_contacts(metacool, chunksize=METACOOL_CHUNKSIZE):
    """Contacts of a .metacool contact map summed by pair of contigs.

    metacooler sums the pixels of the map between the bins of each two
    contigs (metacooler.network at the contig level, by contig id: the row
    of each contig in the contigs table of the map), whatever the bin size
    of the map, at most chunksize pixels read at a time: the memory used
    grows with the number of contig pairs in contact, not with the pixels.
    The contacts of two contigs are those of the pairs that the map was made
    of, each counted once: a pair with a read on the first or the last base
    of a contig too. The map holds the contacts of all its libraries,
    summed.

    Parameters:
    -----------
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection: ``file.mmetacool::entries/<name>``.
    chunksize : int
        Number of pixels read at a time.

    Returns:
    --------
    pandas.DataFrame:
        One row per pair of contigs in contact, in the order of the map:
        columns id1 and id2, the positions of the two contigs in the contigs
        table of the map, from 0 (id1 <= id2: a contact within a contig has
        id1 == id2), and count.

    Raises:
    -------
    metator.metacool.MetacoolError
        If the map cannot be read (a missing file, a file that is not a
        .metacool one, a .mmetacool without the entry), stores both
        triangles of the matrix (storage-mode square), or holds counts that
        are not integers.
    """
    import metacooler

    # metacooler's warnings go to the log of the run; its count of the edges
    # does not.
    with mtmc._metacooler_logging(level=logging.WARNING), mtmc.read_errors(metacool):
        mc = mtmc.open_map(metacool)
        # Each pixel of the table is read: a map that also stores the lower
        # triangle would count the contacts between two bins twice.
        storage_mode = mc.info.get("storage-mode", "symmetric-upper")
        if storage_mode != "symmetric-upper":
            raise mtmc.MetacoolError(
                f"{metacool} stores its contacts with storage-mode {storage_mode!r}: a network is built from a map "
                "that stores each contact once, in the upper triangle (storage-mode 'symmetric-upper'), as "
                "metacooler writes them."
            )
        with mc.open() as group:
            count_dtype = group["pixels"]["count"].dtype
        if not np.issubdtype(count_dtype, np.integer):
            raise mtmc.MetacoolError(
                f"{metacool} holds contact counts of type {count_dtype}: a network is built from counts of pairs, "
                "which are integers."
            )
        # Without out, metacooler returns the edges rather than write them,
        # and with ids, each contig as its row in /contigs: its position in
        # the map, from 0, which no name is read for.
        edges = cast(pd.DataFrame, metacooler.network(metacool, level="contig", out=None, chunksize=chunksize, ids=True))
    return pd.DataFrame(
        {
            "id1": edges["node1"].to_numpy(dtype=np.int64),
            "id2": edges["node2"].to_numpy(dtype=np.int64),
            "count": edges["value"].to_numpy(dtype=np.int64),
        }
    )


def contig_hits(contacts, n_contigs):
    """Hits of each contig of a map, by its position in the map: the sum of
    its contacts (see contig_contacts), a contact within the contig counted
    twice, as a pair with both reads on the contig counts twice in the pairs
    file that the map was made of. A contig in no contact has 0."""
    count = contacts["count"].to_numpy(dtype=np.float64)
    hits = np.bincount(contacts["id1"], weights=count, minlength=n_contigs)
    hits += np.bincount(contacts["id2"], weights=count, minlength=n_contigs)
    # Exact in float64: far fewer than 2**53 contacts.
    return hits.astype(np.int64)


def log_contacts(source, contacts):
    """Count the contacts of a map (see contig_contacts), within and between
    contigs: a map without any is refused (MetacoolError)."""
    all_contacts = int(contacts["count"].sum())
    inter_contacts = int(contacts["count"][contacts["id1"] != contacts["id2"]].sum())
    if not all_contacts:
        raise mtmc.MetacoolError(f"{source} holds no contact: there is no network to build.")
    logger.info(f"Information of {source}:")
    logger.info(f"{all_contacts} contacts in the library.")
    logger.info(f"{inter_contacts} contacts inter-contigs in the library.")
    logger.info(f"3D ratio : {inter_contacts / all_contacts}\n")


def map_network(metacool, weights="stored", contacts=None):
    """The network of the contigs of a .metacool contact map: one edge
    between each two contigs of the map in contact, weighed with the
    normalization of the network.

    The weight of an edge is the count of the contacts of its two contigs
    (see contig_contacts) times the weight of each: ``count * w[i] * w[j]``, as
    cooler applies balancing weights. The weights are those that metator
    network stored in the map (see stored_weights), by default, or those
    given (one per contig of the map, in its order: those of a pipeline run,
    see network_weights); without weights (the normalization None), the
    edges are the counts. Weighing the contacts of two contigs, rather than
    each of their pixels, gives the same network whatever the bin size of
    the map.

    The contacts within a contig are no edge, nor is a pair of contigs whose
    weight is 0, or that has none: a contig that the normalization gives no
    weight (NaN: one without a hit, or without shotgun coverage).

    Parameters:
    -----------
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection: ``file.mmetacool::entries/<name>``.
    weights : "stored", numpy.ndarray or None
        "stored" (the default): the weights that metator network stored in
        the map. Else one weight per contig of the map, in its order, or
        None for the counts. Any other string is refused (ValueError): the
        name of a normalization is no weights.
    contacts : pandas.DataFrame, optional
        The contacts of the map (see contig_contacts), when they are read
        already: the map is not read again.

    Returns:
    --------
    pandas.DataFrame:
        One row per edge, in the order of the map: columns contig1 and
        contig2, the IDs of the two contigs (their positions in the map as
        it is read, from 1: the ID column of the contig data that a run
        writes from it, see contig_table; contig1 < contig2), and weight. A
        run that regroups the map (validation, pipeline) moves its contigs:
        join the network of the map it wrote to its tables by contig name
        (see relabel_network).

    Raises:
    -------
    ValueError
        If weights is a string other than "stored".
    metator.metacool.MetacoolError
        If the map cannot be read or give its contacts (see
        contig_contacts), or, for its stored weights, holds no network, or
        not the weights of its normalization (see stored_weights).
    """
    if isinstance(weights, str):
        if weights != "stored":
            raise ValueError(
                f"map_network weighs the contacts with weights='stored' (the weights stored in the map), one weight "
                f"per contig of the map, or None (the counts): not {weights!r}."
            )
        weights = stored_weights(metacool)
    if contacts is None:
        contacts = contig_contacts(metacool)
    between = contacts[contacts["id1"] != contacts["id2"]]
    id1, id2 = between["id1"].to_numpy(), between["id2"].to_numpy()
    weight = between["count"].to_numpy()
    if weights is not None:
        contig_weights = np.asarray(weights, dtype=np.float64)
        weight = weight * contig_weights[id1] * contig_weights[id2]
    keep = np.isfinite(weight) & (weight != 0)
    return pd.DataFrame({"contig1": id1[keep] + 1, "contig2": id2[keep] + 1, "weight": weight[keep]})


def relabel_network(network, names, ids):
    """The network of a map (see map_network) with other IDs for its
    contigs, those of a contig data table of the map: the IDs of the network
    are the positions of the contigs in the map now, and those of the table
    the positions they had when it was written (by partition), which differ
    once a run has written the map of its clusters (grouped by cluster). The
    edges stay those of the network, listed as map_network lists them: each
    from the contig of the lower ID, in the order of their IDs.

    Parameters:
    -----------
    network : pandas.DataFrame
        The network (columns contig1, contig2 and weight).
    names : sequence of str
        The names of the contigs of the map, in its order.
    ids : dict
        The ID of each contig of the map in the table, by name.

    Returns:
    --------
    pandas.DataFrame:
        The network, with the IDs of the table.
    """
    new_ids = np.array([ids[str(name)] for name in names], dtype=np.int64)
    first, second = new_ids[network["contig1"].to_numpy() - 1], new_ids[network["contig2"].to_numpy() - 1]
    relabeled = pd.DataFrame(
        {"contig1": np.minimum(first, second), "contig2": np.maximum(first, second), "weight": network["weight"].to_numpy()}
    )
    return relabeled.sort_values(["contig1", "contig2"], kind="stable").reset_index(drop=True)


def network_normalization(metacool, mc=None):
    """The normalization of the network that metator network stored in a map
    (see metator.metacool.store_network): that of the entry network of the
    metator key of its metadata, which a pipeline run writes too. "None" for
    the raw counts.

    Parameters:
    -----------
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection, for the message.
    mc : metacooler.MetaCooler, optional
        The map, opened: by default, metacool is opened.

    Raises:
    -------
    metator.metacool.MetacoolError
        If the map holds no network: metator network has not run on it.
    """
    if mc is None:
        mc = mtmc.open_map(metacool)
    metator = mc.metadata.get("metator")
    entry = metator.get("network") if isinstance(metator, dict) else None
    normalization = entry.get("normalization") if isinstance(entry, dict) else None
    if not isinstance(normalization, str):
        raise mtmc.MetacoolError(
            f"{metacool} holds no network of its contigs: run metator network on it first, which stores the network in the map."
        )
    return normalization


class StoredNetwork(NamedTuple):
    """The network that metator network stored in a map, as the steps after
    it read it (see stored_network): the contig data of the contigs of the
    map (table, see map_contig_table), the network (see map_network), its
    normalization, and its identity (digest, see network_digest)."""

    table: pd.DataFrame
    network: pd.DataFrame
    normalization: str
    digest: str


def stored_network(metacool, mc=None):
    """The network that metator network stored in a map, read in full, as
    partition and validation read it: the contig data of its contigs (see
    map_contig_table), the network of its contacts weighed with the stored
    weights (see map_network and stored_weights), its normalization (see
    network_normalization), and its identity (see network_digest), which
    partition records with the partition that it stores, and validation
    checks.

    Given mc, the map opened, every table is read through it, and the
    contacts, which metacooler reads by path, are checked against it: a map
    replaced since mc was opened is refused as replaced (see
    metator.metacool.read_errors), before its Hit column is checked against
    contacts that would be another map's.

    The Hit column of the map must still be the hits of its contacts (see
    contig_hits): an edit that changed the contacts since metator network
    stored them (metacooler merge, which sums the maps of several libraries,
    keeps the contigs table of the first) leaves it stale, and the weights
    of a normalization by the hits with it, but for the normalization None,
    which stores no weights for metacooler to drop.

    Parameters:
    -----------
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection: ``file.mmetacool::entries/<name>``.
    mc : metacooler.MetaCooler, optional
        The map, opened: by default, metacool is opened.

    Returns:
    --------
    StoredNetwork:
        The contig data, the network, its normalization and its digest.

    Raises:
    -------
    metator.metacool.MetacoolError
        If the map cannot be read, holds no network, lacks a column or the
        weights that metator network stores (an edit of the map dropped
        them since), holds a value that is not a number in one of these
        columns, or weights that differ between the bins of a contig, or
        hits that are not those of its contacts: each time, with the step
        to run again.
    """
    with mtmc._metacooler_logging(), mtmc.read_errors(metacool):
        if mc is None:
            mc = mtmc.open_map(metacool)
        normalization = network_normalization(metacool, mc)
    table = map_contig_table(metacool, mc)
    contacts = contig_contacts(metacool)
    # The contacts of the map that mc opened: a read through mc checks its
    # file.
    with mtmc._metacooler_logging(), mtmc.read_errors(metacool):
        with mc.open():
            pass
    hits = contig_hits(contacts, len(table))
    stale = table["Name"][table["Hit"].to_numpy() != hits].tolist()
    if stale:
        raise mtmc.MetacoolError(
            f"The Hit column of {metacool} is not the hits of its contacts for {len(stale)} contig(s): "
            f"{mtmc._shown(stale)}. Its contacts changed since metator network stored the network (a merge of maps "
            "keeps the contigs table of the first): run metator network on the map again."
        )
    weights = stored_weights(metacool, mc)
    return StoredNetwork(
        table,
        map_network(metacool, weights, contacts),
        normalization,
        network_digest(contacts, weights, table["Name"], normalization),
    )


def _check_network_columns(metacool, columns):
    """Refuse a map whose contigs table (its columns) lacks the GC_content or
    Hit column, which metator network stores."""
    missing = [column for column in ("GC_content", "Hit") if column not in columns]
    if missing:
        raise mtmc.MetacoolError(
            f"{metacool} has no {' or '.join(missing)} column in its contigs table, which metator network stores (an "
            "edit of the map dropped it since): run metator network on the map again."
        )


def stored_weights(metacool, mc=None):
    """The weights of the contigs of a map for the normalization of its
    network, which metator network stored in it as the weights
    metator_<normalization> of its bins (see network_normalization and
    metator.metacool.store_network): one per contig, in the order of the
    map, the weight of each of its bins. None for the normalization None,
    whose network is the counts. Given mc, the map opened, it is read
    through it.

    Raises:
    -------
    metator.metacool.MetacoolError
        If the map cannot be read, holds no network, lacks the weights of
        its normalization (an edit of the map dropped them since metator
        network stored them, as metacooler drops weights whose columns it
        changes), or holds weights that differ between the bins of a contig.
    """
    import metacooler

    with mtmc._metacooler_logging(), mtmc.read_errors(metacool):
        if mc is None:
            mc = mtmc.open_map(metacool)
        normalization = network_normalization(metacool, mc)
        if normalization == "None":
            return None
        name = mtmc.normalization_weight_name(normalization)
        with mc.open() as group:
            stored = name in group["bins"]
            offsets = group["indexes/chrom_offset"][()]
        if not stored:
            raise mtmc.MetacoolError(
                f"{metacool} lacks the weights {name} of the {normalization} normalization of its network, which an "
                "edit of the map dropped since metator network stored them: run metator network on the map again."
            )
        weights = metacooler.read_weights(mc, name)
    contig_weights = weights[offsets[:-1]]
    if not np.array_equal(np.repeat(contig_weights, np.diff(offsets)), weights, equal_nan=True):
        raise mtmc.MetacoolError(
            f"The weights {name} of {metacool} differ between the bins of a contig: they are not the weights of the "
            "contigs that metator network stores. Run metator network on the map again."
        )
    return contig_weights


def contig_table(names, lengths, gc, hits, coverage=None, sites=None):
    """The contig data of the contigs of a map, in its order: one row per
    contig, with the columns of CONTIG_TABLE_COLUMNS. ID is the position of
    the contig in the map, from 1. The shotgun coverage and the restriction
    sites are ``-`` when the run has none (no --depth, no --enzyme, nor a
    column of the map that gives them).

    Parameters:
    -----------
    names, lengths, gc, hits : sequences
        The name, length, GC content and Hi-C hits of each contig.
    coverage, sites : sequences, optional
        The shotgun coverage and the restriction sites (plus one, as
        create_contig_data counts them) of each contig.

    Returns:
    --------
    pandas.DataFrame:
        The contig data.
    """
    n_contigs = len(names)
    # The text columns of the type that pandas reads text as (that of
    # contig_data_partition.txt, which the validation reads).
    return pd.DataFrame(
        {
            "ID": np.arange(1, n_contigs + 1, dtype=np.int64),
            "Name": [str(name) for name in names],
            "Size": np.asarray(lengths, dtype=np.int64),
            "GC_content": np.asarray(gc, dtype=np.float64),
            "Hit": np.asarray(hits, dtype=np.int64),
            "Shotgun_coverage": "-" if coverage is None else np.asarray(coverage, dtype=np.float64),
            "Restriction_site": "-" if sites is None else np.asarray(sites, dtype=np.int64),
        },
        columns=CONTIG_TABLE_COLUMNS,
    )


def contig_data_table(contig_data, names, hits, metacool=None):
    """The contig data of the contigs of a map, in its order (see
    contig_table), from what a run measured on the assembly and read from
    the map (contig_data, by contig name: see create_contig_data and
    contig_data_from_map) and the hits of the contigs of the map (see
    contig_hits). The contigs of the assembly that the map lacks (a map
    made with a contig filter, for example) are left out: they have no
    contact, and no node in the network. An INFO message says how many, and
    names the first, when metacool, the map, is given.

    Parameters:
    -----------
    contig_data : dict
        The data of each contig of the assembly, by name.
    names : sequence of str
        The names of the contigs of the map, in its order.
    hits : numpy.ndarray
        The hits of the contigs of the map, in its order.
    metacool : str, optional
        The map, for the message.

    Returns:
    --------
    pandas.DataFrame:
        The contig data.
    """
    names = [str(name) for name in names]
    if metacool is not None:
        in_map = set(names)
        missing = [name for name in contig_data if name not in in_map]
        if missing:
            logger.info(
                f"{len(missing)} contig(s) of the assembly are not in {metacool}: they get no contacts, as contigs that "
                f"no read maps to, and are in no table of the run: {mtmc._shown(missing)}."
            )
    data = [contig_data[name] for name in names]

    def measured(key):
        values = [contig[key] for contig in data]
        return None if all(value == "-" for value in values) else values

    return contig_table(
        names,
        [contig["length"] for contig in data],
        [contig["GC"] for contig in data],
        hits,
        measured("coverage"),
        measured("RS"),
    )


def map_contig_table(metacool, mc=None):
    """The contig data of the contigs of a map, in its order (see
    contig_table), as metator network stored them in its contigs table (see
    metator.metacool.store_network): GC_content and Hit, and
    Shotgun_coverage and Restriction_site when the run measured them. Given
    mc, the map opened, it is read through it.

    Raises:
    -------
    metator.metacool.MetacoolError
        If the map cannot be read, holds no network (see
        network_normalization), or lacks the GC_content or Hit column (an
        edit of the map dropped it since metator network stored it).
    """
    with mtmc._metacooler_logging(), mtmc.read_errors(metacool):
        if mc is None:
            mc = mtmc.open_map(metacool)
        network_normalization(metacool, mc)
        contigs = mc.contigs
    _check_network_columns(metacool, contigs.columns)
    try:
        return contig_table(
            contigs["name"].to_numpy(),
            contigs["length"].to_numpy(),
            contigs["GC_content"].to_numpy(),
            contigs["Hit"].to_numpy(),
            contigs["Shotgun_coverage"].to_numpy() if "Shotgun_coverage" in contigs.columns else None,
            contigs["Restriction_site"].to_numpy() if "Restriction_site" in contigs.columns else None,
        )
    except (TypeError, ValueError) as error:
        raise mtmc.MetacoolError(
            f"The contigs table of {metacool} holds a value that is not a number in a column that metator network "
            f"stores ({error}): run metator network on the map again."
        ) from error


def network_weights(table, normalization):
    """The weights of the contigs of a contig data table (see contig_table)
    for a normalization of the network: those of normalization_weights,
    from its sizes, hits, shotgun coverage and restriction sites, which
    metator network computes from the columns it stores in the map. None for
    the normalization None, whose network is the counts."""
    if normalization == "None":
        return None

    def values(column):
        column = table[column]
        return column.to_numpy(dtype=np.float64) if pd.api.types.is_numeric_dtype(column) else None

    return normalization_weights(
        values("Size"), values("Hit"), values("Shotgun_coverage"), values("Restriction_site"), normalization
    )


def check_map_contigs(metacool, contigs, contig_data):
    """Refuse a .metacool map that is not a map of the contigs of the
    assembly: one that holds contigs that the assembly lacks, or that gives
    them other lengths. It may lack contigs of the assembly. Returns the
    names of its contigs, in its order.

    Raises:
    -------
    metator.metacool.MetacoolError
        With the contigs that are not in the assembly, or that have another
        length in it.
    """
    names = [str(name) for name in contigs["name"]]
    extra = [name for name in names if name not in contig_data]
    if extra:
        raise mtmc.MetacoolError(
            f"{metacool} is not a contact map of the contigs of the assembly: {len(extra)} of its contig(s) are not "
            f"in the assembly: {mtmc._shown(extra)}. Build it from pairs of reads aligned on this assembly."
        )
    lengths = {name: data["length"] for name, data in contig_data.items()}
    mtmc.check_map_lengths(metacool, contigs, lengths, "the assembly", "the assembly")
    return names


def read_depth(depth_file):
    """Read the depth file of jgi_summarize_bam_contig_depths (MetaBAT2): the
    length and the mean coverage of each contig, by name.

    Parameters:
    -----------
    depth_file : str
        Path to the depth file: a header line, then one line per contig, with
        its name, length and mean coverage in its first three columns.

    Returns:
    --------
    dict:
        The length (int) and the coverage (float) of each contig, by name.
    """
    depths = {}
    with open(depth_file, "r") as depth:
        next(depth, None)
        for line in depth:
            fields = line.rstrip("\n").split("\t")
            depths[fields[0]] = (int(fields[1]), float(fields[2]))
    return depths


def create_contig_data(assembly, depth_file=None, enzyme=None):
    """Create a dictionnary with data on each Contig.

    Parameters:
    -----------
    assembly : str
        Path to the assembly fasta file, plain or compressed.
    depth_file : str or None
        Path to the depth.txt file from jgi_summarize_bam_contig_depths from
        Metabat2 Software, read by contig name: it must hold every contig of
        the assembly, with its length. [Default: None]
    enzyme : str or None
        String that contains the names of the enzyme separated by a comma.
        [Default: None]

    Returns:
    --------
    dict:
        The data of each contig of the assembly, by name, in the order of
        the assembly: "length", "GC", "coverage" (``-`` without a depth
        file) and "RS" (the restriction sites plus one; ``-`` without an
        enzyme). The hits of the contigs are counted on the contacts of the
        map (see contig_hits), and their IDs are their positions in the map
        (see contig_table).
    """

    # Create a contig data dictionnary from the assembly file
    contig_data = dict()

    # Extract restriction sites if an enzyme is given.
    if enzyme:
        pattern = mio.get_restriction_site(enzyme)

    # The length and coverage of each contig in the depth file, by name.
    depths = read_depth(depth_file) if depth_file else None
    missing, other_length = [], []

    with mio.read_compressed(assembly) as fasta:
        for contig in SeqIO.parse(fasta, "fasta"):
            gc_content = int(SeqUtils.gc_fraction(contig.seq) * 100 * 10**13) / 10**13
            coverage = "-"
            if depths is not None:
                if contig.id not in depths:
                    missing.append(contig.id)
                elif depths[contig.id][0] != len(contig.seq):
                    other_length.append(contig.id)
                else:
                    coverage = depths[contig.id][1]
            contig_data[contig.id] = {
                "length": len(contig.seq),
                "GC": gc_content,
                "coverage": coverage,
                "RS": (len(re.findall(pattern, str(contig.seq))) + 1) if enzyme else "-",
            }

    # A depth file of another assembly would give the contigs wrong coverages.
    for contigs, problem in [(missing, "are not in"), (other_length, "have another length in")]:
        if contigs:
            logger.error(
                f"{len(contigs)} contigs of the assembly {problem} the depth file {depth_file}, such as {contigs[0]}: "
                "give the depth file of this assembly."
            )
            raise ValueError(f"The depth file {depth_file} is not that of the assembly {assembly}.")

    return contig_data


# The columns of the contigs table of a .metacool map that a run from the map
# reads where it measures nothing: the restriction sites and the shotgun
# coverage of each contig, which metaTOR stores in the maps it writes (see
# metator.metacool.ASSEMBLY_COLUMNS). By column: its key in the contig data
# (see create_contig_data), the option that measures it on the assembly
# instead, and what it holds.
MAP_COLUMNS = {
    "Restriction_site": ("RS", "--enzyme", "restriction sites"),
    "Shotgun_coverage": ("coverage", "--depth", "shotgun coverage"),
}


def map_columns(metacool):
    """The columns of MAP_COLUMNS that the contigs table of a .metacool map
    gives values in, in that order (see _map_values: a column of ``-``, as
    metaTOR writes what it did not measure, gives none).

    Parameters:
    -----------
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection: ``file.mmetacool::entries/<name>``.

    Returns:
    --------
    list of str:
        The columns of MAP_COLUMNS that give values in the contigs table of
        the map.

    Raises:
    -------
    metator.metacool.MetacoolError
        If the map cannot be read, or is replaced or written as it is read.
    """
    with mtmc.read_errors(metacool):
        contigs = mtmc.open_map(metacool).contigs
    return [column for column in MAP_COLUMNS if column in contigs.columns and _map_values(contigs, column)]


def _map_values(contigs, column):
    """The values of a column of MAP_COLUMNS in the contigs table of a map, by
    contig name, as create_contig_data measures them: an int of at least 1
    for the restriction sites (the sites of a contig plus one), a float of
    at least 0 for the coverage. A missing value, or one that is not such a
    number, is left out."""
    values = pd.to_numeric(contigs[column], errors="coerce").to_numpy(dtype=float)
    names = contigs["name"].astype(str).to_numpy()
    if column == "Restriction_site":
        valid = np.isfinite(values) & (values >= 1) & (values == np.floor(values))
        return {name: int(value) for name, value in zip(names[valid], values[valid], strict=True)}
    valid = np.isfinite(values) & (values >= 0)
    return {name: float(value) for name, value in zip(names[valid], values[valid], strict=True)}


def contig_data_from_map(metacool, contig_data, enzyme=None, depth=None):
    """Read the restriction sites and the shotgun coverage of the contigs from
    a .metacool map that a run starts from, where it measures none: the
    Restriction_site and Shotgun_coverage columns of its contigs table
    (MAP_COLUMNS), which metaTOR stores in the maps it writes. --enzyme and
    --depth measure them on the assembly instead (create_contig_data): they
    fill in a column that the map lacks, and replace one that it has, with
    a warning that counts the contigs whose value is not the map's.

    A map of other contigs than those of the assembly is refused first,
    for what it is (see check_map_contigs). A column is read only if it
    gives every contig of the assembly a value: an integer of at least 1
    for the restriction sites (metaTOR counts the sites of a contig plus
    one), a number of at least 0 for the coverage; but a contig of no
    length, which no map holds (metacooler leaves it out), gets the value
    that the option would give it (no site, so a Restriction_site of 1,
    and a coverage of 0). A column that
    gives no valid value (``-``, as metaTOR writes what it did not
    measure) is as no column. One that gives some but not all is not read
    either, and a warning says which contigs it lacks. The run then does
    without it, as without the option: a normalization that needs it is
    refused (see metator.cli.check_normalization), and the map of the bins
    carries the column no more.

    Parameters:
    -----------
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection: ``file.mmetacool::entries/<name>``.
    contig_data : dict
        Dictionary of the all the contigs from the assembly (see
        create_contig_data): the "RS" and "coverage" of each contig are set
        from the columns of the map that are read.
    enzyme : str or None
        The enzymes given to --enzyme, if any.
    depth : str or None
        The depth file given to --depth, if any.

    Returns:
    --------
    list of str:
        The columns of MAP_COLUMNS read from the map, in that order.

    Raises:
    -------
    metator.metacool.MetacoolError
        If the map cannot be read, or is not a map of the contigs of the
        assembly.
    """
    with mtmc.read_errors(metacool):
        contigs = mtmc.open_map(metacool).contigs
    check_map_contigs(metacool, contigs, contig_data)
    given = {"Restriction_site": enzyme, "Shotgun_coverage": depth}
    read = []
    for column, (key, option, what) in MAP_COLUMNS.items():
        if column not in contigs.columns:
            continue
        values = _map_values(contigs, column)
        if not values:
            continue
        if given[column]:
            other = [name for name, value in values.items() if name in contig_data and contig_data[name][key] != value]
            if other:
                logger.warning(
                    f"{option} replaces the {column} column of {metacool}: {len(other)} contig(s) have other {what} "
                    f"in it: {mtmc._shown(other)}."
                )
            continue
        lacking = [name for name, data in contig_data.items() if name not in values and data["length"] > 0]
        if lacking:
            valid = "an integer of at least 1" if key == "RS" else "a number of at least 0"
            logger.warning(
                f"The {column} column of {metacool} gives no {what} to {len(lacking)} contig(s) of the assembly "
                f"(not in the map, or a value that is missing or not {valid}): {mtmc._shown(lacking)}. The run does "
                f"not read it, as without {option}."
            )
            continue
        # A contig of no length: no site (plus one), no coverage.
        empty = 1 if key == "RS" else 0.0
        for name, data in contig_data.items():
            data[key] = values.get(name, empty)
        logger.info(f"The {what} of the contigs: the {column} column of {metacool}.")
        read.append(column)
    return read


def normalization_weights(lengths, hits, coverage, sites, normalization):
    """One weight per contig whose products give the normalization of the
    network.

    Each normalization divides the contacts of two contigs by the geometric
    mean of one value per contig: their shotgun coverage (abundance), their
    hits per bp (length), their hits per restriction site, times 0.01 (RS),
    their hits (empirical_hit), or ten times their coverage times the
    geometric mean of their length and their restriction sites
    (theoritical_hit). So the normalized contact of two contigs is their
    contact count times one weight per contig, ``n * w[i] * w[j]``: that is
    how a .metacool map stores the normalization of the network, the way
    cooler stores its balancing weights (see
    metator.metacool.store_normalization). A contig whose value is 0 has no
    weight (NaN): a contig with no hit, or no shotgun coverage, whose
    contacts are no edge of the network (see map_network).

    Parameters:
    -----------
    lengths : numpy.ndarray
        Length of each contig.
    hits : numpy.ndarray or None
        Hi-C hits of each contig (contig data "hit"), for length, RS and
        empirical_hit.
    coverage : numpy.ndarray or None
        Shotgun coverage of each contig (contig data "coverage"), for
        abundance and theoritical_hit.
    sites : numpy.ndarray or None
        Restriction sites of each contig, plus one, as create_contig_data
        counts them (contig data "RS"), for RS and theoritical_hit.
    normalization : str
        Mode of normalization: abundance, length, RS, empirical_hit or
        theoritical_hit.

    Returns:
    --------
    numpy.ndarray:
        One weight per contig, NaN for a contig that has none.

    Raises:
    -------
    ValueError
        If the normalization is unknown ("None" among them: it has no
        weight), or a value it uses is not given.
    """

    def need(values, what):
        if values is None:
            raise ValueError(f'The "{normalization}" normalization uses the {what} of the contigs.')
        return np.asarray(values, dtype=float)

    with np.errstate(divide="ignore", invalid="ignore"):
        if normalization == "abundance":
            weights = 1 / np.sqrt(need(coverage, "shotgun coverage"))
        elif normalization == "length":
            weights = np.sqrt(need(lengths, "lengths") / need(hits, "Hi-C hits"))
        elif normalization == "RS":
            weights = 10 * np.sqrt(need(sites, "restriction sites") / need(hits, "Hi-C hits"))
        elif normalization == "empirical_hit":
            weights = 1 / np.sqrt(need(hits, "Hi-C hits"))
        elif normalization == "theoritical_hit":
            weights = 1 / np.sqrt(
                10 * need(coverage, "shotgun coverage") * np.sqrt(need(lengths, "lengths") * need(sites, "restriction sites"))
            )
        else:
            raise ValueError(f"No weights for the normalization {normalization!r}.")
    weights = np.array(weights, dtype=float)
    weights[~np.isfinite(weights)] = np.nan
    return weights


def network_digest(contacts, weights, names, normalization):
    """The identity of a network: the SHA-256 of what map_network computes it
    from, by the names of the contigs. It is that of the network, not of a
    run: two runs that store the same network give the same digest, and a run
    that regroups the map (validation, pipeline), which moves the contigs but
    keeps their names, their contacts and their weights, leaves it as it was.
    metator partition records the digest of the network that it partitioned,
    which validation checks against the network stored in the map (see
    metator.metacool.read_partition), and metator network against the
    network that it stores.

    The digest covers the normalization, the names of the contigs of the map,
    the weight of each (none for the normalization None; NaN for a contig
    without one), and the contact count of every two contigs in contact (see
    contig_contacts), each by the names of its two contigs, in the order of
    the names. The weights of the edges are not hashed as map_network
    computes them: it multiplies the count by the weight of the contig that
    comes first in the map, then by the other, and a regroup that swaps the
    two contigs can change the last bit of the product.

    Parameters:
    -----------
    contacts : pandas.DataFrame
        The contacts of the map (see contig_contacts).
    weights : numpy.ndarray or None
        The weight of each contig of the map, in its order (see
        stored_weights and network_weights), or None for the normalization
        None.
    names : sequence of str
        The names of the contigs of the map, in its order.
    normalization : str
        The normalization of the network ("None" for the counts).

    Returns:
    --------
    str:
        The SHA-256, in hexadecimal.
    """
    import hashlib

    names = np.asarray(names, dtype=str)
    order = np.argsort(names, kind="stable")
    rank = np.empty(len(names), dtype=np.int64)
    rank[order] = np.arange(len(names), dtype=np.int64)
    between = contacts[contacts["id1"] != contacts["id2"]]
    first = rank[between["id1"].to_numpy(dtype=np.int64)]
    second = rank[between["id2"].to_numpy(dtype=np.int64)]
    # One key per pair of contigs, by the ranks of their names: one row each
    # in the contacts (see contig_contacts), so that any sort orders them so.
    key = np.minimum(first, second) * len(names) + np.maximum(first, second)
    edges = np.argsort(key)
    # Each field ended by a NUL, which no name holds.
    digest = hashlib.sha256(f"{normalization}\0".encode())
    digest.update("".join(f"{name}\0" for name in names[order].tolist()).encode())
    if weights is not None:
        ordered = np.asarray(weights, dtype=np.float64)[order]
        # One NaN, whatever the bits that the map stored for it.
        digest.update(np.where(np.isnan(ordered), np.nan, ordered).tobytes())
    for values in (key[edges], between["count"].to_numpy(dtype=np.int64)[edges]):
        digest.update(np.ascontiguousarray(values).tobytes())
    return digest.hexdigest()
