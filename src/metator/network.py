#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Generate the metaHiC network of the contigs from a .metacool contact map.

The network sums the contacts of the map by pair of contigs (see
metacool_to_contacts), normalizes them, and writes them with the data of the
contigs measured on the assembly (see create_contig_data).

Core function to build the network are:
    - check_map_contigs
    - contig_data_from_map
    - create_contig_data
    - map_columns
    - metacool_to_contacts
    - normalization_weights
    - normalize_pair
    - read_metacool_contacts
    - write_contig_data
"""

import csv
import numpy as np
import pandas as pd
import re
from Bio import SeqIO
from Bio import SeqUtils
from os.path import join
from metator.log import logger
import metator.io as mio
import metator.metacool as mtmc

# Pixels of a .metacool file read at a time (about 20 bytes each): the
# network is summed chunk by chunk, never from the whole pixel table at once.
METACOOL_CHUNKSIZE = 1_000_000


def metacool_to_contacts(
    metacool,
    contig_data,
    out_dir,
    output_file_network,
    output_file_contig_data,
    tmp_dir,
    n_cpus,
    normalization,
    self_contacts=False,
    chunksize=METACOOL_CHUNKSIZE,
):
    """Generates the network file (in edgelist form) from a .metacool contact
    map: the network of the pairs the map was made of.

    The contacts between two contigs are the sum of the pixels between their
    bins, whatever the bin size of the map. The hits of a contig are the sum
    of the pixels of its bins, a pixel between two of its bins counting
    twice, as a pair with both reads on the contig does in the pairs files
    the map was made of. The intra-contig contacts are left out of the
    network unless self_contacts is set. The network is then normalized
    (see normalize_pair) and written with the contig data file.

    The map may lack contigs of the assembly (a map made with a contig
    filter, for example, which leaves out the pairs with a read on the
    contigs it drops): those have no contacts and 0 hits, as a contig that
    no read maps to. The network is then that of the pairs the map holds.

    It is the network of the pairs the map holds, each counted once: a pair
    with a read on the first or the last base of a contig counts in both
    (and at position 0, which a 1-based pairs file should not hold but some
    do, where metacooler puts the read on the first base), a pair with a
    read past the end of its contig in neither. The map holds the contacts
    of all its libraries, summed.

    Parameters:
    -----------
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection: ``file.mmetacool::entries/<name>``.
    contig_data : dict
        Dictionary of the all the contigs from the assembly (see
        create_contig_data), the contigs names are the keys to the data of the
        contig available with the following keys: "id", "length", "GC",
        "hit", "coverage", "RS". The contigs of the map must be contigs of
        the assembly, with the same lengths, but the map may lack some. Their
        hits are set from the map: 0 for a contig that it lacks.
    out_dir : str
        The output directory to write the network and contig data into.
    output_file_network : str
        The specific file name for the output network file.
    output_file_contig_data : str
        The specific file name for the output contig data file.
    tmp_dir : str
        Path to the temporary directory.
    n_cpus : int
        Number of cpus used to sort the edges.
    normalization : str
        Mode of normalization of the contact counts (see normalize_pair), or
        "None" for the raw counts.
    self_contacts : bool
        Whether to keep the intra-contig contacts in the network. Default is
        False.
    chunksize : int
        Number of pixels read at a time.

    Returns:
    --------
    str:
        Path to the network file.
    str:
        Path to the verbose contig data file.

    Raises:
    -------
    metator.metacool.MetacoolError
        If a map cannot be read, stores both triangles of the matrix, holds
        contigs that are not in the assembly or that have another length in
        it, holds no contact, or holds counts that are not integers (see
        read_metacool_contacts).
    """
    network_file = join(out_dir, output_file_network)
    contig_data_file = join(out_dir, output_file_contig_data)
    pre_network_file = join(tmp_dir, "precompute_network_file.txt")
    pre_network_sorted_file = join(tmp_dir, "tmp_network_sorted.txt")

    hits, contacts = read_metacool_contacts(metacool, contig_data, chunksize)
    names = np.empty(len(hits), dtype=object)
    for name, data in contig_data.items():
        data["hit"] = int(hits[data["id"]])
        names[data["id"]] = name
    _log_contacts(metacool, contacts)
    _write_map_network(
        contacts,
        names,
        contig_data,
        network_file,
        pre_network_file,
        pre_network_sorted_file,
        tmp_dir,
        n_cpus,
        normalization,
        self_contacts,
        chunksize,
    )
    write_contig_data(contig_data, contig_data_file)
    return network_file, contig_data_file


def _log_contacts(source, contacts):
    """Count the contacts of a map, within and between contigs: a map
    without any is refused."""
    all_contacts = int(contacts["count"].sum())
    inter_contacts = int(contacts["count"][contacts["id1"] != contacts["id2"]].sum())
    if not all_contacts:
        raise mtmc.MetacoolError(f"{source} holds no contact: there is no network to build.")
    logger.info(f"Information of {source}:")
    logger.info(f"{all_contacts} contacts in the library.")
    logger.info(f"{inter_contacts} contacts inter-contigs in the library.")
    logger.info(f"3D ratio : {inter_contacts / all_contacts}\n")


def _write_map_network(
    contacts,
    names,
    contig_data,
    network_file,
    pre_network_file,
    pre_network_sorted_file,
    tmp_dir,
    n_cpus,
    normalization,
    self_contacts,
    chunksize,
):
    """Write the network of the contacts of a map: the edges by contig names,
    a block at a time, sorted by the names of their contigs (with sort -V,
    as metaTOR has always listed its edges), normalized with the hits of
    contig_data."""
    with open(pre_network_file, "w") as pre_network:
        for start in range(0, len(contacts), chunksize):
            edges = contacts.iloc[start : start + chunksize]
            if not self_contacts:
                edges = edges[edges["id1"] != edges["id2"]]
            pd.DataFrame(
                {"contig1": names[edges["id1"].to_numpy()], "contig2": names[edges["id2"].to_numpy()], "count": edges["count"]}
            ).to_csv(pre_network, sep="\t", header=False, index=False)
    mio.sort_pairs(pre_network_file, pre_network_sorted_file, tmp_dir=tmp_dir, threads=n_cpus)
    with open(pre_network_sorted_file, "r") as sorted_edges, open(network_file, "w") as net:
        for contig1, contig2, count in csv.reader(sorted_edges, delimiter="\t"):
            write_edge(net, contig_data, (contig1, contig2), int(count), normalization)


def read_metacool_contacts(metacool, contig_data, chunksize=METACOOL_CHUNKSIZE):
    """Contacts between the contigs of the assembly, summed from the pixels
    of a .metacool contact map.

    The pixels are read through metacooler, at most chunksize at a time
    (MetaCooler.iter_pixels), and summed by contig pair: the memory used
    grows with the number of contig pairs in contact, not with the pixels.

    Parameters:
    -----------
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection: ``file.mmetacool::entries/<name>``.
    contig_data : dict
        Dictionary of the all the contigs from the assembly (see
        create_contig_data). The contigs of the map must be contigs of the
        assembly, with the same lengths, but the map may lack some: they have
        no contacts.
    chunksize : int
        Number of pixels read at a time.

    Returns:
    --------
    numpy.ndarray:
        Hits of each contig, indexed by contig id (index 0 is unused): the
        sum of the pixels of its bins, a pixel between two of its bins
        counting twice, and 0 for a contig that the map lacks.
    pandas.DataFrame:
        Contacts between two contigs, one row per contig pair in contact,
        sorted: columns id1, id2 (id1 <= id2, an intra-contig contact has
        id1 == id2) and count.

    Raises:
    -------
    metator.metacool.MetacoolError
        If the map cannot be read (a missing file, a file that is not a
        .metacool one, a .mmetacool without the entry), stores both
        triangles of the matrix (storage-mode square), holds contigs that are
        not in the assembly or that have another length in it, or holds
        counts that are not integers.
    """
    with mtmc._metacooler_logging(), mtmc.unchanged_while_read(metacool):
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
        # Id of each contig of the map, in the order of the map.
        ids = metacool_contig_ids(metacool, mc.contigs, contig_data)
        with mc.open() as group:
            count_dtype = group["pixels"]["count"].dtype
        if not np.issubdtype(count_dtype, np.integer):
            raise mtmc.MetacoolError(
                f"{metacool} holds contact counts of type {count_dtype}: a network is built from counts of pairs, "
                "which are integers."
            )
        size = len(contig_data) + 1

        hits = np.zeros(size, dtype=np.int64)
        summed = (np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64))
        pending, n_pending = [], 0
        # The pixel table as stored, without balancing: the upper triangle.
        for chunk in mc.iter_pixels(chunksize=chunksize, contig_ids=True):
            id1 = ids[chunk["contig1_id"].to_numpy()]
            id2 = ids[chunk["contig2_id"].to_numpy()]
            count = chunk["count"].to_numpy().astype(np.int64)
            # Exact in float64: a chunk sums far fewer than 2**53 contacts.
            hits += np.bincount(id1, count, size).astype(np.int64)
            hits += np.bincount(id2, count, size).astype(np.int64)
            pending.append(_sum_by_key(np.minimum(id1, id2) * size + np.maximum(id1, id2), count))
            n_pending += len(pending[-1][0])
            # Fold the sums of the chunks into the sum of the contig pairs
            # once they outgrow it: each fold merges sorted runs.
            if n_pending >= max(chunksize, len(summed[0])):
                summed = _sum_by_key(*(np.concatenate(arrays) for arrays in zip(summed, *pending, strict=True)))
                pending, n_pending = [], 0
        keys, counts = _sum_by_key(*(np.concatenate(arrays) for arrays in zip(summed, *pending, strict=True)))
    return hits, pd.DataFrame({"id1": keys // size, "id2": keys % size, "count": counts})


def metacool_contig_ids(metacool, contigs, contig_data):
    """Id of each contig of a .metacool contact map in the contig data.

    The map may lack contigs of the assembly (a map made with a contig
    filter, for example): they get no contacts, which an INFO message says,
    with their number and the first names.

    Parameters:
    -----------
    metacool : str
        The map, for the messages.
    contigs : pandas.DataFrame
        The contigs table of the map: columns name and length.
    contig_data : dict
        Dictionary of the all the contigs from the assembly (see
        create_contig_data).

    Returns:
    --------
    numpy.ndarray:
        Id of the contigs of the map, in their order in the map.

    Raises:
    -------
    metator.metacool.MetacoolError
        If the map holds contigs that are not in the assembly (with their
        names), or contigs that have another length in it.
    """
    names = check_map_contigs(metacool, contigs, contig_data)
    in_map = set(names)
    missing = [name for name in contig_data if name not in in_map]
    if missing:
        logger.info(
            f"{len(missing)} contig(s) of the assembly are not in {metacool}: they get no contacts, as contigs that no "
            f"read maps to: {mtmc._shown(missing)}."
        )
    return np.array([contig_data[name]["id"] for name in names], dtype=np.int64)


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


def _sum_by_key(keys, counts):
    """The distinct keys, sorted, and the sum of the counts of each.

    A stable sort merges sorted runs (timsort): summing sorted sums again
    costs a merge, not a sort.
    """
    if not len(keys):
        return keys, counts
    order = np.argsort(keys, kind="stable")
    keys, counts = keys[order], counts[order]
    starts = np.flatnonzero(np.r_[True, keys[1:] != keys[:-1]])
    return keys[starts], np.add.reduceat(counts, starts)


def write_edge(net, contig_data, pair, n_occ, normalization):
    """Write one network edge: the two contig ids and the contact weight.

    Parameters:
    -----------
    net : file object
        Open network file.
    contig_data : dict
        Dictionary of the all the contigs from the assembly, the contigs names
        are the keys to the data of the contig available with the following
        keys: "id", "length", "GC", "hit", "coverage", "RS".
    pair : list of str
        Names of the two contigs.
    n_occ : int
        Contact count between the two contigs.
    normalization : str
        Mode of normalization to use, or "None".
    """
    if normalization != "None":
        weight = normalize_pair(contig_data, pair, n_occ, normalization)
    else:
        weight = n_occ
    net.write(f"{contig_data[pair[0]]['id']}\t{contig_data[pair[1]]['id']}\t{weight}\n")


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
        Dictionary of the all the contigs from the assembly, the contigs names
        are the keys to the data of the contig available with the following
        keys: "id", "length", "GC", "hit", "coverage", "RS". Hit is set to 0:
        the network step counts it on the contacts of the map.
    """

    # Create a contig data dictionnary from the assembly file
    contig_data = dict()

    # Extract restriction sites if an enzyme is given.
    if enzyme:
        pattern = mio.get_restriction_site(enzyme)

    # The length and coverage of each contig in the depth file, by name.
    depths = read_depth(depth_file) if depth_file else None
    missing, other_length = [], []

    global_id = 1
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
                "id": global_id,
                "length": len(contig.seq),
                "GC": gc_content,
                "hit": 0,
                "coverage": coverage,
                "RS": (len(re.findall(pattern, str(contig.seq))) + 1) if enzyme else "-",
            }
            global_id += 1

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
        If the map cannot be read.
    """
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


def normalize_pair(contig_data, pair, n_occ, normalization):
    """Function to do the normalization of an inter contig contact depending on
    the mode of normalisation given.

    Parameters:
    -----------
    contig_data : dict
        Dictionnary of the all the contigs from the assembly, the contigs names
        are the keys to the data of the contig available with the following
        keys: "id", "length", "GC", "hit", "coverage", "RS".
    pair : tuple
        Tuple of the id of the contigs contributing to the contact.
    n_occ : int
        Contact count between the two contigs.
    normalization : str
        Mode of normalization to use.

    Returns:
    --------
    float:
        Normalized contact count.
    """
    # The four first normalizations are normalization by the geometric mean of
    # "coverage".
    if normalization == "abundance":
        factor = np.sqrt(contig_data[pair[0]]["coverage"] * contig_data[pair[1]]["coverage"])
        # If no read mapped from Shotgun libraries could be equal to zero.
        if factor == 0:
            return 0

    elif normalization == "length":
        factor = np.sqrt(
            contig_data[pair[0]]["hit"]
            / contig_data[pair[0]]["length"]
            * contig_data[pair[1]]["hit"]
            / contig_data[pair[1]]["length"]
        )

    elif normalization == "RS":
        factor = 0.01 * np.sqrt(
            contig_data[pair[0]]["hit"] / contig_data[pair[0]]["RS"] * contig_data[pair[1]]["hit"] / contig_data[pair[1]]["RS"]
        )

    # The last two normalizations are normalization by geometric means of hits.
    elif normalization == "empirical_hit":
        factor = np.sqrt(contig_data[pair[0]]["hit"] * contig_data[pair[1]]["hit"])

    elif normalization == "theoritical_hit":
        factor = np.sqrt(
            contig_data[pair[0]]["coverage"]
            * 10
            * np.sqrt(contig_data[pair[0]]["length"] * contig_data[pair[0]]["RS"])
            * contig_data[pair[1]]["coverage"]
            * 10
            * np.sqrt(contig_data[pair[1]]["length"] * contig_data[pair[1]]["RS"])
        )
    # If no read mapped from Shotgun libraries could be equal to zero.
    if factor == 0:
        return 0

    return n_occ / factor


def normalization_weights(lengths, hits, coverage, sites, normalization):
    """One weight per contig whose products give the normalization of the
    network: the factorised form of normalize_pair.

    Each normalization divides the contacts of two contigs by the geometric
    mean of one value per contig (times 0.01 for RS), so the normalized
    contact of two contigs is their contact count times one weight per
    contig: ``n_occ * w[i] * w[j] == normalize_pair(contig_data, (i, j),
    n_occ, normalization)``. That is how a .metacool map stores a
    normalization, the way cooler stores its balancing weights (see
    metator.metacool.store_normalization). A contig whose value gives
    normalize_pair a factor of 0, where it returns 0, has no weight (NaN):
    a contig with no hit, or no shotgun coverage.

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
        Mode of normalization (see normalize_pair).

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


def write_contig_data(contig_data, output_path):
    """Function to write the contig data file at the output path given. The file
    will contains 7 columns separated by a tabulation: id, name, length,
    GC_content, hit, coverage, restriction site for each contig with an header.

    Parameters:
    -----------
    contig_data : dict
        Dictionary of the all the contigs from the assembly, the contigs names
        are the keys to the data of the contig available with the following
        keys: "id", "length", "GC", "hit", "coverage", "RS".
    output_path : str
        Path to the output file where the data from the dictionary will be
        written
    """

    # For each contig extract the data and write them in the file.
    with open(output_path, "w") as contig_data_file_handle:
        line = "ID\tName\tSize\tGC_content\tHit\tShotgun_coverage\tRestriction_site\n"
        contig_data_file_handle.write(line)
        for name in contig_data:
            length = contig_data[name]["length"]
            hit = contig_data[name]["hit"]
            coverage = contig_data[name]["coverage"]
            GC_content = contig_data[name]["GC"]
            idx = contig_data[name]["id"]
            RS = contig_data[name]["RS"]
            line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(idx, name, length, GC_content, hit, coverage, RS)
            contig_data_file_handle.write(line)
