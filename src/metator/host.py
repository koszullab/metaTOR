#!/usr/bin/env python3
# coding: utf-8

"""Find the MAGs that host the MGE clusters of a .metacool map: metator host.

The MGE clusters are those of the round ``mge`` of the map, which ``metator
mge`` stores, and the MAGs are the genomes of the map, whatever made them:
the final clusters of a metaTOR run, or the bins of another binner that
``metacooler reassign`` made the genomes of the map. ``__unbinned__`` is no
MAG, and neither is a genome named after an MGE cluster (an MGE cluster made
a genome of the map, by ``metator mge --reassign`` or ``metacooler reassign
--round mge``). A contig of an MGE cluster counts for its MGE cluster only:
a MAG that holds it too (a validation clusters the MGE contigs with the
other contigs) is taken without it, or the contacts of the MGE cluster
within itself would count as contacts with that MAG.

The signal of an MGE cluster with a MAG is the sum of the weights of the
edges of the network between their contigs: the network that ``metator
network`` stored in the map, weighed with its normalization (see
metator.network.stored_network). A MAG hosts the MGE cluster when that
signal is at least ``min_share`` percent (10 by default) of the signal of the
MGE cluster with all the MAGs, and when the MGE cluster is in contact with at
least ``min_contigs`` of its contigs (5 by default). The contacts of the MGE
cluster with the contigs in no MAG and with the other MGE clusters do not
count in that share: an MGE cluster with few contacts with the MAGs, as a
phage that the library caught in its free particles, can be given a host
from a small part of its contacts. ``share_of_all``, the share of all the
contacts of the MGE cluster with the other contigs, says how small.

``metator host`` stores what it finds in the map (see
metator.metacool.write_map): the hosts of each MGE cluster in the round
``mge`` (HOST_CLUSTER_COLUMNS), and the MGE clusters of each MAG in the
genomes table (HOST_GENOME_COLUMNS); and writes every MGE cluster and MAG in
contact in ``tables/mge_hosts.txt`` (see write_hosts), and the signal and share
of the hosts in ``plots/mge_hosts.pdf`` (see plot_hosts), in its output
directory.
"""

import os
from dataclasses import dataclass

import numpy as np
import pandas as pd
from matplotlib.figure import Figure

import metator.clusters as mtcl
import metator.figures as mtf
import metator.metacool as mtmc
from metator.log import logger

#: The share of the signal of an MGE cluster with the MAGs, in percent, from
#: which a MAG is a host (--min-share).
MIN_SHARE = 10.0

#: The number of contigs of a MAG in contact with an MGE cluster from which
#: the MAG is a host (--min-contigs).
MIN_CONTIGS = 5

#: The table of the MGE clusters and MAGs in contact, and the plot of the
#: hosts, in tables/ and plots/ of the output directory of a run (see
#: metator.outputs.TABLES_DIR and metator.figures.PLOT_DIR): a PDF file.
HOSTS_FILE = "mge_hosts.txt"
PLOT_FILE = "mge_hosts.pdf"

#: The columns of HOSTS_FILE, in order.
PAIRS_COLUMNS = ["mge", "mag", "signal", "contigs", "share", "share_of_all", "host"]

#: What separates the names (and the shares) of a list in a column of the map.
SEPARATOR = ","


@dataclass
class Hosts:
    """The hosts of the MGE clusters of a map (see find_hosts).

    Attributes
    ----------
    pairs : pandas.DataFrame
        One row per MGE cluster and MAG in contact (PAIRS_COLUMNS), the MGE
        clusters in the order of the round mge, the MAGs of each by
        decreasing signal (ties in the order of the genomes table): the
        signal of the MGE cluster with the MAG (signal), the number of
        contigs of the MAG in contact with it (contigs), the signal as a
        percentage of that of the MGE cluster with all the MAGs (share) and
        of all its contacts with the other contigs (share_of_all), and
        whether the MAG hosts it (host).
    clusters : pandas.DataFrame
        One row per MGE cluster, in the order of the round mge: its name, its
        number of hosts (n_hosts), its hosts by decreasing share (hosts) and
        their shares, with two decimals (host_shares), each list separated by
        commas, missing for none.
    genomes : pandas.DataFrame
        One row per MAG, in the order of the genomes table: its name, the
        number of MGE clusters that it hosts (n_mges), and their names, in the
        order of the round mge, separated by commas (mges, an empty text for
        none).
    """

    pairs: pd.DataFrame
    clusters: pd.DataFrame
    genomes: pd.DataFrame


def map_bins(metacool, mc):
    """The MAGs and the MGE clusters of a map, as find_hosts reads them (see the
    module's docstring).

    Parameters
    ----------
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection, for the messages.
    mc : metacooler.MetaCooler
        The map, opened.

    Returns
    -------
    numpy.ndarray
        The MAG of each contig of the map, in its order: an empty text for a
        contig in no MAG (in __unbinned__, in an MGE genome, or in an MGE
        cluster).
    pandas.DataFrame
        The MGE clusters and their contigs: one row per contig of an MGE
        cluster, its cluster (mge) and its position in the map, from 0
        (position).
    list of str
        The MGE clusters, in the order of the round mge.
    list of str
        The MAGs, in the order of the genomes table.

    Raises
    ------
    metator.metacool.MetacoolError
        If the map holds no round mge, or no MAG."""
    with mtmc.read_errors(metacool):
        rounds = [str(entry.get("name")) for entry in mc.cluster_rounds]
        if mtcl.MGE not in rounds:
            raise mtmc.MetacoolError(
                f"{metacool} holds no MGE clusters (no clustering round {mtcl.MGE}): run metator mge on it first, which "
                "stores them in the map."
            )
        clusters = mc.clusters
        mges = [str(name) for name in clusters.loc[(clusters["round"] == mtcl.MGE).to_numpy(), "name"]]
        members = mc.cluster_members(mtcl.MGE)
        contigs = mc.contigs
        genome_names = mc.genomes["name"].astype(str).to_numpy()
    names = contigs["name"].astype(str).to_numpy()
    position = pd.Series(np.arange(len(names)), index=names)
    members = pd.DataFrame(
        {
            "mge": members["cluster"].astype(str).to_numpy(),
            "position": position.reindex(members["contig"].astype(str)).to_numpy(dtype=np.int64),
        }
    )
    in_mge = np.zeros(len(names), dtype=bool)
    in_mge[members["position"].to_numpy()] = True
    # Genomes named after an MGE cluster are MGE clusters made genomes, no
    # MAGs.
    mge_genomes = set(mges)
    mags = [name for name in genome_names if name != mtmc.UNBINNED_GENOME and name not in mge_genomes]
    if not mags:
        raise mtmc.MetacoolError(
            f"{metacool} holds no MAG to host its MGE clusters: none of its genomes is a MAG (but __unbinned__ and the "
            "MGE clusters made genomes). Run metator validation (or pipeline) on it first, which makes its final "
            "clusters its genomes, or make the bins of another binner its genomes with metacooler reassign."
        )
    genome = genome_names[contigs["genome_id"].to_numpy()]
    mag = np.where(~in_mge & np.isin(genome, mags), genome, "")
    return mag, members, mges, mags


def contact_pairs(network, mag, members, mges, mags):
    """The signal of each MGE cluster with each MAG it is in contact with (see the
    module's docstring), from the network of the map.

    Parameters
    ----------
    network : pandas.DataFrame
        The network of the map, by the positions of its contigs, from 1 (see
        metator.network.map_network): contig1, contig2 and weight.
    mag, members, mges, mags :
        The MAGs and the MGE clusters of the map (see map_bins).

    Returns
    -------
    pandas.DataFrame
        One row per MGE cluster and MAG in contact: mge, mag, signal,
        contigs, share and share_of_all (see Hosts.pairs), the MGE clusters
        in the order of mges, the MAGs of each by decreasing signal, ties in
        the order of mags."""
    in_mge = np.zeros(len(mag), dtype=bool)
    in_mge[members["position"].to_numpy()] = True
    first = network["contig1"].to_numpy() - 1
    second = network["contig2"].to_numpy() - 1
    weight = network["weight"].to_numpy(dtype=np.float64)
    # The edges of the MGE contigs, from each of their two contigs that is
    # one: the contacts of the MGE clusters.
    touched = in_mge[first] | in_mge[second]
    edges = pd.DataFrame(
        {
            "position": np.concatenate([first[touched], second[touched]]),
            "other": np.concatenate([second[touched], first[touched]]),
            "weight": np.concatenate([weight[touched], weight[touched]]),
        }
    )
    edges = edges.merge(members, on="position")
    # The contacts of an MGE cluster with the other contigs: not those within
    # it.
    within = members.rename(columns={"position": "other"}).assign(within=True)
    edges = edges.merge(within, on=["mge", "other"], how="left")
    edges = edges[edges["within"].isna()]
    total = edges.groupby("mge", sort=False)["weight"].sum()
    edges = edges.assign(mag=mag[edges["other"].to_numpy()])
    edges = edges[(edges["mag"] != "").to_numpy()]
    pairs = edges.groupby(["mge", "mag"], sort=False).agg(signal=("weight", "sum"), contigs=("other", "nunique")).reset_index()
    pairs["contigs"] = pairs["contigs"].astype(np.int64)
    # As numbers, whatever the types of an empty table (no MGE cluster in
    # contact with a MAG).
    signal = pairs["signal"].to_numpy(dtype=np.float64)
    pairs["signal"] = signal
    pairs["share"] = 100 * signal / pairs.groupby("mge", sort=False)["signal"].transform("sum").to_numpy(dtype=np.float64)
    pairs["share_of_all"] = 100 * signal / pairs["mge"].map(total).to_numpy(dtype=np.float64)
    order = pd.DataFrame(
        {
            "mge_rank": pairs["mge"].map({name: rank for rank, name in enumerate(mges)}),
            "mag_rank": pairs["mag"].map({name: rank for rank, name in enumerate(mags)}),
            "signal": -pairs["signal"],
        }
    )
    ranked = order.sort_values(["mge_rank", "signal", "mag_rank"], kind="stable").index
    return pairs.loc[ranked].reset_index(drop=True)


def find_hosts(metacool, mc, network, min_share=MIN_SHARE, min_contigs=MIN_CONTIGS):
    """The hosts of the MGE clusters of a map (see the module's docstring).

    Parameters
    ----------
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection, for the messages.
    mc : metacooler.MetaCooler
        The map, opened: it gives the MAGs and the MGE clusters (see
        map_bins).
    network : pandas.DataFrame
        The network of the map, by the positions of its contigs, from 1 (see
        metator.network.stored_network).
    min_share : float
        The share of the signal of an MGE cluster with the MAGs, in percent,
        from which a MAG is a host.
    min_contigs : int
        The number of contigs of a MAG in contact with an MGE cluster from
        which the MAG is a host.

    Returns
    -------
    Hosts
        The MGE clusters and MAGs in contact, the hosts of each MGE cluster,
        and the MGE clusters of each MAG.

    Raises
    ------
    metator.metacool.MetacoolError
        If the map holds no round mge, or no MAG (see map_bins)."""
    mag, members, mges, mags = map_bins(metacool, mc)
    pairs = contact_pairs(network, mag, members, mges, mags)
    pairs["host"] = (pairs["share"] >= min_share) & (pairs["contigs"] >= min_contigs)
    hosts = pairs[pairs["host"].to_numpy()]
    by_mge = {name: group for name, group in hosts.groupby("mge", sort=False)}
    clusters = pd.DataFrame(
        [
            (
                name,
                len(by_mge[name]) if name in by_mge else 0,
                SEPARATOR.join(by_mge[name]["mag"]) if name in by_mge else None,
                SEPARATOR.join(f"{share:.2f}" for share in by_mge[name]["share"]) if name in by_mge else None,
            )
            for name in mges
        ],
        columns=["name", *mtmc.HOST_CLUSTER_COLUMNS],
    )
    rank = {name: position for position, name in enumerate(mges)}
    hosted = {name: sorted(group["mge"], key=rank.__getitem__) for name, group in hosts.groupby("mag", sort=False)}
    genomes = pd.DataFrame(
        [(name, len(hosted.get(name, [])), SEPARATOR.join(hosted.get(name, []))) for name in mags],
        columns=["name", *mtmc.HOST_GENOME_COLUMNS],
    )
    n_hosts = clusters["n_hosts"].to_numpy()
    logger.info(
        f"{len(mges)} MGE clusters and {len(mags)} MAGs: {int((n_hosts == 1).sum())} MGE clusters with one host, "
        f"{int((n_hosts > 1).sum())} with several, {int((n_hosts == 0).sum())} with none; "
        f"{int((genomes['n_mges'] > 0).sum())} MAGs host at least one."
    )
    return Hosts(pairs[PAIRS_COLUMNS], clusters, genomes)


def write_hosts(pairs, path):
    """Write the MGE clusters and MAGs in contact (pairs, see Hosts.pairs) at path,
    tab-separated, with a header: PAIRS_COLUMNS, host True or False. The file is
    replaced if it exists (a symbolic link there, not its target)."""
    if os.path.islink(path):
        os.remove(path)
    pairs[PAIRS_COLUMNS].to_csv(path, sep="\t", index=False)


def plot_hosts(pairs, path):
    """Plot the hosts of the MGE clusters (the rows of pairs that are hosts, see
    Hosts.pairs) at path: the distribution of their signal, on a log scale,
    and of their share. The axes say so when there is no host."""
    hosts = pairs[pairs["host"].to_numpy()]
    fig = Figure(figsize=(12, 5))
    signal, share = fig.subplots(1, 2)
    if len(hosts):
        values = hosts["signal"].to_numpy(dtype=np.float64)
        low, high = values.min(), values.max()
        # geomspace, whose last edge is the highest signal exactly: a host
        # there is counted.
        signal.hist(values, bins=np.geomspace(low, high if high > low else 10 * low, 21), edgecolor="black")
        signal.set_xscale("log")
        share.hist(hosts["share"].to_numpy(dtype=np.float64), bins=np.linspace(0, 100, 21), edgecolor="black")
    else:
        for ax in (signal, share):
            ax.text(0.5, 0.5, "No host", ha="center", va="center", transform=ax.transAxes)
    signal.set_xlabel("Signal of the MGE cluster with its host (log scale)")
    share.set_xlabel("Share of the signal of the MGE cluster with the MAGs (%)")
    share.set_xlim(0, 100)
    for ax in (signal, share):
        ax.set_ylabel("Hosts")
    fig.suptitle(f"The {len(hosts)} hosts of {hosts['mge'].nunique()} MGE clusters")
    mtf._saved(fig, path)
