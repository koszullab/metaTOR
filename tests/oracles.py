"""Oracles of the tests: what metaTOR computed before it computed it
otherwise."""

import numpy as np


def normalize_pair(contig_data, pair, n_occ, normalization):
    """The normalization of a contact of two contigs as metaTOR 1.x defined
    it, edge by edge (its metator.network.normalize_pair): the count divided
    by the geometric mean of one value of each contig, 0 when that mean is
    0. metaTOR 2.0 stores one weight per contig instead (see
    metator.network.normalization_weights), whose products give the same
    normalized contacts.

    Parameters
    ----------
    contig_data : dict
        The data of the contigs by name, as metator.network.create_contig_data
        gives them, with their hits ("hit").
    pair : tuple of str
        The names of the two contigs.
    n_occ : int
        Their contact count.
    normalization : str
        abundance, length, RS, empirical_hit or theoritical_hit.
    """
    first, second = contig_data[pair[0]], contig_data[pair[1]]

    def mean(value):
        return np.sqrt(value(first) * value(second))

    factor = {
        "abundance": lambda: mean(lambda c: c["coverage"]),
        "length": lambda: mean(lambda c: c["hit"] / c["length"]),
        "RS": lambda: 0.01 * mean(lambda c: c["hit"] / c["RS"]),
        "empirical_hit": lambda: mean(lambda c: c["hit"]),
        "theoritical_hit": lambda: mean(lambda c: c["coverage"] * 10 * np.sqrt(c["length"] * c["RS"])),
    }[normalization]()
    return 0 if factor == 0 else n_occ / factor


def partition_table(contig_data, partition):
    """The contig data table of a partition as metaTOR wrote it before it
    vectorised it (metator.partition.partition_table, cluster by cluster
    with pandas' .loc): the reference of its output, byte for byte."""
    import pandas as pd

    contigs_data = contig_data.reset_index(drop=True).copy()
    for kind in ("Core", "Overlapping"):
        for column in ("ID", "contigs", "size"):
            contigs_data[f"{kind}_cluster_{column}"] = pd.Series("-", index=contigs_data.index, dtype=object)
    for kind, clusters in (("Core", partition.core), ("Overlapping", partition.overlapping)):
        for number, contigs in clusters.items():
            rows = [contig - 1 for contig in contigs]
            length = sum(contigs_data.iloc[rows].Size)
            contigs_data.loc[rows, f"{kind}_cluster_ID"] = f"{number:05d}"
            contigs_data.loc[rows, f"{kind}_cluster_contigs"] = len(rows)
            contigs_data.loc[rows, f"{kind}_cluster_size"] = length
    return contigs_data
