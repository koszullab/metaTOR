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
