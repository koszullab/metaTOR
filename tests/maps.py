"""The map of a binning written from its pairs, the oracle of the maps that
metaTOR writes from a map of the contigs."""

import metacooler

import metator.metacool as mtmc
from metator.version import __version__


def pairs_map(pairs, contig_data, outfile, bin_summary=None, assembly_name=None, metadata=None):
    """The map of a binning written from its pairs, as metator metacool of
    metaTOR 1.4.0 wrote it (metator.metacool.write_metacool): metacooler create pairs
    with the final clusters, then the contig columns stored as metaTOR
    stores them. The oracle of the map that metaTOR writes from a map of the
    contigs (metator.metacool.reassign_metacool): the same bins, contacts
    and tables."""
    bins, contig_bin, contig_columns, header = mtmc.read_contig_data(str(contig_data))
    metacooler.from_pairs(
        str(outfile),
        bins=bins,
        pairs=[str(pairs)],
        contig_genome_map={name: contig_bin.get(name, "ND") for name in bins["chrom"]},
        assembly_name=assembly_name,
        genomes_table=str(bin_summary) if bin_summary and contig_bin else None,
        genomes_table_format="1",
        metadata={"metator": {"version": __version__, **(metadata or {})}},
        ps=False,
        keep_partial=False,
    )
    mtmc._write_back_columns(str(outfile), str(contig_data), contig_columns, header)
    mtmc._rewrite_in_a_fresh_file(str(outfile))
    return outfile
