"""Maps that metaTOR reads and writes, made from pairs files in the tests: the
map of the contigs that a run starts from, and the map of a binning written
from its pairs, the oracle of the maps that metaTOR writes from a map of the
contigs."""

import os
import tempfile

import metacooler
import pandas as pd

import metator.metacool as mtmc
from metator.version import __version__

# What metaTOR measures on the assembly, by column of the contigs table of a
# map and by key of the contig data of metator.network.create_contig_data: a
# run reads the restriction sites and the coverage from the map, where
# --enzyme and --depth measure none.
ASSEMBLY_COLUMNS = {"GC_content": "GC", "Shotgun_coverage": "coverage", "Restriction_site": "RS"}


def contig_map(pairs, contig_data, outfile, assembly_name=None, columns=True, metadata=None):
    """The map of the contigs that a run starts from, made from its pairs
    files as metacooler create pairs makes it: one bin per contig (hicstuff
    pipeline -M metacool makes one bin per restriction fragment, which gives
    the same network), every contig in __unbinned__, the longest first, with
    the contact statistics of the reads (the QC block). With
    columns, its contigs table also holds what metaTOR measures on the
    assembly (ASSEMBLY_COLUMNS: GC_content, and Shotgun_coverage and
    Restriction_site where measured), added with metacooler annotate, as
    metaTOR stored them in the maps it built before it read only maps.

    Parameters
    ----------
    pairs : str or list of str
        Pairs file(s) of the libraries, summed.
    contig_data : dict
        The contigs of the assembly, as metator.network.create_contig_data
        returns them.
    outfile : str or pathlib.Path
        Path to the .metacool file to write.
    metadata : dict, optional
        The metadata of the map.
    """
    pairs = [str(pairs)] if isinstance(pairs, (str, os.PathLike)) else [str(path) for path in pairs]
    names = list(contig_data)
    metacooler.from_pairs(
        str(outfile),
        bins=pd.DataFrame({"chrom": names, "start": 0, "end": [contig_data[name]["length"] for name in names]}),
        pairs=pairs,
        contig_genome_map=None,
        assembly_name=assembly_name,
        metadata=metadata,
        qc=True,
        ps=False,
        keep_partial=False,
    )
    if columns:
        store_assembly_columns(outfile, contig_data)
    return outfile


def pairs_map(pairs, contig_data, outfile, bin_summary=None, assembly_name=None, metadata=None):
    """The map of a binning written from its pairs, as metator metacool of
    metaTOR 1.4.0 wrote it (metator.metacool.write_metacool): metacooler create pairs
    with the final clusters (Final_cluster of the contig data file
    contig_data, ND for none), then the contig columns that metaTOR measured
    (those of CONTIG_COLUMNS that the table fills) stored as metaTOR stores
    them. The oracle of the map that metaTOR writes from a map of the
    contigs (metator.metacool.write_map): the same bins, contacts and
    tables."""
    table = pd.read_csv(str(contig_data), sep="\t", dtype=str, keep_default_na=False)
    bins = pd.DataFrame({"chrom": table["Name"], "start": 0, "end": table["Size"].astype("int64")})
    final = table["Final_cluster"].str.strip()
    binned = ~final.isin(mtmc.UNBINNED)
    contig_bin = dict(zip(table["Name"][binned], final[binned], strict=True))
    contig_columns = [column for column in mtmc.CONTIG_COLUMNS if column in table and not table[column].str.strip().eq("-").all()]
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
    # In an edit session, as metaTOR stores its columns: a fresh copy of the
    # map, as metacooler stages every map it writes.
    with metacooler.staged_edit(str(outfile)) as staged:
        mtmc._write_back_columns(staged, str(contig_data), contig_columns)
    return outfile


def store_assembly_columns(metacool, contig_data):
    """Add to the contigs table of a map what metaTOR measures on the
    assembly (ASSEMBLY_COLUMNS: GC_content, and Shotgun_coverage and
    Restriction_site where measured, as contig_data gives them), with
    metacooler annotate, as metaTOR stored them in the maps it built from
    the pairs before 2.0."""
    measured = [column for column, key in ASSEMBLY_COLUMNS.items() if any(str(data[key]) != "-" for data in contig_data.values())]
    if not measured:
        return
    with tempfile.TemporaryDirectory() as scratch:
        table = os.path.join(scratch, "contig_columns.txt")
        with open(table, "w") as handle:
            handle.write("\t".join(["Name", *measured]) + "\n")
            for name, data in contig_data.items():
                handle.write("\t".join([name, *(str(data[ASSEMBLY_COLUMNS[column]]) for column in measured)]) + "\n")
        # A contig of no length, which metacooler leaves out of the map, has a
        # row in the table: no warning about it.
        metacooler.annotate(str(metacool), contigs_table=table, contigs_table_format="Name", unmatched="ignore", inplace=True)
