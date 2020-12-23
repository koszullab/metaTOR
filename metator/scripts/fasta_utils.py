#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A few handy functions for processing fasta files.
"""

from Bio import SeqIO
import argparse
import re

try:
    from metator.scripts.network import DEFAULT_CHUNK_SIZE
except ImportError:
    DEFAULT_CHUNK_SIZE = 1000


def rename_genome(genome_in, genome_out=None):
    """Rename genome and slugify headers

    Rename genomes according to a simple naming scheme; this
    is mainly done to avoid special character weirdness.

    Parameters
    ----------
    genome_in : file, str or pathlib.Path
        The input genome to be renamed and slugify.
    genome_out : file, str or pathlib.Path
        The output genome to be written into. Defaults is <base>_renamed.fa,
        where <base> is genome_in without its extension.
    """

    if genome_out is None:
        genome_out = "{}_renamed.fa".format(genome_in.split(".")[0])

    with open(genome_out, "w") as output_handle:
        for record in SeqIO.parse(genome_in, "fasta"):

            # Replace hyphens, tabs and whitespace with underscores
            new_record_id = record.id.replace(" ", "_")
            new_record_id = new_record_id.replace("-", "_")
            new_record_id = new_record_id.replace("\t", "_")

            # Remove anything that's weird, i.e. not alphanumeric
            # or an underscore
            new_record_id = re.sub("[^_A-Za-z0-9]+", "", new_record_id)
            header = ">{}\n".format(new_record_id)

            output_handle.write(header)
            output_handle.write("{}\n".format(str(record.seq)))


def filter_genome(genome_in, threshold=500, list_records=None):
    """Filter fasta file according to various parameters.

    Filter a fasta file according to size and/or an explicit list of records
    to keep.

    Parameters
    ----------
    genome_in: file, str or pathlib.Path
        The input genome in FASTA format.
    threshold: int, optional
        The size below which genome records are discarded. Default is the
        default minimum chunk size, i.e. 500.
    list_records: array_like, optional
        A list of record ids to keep. If not None, records that don't belong
        to that list are discarded. Default is None, i.e. all records are
        kept.

    Returns
    -------
    records_to_write: generator
        Filtered records that were kept.
    """

    if list_records is None:

        def truth(*args):
            del args
            return True

        is_a_record_to_keep = truth

    else:
        try:
            with open(list_records) as records_handle:
                records_to_keep = records_handle.readlines()
        except OSError:
            if not hasattr(list_records, "__contains__"):
                raise
            else:
                records_to_keep = list_records

        is_a_record_to_keep = records_to_keep.__contains__

    records_to_write = (
        record
        for record in SeqIO.parse(genome_in, "fasta")
        if (len(record.seq) >= threshold and is_a_record_to_keep(record.id))
    )

    return records_to_write


def rename_proteins(prot_in, prot_out=None, chunk_size=DEFAULT_CHUNK_SIZE):
    """Rename prodigal output files

    Rename output files from prodigal according to the following naming
    scheme: >contigX_chunkY__geneZ

    Chunk numbering starts at 0 and gene identification is taken from prodigal.

    Parameters
    ----------
    prot_in : file, str or pathlib.Path
        The input protein file in FASTA format to be renamed.
    prot_out : file, str or pathlib.Path
        The output protein file to be renamed into.
    chunk_size : int, optional
        The size of the chunks (in bp) used in the pipeline. Default is 1000.
    """

    if prot_out is None:
        prot_out = "{}_renamed.fa".format(prot_in.split(".")[0])

    with open(prot_out, "w") as prot_out_handle:

        for record in SeqIO.parse(prot_in, "fasta"):
            header = record.description
            name, pos_start, _, _, _ = header.split("#")

            chunk_start = int(pos_start) // chunk_size

            name_split = name.split("_")
            contig_name = "_".join(name_split[:-1])
            gene_id = name_split[-1]

            new_record_id = "{}_{}__gene{}".format(
                contig_name, chunk_start, gene_id
            )

            prot_out_handle.write(">{}\n".format(new_record_id))
            prot_out_handle.write("{}\n".format(str(record.seq)))


def write_records(records, output_file, split=False):

    """Write FASTA records

    Write a FASTA file from an iterable of records.

    Parameters
    ----------
    records : iterable
        Input records to write.
    output_file : file, str or pathlib.Path
        Output FASTA file to be written into.
    split : bool, optional
        If True, each record is written into its own separate file. Default is
        False.
    """

    if split:
        for record in records:
            with open(
                "{}{}.fa".format(output_file, record.id), "w"
            ) as record_handle:
                SeqIO.write(record, record_handle, "fasta")
    else:
        SeqIO.write(records, output_file, "fasta")


def main():

    parser = argparse.ArgumentParser(
        description="A few utilities"
        " to manipulate "
        "and filter fasta files."
    )
    parser.add_argument("-i", "--input", help="Input fasta file", required=True)
    parser.add_argument("-o", "--output", help="Output file", required=True)

    parser.add_argument(
        "-t",
        "--threshold",
        type=int,
        help="Threshold size of records to keep",
        default=0,
    )

    parser.add_argument(
        "-r", "--records", help="List of records to keep", default=None
    )

    parser.add_argument(
        "-s",
        "--split",
        help="Split input file into single record fasta files",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "-n",
        "--rename",
        help="Rename fasta file ids " "and remove special characters",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-p",
        "--proteins",
        help="Rename prodigal output file",
        action="store_true",
        default=False,
    )

    args = parser.parse_args()

    fasta_file = args.input
    output_file = args.output
    threshold = args.threshold
    records = args.records
    split = args.split
    rename = args.rename
    proteins = args.proteins

    if proteins:
        rename_proteins(prot_in=fasta_file, prot_out=output_file)

    elif rename:
        rename_genome(genome_in=fasta_file, genome_out=output_file)

    else:
        records_to_write = filter_genome(
            genome_in=fasta_file, threshold=threshold, list_records=records
        )

        write_records(
            records=records_to_write, output_file=output_file, split=split
        )


if __name__ == "__main__":

    main()
