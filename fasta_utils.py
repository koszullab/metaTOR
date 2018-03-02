#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A few handy functions for processing fasta files.
"""

from Bio import SeqIO
import argparse
import re

try:
    from network import DEFAULT_CHUNK_SIZE
except ImportError:
    DEFAULT_CHUNK_SIZE = 1000


def rename_genome(genome_in, genome_out=None):
    """Rename genomes according to a simple naming scheme; this
    is mainly done to avoid special character weirdness
    """

    if genome_out is None:
        genome_out = "{}_renamed.fa".format(genome_in.split('.')[0])

    with open(genome_out, "w") as output_handle:
        for record in SeqIO.parse(genome_in, "fasta"):

            # Replace hyphens, tabs and whitespace with underscores
            new_record_id = record.id.replace(" ", "_")
            new_record_id = new_record_id.replace("-", "_")
            new_record_id = new_record_id.replace("\t", "_")

            # Remove anything that's weird, i.e. not alphanumeric
            # or an underscore
            new_record_id = re.sub('[^_A-Za-z0-9]+', '', new_record_id)
            header = ">{}\n".format(new_record_id)

            output_handle.write(header)
            output_handle.write("{}\n".format(str(record.seq)))


def filter_genome(genome_in, threshold=500, list_records=None):
    """Filter fasta file according to various parameters.
    """

    if list_records is None:

        def truth(*x):
            return True
        is_a_record_to_keep = truth

    else:
        try:
            with open(list_records) as records_handle:
                records_to_keep = records_handle.readlines()
        except OSError:
            if not hasattr(list_records, '__contains__'):
                raise
            else:
                records_to_keep = list_records

        is_a_record_to_keep = records_to_keep.__contains__

    records_to_write = (record for record in SeqIO.parse(genome_in, "fasta")
                        if (len(record.seq) >= threshold and
                            is_a_record_to_keep(record.id)))

    return records_to_write


def rename_proteins(prot_in, prot_out=None, chunk_size=DEFAULT_CHUNK_SIZE):
    """Rename prodigal output files
    """

    if prot_out is None:
        prot_out = "{}_renamed.fa".format(prot_in.split('.')[0])

    with open(prot_out, "w") as prot_out_handle:

        for record in SeqIO.parse(prot_in, "fasta"):
            header = record.description
            name, pos_start, pos_end, ori, orf_data = header.split('#')

            chunk_start = int(pos_start) // chunk_size

            name_split = name.split('_')
            contig_name = '_'.join(name_split[:-1])
            gene_id = name_split[-1]

            new_record_id = "{}_{}__gene{}".format(contig_name,
                                                   chunk_start,
                                                   gene_id)

            prot_out_handle.write(">{}\n".format(new_record_id))
            prot_out_handle.write("{}\n".format(str(record.seq)))


def write_records(records, output_file, split=False):

    if split:
        for record in records:
            with open("{}{}.fa".format(output_file,
                                       record.id), "w") as record_handle:
                SeqIO.write(record, record_handle, "fasta")
    else:
        SeqIO.write(records, output_file, "fasta")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='A few utilities'
                                                 ' to manipulate '
                                                 'and filter fasta files.')
    parser.add_argument('-i', '--input',
                        help='Input fasta file', required=True)
    parser.add_argument('-o', '--output',
                        help='Output file', required=True)

    parser.add_argument('-t', '--threshold', type=int,
                        help='Threshold size of records to keep', default=0)

    parser.add_argument('-r', '--records',
                        help='List of records to keep', default=None)

    parser.add_argument('-s', '--split',
                        help='Split input file into single record fasta files',
                        action='store_true', default=False)

    parser.add_argument('-n', '--rename',
                        help='Rename fasta file ids '
                             'and remove special characters',
                        action='store_true', default=False)
    parser.add_argument('-p', '--proteins', help='Rename prodigal output file',
                        action='store_true', default=False)

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
        records_to_write = filter_genome(genome_in=fasta_file,
                                         threshold=threshold,
                                         list_records=records)

        write_records(records=records_to_write,
                      output_file=output_file,
                      split=split)
