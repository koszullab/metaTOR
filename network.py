#!/usr/bin/env python
# coding:utf-8

"""
General utility functions for handling BAM files and generating 3C networks.
"""


import argparse
import collections
import copy
import itertools
import os
import numpy as np
import pysam
import gzip
import operator
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq

DEFAULT_SIZE_CHUNK_THRESHOLD = 500
DEFAULT_MAPQ_THRESHOLD = 10
DEFAULT_CHUNK_SIZE = 1000
DEFAULT_SELF_CONTACTS = False
DEFAULT_NORMALIZED = False
DEFAULT_READ_SIZE = 65
DEFAULT_NETWORK_FILE_NAME = "network.txt"
DEFAULT_CHUNK_DATA_FILE_NAME = "idx_contig_hit_size_cov.txt"

DEFAULT_PARAMETERS = {
    "size_chunk_threshold": DEFAULT_SIZE_CHUNK_THRESHOLD,
    "mapq_threshold": DEFAULT_MAPQ_THRESHOLD,
    "chunk_size": DEFAULT_CHUNK_SIZE,
    "read_size": DEFAULT_READ_SIZE,
    "self_contacts": DEFAULT_SELF_CONTACTS,
    "normalized": DEFAULT_NORMALIZED,
}


def alignment_to_contacts(
    sam_merged,
    assembly,
    output_dir,
    output_file_network=DEFAULT_NETWORK_FILE_NAME,
    output_file_chunk_data=DEFAULT_CHUNK_DATA_FILE_NAME,
    parameters=DEFAULT_PARAMETERS,
):
    """Generates a network file (in edgelist form) from an
    alignment in sam or bam format. Contigs are virtually split into
    'chunks' of nearly fixed size (by default between 500 and 1000 bp)
    to reduce size bias. The chunks are the network nodes and the edges
    are the contact counts.

    The network is in a strict barebone form so that it can be reused and
    imported quickly into other applications etc. Verbose information about
    every single node in the network is written on a 'chunk data' file,
    by default called 'idx_contig_hit_size_cov.txt'
    """

    all_contacts = collections.Counter()
    all_chunks = collections.Counter()

    #   Initialize parameters
    chunk_size = int(parameters["chunk_size"])
    mapq_threshold = int(parameters["mapq_threshold"])
    size_chunk_threshold = int(parameters["size_chunk_threshold"])
    read_size = int(parameters["read_size"])
    self_contacts = parameters["self_contacts"]
    normalized = parameters["normalized"]

    print("Establishing chunk list...")
    chunk_complete_data = dict()

    #   Get all information about all chunks from all contigs
    #   (this gets updated at the end)
    global_id = 1
    for record in SeqIO.parse(assembly, "fasta"):
        length = len(record.seq)

        n_chunks = length // chunk_size
        n_chunks += (length % chunk_size) >= size_chunk_threshold

        for i in range(n_chunks):

            if (i + 1) * chunk_size <= length:
                size = chunk_size
            else:
                size = length % chunk_size

            chunk_name = "{}_{}".format(record.id, i)
            chunk_complete_data[chunk_name] = {
                "id": global_id,
                "hit": 0,
                "size": size,
                "coverage": 0,
            }
            global_id += 1

    print("Opening alignment files...")

    current_read = None

    # Read the BAM file to detect contacts.
    with pysam.AlignmentFile(sam_merged, "rb") as alignment_merged_handle:

        names = alignment_merged_handle.references
        lengths = alignment_merged_handle.lengths
        names_and_lengths = {
            name: length for name, length in itertools.izip(names, lengths)
        }

        print("Reading contacts...")

        # Since the BAM file is supposed to be sorted and interleaved,
        # pairs should be always grouped with one below the other (the exact
        # order doesn't matter since the network is symmetric, so we simply
        # treat the first one as 'forward' and the second one as 'reverse')

        # We keep iterating until two consecutive reads have the same name,
        # discarding ones that don't.

        while "Reading forward and reverse alignments alternatively":
            try:
                my_read = next(alignment_merged_handle)
                if current_read is None:
                    # First read
                    current_read = my_read
                    continue

                elif current_read.query_name != my_read.query_name:

                    # print("{}_{}".format(current_read, my_read))
                    current_read = my_read
                    continue

                read_forward, read_reverse = current_read, my_read

            except StopIteration:
                break

            # Get a bunch of info about the alignments to pass the tests below
            read_name_forward = read_forward.query_name
            read_name_reverse = read_reverse.query_name

            flag_forward, flag_reverse = read_forward.flag, read_reverse.flag

            try:
                assert read_name_forward == read_name_reverse
            except AssertionError:
                print(
                    "Reads don't have the same name: "
                    "{} and {}".format(read_name_forward, read_name_reverse)
                )
                raise

            # To check if a flag contains 4
            # (digit on the third position from the right in base 2),
            # 4 = unmapped in SAM spec
            def is_unmapped(flag):
                return np.base_repr(flag, padding=3)[-3] == "1"

            if is_unmapped(flag_forward) or is_unmapped(flag_reverse):
                # print("Detected unmapped read on one end, skipping")
                continue

            contig_name_forward = read_forward.reference_name
            contig_name_reverse = read_reverse.reference_name

            len_contig_for = names_and_lengths[contig_name_forward]
            len_contig_rev = names_and_lengths[contig_name_reverse]

            position_forward = read_forward.reference_start
            position_reverse = read_reverse.reference_start

            mapq_forward = read_forward.mapping_quality
            mapq_reverse = read_reverse.mapping_quality

            # Some more tests: checking for size, map quality, map status etc.
            mapq_test = min(mapq_forward, mapq_reverse) > mapq_threshold

            min_length = min(len_contig_for, len_contig_rev)
            length_test = min_length > size_chunk_threshold

            # Trickest test:
            #
            #
            #                contig
            #    pos1                          pos2
            #     ^                             ^
            # |-------|-------|-------|-------|---|
            # <-------><------><------><------><-->            <->
            #   chunk   chunk                  tail   size_chunk_threshold
            #
            # Test is passed if tail >= size_chunk_threshold (pos2)
            # or if the position is a non-tail chunk (pos1)

            if position_forward < chunk_size * (len_contig_for // chunk_size):
                current_chunk_forward_size = chunk_size
            else:
                current_chunk_forward_size = len_contig_for % chunk_size

            if position_reverse < chunk_size * (len_contig_rev // chunk_size):
                current_chunk_reverse_size = chunk_size
            else:
                current_chunk_reverse_size = len_contig_rev % chunk_size

            min_chunk_size = min(
                current_chunk_forward_size, current_chunk_reverse_size
            )

            chunk_test = min_chunk_size >= size_chunk_threshold

            if mapq_test and length_test and chunk_test:

                chunk_forward = position_forward // chunk_size
                chunk_reverse = position_reverse // chunk_size

                chunk_name_forward = "{}_{}".format(
                    contig_name_forward, chunk_forward
                )
                chunk_name_reverse = "{}_{}".format(
                    contig_name_reverse, chunk_reverse
                )

                # print("Detected contact between "
                #       "{} and {}".format(chunk_name_forward,
                #                          chunk_name_reverse))

                if self_contacts or chunk_name_forward != chunk_name_reverse:

                    contact = tuple(
                        sorted((chunk_name_forward, chunk_name_reverse))
                    )

                    all_contacts[contact] += 1

                    # print("I add {} of"
                    #       " size {} to"
                    #       " my chunk "
                    #       "dictionary".format(chunk_name_forward,
                    #                           current_chunk_forward_size))
                    chunk_key_forward = (
                        chunk_name_forward,
                        current_chunk_forward_size,
                    )
                    all_chunks[chunk_key_forward] += 1

                    # print("I add {} of"
                    #       " size {} to"
                    #       " my chunk "
                    #       "dictionary".format(chunk_name_reverse,
                    #                           current_chunk_reverse_size))
                    chunk_key_reverse = (
                        chunk_name_reverse,
                        current_chunk_reverse_size,
                    )
                    all_chunks[chunk_key_reverse] += 1

    print("Writing chunk data...")

    # Now we can update the chunk dictionary
    # with the info we gathered from the BAM file

    output_chunk_data_path = os.path.join(output_dir, output_file_chunk_data)

    with open(output_chunk_data_path, "w") as chunk_data_file_handle:

        for name in sorted(chunk_complete_data.keys()):

            chunk_data = chunk_complete_data[name]
            size = chunk_data["size"]
            chunk = (name, chunk_data["size"])
            hit = all_chunks[chunk]
            coverage = hit * read_size * 1.0 / size
            try:
                chunk_complete_data[name]["hit"] = hit
                chunk_complete_data[name]["coverage"] = coverage
            except KeyError:
                print(
                    "A mismatch was detected between the reference "
                    "genome and the genome used for the alignment "
                    "file, some sequence names were not found"
                )
                raise

            idx = chunk_complete_data[name]["id"]
            line = "{}\t{}\t{}\t{}\t{}\n".format(
                idx, name, hit, size, coverage
            )
            chunk_data_file_handle.write(line)

    # Lastly, generate the network proper

    print("Writing network...")

    output_network_path = os.path.join(output_dir, output_file_network)

    with open(output_network_path, "w") as network_file_handle:

        for chunks in sorted(all_contacts.keys()):

            chunk_name1, chunk_name2 = chunks
            contact_count = all_contacts[chunks]

            if normalized:
                coverage1 = chunk_complete_data[chunk_name1]["coverage"]
                coverage2 = chunk_complete_data[chunk_name2]["coverage"]
                mean_coverage = np.sqrt(coverage1 * coverage2)
                effective_count = contact_count * 1.0 / mean_coverage
            else:
                effective_count = contact_count

            try:
                idx1 = chunk_complete_data[chunk_name1]["id"]
                idx2 = chunk_complete_data[chunk_name2]["id"]
                line = "{}\t{}\t{}\n".format(idx1, idx2, effective_count)
                network_file_handle.write(line)
            except KeyError as e:
                print("Mismatch detected: {}".format(e))

    return chunk_complete_data, all_contacts


def merge_networks(output_file="merged_network.txt", *files):
    """A naive implementation for merging two edgelists.

    :note: The partitioning step doesn't mind redundant
    edges and handles them pretty well, so if you are not
    using the merged edgelist for anything else you can just
    concatenate the edgelists without using this function.
    """

    contacts = dict()
    for network_file in files:
        with open(network_file) as network_file_handle:
            for line in network_file_handle:
                id_a, id_b, n_contacts = line.split("\t")
                pair = sorted((id_a, id_b))
                try:
                    contacts[pair] += n_contacts
                except KeyError:
                    contacts[pair] = n_contacts

    sorted_contacts = sorted(contacts)
    with open(output_file, "w") as output_handle:
        for index_pair in sorted_contacts:
            id_a, id_b = index_pair
            n_contacts = contacts[index_pair]
            output_handle.write("{}\t{}\t{}\n".format(id_a, id_b, n_contacts))


def merge_chunk_data(output_file="merged_idx_contig_hit_size_cov.txt", *files):
    """Merge any number of chunk data files.
    """

    chunks = dict()
    for chunk_file in files:
        with open(chunk_file) as chunk_file_handle:
            for line in chunk_file_handle:
                chunk_id, chunk_name, hit, size, cov = line.split("\t")
                try:
                    chunks[chunk_id]["hit"] += hit
                    chunks[chunk_id]["cov"] += cov
                except KeyError:
                    chunks[chunk_id] = {
                        "name": chunk_name,
                        "hit": hit,
                        "size": size,
                        "cov": cov,
                    }

    sorted_chunks = sorted(chunks)
    with open(output_file, "w") as output_handle:
        for chunk_id in sorted_chunks:
            my_chunk = chunks[chunk_id]
            name, hit, size, cov = (
                my_chunk["name"],
                my_chunk["hit"],
                my_chunk["size"],
                my_chunk["cov"],
            )

            my_line = "{}\t{}\t{}\t{}\t{}".format(
                chunk_id, name, hit, size, cov
            )
            output_handle.write(my_line)


def alignment_to_reads(
    sam_merged,
    output_dir,
    parameters=DEFAULT_PARAMETERS,
    save_memory=True,
    *bin_fasta
):
    """Extract reads found to be mapping an input FASTA bin.
    If one read maps, the whole pair is extracted and written
    to the output paired-end FASTQ files. Reads that mapped
    and weren't part of a pair are kept in a third 'single'
    file for people who need it (e.g. to get extra paired reads
    by fetching the opposite one from the original FASTQ library).

    :note: will throw an IOError ('close failed in file
    object destructor') on exit with older versions
    of pysam for some reason. It's harmless but annoying
    so consider upgrading that if it comes up in a pipeline.
    """

    #   Just in case file objects are sent as input
    def get_file_string(file_thing):
        try:
            file_string = file_thing.name
        except AttributeError:
            file_string = str(file_thing)
        return file_string

    #   Global set of chunks against which reads are required to
    #   map - we store them in a tuple that keeps track of the
    #   original bin each chunk came from so we can reattribute the reads later

    bin_chunks = set()
    for bin_file in bin_fasta:
        for record in SeqIO.parse(bin_file, "fasta"):
            bin_chunks.add((get_file_string(bin_file), record.id))

    chunk_size = int(parameters["chunk_size"])

    mapq_threshold = int(parameters["mapq_threshold"])

    def read_name(read):
        return read.query_name.split()[0]

    #   Since reading a huge BAM file can take up a
    #   lot of time and resources, we only do it once
    #   but that requires opening fastq files for writing
    #   as matching reads get detected along the
    #   bam and keeping track of which ones are
    #   currently open.

    def get_base_name(bin_file):

        base_name = ".".join(os.path.basename(bin_file).split(".")[:-1])

        output_path = os.path.join(
            output_dir, "{}.readnames".format(base_name)
        )

        return output_path

    if save_memory:
        opened_files = dict()
    else:
        read_names = dict()

    with pysam.AlignmentFile(sam_merged, "rb") as alignment_merged_handle:

        for (my_read_name, alignment_pool) in itertools.groupby(
            alignment_merged_handle, read_name
        ):

            for my_alignment in alignment_pool:

                relative_position = my_alignment.reference_start
                contig_name = my_alignment.reference_name

                chunk_position = relative_position // chunk_size

                # The 'chunk name' is used to detect macthing positions
                chunk_name = "{}_{}".format(contig_name, chunk_position)

                # But such matching positions have to map acceptably
                quality_test = my_alignment.mapping_quality > mapq_threshold

                for bin_file in bin_fasta:
                    chunk_tuple = (bin_file, chunk_name)
                    if chunk_tuple in bin_chunks and quality_test:
                        if save_memory:
                            output_path = get_base_name(bin_file)
                            try:
                                output_handle = opened_files[bin_file]
                            except KeyError:
                                output_handle = open(output_path, "w")
                                opened_files[bin_file] = output_handle

                            output_handle.write("@{}\n".format(my_read_name))

                        else:
                            try:
                                read_names[my_read_name].append(bin_file)
                            except KeyError:
                                read_names[my_read_name] = [bin_file]

    for file_handle in opened_files.values():
        file_handle.close()
    #   Return unpaired file names for pair_unpaired_reads() to process
    if save_memory:
        return opened_files.keys()
    else:
        return read_names


def retrieve_reads_from_fastq(
    fastq_forward, fastq_reverse, read_names, output_dir
):

    opened_files = dict()

    read_set = set(read_names.keys())

    def get_base_names(bin_file):

        base_name = ".".join(os.path.basename(bin_file).split(".")[:-1])

        for_fastq_path = os.path.join(
            output_dir, "{}_for.fastq".format(base_name)
        )
        rev_fastq_path = os.path.join(
            output_dir, "{}_rev.fastq".format(base_name)
        )

        return for_fastq_path, rev_fastq_path

    def seamless_open(my_file):

        try:
            return gzip.open(my_file)
        except IOError:
            return open(my_file)

    with seamless_open(fastq_forward) as forward_handle:
        with seamless_open(fastq_reverse) as reverse_handle:

            forward_fastq_iterator = FastqGeneralIterator(forward_handle)
            reverse_fastq_iterator = FastqGeneralIterator(reverse_handle)

            for read_for, read_rev in itertools.izip(
                forward_fastq_iterator, reverse_fastq_iterator
            ):

                name_for, sequence_for, quality_for = read_for
                name_rev, sequence_rev, quality_rev = read_rev
                short_name_for = name_for.split()[0]
                short_name_rev = name_rev.split()[0]

                if short_name_for in read_set or short_name_rev in read_set:
                    try:
                        bin_files = read_names[short_name_for]
                    except KeyError:
                        bin_files = read_names[short_name_rev]

                    for bin_file in bin_files:
                        (for_fastq_path, rev_fastq_path) = get_base_names(
                            bin_file
                        )

                        # Files are opened 'lazily' so as to not end up
                        # with a bunch of empty FASTQ files
                        try:
                            output_for_handle = opened_files[for_fastq_path]
                        except KeyError:
                            output_for_handle = open(for_fastq_path, "w")
                            opened_files[for_fastq_path] = output_for_handle

                        line_for = "@{}\n{}\n+\n{}\n".format(
                            name_for, sequence_for, quality_for
                        )
                        output_for_handle.write(line_for)

                        try:
                            output_rev_handle = opened_files[rev_fastq_path]
                        except KeyError:
                            output_rev_handle = open(rev_fastq_path, "w")
                            opened_files[rev_fastq_path] = output_rev_handle

                        line_rev = "@{}\n{}\n+\n{}\n".format(
                            name_rev, sequence_rev, quality_rev
                        )
                        output_rev_handle.write(line_rev)

    #   Close all the files we've opened for writing
    for file_handle in opened_files.values():
        file_handle.close()


def retrieve_reads_contig_wise(sam_merged, contig_data, output_dir):

    contig_dict = dict()
    opened_files = dict()

    def close_all_files():
        for my_file in opened_files.values():
            my_file.close()

    with open("contig_data.txt") as contig_data:
        for line in contig_data:
            fields = line.split()
            node = fields[0]
            core = fields[-3]
            contig_dict[node] = core

    query_getter = operator.attrgetter("query_name")

    with pysam.AlignmentFile(sam_merged, "rb") as sam_handle:
        for (my_read_name, alignment_pool) in itertools.groupby(
            sam_handle, query_getter
        ):
            my_read_set = dict()
            my_core_set = set()
            while "Reading alignments from alignment pool":
                try:
                    my_alignment = next(alignment_pool)
                    # print(contig_dict[my_alignment.reference_name])

                    is_reverse = my_alignment.is_reverse

                    my_seq = my_alignment.query_sequence

                    my_qual = my_alignment.query_qualities
                    my_qual_string = pysam.array_to_qualitystring(my_qual)

                    if is_reverse:
                        my_seq_string = str(Seq(my_seq).reverse_complement())
                        my_qual_string = my_qual_string[::-1]
                    else:
                        my_seq_string = my_seq

                    my_seq_tuple = (my_seq_string, my_qual_string)

                    if len(my_read_set) > 2:
                        print("We've got a serious problem")
                        print(my_read_set)
                    elif len(my_read_set) == 0:
                        my_read_set[my_seq_tuple] = "forward"
                    elif my_seq_tuple not in my_read_set.keys():
                        my_read_set[my_seq_tuple] = "reverse"
                    try:
                        ref = contig_dict[my_alignment.reference_name]
                        my_core_set.add(ref)
                    except KeyError:
                        my_core_set.add("unknown")
                except StopIteration:
                    if len(my_read_set) == 2:
                        for core_name in my_core_set:
                            for my_tuple, file_id in my_read_set.iteritems():
                                if file_id == "forward":
                                    file_end = ".end1"
                                elif file_id == "reverse":
                                    file_end = ".end2"

                                basename = "{}{}".format(core_name, file_end)
                                filename = os.path.join(output_dir, basename)
                                try:
                                    file_to_write = opened_files[filename]
                                except KeyError:
                                    file_to_write = open(filename, "w")
                                    opened_files[filename] = file_to_write
                                except IOError:
                                    print(len(opened_files))
                                    raise

                                seq, qual = my_tuple
                                line = "@{}\n{}\n+\n{}\n".format(
                                    my_read_name, seq, qual
                                )
                                file_to_write.write(line)
                    elif len(my_read_set) == 1:
                        for core_name in my_core_set:
                            file_end = ".end"
                            basename = "{}{}".format(core_name, file_end)
                            filename = os.path.join(output_dir, basename)
                            try:
                                file_to_write = opened_files[filename]
                            except KeyError:
                                file_to_write = open(filename, "w")
                                opened_files[filename] = file_to_write
                            except IOError:
                                print(len(opened_files))
                                raise
                            seq, qual = my_read_set.keys()[0]
                            line = "@{}\n{}\n+\n{}\n".format(
                                my_read_name, seq, qual
                            )
                            file_to_write.write(line)
                    else:
                        print("We've got a very serious problem")
                    break

    close_all_files()


def main():

    parser = argparse.ArgumentParser(description="Process alignment files.")

    parser.add_argument(
        "-i",
        "--input",
        type=argparse.FileType("r"),
        help="Merged (interlaced) alignment file",
        required=True,
    )

    parser.add_argument(
        "-f", "--reference", help="Reference fasta file", nargs="+"
    )

    parser.add_argument(
        "-o", "--output", help="Output directory", required=True
    )

    parser.add_argument(
        "-q",
        "--map-quality",
        type=int,
        help="Minimum mapping quality threshold",
        default=DEFAULT_PARAMETERS["mapq_threshold"],
    )

    parser.add_argument(
        "-c",
        "--chunk-size",
        type=int,
        help="Standard chunk size for nodes in the network",
        default=DEFAULT_PARAMETERS["chunk_size"],
    )

    parser.add_argument(
        "-s",
        "--size-chunk-threshold",
        type=int,
        help="Minimum size for tail ends to be integrated"
        " as chunks in the network",
        default=DEFAULT_PARAMETERS["size_chunk_threshold"],
    )

    parser.add_argument(
        "-a",
        "--self-contacts",
        action="store_true",
        help="Do not discard self contacts",
        default=DEFAULT_PARAMETERS["self_contacts"],
    )

    parser.add_argument(
        "-n",
        "--normalize",
        action="store_true",
        help="Normalize contacts by the geometric mean"
        " of both coverages of the chunks",
        default=DEFAULT_PARAMETERS["normalized"],
    )

    parser.add_argument(
        "-r",
        "--read-size",
        help="Read size",
        default=DEFAULT_PARAMETERS["read_size"],
    )

    parser.add_argument(
        "-F",
        "--fastq",
        help="Reconstruct FASTQ from "
        "reference and interleaved alignment file",
        nargs=2,
    )

    parser.add_argument(
        "-D",
        "--fastq-contig",
        help="Reconstruct FASTQ from "
        "contig data and interleaved alignment file",
        nargs=1,
    )

    parser.add_argument(
        "-m",
        "--mem",
        help="Save memory by only writing " " read names in FASTQ mode",
        action="store_true",
    )

    args = parser.parse_args()

    merged_file = args.input
    reference_file = args.reference
    output_dir = args.output
    fastq_files = args.fastq
    contig_data = args.fastq_contig
    save_memory = args.mem

    parameters = copy.deepcopy(DEFAULT_PARAMETERS)

    parameters["mapq_threshold"] = args.map_quality
    parameters["chunk_size"] = args.chunk_size
    parameters["read_size"] = args.read_size
    parameters["size_chunk_threshold"] = args.size_chunk_threshold
    parameters["self_contacts"] = args.self_contacts
    parameters["normalized"] = args.normalize

    if fastq_files:

        fastq_forward, fastq_reverse = fastq_files
        read_names = alignment_to_reads(
            merged_file, output_dir, parameters, save_memory, *reference_file
        )

        if save_memory:
            retrieve_reads_from_fastq(
                fastq_forward, fastq_reverse, read_names, output_dir
            )

    elif contig_data:
        retrieve_reads_contig_wise(merged_file, contig_data, output_dir)

    else:

        my_assembly, = reference_file
        alignment_to_contacts(
            sam_merged=merged_file,
            assembly=my_assembly,
            output_dir=output_dir,
            parameters=parameters,
        )


if __name__ == "__main__":

    main()
