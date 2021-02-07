#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Generate 3C networks of teh contig.
General utility functions for handling BAM files and generating 3C networks.

Core function to build the network are:
    -alignement_to_contacts
    -compute_contig_coverage
    -compute_network
    -create_contig_data
    -precompute_network
    -write_contig_data

Some old functions have been kept but won't work with the new files:
    -merge_networks
    -merge_chunk_data
    -alignment_to_reads
    -retrieve_reads_from_fastq
    -retrieve_reads_contig_wise
"""

import csv
import numpy as np
from Bio import SeqIO
from Bio import SeqUtils
from os.path import join
from metator.log import logger
import metator.io as mio


def alignment_to_contacts(
    bed2D_file,
    genome,
    output_dir,
    output_file_network="network.txt",
    output_file_contig_data="contig_data_network.txt",
    tmpdir=".",
    n_cpus=1,
    normalized=True,
    self_contacts=False,
):
    """Generates a network file (in edgelist form) from an alignment in bed2D
    format. Contigs are the network nodes and the edges are the contact counts.

    The network is in a strict barebone form so that it can be reused and
    imported quickly into other applications etc. Verbose information about
    every single node in the network is written on a 'contig data' file, by
    default called 'idx_contig_length_GC_hit_cov.txt'

    Parameters
    ----------
    bed2D_file : str
        Path to the bed2D file used as input
    genome : str
        The initial assembly path acting as the alignment file's reference
        genome.
    output_dir : str
        The output directory to write the network and chunk data into.
    output_file_network : str, optional
        The specific file name for the output network file. Default is
        'network.txt'
    output_file_contig_data : str, optional
        The specific file name for the output chunk data file. Default is
        'idx_contig_length_GC_hit_cov.txt'
    tmpdir : str
        Path to th temporary directory. Default in the working directory
    self_contacts : bool
        Whether to count alignments between a contig and
        itself. Default is False.
    normalized : bool
        Whether to normalize contacts between contigs by their geometric mean
        coverage. Default is True.

    Returns
    -------
    dict
        A dictionary where the keys are chunks in (contig, position) form and
        the values are their id, name, total contact count, size and coverage.
    """

    # Create temporary and output file which will be necessary
    precompute_network_file = join(tmpdir, "precompute_network_file.txt")
    network_file = join(output_dir, output_file_network)
    contig_data_file = join(output_dir, output_file_contig_data)

    # Create the contig data dictionnary
    contig_data = create_contig_data(genome)

    # Create a contact file easily readable for counting the contacts
    contig_data = precompute_network(
        bed2D_file, contig_data, precompute_network_file, self_contacts
    )

    # Compute the coverage of the contigs
    contig_data = compute_contig_coverage(contig_data=contig_data)

    # Compute network
    compute_network(
        precompute_network_file,
        network_file,
        contig_data,
        tmpdir,
        n_cpus,
        normalized,
    )

    # Write the data from the contigs
    write_contig_data(contig_data, contig_data_file)

    return contig_data


def compute_contig_coverage(contig_data):
    """Compute the coverage of the contigs from the hit information compute
    previously

    Parameters
    ----------
    contig_data : dict
        Dictionnary of the all the contigs from the assembly, the contigs names
        are the keys to the data of the contig available with the following
        keys: "id", "length", "GC", "hit", "coverage".

    Return
    ------
    dict:
        Dictionnary of the all the contigs from the assembly, the contigs names
        are the keys to the data of the contig available with the following
        keys: "id", "length", "GC", "hit", "coverage".
    """

    # Loop on the dictionnary and compute the coverage
    for name in contig_data:
        contig_data[name]["coverage"] = (
            contig_data[name]["hit"] * 1.0 / contig_data[name]["length"]
        )
    return contig_data


def compute_network(
    pre_network_file,
    network_file,
    contig_data,
    tmp_dir,
    n_cpus,
    normalized=True,
):
    """Generate the network file from the prenetwork file.

    Parameters
    ----------
    pre_network_file  : str
        Path of the input file (prenetwork file, output from teh precompute
        network function)
    network_file : str
        Path of the output file (network file).
    contig_data : dict
        Dictionnary of the all the contigs from the assembly, the contigs names
        are the keys to the data of the contig available with the following
        keys: "id", "length", "GC", "hit", "coverage".
    tmp_dir : str
        Path of the temporary directory to write files.
    n_cpus : int
        Number of cpus used to sort the prenetwork.
    normalized : bool
        If true, normalized the count of a contact by the geometric mean of the
        coverage od the contigs.
    """

    # Create a temporary directory for the sorted pre-network
    pre_network_sorted = join(tmp_dir, "tmp_network_sorted.txt")

    # Sort the the pre-network file
    mio.sort_pairs(
        pre_network_file,
        pre_network_sorted,
        tmp_dir=None,
        threads=8,
        buffer="2G",
    )

    # Set the variables used in the loop
    n_occ = 1  # Number of occurences of each frag combination
    n_nonzero = 1  # Total number of nonzero matrix entries
    n_pairs = 0  # Total number of pairs entered in the matrix

    # Read the sorted paors
    with open(pre_network_sorted, "r") as pairs, open(network_file, "w") as net:
        pairs_reader = csv.reader(pairs, delimiter="\t")
        prev_pair = next(pairs_reader)
        for pair in pairs_reader:
            # Extract the contig names and their id
            curr_pair = pair
            # Increment number of occurences for fragment pair
            if prev_pair == curr_pair:
                n_occ += 1
            # Write previous pair and start a new one
            else:
                if n_occ > 0:

                    # Normalisation using the geometric mean of the coverage
                    if normalized:
                        mean_coverage = np.sqrt(
                            contig_data[prev_pair[0]]["coverage"]
                            * contig_data[prev_pair[1]]["coverage"]
                        )
                        effective_count = n_occ / mean_coverage
                    else:
                        effective_count = n_occ

                    # Write the line in the file
                    net.write(
                        "\t".join(
                            map(
                                str,
                                [
                                    contig_data[prev_pair[0]]["id"],
                                    contig_data[prev_pair[1]]["id"],
                                    effective_count,
                                ],
                            )
                        )
                        + "\n"
                    )
                # Reset the values for the next pair
                prev_pair = curr_pair
                n_pairs += n_occ
                n_occ = 1
                n_nonzero += 1

        # Write the last value
        if normalized:
            mean_coverage = np.sqrt(
                contig_data[curr_pair[0]]["coverage"]
                * contig_data[curr_pair[1]]["coverage"]
            )
            effective_count = n_occ / mean_coverage
        else:
            effective_count = n_occ

        net.write(
            "\t".join(
                map(
                    str,
                    [
                        contig_data[curr_pair[0]]["id"],
                        contig_data[curr_pair[1]]["id"],
                        effective_count,
                    ],
                )
            )
            + "\n"
        )
        n_nonzero += 1
        n_pairs += n_occ


def create_contig_data(assembly):
    """Create a dictionnary with data on each Contig.

    Paramaters
    ----------
    assembly : str
        Path to the assembly fasta file

    Return
    ------
    dict():
        Dictionnary of the all the contigs from the assembly, the contigs names
        are the keys to the data of the contig available with the following
        keys: "id", "length", "GC", "hit", "coverage". Hit and coverage are set
        to 0 and need to be updated later.
    """
    # TODO: Avoid using the assembly file or make something with the choice. Now
    # we are avoiding to build index each time.
    # # Throw error if index does not exist
    # index = mio.check_fasta_index(ref, mode="bowtie2")
    # if index is None:
    #     logger.error(
    #         "Reference index is missing, please build the bowtie2 index first."
    #     )
    #     sys.exit(1)

    # # Extract the name of the contigs and their length from the bowtie2 index
    # cmd = (
    #     "bowtie2-inspect -s {}"
    # ).format(index)

    # Create a contig data dictionnary from the assembly file
    contig_data = dict()
    global_id = 1
    for contig in SeqIO.parse(assembly, "fasta"):
        contig_name = contig.id
        contig_length = len(contig.seq)
        contig_GC = SeqUtils.GC(contig.seq)
        contig_data[contig_name] = {
            "id": global_id,
            "length": contig_length,
            "GC": contig_GC,
            "hit": 0,
            "coverage": 0,
        }
        global_id += 1
    return contig_data


def precompute_network(bed2D_file, contig_data, out_file, self_contacts=False):
    """Write a file with only the contig id separated by a tabulation and count
    the contacts by contigs to be able to compute directlty the normalized
    network.

    Parameters
    ----------
    bed2D_file : str
        Path to the bed2D file.
    contig_data : dict
        Dictionnary of the all the contigs from the assembly, the contigs names
        are the keys to the data of the contig available with the following
        keys: "id", "length", "GC", "hit", "coverage". Coverage still at 0 and
        need to be updated later.
    out_file : str
        Path to the write the output_file which will be necessary to compute the
        network.
    self_contacts : bool
        If True, the contacts on the same contigs will be displayed. Otherwise
        only displays the inter contigs contacts.

    Return
    ------
    dict:
        Dictionnary of the all the contigs from the assembly, the contigs names
        are the keys to the data of the contig available with the following
        keys: "id", "length", "GC", "hit", "coverage". Coverage still at 0 and
        need to be updated later.
    """
    # Initiate value to compute 3D ratio
    all_contacts = 0
    inter_contacts = 0

    # Prepare a file to save contact with their global ID
    with open(out_file, "w") as pre_net:

        # Read the 2D bed file and build pairs for the network
        with open(bed2D_file, "r") as pairs:
            for pair in pairs:

                # Split the line on the tabulation
                p = pair.split("\t")

                # Extract the contig names which are at the position 2 and 6.
                contig1, contig2 = p[1], p[5]

                id1 = contig_data[contig1]["id"]
                id2 = contig_data[contig2]["id"]

                # Count the contact
                all_contacts += 1
                contig_data[contig1]["hit"] += 1
                contig_data[contig2]["hit"] += 1

                # Write the file used for the computation of the network.
                if self_contacts and id1 == id2:
                    pre_net.write(
                        "\t".join(map(str, [contig1, contig2])) + "\n"
                    )
                    # pre_net[pair] += 1
                elif id1 < id2:
                    inter_contacts += 1
                    pre_net.write(
                        "\t".join(map(str, [contig1, contig2])) + "\n"
                    )
                    # pre_net[pair] += 1
                elif id1 > id2:
                    inter_contacts += 1
                    pre_net.write(
                        "\t".join(map(str, [contig2, contig1])) + "\n"
                    )
                    # pre_net[pair] += 1
        pairs.close()
    pre_net.close()

    # Return information about the network
    logger.info("{0} contacts in the library.".format(all_contacts))
    logger.info(
        "{0} contacts inter-contigs in the library.".format(inter_contacts)
    )
    logger.info("3D ratio : {0}".format(inter_contacts / all_contacts))

    return contig_data


def write_contig_data(contig_data, output_path):
    """Function to write the contig data file at the output path given. The file
    will contains 6 columns separated by a tabulation: id, name, length,
    GC_content, hit, coverage for each contig.

    Parameters
    ----------
    contig_data : dict
        Dictionnary of the all the contigs from the assembly, the contigs names
        are the keys to the data of the contig available with the following
        keys: "id", "length", "GC", "hit", "coverage".
    output_path : str
        Path to the output file where the data from the dictionnary will be
        written
    """

    # For each contig extract the data and write them in the file.
    with open(output_path, "w") as contig_data_file_handle:
        for name in contig_data:
            length = contig_data[name]["length"]
            hit = contig_data[name]["hit"]
            coverage = contig_data[name]["coverage"]
            GC_content = contig_data[name]["GC"]
            idx = contig_data[name]["id"]
            line = "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                idx, name, length, GC_content, hit, coverage
            )
            contig_data_file_handle.write(line)


# TO REMOVE ?:


def merge_networks(output_file="merged_network.txt", *files):
    """Merge networks into a larger network.

    A naive implementation for merging two networks in edgelist format.

    Parameters
    ---------
    output_file : file, str, or pathlib.Path, optional
        The output file to write the merged network into. Default is
        merged_network.txt
    `*files` : file, str or pathlib.Path
        The network files to merge.

    Note
    ----
    The partitioning step doesn't mind redundant edges and handles them pretty
    well, so if you are not using the merged edgelist for anything else you can
    just concatenate the edgelists without using this function.
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
    """Merge chunk data from different networks

    Similarly to merge_network, this merges any number of chunk data files.

    Parameters
    ---------
    output_file : file, str, or pathlib.Path, optional
        The output file to write the merged chunk data files into. Default is
        merged_idx_contig_hit_size_cov.txt
    `*files` : file, str or pathlib.Path
        The chunk data files to merge.
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
    sam_merged, output_dir, parameters, save_memory=True, *bin_fasta
):
    """Generate reads from ambiguous alignment file

    Extract reads found to be mapping an input FASTA bin. If one read maps, the
    whole pair is extracted and written to the output paired-end FASTQ files.
    Reads that mapped and weren't part of a pair are kept in a third 'single'
    file for people who need it (e.g. to get extra paired reads by fetching the
    opposite one from the original FASTQ library).

    Parameters
    ----------
    sam_merged : file, str or pathlib.Path
        The input alignment file in SAM/BAM format to be processed.
    output_dir : str or pathlib.Path
        The output directory to write the network and chunk data into.
    parameters : dict, optional
        Parameters for the network to read conversion, similar to
        alignment_to_network.
    save_memory : bool, optional
        Whether to keep the read names into memory or write them in different
        files, which takes longer but may prevent out-of-memory crashes. Default
        is True.
    `*bin_fasta` : file, str or pathlib.Path
        The bin FASTA files with appropriately named records.

    Returns
    -------
    A dictionary of files with read names for each bin if save_memory is True,
    and a dictionary of the read names lists themselves otherwise.


    Note
    ----
    This will throw an IOError ('close failed in file object destructor') on
    exit with older versions of pysam for some reason. It's harmless but you may
    consider upgrading to a later version of pysam if it comes up in a pipeline.
    """

    #   Just in case file objects are sent as input
    def get_file_string(file_thing):
        try:
            file_string = file_thing.name
        except AttributeError:
            file_string = str(file_thing)
        return file_string

    #   Global set of chunks against which reads are required to map - we store
    #   them in a tuple that keeps track of the original bin each chunk came
    #   from so we can reattribute the reads later

    bin_chunks = set()
    for bin_file in bin_fasta:
        for record in SeqIO.parse(bin_file, "fasta"):
            bin_chunks.add((get_file_string(bin_file), record.id))

    chunk_size = int(parameters["chunk_size"])

    mapq_threshold = int(parameters["mapq_threshold"])

    def read_name(read):
        return read.query_name.split()[0]

    #   Since reading a huge BAM file can take up a lot of time and resources,
    #   we only do it once but that requires opening fastq files for writing as
    #   matching reads get detected along the bam and keeping track of which
    #   ones are currently open.

    def get_base_name(bin_file):

        base_name = ".".join(os.path.basename(bin_file).split(".")[:-1])

        output_path = os.path.join(output_dir, "{}.readnames".format(base_name))

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

            for read_for, read_rev in zip(
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
                        logger.warning(
                            "Something's gone wrong with read set %s, as "
                            "there are %s of them",
                            my_read_name,
                            len(my_read_set),
                        )
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
                            for my_tuple, file_id in my_read_set.items():
                                if file_id == "forward":
                                    file_end = ".end1"
                                elif file_id == "reverse":
                                    file_end = ".end2"

                                basename = "{core_name}{file_end}".format(
                                    core_name, file_end
                                )
                                filename = os.path.join(output_dir, basename)
                                try:
                                    file_to_write = opened_files[filename]
                                except KeyError:
                                    file_to_write = open(filename, "w")
                                    opened_files[filename] = file_to_write
                                except IOError:
                                    logger.error(
                                        "Error when trying to handle"
                                        "%s. Maybe there are too many opened"
                                        "files at once: %s",
                                        filename,
                                        len(opened_files),
                                    )
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
                                logger.error(
                                    "Error when trying to handle"
                                    "%s. Maybe there are too many opened"
                                    "files at once: %s",
                                    filename,
                                    len(opened_files),
                                )
                                raise
                            seq, qual = my_read_set.keys()[0]
                            line = "@{}\n{}\n+\n{}\n".format(
                                my_read_name, seq, qual
                            )
                            file_to_write.write(line)
                    else:
                        logger.warning(
                            "Something's gone wrong with read set %s, as "
                            "there are %s of them",
                            my_read_name,
                            len(my_read_set),
                        )
                    break

    close_all_files()