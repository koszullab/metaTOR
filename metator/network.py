#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Generate metaHiC contigs network from fastq reads or bam files.

General utility functions for handling aligment files from align module and
generating metaHiC networks.

Core function to build the network are:
    - alignment_to_contacts
    - compute_network
    - create_contig_data
    - normalize_pair
    - precompute_network
    - write_contig_data
    - write_hit_data
"""

import csv
import numpy as np
import re
from Bio import SeqIO
from Bio import SeqUtils
from os.path import join, basename
from metator.log import logger
import metator.io as mio


def alignment_to_contacts(
    alignment_files,
    contig_data,
    hit_data,
    out_dir,
    output_file_network,
    output_file_contig_data,
    tmp_dir,
    n_cpus,
    normalization,
    self_contacts,
):
    """Generates a network file (in edgelist form) from an alignment. Contigs
    are the network nodes and the edges are the contact counts.

    The network is in a strict barebone form so that it can be reused and
    imported quickly into other applications etc. Verbose information about
    every single node in the network is written on a 'contig data' file.

    Parameters:
    -----------
    alignment_files : list of str
        List of path to the alignment file(s) used as input.
    contig_data : dict
        Dictionnary of the all the contigs from the assembly, the contigs names
        are the keys to the data of the contig available with the following
        keys: "id", "length", "GC", "hit", "coverage". Coverage still at 0 and
        need to be updated later.
    hit_data : dict:
        Dictionnary for hit information on each contigs.
    out_dir : str
        The output directory to write the network and chunk data into.
    output_file_network : str, optional
        The specific file name for the output network file. Default is
        'network.txt'
    output_file_contig_data : str, optional
        The specific file name for the output chunk data file. Default is
        'idx_contig_length_GC_hit_cov.txt'
    tmp_dir : str
        Path to th temporary directory. Default in the working directory
    normalization : str
        If None, do not normalized the count of a contact by the geometric mean
        of the coverage of the contigs. Otherwise it's the type of
        normalization.
    self_contacts : bool
        Whether to return network with self contact. Default is False.

    Returns:
    --------
    str:
        Path to the network file.
    str:
        Path to the verbose contig data file.
    """

    # Create temporary and output file which will be necessary
    precompute_network_file = join(tmp_dir, "precompute_network_file.txt")
    pre_network_sorted_file = join(tmp_dir, "tmp_network_sorted.txt")
    network_file = join(out_dir, output_file_network)
    contig_data_file = join(out_dir, output_file_contig_data)
    hit_data_file = join(out_dir, "hit_data_alignment.txt")
    nb_alignment = len(alignment_files)
    logger.info("New time course network")

    # Create a contact file easily readable for counting the contacts.
    contig_data, out_files_list = precompute_network(
        alignment_files,
        contig_data,
        hit_data,
        precompute_network_file,
        tmp_dir,
        self_contacts,
    )

    # Compute network
    compute_network(
        precompute_network_file,
        network_file,
        contig_data,
        tmp_dir,
        pre_network_sorted_file,
        n_cpus,
        normalization,
    )

    # Compute sample network
    for i, precompute_network_file_sample in enumerate(out_files_list):
        network_file_sample = join(out_dir, "network_{0}.txt".format(i))
        pre_network_sorted_file = join(
            tmp_dir, "tmp_network_sorted_{0}.txt".format(i)
        )
        compute_network(
            precompute_network_file_sample,
            network_file_sample,
            contig_data,
            tmp_dir,
            pre_network_sorted_file,
            n_cpus,
            normalization,
        )

    # Write the data from the contigs
    write_contig_data(contig_data, contig_data_file)
    if nb_alignment > 1:
        write_hit_data(hit_data, hit_data_file)

    return network_file, contig_data_file


def compute_network(
    pre_network_file,
    network_file,
    contig_data,
    tmp_dir,
    tmp_file,
    n_cpus,
    normalization,
):
    """Generate the network file from the prenetwork file.

    Parameters:
    -----------
    pre_network_file  : str
        Path of the input file (prenetwork file, output from teh precompute
        network function)
    network_file : str
        Path of the output file (network file).
    contig_data : dict
        Dictionnary of the all the contigs from the assembly, the contigs names
        are the keys to the data of the contig available with the following
        keys: "id", "length", "GC", "hit", "coverage", "RS".
    tmp_dir : str
        Path to the temporary directory.
    tmp_file : str
        Path of the temporary file to write the sorted precompute network.
    n_cpus : int
        Number of cpus used to sort the prenetwork.
    normalization : str
        If None, do not normalized the count of a contact by the geometric mean
        of the coverage of the contigs. Otherwise it's the type of
        normalization.
    """

    # Sort the the pre-network file
    mio.sort_pairs(
        pre_network_file,
        tmp_file,
        tmp_dir=tmp_dir,
        threads=n_cpus,
    )

    # Set the variables used in the loop
    n_occ = 1  # Number of occurences of each frag combination
    n_nonzero = 1  # Total number of nonzero matrix entries
    n_pairs = 0  # Total number of pairs entered in the matrix

    # Read the sorted paors
    with open(tmp_file, "r") as pairs, open(network_file, "w") as net:
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
                    if normalization != "None":
                        effective_count = normalize_pair(
                            contig_data, prev_pair, n_occ, normalization
                        )
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
        if normalization != "None":
            effective_count = normalize_pair(
                contig_data, prev_pair, n_occ, normalization
            )
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


def create_contig_data(assembly, nb_alignment, depth_file, enzyme):
    """Create a dictionnary with data on each Contig.

    Parameters:
    -----------
    assembly : str
        Path to the assembly fasta file
    nb_alignement : int
        Numbers of alignment files.
    depth_file : str or None
        Path to the depth.txt file from jgi_summarize_bam_contig_depths from
        Metabat2 Software.
    enzyme : str or None
        String that contains the names of the enzyme separated by a comma.

    Returns:
    --------
    dict:
        Dictionnary of the all the contigs from the assembly, the contigs names
        are the keys to the data of the contig available with the following
        keys: "id", "length", "GC", "hit", "coverage", "RS". Hit is set to 0 and
        need to be updated later.
    dict:
        Dictionnary for hit information on each contigs.
    """

    # Create a contig data dictionnary from the assembly file
    contig_data = dict()
    if nb_alignment > 1:
        hit_data = dict()
    else:
        hit_data = None

    # Extract restriction sites if an enzyme is given.
    if enzyme:
        pattern = mio.get_restriction_site(enzyme)

    global_id = 1
    if depth_file:
        with open(depth_file, "r") as depth:
            line = depth.readline()
            for contig in SeqIO.parse(assembly, "fasta"):
                line = depth.readline().split("\t")
                contig_data[contig.id] = {
                    "id": global_id,
                    "length": int(line[1]),
                    "GC": SeqUtils.GC(contig.seq),
                    "hit": 0,
                    "coverage": float(line[2]),
                    "RS": (len(re.findall(pattern, str(contig.seq))) + 1)
                    if enzyme
                    else "-",
                }
                if nb_alignment > 1:
                    hit_data[contig.id] = {
                        "id": global_id,
                        "hit": [0] * nb_alignment,
                    }
                global_id += 1
    else:
        for contig in SeqIO.parse(assembly, "fasta"):
            contig_data[contig.id] = {
                "id": global_id,
                "length": len(contig.seq),
                "GC": SeqUtils.GC(contig.seq),
                "hit": 0,
                "coverage": "-",
                "RS": (len(re.findall(pattern, str(contig.seq))) + 1)
                if enzyme
                else "-",
            }
            if nb_alignment > 1:
                hit_data[contig.id] = {
                    "id": global_id,
                    "hit": [0] * nb_alignment,
                }
            global_id += 1

    return contig_data, hit_data


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
        Contac count between the two contigs.
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
        factor = np.sqrt(
            contig_data[pair[0]]["coverage"] * contig_data[pair[1]]["coverage"]
        )
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
            contig_data[pair[0]]["hit"]
            / contig_data[pair[0]]["RS"]
            * contig_data[pair[1]]["hit"]
            / contig_data[pair[1]]["RS"]
        )

    # The last two normalizations are normalization by geometric means of hits.
    elif normalization == "empirical_hit":
        factor = np.sqrt(
            contig_data[pair[0]]["hit"] * contig_data[pair[1]]["hit"]
        )

    elif normalization == "theoritical_hit":
        factor = np.sqrt(
            contig_data[pair[0]]["coverage"]
            * 10
            * np.sqrt(
                contig_data[pair[0]]["length"] * contig_data[pair[0]]["RS"]
            )
            * contig_data[pair[1]]["coverage"]
            * 10
            * np.sqrt(
                contig_data[pair[1]]["length"] * contig_data[pair[1]]["RS"]
            )
        )
        # If no read mapped from Shotgun libraries could be equal to zero.
        if factor == 0:
            return 0

    return n_occ / factor


def precompute_network(
    alignment_files,
    contig_data,
    hit_data,
    out_file,
    tmp_dir,
    self_contacts=False,
):
    """Write a file with only the contig id separated by a tabulation and count
    the contacts by contigs to be able to compute directlty the normalized
    network.

    Parameters:
    -----------
    alignment_files : list of str
        List of path to the alignment file(s).
    contig_data : dict
        Dictionnary of the all the contigs from the assembly, the contigs names
        are the keys to the data of the contig available with the following
        keys: "id", "length", "GC", "hit", "coverage". Coverage still at 0 and
        need to be updated later.
    hit_data : dict
        Dictionnary with the count of hits for each aligment file.
    out_file : str
        Path to the write the output_file which will be necessary to compute the
        network.
    self_contacts : bool
        If True, the contacts on the same contigs will be kept. Otherwise only
        displays the inter contigs contacts. [Default False]

    Return:
    -------
    dict:
        Dictionnary of the all the contigs from the assembly, the contigs names
        are the keys to the data of the contig available with the following
        keys: "id", "length", "GC", "hit", "coverage" "RS". Coverage still at 0
        and need to be updated later.
    """
    # Initiate value to compute 3D ratio
    all_contacts = 0
    inter_contacts = 0
    out_files_list = []

    # Prepare a file to save contact with their global ID
    with open(out_file, "w") as pre_net:

        # Iterates on the alignment files
        for i, aligment_file in enumerate(alignment_files):

            all_contacts_temp = 0
            inter_contacts_temp = 0
            out_file_sample = join(tmp_dir, "prenetwork" + str(i) + ".txt")
            out_files_list.append(out_file_sample)

            # Read the alignment_file and build pairs for the network
            with open(aligment_file, "r") as pairs, open(
                out_file_sample, "w"
            ) as pre_net_sample:
                for pair in pairs:
                    # Ignore header lines
                    if pair.startswith("#"):
                        continue

                    # Split the line on the tabulation
                    p = pair.split("\t")

                    # Extract the contig names which are at the position 2 and
                    # 4.
                    contig1, contig2 = p[1], p[3]
                    id1 = contig_data[contig1]["id"]
                    id2 = contig_data[contig2]["id"]

                    # Count the contact
                    all_contacts_temp += 1
                    contig_data[contig1]["hit"] += 1
                    contig_data[contig2]["hit"] += 1
                    if len(alignment_files) > 1:
                        hit_data[contig1]["hit"][i] += 1
                        hit_data[contig2]["hit"][i] += 1

                    # Write the file used for the computation of the network.
                    if self_contacts and id1 == id2:
                        pre_net.write(
                            "\t".join(map(str, [contig1, contig2])) + "\n"
                        )
                        pre_net_sample.write(
                            "\t".join(map(str, [contig1, contig2])) + "\n"
                        )
                    elif id1 < id2:
                        inter_contacts_temp += 1
                        pre_net.write(
                            "\t".join(map(str, [contig1, contig2])) + "\n"
                        )
                        pre_net_sample.write(
                            "\t".join(map(str, [contig1, contig2])) + "\n"
                        )
                    elif id1 > id2:
                        inter_contacts_temp += 1
                        pre_net.write(
                            "\t".join(map(str, [contig2, contig1])) + "\n"
                        )
                        pre_net_sample.write(
                            "\t".join(map(str, [contig2, contig1])) + "\n"
                        )

            # Count contacts and return sample informations.
            all_contacts += all_contacts_temp
            inter_contacts += inter_contacts_temp
            logger.info("Information of {0}:".format(basename(aligment_file)))
            logger.info(
                "{0} contacts in the library.".format(all_contacts_temp)
            )
            logger.info(
                "{0} contacts inter-contigs in the library.".format(
                    inter_contacts_temp
                )
            )
            logger.info(
                "3D ratio : {0}\n".format(
                    inter_contacts_temp / all_contacts_temp
                )
            )

    # Return information about the network
    if len(alignment_files) > 1:
        logger.info("General information:")
        logger.info("{0} contacts in the library.".format(all_contacts))
        logger.info(
            "{0} contacts inter-contigs in the library.".format(inter_contacts)
        )
        logger.info("3D ratio : {0}\n".format(inter_contacts / all_contacts))

    return contig_data, out_files_list


def write_contig_data(contig_data, output_path):
    """Function to write the contig data file at the output path given. The file
    will contains 7 columns separated by a tabulation: id, name, length,
    GC_content, hit, coverage, restriction site for each contig with an header.

    Parameters:
    -----------
    contig_data : dict
        Dictionnary of the all the contigs from the assembly, the contigs names
        are the keys to the data of the contig available with the following
        keys: "id", "length", "GC", "hit", "coverage", "RS".
    output_path : str
        Path to the output file where the data from the dictionnary will be
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
            line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                idx, name, length, GC_content, hit, coverage, RS
            )
            contig_data_file_handle.write(line)


def write_hit_data(hit_data, output_path):
    """Function to write the contig data file at the output path given. The file
    will contains 6 columns separated by a tabulation: id, name, length,
    GC_content, hit, coverage for each contig.

    Parameters:
    -----------
    hit_data : dict
        Dictionnary of the all the contigs from the assembly, the contigs names
        are the keys to a list of hits from each alignment files separately.
    output_path : str
        Path to the output file where the data from the dictionnary will be
        written.
    """

    # For each contig extract the data and write them in the file.
    with open(output_path, "w") as hit_data_file_handle:
        for name in hit_data:
            hit_list = map(str, hit_data[name]["hit"])
            hit_str = "\t".join(hit_list)
            idx = hit_data[name]["id"]
            line = "{0}\t{1}\t{2}\n".format(idx, name, hit_str)
            hit_data_file_handle.write(line)
