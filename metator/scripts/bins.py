#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Process metagenomic bins and their networks

Extract bin data and figures from partitions, namely bin subnetworks
(notably for recursion purposes) and bin fasta files. Also draws the subnetwork
in adjacency matrix form.

"""

from __future__ import absolute_import

import numpy as np
from scipy import sparse
from scipy import stats
from Bio import SeqIO
import os
import argparse
import itertools
import operator
from metator.scripts import hicstuff as hcs

from metator.scripts.log import logger

from metator.scripts.figures import spaceless_pdf_plot_maker

from metator.scripts.network import DEFAULT_CHUNK_SIZE
from metator.scripts.figures import DEFAULT_MAX_SIZE_MATRIX
from metator.scripts.figures import DEFAULT_SATURATION_THRESHOLD

DEFAULT_MAX_CORES = 100
DEFAULT_MAX_CONTAMINATION = 50
DEFAULT_MIN_COMPLETENESS = 80


def extract_subnetworks(
    partition_file,
    network_file,
    output_dir,
    max_cores=DEFAULT_MAX_CORES,
    max_size_matrix=DEFAULT_MAX_SIZE_MATRIX,
    saturation_threshold=DEFAULT_SATURATION_THRESHOLD,
):
    """Extract bin subnetworks from the main network

    Identify bins, extract subnets, draws the adjacency matrices,
    saves it all in a specified output directory.

    Parameters
    ----------
    partition_file : file, str or pathlib.Path
        The file containing, for each chunk, the communities it was
        assigned to at each iteration.
    network_file : file, str or pathlib.Path
        The file containing the network in sparse (edge list) format
    output_dir : str or pathlib.Path
        The output directory to write the subnetworks into.
    max_cores : int, optional
        The maximum number of bins to extract. Default is 100.
    max_size_matrix : int, optional
        When rendering contact maps for each bin, the maximum size for the
        matrix. Default is 2000.
    saturation_threshold : float, optional
        When rendering contact maps for each bin, the percentile value over
        which the color map should be saturated. Default is 80.
    """

    logger.info("Loading partition...")
    data_chunks = np.loadtxt(partition_file, usecols=(1,), dtype=np.int32)

    logger.info("Loading network...")
    network = np.loadtxt(network_file, dtype=np.int32)
    cores = data_chunks

    core_network = np.copy(network)

    core_network[:, 0] = cores[network[:, 0]]
    core_network[:, 1] = cores[network[:, 1]]

    n = np.amax(cores) + 1

    def extract(network_to_keep, filename):

        subnetwork = np.copy(network[network_to_keep])
        subnetwork[:, 0] -= 1
        subnetwork[:, 1] -= 1

        np.savetxt(filename, subnetwork, fmt="%i")

        return subnetwork

    def draw(subnetwork, filename):

        try:
            # Numpy array format
            row = subnetwork[:, 0]
            col = subnetwork[:, 1]
            data = subnetwork[:, 2]
        except TypeError:
            # Scipy sparse format
            row = subnetwork.row
            col = subnetwork.col
            data = subnetwork.data

        row_indices = stats.rankdata(
            np.concatenate((row, col)), method="dense"
        )
        col_indices = stats.rankdata(
            np.concatenate((col, row)), method="dense"
        )
        data = np.concatenate((data, data))

        # print("Row length: {}, col length: {}, data length: {}"
        #       "".format(len(row_indices), len(col_indices), len(data)))

        unique_row = np.unique(row)
        unique_col = np.unique(col)

        # print("Network shape: {},{}".format(len(unique_row),
        #                                     len(unique_col)))

        size = len(np.unique(np.concatenate((unique_row, unique_col)))) + 1
        # print("Size of matrix to draw: {}".format(size))

        try:
            sparse_subnet = sparse.coo_matrix(
                (data, (row_indices, col_indices)), shape=(size, size)
            )
            binning_factor = (size // max_size_matrix) + 1
            binned_subnet = hcs.bin_sparse(
                sparse_subnet, subsampling_factor=binning_factor
            )
            dense_subnet = binned_subnet.todense()

            diagonal = np.diag(np.diag(dense_subnet))
            normed_subnet = hcs.normalize_dense(dense_subnet - diagonal)

            vmax = np.percentile(normed_subnet, saturation_threshold)

            spaceless_pdf_plot_maker(normed_subnet, filename, vmax=vmax)

        except MemoryError:
            logger.warning(
                "Warning, couldn't save matrix due to memory issues"
            )

    def extract_and_draw(network_to_keep, filename_text, filename_image):

        subnetwork = extract(network_to_keep, filename=filename_text)
        draw(subnetwork, filename=filename_image)

    #   Extract and draw subnetworks for chosen cores and draw 2D arrays
    global_network_indices_list = []

    for i in range(1, n):

        if i > max_cores:
            break

        # print("Bin {}:".format(i))
        network_to_keep_1 = core_network[:, 0] == i
        network_to_keep_2 = core_network[:, 1] == i

        network_to_keep = network_to_keep_1 * network_to_keep_2

        nonzero_indices, = np.nonzero(network_to_keep)
        global_network_indices_list += nonzero_indices.tolist()

        subnetwork_file = os.path.join(
            output_dir, "subnetwork_core_{}.dat".format(i)
        )

        image_name = os.path.join(output_dir, "core_{}.eps".format(i))

        extract_and_draw(
            network_to_keep=network_to_keep,
            filename_text=subnetwork_file,
            filename_image=image_name,
        )


def extract_fasta(
    partition_file,
    fasta_file,
    output_dir,
    chunk_size=DEFAULT_CHUNK_SIZE,
    max_cores=DEFAULT_MAX_CORES,
):

    """Extract sequences from bins

    Identify bins, extract chunks belonging to each bins and gather them
    in a single FASTA file.

    Parameters
    ----------
    partition_file : file, str or pathlib.Path
        The file containing, for each chunk, the communities it was
        assigned to at each iteration.
    fasta_file : file, str or pathlib.Path
        The initial assembly from which chunks were initialized.
    output_dir : str or pathlib.Path
        The output directory to write the FASTA chunks into.
    chunk_size : int, optional
        The size of the chunks (in bp) used in the pipeline. Default is 1000.
    max_cores : int, optional
        How many bins to extract FASTA sequences from. Default is 100.
    """

    genome = {
        record.id: record.seq for record in SeqIO.parse(fasta_file, "fasta")
    }

    data_chunks = list(
        zip(*np.genfromtxt(partition_file, usecols=(0, 1), dtype=None))
    )

    chunk_names = np.array(data_chunks[0], dtype=object)
    cores = np.array(data_chunks[1])
    for core in set(cores):
        if core > max_cores:
            continue
        chunks_to_keep = chunk_names[cores == core]
        core_name = "core_{}.fa".format(core)

        core_file = os.path.join(output_dir, core_name)

        with open(core_file, "w") as core_handle:
            for name in chunks_to_keep:
                fields = name.split("_")
                header_name = "_".join(fields[:-1])
                chunk = int(fields[-1])

                pos_start = chunk * chunk_size
                pos_end = min(
                    (chunk + 1) * chunk_size, len(genome[header_name])
                )

                sequence = str(genome[header_name][pos_start:pos_end])

                core_handle.write(">{}\n".format(name))
                core_handle.write("{}\n".format(sequence))


def merge_fasta(fasta_file, output_dir):

    """Merge chunks into complete FASTA bins

    Merge bin chunks by appending consecutive chunks to one another.

    Parameters
    ----------
    fasta_file : file, str or pathlib.Path
        The FASTA file containing the chunks to merge.
    output_dir : str or pathlib.Path
        The output directory to write the merged FASTA bin into.
    """

    #   First, define some functions for ordering chunks and detecting
    #   consecutive chunk sequences

    def chunk_lexicographic_order(chunk):
        """A quick callback to sort chunk ids lexicographically
        (first on original names alphabetically, then on relative
        position on the original contig)
        """
        chunk_fields = chunk.split("_")
        chunk_name = chunk_fields[:-1]
        chunk_id = chunk_fields[-1]
        return (chunk_name, int(chunk_id))

    def are_consecutive(chunk1, chunk2):
        if None in {chunk1, chunk2}:
            return False
        else:
            ord1 = chunk_lexicographic_order(chunk1)
            ord2 = chunk_lexicographic_order(chunk2)
            return (ord1[0] == ord2[0]) and (ord1[1] == ord2[1] + 1)

    def consecutiveness(key_chunk_pair):
        """A callback for the groupby magic below
        """
        key, chunk = key_chunk_pair
        chunk_name, chunk_id = chunk_lexicographic_order(chunk)
        return (chunk_name, chunk_id - key)

    #   Read chunks and sort them
    genome = {
        record.id: record.seq for record in SeqIO.parse(fasta_file, "fasta")
    }

    sorted_ids = sorted(genome, key=chunk_lexicographic_order)

    #   Identify consecutive ranges and merge them
    new_genome = dict()
    for _, g in itertools.groupby(enumerate(sorted_ids), consecutiveness):
        chunk_range = map(operator.itemgetter(1), g)
        first_chunk = next(chunk_range)
        my_sequence = genome[first_chunk]
        my_chunk = None
        while "Reading chunk range":
            try:
                my_chunk = next(chunk_range)
                my_sequence += genome[my_chunk]
            except StopIteration:
                break

        try:
            last_chunk_id = my_chunk.split("_")[-1]
        except AttributeError:
            last_chunk_id = ""

        if last_chunk_id:
            new_chunk_id = "{}_{}".format(first_chunk, last_chunk_id)
        else:
            new_chunk_id = first_chunk
        new_genome[new_chunk_id] = my_sequence

    #   Write the result
    base_name = ".".join(os.path.basename(fasta_file).split(".")[:-1])
    output_name = "{}_merged.fa".format(base_name)
    merged_core_file = os.path.join(output_dir, output_name)
    with open(merged_core_file, "w") as output_handle:
        for my_id in sorted(new_genome, key=chunk_lexicographic_order):
            output_handle.write(">{}\n".format(my_id))
            output_handle.write("{}\n".format(new_genome[my_id]))


def pick_contaminated_bins(
    checkm_file,
    fasta_dir,
    output_dir,
    min_completeness=DEFAULT_MIN_COMPLETENESS,
    max_contamination=DEFAULT_MAX_CONTAMINATION,
    verbose=True,
):
    contaminated_and_complete_bins = set()

    with open(checkm_file) as checkm_handle:
        first_line = checkm_handle.readline()
        second_line = checkm_handle.readline()
        third_line = checkm_handle.readline()
        assert all(
            first_line.startswith("----"),
            "Completeness" in second_line,
            third_line.startswith("----"),
        ), "{} does not look like a proper checkM output file.".format(
            checkm_file
        )
        for line in checkm_handle:
            bin_name, *_, completeness, contamination, _ = line.split()

            if int(min_completeness) > completeness:
                break

            if contamination < max_contamination:
                continue

            contaminated_and_complete_bins.add(bin_name)

    if verbose:
        for my_bin in contaminated_and_complete_bins:
            print(my_bin)

    return contaminated_and_complete_bins


def main():

    parser = argparse.ArgumentParser(
        description="Extract bin data," " subnets and figures"
    )

    parser.add_argument(
        "-f", "--fasta", help="Reference FASTA file for extracting bins"
    )

    parser.add_argument(
        "-i", "--input", help="Input partition file for bin extraction"
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output directory for bin data and figures",
    )

    parser.add_argument(
        "-n", "--network", help="Reference network file for extracting subnets"
    )

    parser.add_argument(
        "-c",
        "--chunk-size",
        type=int,
        help="Chunk size for bin reconstruction",
        default=DEFAULT_CHUNK_SIZE,
    )

    parser.add_argument(
        "-m", "--merge", help="Merge consecutive chunks in a FASTA file"
    )

    parser.add_argument(
        "-N",
        "--n-bins",
        help="Number of bins to extract in size order",
        type=int,
        default=DEFAULT_MAX_CORES,
    )

    args = parser.parse_args()

    fasta_file = args.fasta
    partition_file = args.input
    output_dir = args.output
    network_file = args.network
    chunk_size = args.chunk_size
    fasta_to_merge = args.merge
    max_cores = args.n_bins

    if network_file is not None:
        extract_subnetworks(
            partition_file=partition_file,
            network_file=network_file,
            max_cores=max_cores,
            output_dir=output_dir,
        )

    elif fasta_file is not None:
        extract_fasta(
            partition_file=partition_file,
            fasta_file=fasta_file,
            output_dir=output_dir,
            max_cores=max_cores,
            chunk_size=chunk_size,
        )

    elif fasta_to_merge is not None:
        merge_fasta(fasta_file=fasta_to_merge, output_dir=output_dir)


if __name__ == "__main__":

    main()
