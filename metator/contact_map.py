#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Generates a HiC contact map of a MetaTOR object.

General utility functions and class for handling one MetaTOR object and for
generating its HiC contact map using hicstuff pipeline.

Core class to handle bin object easily:
    - MetatorObject

Core function to build the network are:
    - extract_pairs
    - generate_contact_map
"""


import hicstuff.pipeline as hcp
import os
import pandas as pd
import subprocess as sp
from metator.log import logger
from os.path import join


class MetatorObject:
    """
    Class to handle MetaTOR object informations and to extract easily aligned
    pairs from the object given.

    The class is initiate by the MetaTOR object type, its name, the assembly
    file, the contig data file from MetaTOR and the pairs files.

    It's possible to add a threshold on the contigs size.
    """

    def __init__(
        self,
        metator_object,
        name,
        assembly,
        contig_data_file,
        pairs,
        min_size=5000,
    ):
        """Initiates the MetaTOR objects and set useful informations.

        Parameters:
        -----------
        metator_object : str
            Object to extract contigs to build the matrix. Either "contig",
            "core_bin", "overlapping_bin", "recursive_bin", "final_bin" or
            "other".
        name : str
            Name of the object. Could be the name of a contig, an id of a bin or
            the name of the bin. Example: "NODE_1" or "MetaTOR_1_0".
        assembly : str
            Path to the fasta file containing the contigs of interest. Could be
            the whole or the extracted contigs of one bin.
        contig_data_file : str
            Path to the contig_data_final.txt file form MetaTOR output.
        pairs : str
            Path of the ".pairs" file. If more than one is given, files should
            be separated by a comma.
        min_size : int
            Size threshold used to filter contigs. Default: 5000.
        """
        self.assembly = assembly
        self.contigs_data = contig_data_file
        self.pairs_files = pairs.split(",")
        self.min_size = min_size
        self.set_metator_object(metator_object, name)

    def get_contigs_data(self):
        """Method to get contigs data Dataframe

        Returns:
        --------
        pandas.core.frame.DataFrame
            Table with informations from metaTOR binning for each contigs.
        """
        contigs_data = pd.read_csv(self.contigs_data, sep="\t")
        return contigs_data

    def set_contigs(self):
        """Method to extract the list of contigs of the class from the contig
        data file and their size.
        """
        contigs_data = self.get_contigs_data()
        self.contigs = list(
            contigs_data["Name"][contigs_data[self.object] == self.name]
        )
        self.contigs_size = list(
            contigs_data["Size"][contigs_data[self.object] == self.name]
        )

    def set_large_contigs(self):
        """Method to keep only contigs bigger than the threshold given to remove
        small bins which will be smaller than the binning value and which will
        prevent a good scaffolding with Instagraal.
        """
        self.contigs_size = [x for x in self.contigs_size if x >= self.min_size]
        self.contigs = self.contigs[: len(self.contigs_size)]

    def set_metator_object(self, metator_object, name):
        """Method to get the metator object and name of the object usable for
        the algorithm.

        Parameters:
        -----------
        metator_object : str Object to extract contigs to build the matrix.
            Either "contig", "core_bin", "overlapping_bin", "recursive_bin",
            "final_bin" or "other".
        name : str
            Name of the object. Could be the name of a contig, an id of a bin or
            the name of the bin. Example: "NODE_1" or "MetaTOR_1_0".
        """
        if metator_object == "contig":
            self.object = "Name"
            self.name = name
        elif metator_object == "core_bin":
            self.object = "Core_bin_ID"
            try:
                int(name)
                self.name = str(name)
            except ValueError:
                logger.error(
                    "With core bin object, the name should be the numeric ID of the core bin."
                )
                raise ValueError
        elif metator_object == "overlapping_bin":
            self.object = "Overlapping_bin_ID"
            try:
                int(name)
                self.name = str(name)
            except ValueError:
                self.name = str(name.split("_")[1])
                try:
                    int(self.name)
                except ValueError:
                    logger.error(
                        "Overlapping bin name should be either numerical ID or a name like 'MetaTOR_1' or 'MetaTOR_1_0'."
                    )
                    raise ValueError
        elif metator_object == "recursive_bin":
            self.object = "Recursive_bin_ID"
            try:
                int(name)
                self.name = str(name)
            except ValueError:
                self.name = str(name.split("_")[2])
                try:
                    int(self.name)
                except ValueError:
                    logger.error(
                        "Overlapping bin name should be either numerical ID or a name like 'MetaTOR_1_1'."
                    )
                    raise ValueError
                if int(self.name) <= 0:
                    logger.error(
                        "A recursive bin should have an id bigger than 0."
                    )
        elif metator_object == "final_bin":
            self.object = "Final_bin"
            self.name = name
        elif metator_object == "other":
            self.object = "Other"
            self.name = name
        else:
            logger.error(
                'Metator object should be one of these value: "contig", "core_bin", "overlapping_bin", "recursive_bin", "final_bin", "other"'
            )
            raise ValueError

    def write_fasta(self, tmp_dir, out_dir):
        """Method to write the new fasta with only the contigs of interest in
        the outdir directory.

        Parameters:
        -----------
        tmp_dir : str
            Path to the temporary directory to write the list of contigs of the
            object.
        out_dir : str
            Path to the output directory where to write the new fasta file.
        """
        # Create a fasta ouput file.
        self.fasta = join(out_dir, self.name + ".fa")
        # Write a contigs used by pyfastx extract.
        contigs_list = join(tmp_dir, self.name + ".txt")
        with open(contigs_list, "w") as f:
            for contig_name in self.contigs:
                f.write("%s\n" % contig_name)
        # Extract contigs from the fastq.
        cmd = "pyfastx extract {0} -l {1} > {2}".format(
            self.assembly, contigs_list, self.fasta
        )
        process = sp.Popen(cmd, shell=True)
        process.communicate()


def extract_pairs(metator_data):
    """Function to extract the pairs of a MetaTOR object from pairs files.

    Parameters:
    -----------
    metator_data : object from class MetatorObject
        Object of the class MetatorObject. It should have the parameters contigs
        and the output pairs file of the metaTOR object set.

    Returns:
    --------
    int:
        Number of pairs from the MetaTOR object.
    """
    # Initiation
    n_pairs = 0
    output_file = metator_data.pairs

    # Write one pair file for all the ones given.
    with open(output_file, "w") as output_pairs:
        # Write the header of the output pairs
        output_pairs.write("## pairs format v1.0\n")
        output_pairs.write(
            "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n"
        )
        for i in range(len(metator_data.contigs)):
            output_pairs.write(
                "#chromsize: {0} {1}\n".format(
                    metator_data.contigs[i],
                    metator_data.contigs_size[i],
                )
            )
        for pair_file in metator_data.pairs_files:
            # Iterates on the input pairs file
            with open(pair_file, "r") as input_pairs:
                for pair in input_pairs:
                    # Ignore header lines.
                    if pair.startswith("#"):
                        continue
                    # Split the line on the tabulation and check if both contigs
                    # are in the bin.
                    p = pair.split("\t")
                    if (
                        p[1] in metator_data.contigs
                        and p[3] in metator_data.contigs
                    ):
                        n_pairs += 1
                        output_pairs.write(pair)
    return n_pairs


def generate_contact_map(
    assembly,
    contig_data_file,
    enzyme,
    name,
    pairs,
    out_dir,
    tmp_dir,
    filter_events=False,
    force=False,
    mat_fmt="graal",
    metator_object="final_bin",
    min_size=5000,
    no_cleanup=False,
    pcr_duplicates=False,
    threads=1,
):
    """General function to extract pairs of the MetaTOR object and generate its
    the contact map.

    Parameters:
    -----------
    assembly : str
        Path to the fasta file containing the contigs of interest. Could be the
        whole or the extracted contigs of one bin.
    contig_data_file : str
        Path to the contig_data_final.txt file form MetaTOR output.
    enzyme : str
        Enzyme used to digest the genome in the HiC experiment. Example:
        HpaII,MluCI.
    name : str
        Name of the object. Could be the name of a contig, an id of a bin or the
        name of the bin. Example: "NODE_1" or "MetaTOR_1_0".
    pairs : str
        Path of the ".pairs" file. If more than one is given, files should be
        separated by a comma.
    out_dir : str
        Path where output files should be written. Current directory by default.
    tmp_dir : str
        Path where temporary files will be written.
    filter_events : bool
        Filter spurious or uninformative 3C events. Requires a restriction
        enzyme. Default: False.
    force : bool
        If True, overwrite existing files with the same name as output.
        Default: False.
    mat_fmt : str
        Select the output matrix format. Can be either "bg2" for the bedgraph2
        format, "cool" for Mirnylab's cool format, or graal for a plain text COO
        format compatible with Koszullab's instagraal software.
        Default: "graal".
    metator_object : str
        Object to extract contigs to build the matrix. Either "contig",
        "core_bin", "overlapping_bin", "recursive_bin", "final_bin" or "other".
    min_size : int
        Minimum contig size required to keep it.
    no_cleanup : bool
        Whether temporary files should be deleted at the end of the pipeline.
        Default: False.
    pcr_duplicates : bool
        If True, PCR duplicates will be filtered based on genomic positions.
        Pairs where both reads have exactly the same coordinates are considered
        duplicates and only one of those will be conserved. Default: False.
    threads : int
        Numbers of threads to use. Default: 1.
    """

    # Extract bin information from metaTOR outdir.
    logger.info("Generate HiC contact map for {0}".format(name))
    metator_data = MetatorObject(
        metator_object, name, assembly, contig_data_file, pairs, min_size
    )
    metator_data.set_contigs()
    if min_size > 0:
        metator_data.set_large_contigs()
    metator_data.write_fasta(tmp_dir, out_dir)
    metator_data.pairs = join(tmp_dir, name + ".pairs")

    # Extract pairs of the bin.
    n_pairs = extract_pairs(metator_data)
    logger.info("{0} pairs have been extracted.".format(n_pairs))

    # Launch hicstuff pipeline.
    hcp.full_pipeline(
        genome=metator_data.fasta,
        input1=metator_data.pairs,
        distance_law=False,
        enzyme=enzyme,
        filter_events=filter_events,
        force=force,
        out_dir=out_dir,
        pcr_duplicates=pcr_duplicates,
        plot=False,
        prefix=metator_data.name,
        start_stage="pairs",
        threads=threads,
        tmp_dir=tmp_dir,
    )
