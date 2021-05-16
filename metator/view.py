#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Generate HiC contact map of one bin from metaTOR outdir.

General utility functions and class for handling one bin from metaTOR output and
generating HiC contact map using hicstuff.

Core class to handle bin object easily:
    - Bin

Core function to build the network are:
    - extract_pairs_bins
    - generate_contact_map_bin
"""


import fnmatch
import hicstuff.pipeline as hcp
import os
import pandas as pd
from metator.log import logger
from os.path import join


class Bin:
    """
    Class to handle bin informations to extract easily aligned pairs from the
    bin.

    The class is initiate by the name of the bin and the "project", i.e. the
    path of the output directory from metator validation or metator pipeline
    with the validation step done.
    """

    def __init__(self, name, project, min_size=5000):
        """Initiate the bin objects and set the path to useful files related to
        the bin.

        Parameters:
        -----------
        name : str
            Name of the bin. Example: "MetaTOR_1_0".
        project :
            Path of the output directory from metator validation or metator
            pipeline with the validation step done.
        min_size : int
            Size threshold used to filter contigs. Default: 5000.
        """
        self.name = name
        self.project = project
        self.genome = join(project, "final_bin", name + ".fa")
        self.contigs_data = join(project, "contig_data_final.txt")
        self.min_size = min_size

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
            contigs_data["Name"][contigs_data["Final_bin"] == self.name]
        )
        self.contigs_size = list(
            contigs_data["Size"][contigs_data["Final_bin"] == self.name]
        )

    def set_large_contigs(self):
        """Method to keep only contigs bigger than the threshold given to remove
        small bins which will be smaller than the binning value and which will
        prevent a good scaffolding with Instagraal.
        """
        self.large_contigs_size = [
            x for x in self.contigs_size if x >= self.min_size
        ]
        self.large_contigs = self.contigs[: len(self.large_contigs_size)]

    def set_pairs_files_list(self):
        """Method to set the pairs file list from the output of MetaTOR."""
        files = os.listdir(self.project)
        pattern = "*.pairs"
        self.full_pairs = []
        for file in files:
            if fnmatch.fnmatch(file, pattern):
                self.full_pairs.append(file)


def extract_pairs_bin(bin_data):
    """Function to extract the pairs from one bin from a pair file.

    Parameters:
    -----------
    bin_data : object from class bin
        Object of the class Bin. It should have the parameters large_contigs set
        and the pairs file of the bin defined.

    Returns:
    --------
    int:
        Number of pairs from the bin.
    """
    # Initiation
    bin_data.set_pairs_files_list()
    n_pairs = 0
    output_file = bin_data.pairs

    # Write one pair file for all the ones given.
    with open(output_file, "w") as output_pairs:
        # Write the header of the output pairs
        output_pairs.write("## pairs format v1.0\n")
        output_pairs.write(
            "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n"
        )
        for i in range(len(bin_data.large_contigs)):
            output_pairs.write(
                "#chromsize: {0} {1}\n".format(
                    bin_data.large_contigs[i],
                    bin_data.large_contigs_size[i],
                )
            )
        for pair_file in bin_data.full_pairs:
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
                        p[1] in bin_data.large_contigs
                        and p[3] in bin_data.large_contigs
                    ):
                        n_pairs += 1
                        output_pairs.write(pair)
    return n_pairs


def generate_contact_map_bin(
    bin_name,
    project,
    out_dir,
    tmp_dir,
    enzyme,
    filter_events=False,
    force=False,
    mat_fmt="graal",
    min_size=5000,
    no_cleanup=False,
    pcr_duplicates=False,
    threads=1,
):
    """General function to extract pairs of the bin and generate the contact map
    of the bin.

    Parameters:
    -----------
    bin_name : str
        Name of the bin. Example: "MetaTOR_1_0".
    project : str
        Path of the output directory from metator validation or metator
        pipeline with the validation step done.
    out_dir : str
        Path where output files should be written. Current directory by default.
    tmp_dir : str
        Path where temporary files will be written.
    enzyme : str
        Enzyme used to digest the genome in the HiC experiment.
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
    logger.info("Generate HiC contact map for {0}".format(bin_name))
    bin_data = Bin(bin_name, project, min_size)
    bin_data.set_contigs()
    bin_data.set_large_contigs()
    bin_data.pairs = join(tmp_dir, bin_data.name + ".pairs")

    # Extract pairs of the bin.
    n_pairs = extract_pairs_bin(bin_data)
    logger.info("{0} pairs have been extracted.".format(n_pairs))

    # Launch hicstuff pipeline.
    hcp.full_pipeline(
        genome=bin_data.genome,
        input1=bin_data.pairs,
        distance_law=False,
        enzyme=enzyme,
        filter_events=filter_events,
        force=force,
        out_dir=out_dir,
        pcr_duplicates=pcr_duplicates,
        plot=False,
        prefix=bin_data.name,
        start_stage="pairs",
        threads=threads,
        tmp_dir=tmp_dir,
    )
