#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Generates a HiC contact map of a MetaTOR object.

General utility functions and class for handling one MetaTOR object and for
generating its HiC contact map using hicstuff pipeline, from the pairs files,
at the resolution of the restriction fragments, or from a .metacool contact
map, at that of its bins (one per contig in a map of metaTOR), as metacooler
extract writes it.

Core class to handle bin object easily:
    - MetatorObject

Core function to build the network are:
    - extract_pairs
    - generate_contact_map
    - extract_contact_map
"""

import os
import hicstuff.pipeline as hcp
import metator.io as mio
import metator.metacool as mtmc
from metator.log import logger
from os.path import join
import pandas as pd
import pyfastx


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
        pairs : List of str
            List of path of the ".pairs" file.
        min_size : int
            Size threshold used to filter contigs. Default: 5000.
        """
        self.assembly = assembly
        self.contigs_data = contig_data_file
        self.pairs_files = pairs
        self.min_size = min_size
        # The name given, which names the output files.
        self.label = str(name)
        self.set_metator_object(metator_object, name)
        self.contigs = None
        self.large_contigs = None
        self.contigs_size = None
        self.fasta = None
        self.pairs = None

    def get_contigs_data(self):
        """Method to get contigs data Dataframe

        Returns:
        --------
        pandas.core.frame.DataFrame
            Table with informations from metator binning for each contigs.
        """
        return mio.rename_1x_columns(pd.read_csv(self.contigs_data, sep="\t"))

    def set_contigs(self):
        """Method to extract the list of contigs of the class from the contig
        data file and their size.
        """
        contigs_data = self.get_contigs_data()
        # A column named as metaTOR 1.x names it (Final_bin) is read under its
        # 2.0 name (see metator.io.rename_1x_columns), unless the table also
        # has a column of that name.
        if self.object not in contigs_data.columns:
            self.object = mio.CONTIG_DATA_1X_COLUMNS.get(self.object, self.object)
        try:
            contigs_data[self.object]
        except KeyError:
            logger.error("The object that you gave is not in the table.")
            raise ValueError from None
        column = contigs_data[self.object]
        ids = ("Core_cluster_ID", "Overlapping_cluster_ID", "Recursive_cluster_ID")
        # Their 1.x names too, in a table that has both (merged by hand).
        if self.object in ids or mio.CONTIG_DATA_1X_COLUMNS.get(self.object) in ids:
            # Bin ids are written zero-padded ("00008") and read back as strings
            # or integers depending on the column: compare them as numbers.
            selected = pd.to_numeric(column, errors="coerce") == int(self.name)
            # Recursive bins are numbered within their overlapping bin.
            if self.overlapping is not None:
                overlapping = pd.to_numeric(contigs_data["Overlapping_cluster_ID"], errors="coerce")
                selected &= overlapping == int(self.overlapping)
        else:
            selected = column.astype(str) == str(self.name)
        self.contigs = list(contigs_data["Name"][selected])
        self.contigs_size = list(contigs_data["Size"][selected])

    def set_large_contigs(self):
        """Method to keep only contigs bigger than the threshold given to remove
        small bins which will be smaller than the binning value and which will
        prevent a good scaffolding with Instagraal.
        """
        self.large_contigs = []
        for index, value in enumerate(self.contigs_size):
            if value >= self.min_size:
                self.large_contigs.append(self.contigs[index])
        self.contigs = self.large_contigs
        self.contigs_size = [size for size in self.contigs_size if size >= self.min_size]

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
            the name of the bin. Example: "NODE_1" or "MetaTOR_1_0". A
            recursive bin is numbered within its overlapping bin: it needs its
            name, <prefix>_<overlapping bin ID>_<recursive bin ID>, whose prefix
            may hold "_" (--prefix).
        """
        # The overlapping bin of a recursive bin.
        self.overlapping = None
        if metator_object == "contig":
            self.object = "Name"
            self.name = name
        elif metator_object == "core_bin":
            self.object = "Core_cluster_ID"
            try:
                int(name)
                self.name = str(name)
            except ValueError as object_no_exist:
                logger.error(
                    "With the core_bin object, the name should be the numeric ID of\
                     the core cluster."
                )
                raise ValueError from object_no_exist
        elif metator_object == "overlapping_bin":
            self.object = "Overlapping_cluster_ID"
            # A numeric ID, or a name: <prefix>_<ID> or
            # <prefix>_<ID>_<recursive bin ID>, whose prefix may hold "_".
            fields = str(name).split("_")
            self.name = fields[-2] if len(fields) >= 3 else fields[-1]
            if not self.name.isdigit():
                logger.error(
                    "Overlapping cluster name should be either numerical ID or a name like 'MetaTOR_1' or 'MetaTOR_1_0'."
                )
                raise ValueError
        elif metator_object == "recursive_bin":
            self.object = "Recursive_cluster_ID"
            fields = str(name).split("_")
            if len(fields) < 3 or not (fields[-2].isdigit() and fields[-1].isdigit()):
                logger.error(
                    "A recursive cluster is numbered within its overlapping cluster: give its name, "
                    "<prefix>_<overlapping cluster ID>_<recursive cluster ID>, such as 'metator_00012_00003'."
                )
                raise ValueError
            self.overlapping, self.name = fields[-2], fields[-1]
            if int(self.name) <= 0:
                logger.error("A recursive cluster should have an id bigger than 0.")
        elif metator_object == "final_bin":
            self.object = "Final_cluster"
            self.name = name
        else:
            # A column of the contig data (see set_contigs for its name of 1.x).
            self.object = str(metator_object)
            try:
                self.name = int(name)
            except ValueError:
                self.name = str(name)

    def write_fasta(self, out_dir):
        """Method to write the new fasta with only the contigs of interest in
        the outdir directory.

        Parameters:
        -----------
        out_dir : str
            Path to the output directory where to write the new fasta file.
        """
        self.fasta = join(out_dir, f"{self.label}.fa")
        mio.write_contigs(pyfastx.Fasta(self.assembly), self.contigs, self.fasta)


def extract_pairs(metator_data):
    """Function to extract the pairs of a MetaTOR object from pairs files.

    Parameters:
    -----------
    metator_data : object from class MetatorObject
        Object of the class MetatorObject. It should have the parameters contigs
        and the output pairs file of the metator object set.

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
        output_pairs.write("#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n")
        for contig_id, contig in enumerate(metator_data.contigs):
            output_pairs.write(
                "#chromsize: {0} {1}\n".format(
                    contig,
                    metator_data.contigs_size[contig_id],
                )
            )
        for pairs_file in metator_data.pairs_files:
            # Check if the pypairix index exist.
            pairs_data = mio.get_pairs_data(pairs_file)
            # Need a sorted (chr1 chr2 pos1 pos2) pair file indexed with pypairix.
            for contig_id1, contig in enumerate(metator_data.contigs):
                # Only need to retrieve the upper triangle.
                for contig_id2 in range(contig_id1, len(metator_data.contigs)):
                    pairs_lines = pairs_data.query2D(
                        contig,
                        0,
                        metator_data.contigs_size[contig_id1],
                        metator_data.contigs[contig_id2],
                        0,
                        metator_data.contigs_size[contig_id2],
                        1,
                    )
                    for pairs_line in pairs_lines:
                        pairs_line = pairs_line[:7]
                        n_pairs += 1
                        output_pairs.write("\t".join(pairs_line) + "\n")
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
        Path of the ".pairs" file or bgzip indexed pair file. If more than one
        is given, files should be separated by a comma.
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
    pcr_duplicates : bool
        If True, PCR duplicates will be filtered based on genomic positions.
        Pairs where both reads have exactly the same coordinates are considered
        duplicates and only one of those will be conserved. Default: False.
    threads : int
        Numbers of threads to use. Default: 1.
    """

    # Extract bin information from metator outdir.
    logger.info("Generate HiC contact map for %s", name)
    metator_data = MetatorObject(metator_object, name, assembly, contig_data_file, pairs, min_size)
    metator_data.set_contigs()
    if min_size > 0:
        metator_data.set_large_contigs()
    metator_data.write_fasta(out_dir)
    metator_data.pairs = join(tmp_dir, name + ".pairs")

    # Extract pairs of the bin.
    n_pairs = extract_pairs(metator_data)

    if n_pairs == 0:
        logger.info("No pairs have been extracted")

    else:
        logger.info("%d pairs have been extracted.", n_pairs)

        # Launch hicstuff pipeline.
        hcp.full_pipeline(
            genome=metator_data.fasta,
            input1=metator_data.pairs,
            distance_law=False,
            enzyme=enzyme,
            filter_events=filter_events,
            force=force,
            mat_fmt=mat_fmt,
            out_dir=out_dir,
            pcr_duplicates=pcr_duplicates,
            prefix=name,
            plot=False,
            start_stage="pairs",
            threads=threads,
            tmp_dir=tmp_dir,
        )
    return n_pairs


def extract_contact_map(
    assembly,
    contig_data_file,
    name,
    metacool,
    out_dir,
    force=False,
    metator_object="final_bin",
    min_size=5000,
):
    """Generate the contact map of a MetaTOR object from a .metacool contact
    map instead of the pairs files: the map of its contigs, as metacooler
    extract writes it, at the bins of the map (one per contig in a map of
    metaTOR), where generate_contact_map builds the map of their restriction
    fragments with hicstuff. Its contigs are those that generate_contact_map
    takes, and their contacts those of the pairs that it extracts from the
    pairs files that the map was made of: each pair whose two reads lie on
    the contigs, within their length.

    It writes <name>.fa, the contigs of the object, as generate_contact_map
    does, and <name>.cool, its contact map: a cooler of its contigs, in the
    order of the map, with their contacts. A contig that the map lacks is
    left out, with a warning; with none of the contigs, no map is written,
    and an earlier <name>.cool is removed (force).

    Parameters:
    -----------
    assembly : str
        Path to the fasta file containing the contigs of interest. Could be the
        whole or the extracted contigs of one bin.
    contig_data_file : str
        Path to the contig_data_final.txt file form MetaTOR output: the
        lengths of the contigs (Size) must be theirs in the map.
    name : str
        Name of the object. Could be the name of a contig, an id of a bin or the
        name of the bin. Example: "NODE_1" or "MetaTOR_1_0".
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection: ``file.mmetacool::entries/<name>``.
    out_dir : str
        Path where output files should be written.
    force : bool
        If True, overwrite an existing <name>.cool. Default: False.
    metator_object : str
        Object to extract contigs to build the matrix. Either "contig",
        "core_bin", "overlapping_bin", "recursive_bin", "final_bin" or "other".
    min_size : int
        Minimum contig size required to keep it.

    Returns:
    --------
    list of str:
        The contigs of the map written, in its order.

    Raises:
    -------
    FileExistsError
        If <name>.cool exists, without force: refused before any work, the
        reason logged as an error.
    metator.metacool.MetacoolError
        If the map cannot be read, stores both triangles of the matrix, or
        gives a contig of the object another length than the contig data:
        refused before <name>.fa is written and <name>.cool replaced. Or if
        metacooler cannot write <name>.cool.
    """
    logger.info("Generate HiC contact map for %s from %s", name, metacool)
    cool = join(out_dir, f"{name}.cool")
    if os.path.lexists(cool) and not force:
        reason = f"{cool} exists: give -F to replace it with the contact map of {name}."
        logger.error(reason)
        raise FileExistsError(reason)
    metator_data = MetatorObject(metator_object, name, assembly, contig_data_file, [], min_size)
    metator_data.set_contigs()
    if min_size > 0:
        metator_data.set_large_contigs()
    names, sizes = metator_data.contigs or [], metator_data.contigs_size or []
    contigs = {str(contig): int(size) for contig, size in zip(names, sizes, strict=True)}
    # The map is checked before any output is written or replaced: a
    # refused run leaves those of an earlier run as they were.
    written = mtmc.held_contigs(metacool, contigs, contig_data_file, f"the contigs of {name}")
    metator_data.write_fasta(out_dir)
    if os.path.lexists(cool):
        os.remove(cool)
    if written:
        mtmc.extract_contigs(metacool, written, cool)
    else:
        logger.info("%s holds no contig of %s of at least %d bp: no contact map is written.", metacool, name, min_size)
    return written
