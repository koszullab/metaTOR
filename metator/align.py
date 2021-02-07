#!/usr/bin/env python3
# coding: utf-8
"""Alignement of reads

Align read files onto the assembly and return a 2D BED file with columns:
readIDA, contigA, posA, strandA, readIDB, contigB, posB, strandB. Reads are
mapped separately, sorted by names, then interleaved (rather than mapped in
paired-end mode) to capture the pairs mapping on two different contigs.

If the ligation sites are given, it will make an digestion at the ligation site
before teh alignement and create new pairs of reads according to the new 
fragments. For example if the forward read have a ligation site, it will cut the
read in two and create two reads the first one with the first part of the
forward read and the reverse read and the second one with one. 

This module contains all these alignment functions:
    - align
    - digest_liagtion_sites
    - digest_pair
    - merge_alignement
    - pairs_alignement
    - process_bamfile 
    - Writer
"""

import csv
import gzip
import multiprocessing
import pyfastx
import os
import re
import sys
import metator.io as mio
import pandas as pd
import pysam as ps
import subprocess as sp
from metator.log import logger
from os.path import join, basename
from pkg_resources import parse_version


def align(fq_in, index, bam_out, n_cpu):
    """Alignment
    Aligns reads of fq_in with bowtie2. Parameters of bowtie2 are set as
    --very-sensitive-local

    Parameters
    ----------
    fq_in : str
        Path to input fastq file to align. If multiple files are given, list of
        path separated by a comma.
    index : str
        Path to the bowtie2 index genome.
    bam_out : str
        Path where the alignment should be written in BAM format.
    n_cpu : int
        The number of CPUs to use for the alignment.
    """

    # Align the reads on the reference genome
    map_args = {"cpus": n_cpu, "fq": fq_in, "idx": index, "bam": bam_out}
    cmd = (
        "bowtie2 -x {idx} -p {cpus} --quiet --very-sensitive-local {fq}"
    ).format(**map_args)

    map_process = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
    sort_process = sp.Popen(
        "samtools sort -n -@ {cpus} -o {bam}".format(**map_args),
        shell=True,
        stdin=map_process.stdout,
    )
    out, err = sort_process.communicate()

    return 0


def digest_ligation_sites(fq_for, fq_rev, enzyme, mode, outdir, n_cpu):
    """Create new reads to manage pairs with a digestion and create multiple
    pairs to take into account all the contact present.

    The function write two files for both the forward and reverse fastq with the
    new reads. The new reads have at the end of the ID ":1" added to
    differentiate the different pairs created from one read.

    To make the function faster, for each reads only the first site of ligation
    is kept and the algorithm stop to search for others sites as the probability
    is very low. We already have very few pairs with ligation sites in both
    reads.

    Parameters
    ----------
    fq_for : str
        Path to all the forward fastq file to digest separated by a comma.
    fq_rev : str
        Path to all the reverse fatsq file to digest separated by a comma..
    enzyme : str
        The list of restriction enzyme used to digest the genome separated by a
        comma. Example: DpnII,HinfI
    mode : str
        Mode to use to make the digestion. Three values possible: "all",
        "for_vs_rev", "pile".
    outpdir : str
        Path for the ouput directory. The fastq will have their initial names
        with a suffix added "_digested".
    n_cpu : int
        Number of CPUs.
    """

    # Process the ligation sites given
    ligation_sites = mio.process_enzyme(enzyme)

    # Defined stop_token and stack_size for processing
    stop_token = "STOP"
    stack_size = 1000

    # Manage multiple fastq
    list_fq_for = fq_for.split(",")
    list_fq_rev = fq_rev.split(",")

    output_list_for = []
    output_list_rev = []

    # Sanity check
    if len(list_fq_for) != len(list_fq_rev):
        logger.error("Number of forward and reverse files are different")
        sys.exit(1)

    for i in range(len(list_fq_for)):
        fq_for = list_fq_for[i]
        fq_rev = list_fq_rev[i]

        logger.info("Reads in progress:\n{0}\n{1}".format(fq_for, fq_rev))

        # Create the two output file and add them to the list
        name_for = basename(fq_for).split(".")[0]
        name_rev = basename(fq_rev).split(".")[0]
        output_for = join(outdir, "{0}_digested.fq.gz".format(name_for))
        output_rev = join(outdir, "{0}_digested.fq.gz".format(name_rev))
        output_list_for.append(output_for)
        output_list_rev.append(output_rev)

        # Create count to have an idea of the digested pairs repartition.
        original_number_of_pairs = 0
        final_number_of_pairs = 0
        new_reads_for = ""
        new_reads_rev = ""
        n = 0

        # Start parrallel threading to compute the
        ctx = multiprocessing.get_context("spawn")
        queue = multiprocessing.Queue(n_cpu - 1)
        writer_process = multiprocessing.Process(
            target=Writer, args=(output_for, output_rev, queue, stop_token)
        )
        writer_process.start()

        # Iterate on all pairs
        for read_for, read_rev in zip(
            pyfastx.Fastq(fq_for, build_index=False),
            pyfastx.Fastq(fq_rev, build_index=False),
        ):

            # Count the numbers of original reads processed.
            original_number_of_pairs += 1

            # Count for stack size.
            n += 1

            # Extract components of the reads.
            for_name, for_seq, for_qual = read_for
            rev_name, rev_seq, rev_qual = read_rev

            # Sanity check to be sure all reads are with their mate.
            if for_name != rev_name:
                logger.error(
                    "The fastq files contains reads not sorted :\n{0}\n{1}".format(
                        read_for.id, read_rev.id
                    )
                )
                sys.exit(1)

            # Cut the forward and reverse reads at the ligation sites.
            for_seq_list, for_qual_list = cutsite_read(
                ligation_sites, for_seq, for_qual
            )
            rev_seq_list, rev_qual_list = cutsite_read(
                ligation_sites, rev_seq, rev_qual
            )

            # Write the new combinations of fragments.
            new_reads_for, new_reads_rev, final_number_of_pairs = write_pair(
                new_reads_for,
                new_reads_rev,
                for_name,
                for_seq_list,
                for_qual_list,
                rev_seq_list,
                rev_qual_list,
                mode,
                final_number_of_pairs,
            )

            # If stack full, add it in the queue.
            if n == stack_size:

                # Add the pair in the queue.
                pairs = (new_reads_for.encode(), new_reads_rev.encode())
                queue.put(pairs)

                # Empty the stack
                n = 0
                new_reads_for = ""
                new_reads_rev = ""

        # End the parallel processing.
        queue.put(stop_token)
        writer_process.join()

        # Return information on the different pairs created
        logger.info(
            "Number of pairs before digestion: {0}".format(
                original_number_of_pairs
            )
        )
        logger.info(
            "Number of pairs after digestion: {0}".format(final_number_of_pairs)
        )

    # Create list of the new localization of the FastQ.
    output_for = ",".join(output_list_for)
    output_rev = ",".join(output_list_rev)
    return output_for, output_rev


def cutsite_read(ligation_sites, seq, qual):
    """Find ligation sites in a given sequence.

    Parameters:
    -----------
    ligation_sites : str
        Regex of all possible ligations according to the given enzymes.
    seq : str
        Sequence where to search for ligation_sites.
    qual : str
        Quality values of the sequence given.

    Returns:
    --------
    list of str
        List of string of the sequences. The split is made at the start of the
        ligation sites.
    list of str
        List of string of the qualities.

    Examples:
    ---------
    >>> find_ligation_site("GA.TA.TC", "AAGAGTATTC", "FFF--FAFAF")
    (['AA', 'GAGTATTC'], ['FF', 'F--FAFAF'])
    """

    # Find the ligation sites.
    ligation_sites_list = []
    if re.search(ligation_sites, seq):
        ligation_sites_list = [
            site.start() for site in re.finditer(ligation_sites, seq)
        ]
    ligation_sites_list.append(len(seq))

    # Split the sequences on the ligation sites.
    seq_list = []
    qual_list = []
    left_site = 0
    for right_site in ligation_sites_list:
        seq_list.append(seq[left_site:right_site])
        qual_list.append(qual[left_site:right_site])
        left_site = right_site

    return seq_list, qual_list


def write_pair(
    new_reads_for,
    new_reads_rev,
    name,
    for_seq_list,
    for_qual_list,
    rev_seq_list,
    rev_qual_list,
    mode,
    final_number_of_pairs,
):
    """Function to write one pair with the combinations of fragment depending on
    the chosen mode.

    Parameters:
    -----------
    new_reads_for : str
        Stack of the new forward reads ready to be written.
    new_reads_rev : str
        Stack of the new reverse reads ready to be written.
    name : str
        Name of the fastq read.
    for_seq : str
        Forward sequence of the fastq read.
    for_qual : str
        Forward quality of the fastq read.
    rev_seq : str
        Reverse sequence of the fastq read.
    rev_qual : str
        Reverse quality of the fastq read.
    mode : str
        Mode to use to make the digestion. Three values possible: "all",
        "for_vs_rev", "pile".
    final_numbers_of_pairs : int
        Count of pairs after cutting.

    Returns:
    --------
    str
        Stack of forward reads ready to be written with the last pairs added.
    str
        Stack of reverse reads ready to be written with the last pairs added.
    int
        Count of pairs after cutting.
    """

    # Mode "for_vs_rev": Make contacts only between fragments from different
    # reads (one fragment from forward and one from reverse).
    if mode == "for_vs_rev":
        for i in range(len(for_seq_list)):
            for j in range(len(rev_seq_list)):
                final_number_of_pairs += 1
                new_reads_for += "@%s\n%s\n+\n%s\n" % (
                    name + ":" + str(i) + str(j),
                    for_seq_list[i],
                    for_qual_list[i],
                )
                new_reads_rev += "@%s\n%s\n+\n%s\n" % (
                    name + ":" + str(i) + str(j),
                    rev_seq_list[j],
                    rev_qual_list[j],
                )

    #  Mode "all": Make all the possible contacts between the fragments.
    elif mode == "all":
        seq_list = for_seq_list + rev_seq_list
        qual_list = for_qual_list + rev_qual_list
        for i in range(len(seq_list)):
            for j in range(i + 1, len(seq_list)):
                final_number_of_pairs += 1
                new_reads_for += "@%s\n%s\n+\n%s\n" % (
                    name + ":" + str(i) + str(j),
                    seq_list[i],
                    qual_list[i],
                )
                new_reads_rev += "@%s\n%s\n+\n%s\n" % (
                    name + ":" + str(i) + str(j),
                    seq_list[j],
                    qual_list[j],
                )

    # Mode "pile": Only make contacts bewteen two adjacent fragments.
    elif mode == "pile":
        seq_list = for_seq_list + rev_seq_list
        qual_list = for_qual_list + rev_qual_list
        for i in range(len(seq_list) - 1):
            final_number_of_pairs += 1
            new_reads_for += "@%s\n%s\n+\n%s\n" % (
                name + ":" + str(i) + str(i + 1),
                seq_list[i],
                qual_list[i],
            )
            new_reads_rev += "@%s\n%s\n+\n%s\n" % (
                name + ":" + str(i) + str(i + 1),
                seq_list[i + 1],
                qual_list[i + 1],
            )

    return new_reads_for, new_reads_rev, final_number_of_pairs


def merge_alignment(forward_aligned, reverse_aligned, out_file):
    """Merge forward and reverse alignment into one file with pairs which both
    reads are aligned on the genome. The final alignment  dataframe is written
    in the output file.

    Parameters
    ----------
    forward_alignement : str
        File with the table containing the data of the forward reads kept after
        the alignment. With five columns: ReadID, Contig, Position_start,
        Position_end, strand.
    reverse_alignement : str
        File with the table containing the data of the forward reads kept after
        the alignment. With five columns: ReadID, Contig, Position_start,
        Position_end, strand.
    out_file : str
        Path to write the output file.

    Returns
    -------
    str:
        File with the Table containing the alignement data of the pairs: ReadID,
        ContigA, Position_startA, Position_endA, StrandA, ContigB,
        Position_startB, Position_endB, StrandB
    """

    # Open files for reading and writing.
    with open(forward_aligned, "r") as for_file, open(
        reverse_aligned, "r"
    ) as rev_file, open(out_file, "w") as merge_bed:
        for_bed = csv.reader(for_file, delimiter="\t")
        rev_bed = csv.reader(rev_file, delimiter="\t")

        # Initialization.
        n_pairs = 0
        for_read = next(for_bed)
        rev_read = next(rev_bed)

        # Loop while at least one end of one fild is reached. It's possible to
        # advance like that as the two bed files are sorted on the id of the
        # reads.
        for i in range(10 ** 12):
            # Case of both reads of the pair map.
            if for_read[0] == rev_read[0]:
                merge_bed.write("\t".join(for_read + rev_read[1:]) + "\n")
                n_pairs += 1
                # print("w")
                try:
                    for_read = next(for_bed)
                except StopIteration:
                    break
                try:
                    rev_read = next(rev_bed)
                except StopIteration:
                    break
            # As the file is version sorted we have to do compare the two names
            # according to the version order.
            else:
                names = [for_read[0], rev_read[0]]
                names_sorted = sorted(names, key=parse_version)
                # Case of the forward read mapped but not the reverse. Indeed,
                # no line would have been added previously if the read didn't
                # map.
                if names == names_sorted:
                    try:
                        for_read = next(for_bed)
                    except StopIteration:
                        break
                # Same but with no mapped forward reads.
                else:
                    try:
                        rev_read = next(rev_bed)
                    except StopIteration:
                        break

    logger.info("{0} pairs aligned.".format(n_pairs))

    return out_file


def pairs_alignment(
    for_fq_in,
    rev_fq_in,
    min_qual,
    tmp_dir,
    ref,
    digestion_only,
    enzyme,
    mode,
    out_dir,
    n_cpu,
):
    """General function to do the whole alignment of both fastq.

    The Function write at the output directory location given as an argument and
    return a 2D bed file of the aligned reads with 9 columns: ReadID, ContigA,
    Position_startA, Position_endA, StrandA, ContigB, Position_startB,
    Position_endB, StrandB. The name of the file will be alignment.bed.

    The function will also digest the Fastq files on the sequences of the
    ligation sites if there are given and write the new fastq in the output
    directory. The name of the fastq will be the same as the previous ones with
    a "_digested" added in their names.

    Parameters
    ----------
    digestion_only : bool
        If True, will only make the digestion with the ligation sites and then
        stop.
    for_fq_in : str
        Path to input forward fastq file to align. If multiple files are given,
        list of path separated by a comma.
    rev_fq_in : str
        Path to input reverse fastq file to align. If multiple files are given,
        list of path separated by a comma.
    min_qual : int
        Minimum mapping quality required to keep Hi-C pairs.
    tmp_dir : str
        Path where temporary files should be written.
    ref : str
        Path to the index genome.
    enzyme : str
        The list of enzymes used to digested separated by a comma. Example:
        DpnII,HinfI.
    mode : str
        Mode to use to make the digestion. Three values possible: "all",
        "for_vs_rev", "pile".
    out_file : str
        Path to directory where to write the output file.
    n_cpu : int
        The number of CPUs to use for the alignment.

    Returns
    -------
    str:
        File with the Table containing the alignement data of the pairs: ReadID,
        ContigA, Position_startA, Position_endA, StrandA, ContigB,
        Position_startB, Position_endB, StrandB
    """

    # Throw error if index does not exist and digestion mode not enabled
    if not digestion_only:
        index = mio.check_fasta_index(ref, mode="bowtie2")
        if index is None:
            logger.error(
                "Reference index is missing, please build the bowtie2 index first."
            )
            sys.exit(1)

    # If asked will digest the reads with the ligtion sites before
    # alignment.
    if isinstance(enzyme, str):

        logger.info("Digestion of the reads:")
        logger.info("Enzyme used: {0}".format(enzyme))
        logger.info("Mode used to cut the reads: {0}".format(mode))
        # Create a temporary fastq with the trimmed reads at the ligation
        # sites.
        (for_fq_in, rev_fq_in,) = digest_ligation_sites(
            for_fq_in,
            rev_fq_in,
            enzyme,
            mode,
            out_dir,
            int(n_cpu),
        )

        if digestion_only:
            return 0

    # Create a temporary file to save the alignment.
    temp_alignment_for = join(tmp_dir, "temp_alignment_for.bam")
    temp_alignment_rev = join(tmp_dir, "temp_alignment_rev.bam")
    filtered_out_for = join(tmp_dir, "for_temp_alignment.bed")
    filtered_out_rev = join(tmp_dir, "rev_temp_alignment.bed")
    out_file = join(out_dir, "alignment.bed")

    # Align the forward reads
    logger.info("Alignement of the forward reads:")
    align(for_fq_in, index, temp_alignment_for, n_cpu)

    # Filters the aligned and non aligned reads
    process_bamfile(temp_alignment_for, min_qual, filtered_out_for)

    # forward_aligned = pd.DataFrame(
    #     csv.reader(open(filtered_out_for), delimiter="\t")
    # )

    # Align the reverse reads
    logger.info("Alignement of the reverse reads:")
    align(rev_fq_in, index, temp_alignment_rev, n_cpu)

    # Filters the aligned and non aligned reads
    process_bamfile()

    # reverse_aligned = pd.DataFrame(
    #     csv.reader(open(filtered_out_rev), delimiter="\t")
    # )

    # Merge alignement to create a pairs file
    logger.info("Merging the pairs:")
    merge_alignment(filtered_out_for, filtered_out_rev, out_file)

    # pairs.to_csv(out_file, sep="\t", index=False, header=False)
    return out_file


def process_bamfile(alignment, min_qual, filtered_out):
    """Filter alignment BAM files

    Reads all the reads in the input BAM alignment file. Keep reads in the
    output if they are aligned with a good quality saving their uniquely ReadID,
    Contig, Position_start, Position_end, strand to save memory.

    Parameters
    ----------
    alignment : str
        Path to the input temporary alignment.
    min_qual : int
        Minimum mapping quality required to keep a Hi-C pair.
    filtered_out : str
        Path to the output temporary bed alignement.

    Returns
    -------
    str:
        Path to the table containing the data of the reads mapping unambiguously
        and with a mapping quality superior to the threshold given. Five
        columns: ReadID, Contig, Position_start, Position_end, strand
    """
    # Check the quality and status of each aligned fragment.
    # Write the ones with good quality in the aligned dataframe.
    # Keep ID of those that do not map unambiguously to be trimmed.

    aligned_reads = 0
    temp_bam = ps.AlignmentFile(alignment, "rb", check_sq=False)
    with open(filtered_out, "a") as f:
        for r in temp_bam:
            if r.mapping_quality >= min_qual:
                if r.flag == 0:
                    aligned_reads += 1
                    read = str(
                        r.query_name
                        + "\t"
                        + r.reference_name
                        + "\t"
                        + str(r.reference_start)
                        + "\t"
                        + str(r.reference_end)
                        + "\t"
                        + "+\n"
                    )
                    f.write(read)
                elif r.flag == 16:
                    aligned_reads += 1
                    read = str(
                        r.query_name
                        + "\t"
                        + r.reference_name
                        + "\t"
                        + str(r.reference_start)
                        + "\t"
                        + str(r.reference_end)
                        + "\t"
                        + "-\n"
                    )
                    f.write(read)

    f.close()
    temp_bam.close()

    # Display alignement informations
    logger.info("{0} reads aligned.".format(aligned_reads))

    return 0


def Writer(output_for, output_rev, queue, stop_token):
    """Function to write the pair throw parallel processing.

    Parameters:
    -----------
    output_for : str
        Path to the output forward compressed fastq file.
    output_rev : str
        Path to the output reverse compressed fastq file.
    queue : multiprocessing.queues.Queue
        Queue for the multiprocesing of the whole file.
    stop_token : str
        Token to signal that the end of the file have been reached.
    """
    with gzip.open(output_for, "wb") as for_fq, gzip.open(
        output_rev, "wb"
    ) as rev_fq:
        while True:
            line = queue.get()
            if line == stop_token:
                return 0
            for_fq.write(line[0])
            rev_fq.write(line[1])