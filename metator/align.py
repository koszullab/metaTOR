#!/usr/bin/env python3
# coding: utf-8

"""Alignement of reads

Align read files onto the assembly and return a tsv file with 9 columns: ReadID,
ContigA, Position_startA, Position_endA, StrandA, ContigB, Position_startB,
Position_endB, StrandB. Reads are mapped separately, sorted by names, then
interleaved (rather than mapped in paired-end mode) to capture the pairs mapping
on two different contigs.

This module contains all these alignment functions:
    - align
    - get_contact_pairs
    - merge_alignement
    - process_bamfile
"""

import csv
import pysam
import subprocess as sp
from metator.log import logger
from os.path import join
from pkg_resources import parse_version


def align(fq_in, index, bam_out, n_cpu):
    """Alignment
    Aligns reads of fq_in with bowtie2. Parameters of bowtie2 are set as
    -D 20 -R 3 -N 0 -L 16 -i S,1,0.50 --local and only the aligned reads are
    kept in the bam file.
    Multiple files could be given, but only one file will be written as output.

    Parameters:
    -----------
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
        "bowtie2 -x {idx} -p {cpus} --quiet -D 20 -R 3 -N 0 -L 16 -i S,1,0.50 --local {fq} --no-unal"
    ).format(**map_args)

    map_process = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
    sort_process = sp.Popen(
        "samtools sort -n -@ {cpus} -o {bam}".format(**map_args),
        shell=True,
        stdin=map_process.stdout,
    )
    out, err = sort_process.communicate()

    return 0


def get_contact_pairs(
    for_in,
    rev_in,
    index,
    min_qual,
    start,
    out_dir,
    tmp_dir,
    n_cpu,
):
    """General function to do the whole alignment of both fastq.

    The Function write at the output directory location given as an argument and
    return a tsv file of the aligned reads with 9 columns: ReadID, ContigA,
    Position_startA, Position_endA, StrandA, ContigB, Position_startB,
    Position_endB, StrandB. The name of the file will be alignment.txt.

    Two start stages are possible, from fastq or bam files.

    Parameters:
    -----------
    for_in : str
        Path to input forward fastq or bam file to align. If multiple files are
        given, list of path separated by a comma.
    rev_in : str
        Path to input reverse fastq or bam file to align. If multiple files are
        given, list of path separated by a comma.
    index : str
        Path to the bowtie2 index of the assembly.
    min_qual : int
        Minimum mapping quality required to keep Hi-C pairs.
    start : str
        Either fastq or bam. Starting point for the pipeline.
    out_dir : str
        Path to directory where to write the output file.
    tmp_dir : str
        Path where temporary files should be written.
    n_cpu : int
        The number of CPUs to use for the alignment.

    Returns:
    --------
    list of str:
        List of path of the Files with the table containing the alignement data
        of the pairs: ReadID, ContigA, Position_startA, Position_endA, StrandA,
        ContigB, Position_startB, Position_endB, StrandB.
    """

    # Iterates on all the input files:
    for_list = for_in.split(",")
    rev_list = rev_in.split(",")
    out_file_list = []

    for i in range(len(for_list)):
        for_in = for_list[i]
        rev_in = rev_list[i]
        name = "alignment_" + str(i)
        out_file = join(tmp_dir, name + ".txt")
        out_file_list.append(out_file)

        # Align if necessary
        if start == "fastq":
            # Create files to save the alignment.
            alignment_for = join(out_dir, name + "_for.bam")
            alignment_rev = join(out_dir, name + "_rev.bam")

            # Align the forward reads
            logger.info("Alignment of {0}:".format(for_in))
            align(for_in, index, alignment_for, n_cpu)

            # Align the reverse reads
            logger.info("Alignment of {0}:".format(rev_in))
            align(rev_in, index, alignment_rev, n_cpu)

        elif start == "bam":
            alignment_for = for_in
            alignment_rev = rev_in

        else:
            logger.error("Start argument should be either 'fastq' or 'bam'.")
            raise ValueError

        # Create files to save the alignment.
        alignment_temp_for = join(tmp_dir, name + "_for_temp.txt")
        alignment_temp_rev = join(tmp_dir, name + "_rev_temp.txt")

        # Filters the aligned and non aligned reads from the forward and
        # reverse bam files.
        aligned_reads_for = process_bamfile(
            alignment_for, min_qual, alignment_temp_for
        )
        aligned_reads_rev = process_bamfile(
            alignment_rev, min_qual, alignment_temp_rev
        )
        logger.info(
            "{0} forward reads aligned and {1} reverse reads aligned".format(
                aligned_reads_for, aligned_reads_rev
            )
        )

        # Merge alignement to create a pairs file
        logger.info("Merging the pairs:")
        merge_alignment(alignment_temp_for, alignment_temp_rev, out_file)

    return out_file_list


def merge_alignment(forward_aligned, reverse_aligned, out_file):
    """Merge forward and reverse alignment into one file with only pairs which
    have both reads are aligned on the genome with 9 columns: ReadID, ContigA,
    Position_startA, Position_endA, StrandA, ContigB, Position_startB,
    Position_endB, StrandB.

    Parameters:
    -----------
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

    Returns:
    --------
    str:
        File with the table containing the alignement data of the pairs: ReadID,
        ContigA, Position_startA, Position_endA, StrandA, ContigB,
        Position_startB, Position_endB, StrandB
    """

    # Open files for reading and writing.
    with open(forward_aligned, "r") as for_file, open(
        reverse_aligned, "r"
    ) as rev_file, open(out_file, "w") as merged:
        for_bam = csv.reader(for_file, delimiter="\t")
        rev_bam = csv.reader(rev_file, delimiter="\t")

        # Initialization.
        n_pairs = 0
        for_read = next(for_bam)
        rev_read = next(rev_bam)

        # Loop while at least one end of one endd is reached. It's possible to
        # advance like that as the two tsv files are sorted on the id of the
        # reads.
        for i in range(10 ** 12):
            # Case of both reads of the pair map.
            if for_read[0] == rev_read[0]:
                merged.write("\t".join(for_read + rev_read[1:]) + "\n")
                n_pairs += 1
                # print("w")
                try:
                    for_read = next(for_bam)
                except StopIteration:
                    break
                try:
                    rev_read = next(rev_bam)
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
                        for_read = next(for_bam)
                    except StopIteration:
                        break
                # Same but with no mapped forward reads.
                else:
                    try:
                        rev_read = next(rev_bam)
                    except StopIteration:
                        break

    logger.info("{0} pairs aligned.".format(n_pairs))

    return out_file


def process_bamfile(alignment, min_qual, filtered_out):
    """Filter alignment BAM files

    Reads all the reads in the input BAM alignment file. Keep reads in the
    output if they are aligned with a good quality (greater than min quality
    threshold given) saving their only some columns: ReadID, Contig,
    Position_start, Position_end, strand to save memory.

    Parameters:
    -----------
    alignment : str
        Path to the input temporary alignment.
    min_qual : int
        Minimum mapping quality required to keep a Hi-C pair.
    filtered_out : str
        Path to the output temporary tsv alignement.

    Returns:
    --------
    str:
        Path to the table containing the data of the reads mapping unambiguously
        and with a mapping quality superior to the threshold given. Five
        columns: ReadID, Contig, Position_start, Position_end, strand.
    """

    # Check the quality and status of each aligned fragment.
    aligned_reads = 0
    temp_bam = pysam.AlignmentFile(alignment, "rb", check_sq=False)
    with open(filtered_out, "a") as f:
        for r in temp_bam:
            # Check mapping quality
            if r.mapping_quality >= min_qual:
                # Check Mapping (0 or 16 flags are kept only)
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
    temp_bam.close()

    return aligned_reads
