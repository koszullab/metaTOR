#!/usr/bin/env python3
# coding: utf-8

"""Alignement of reads

Align read files onto the assembly and return a 2D BED file with columns:
readIDA, contigA, posA, strandA, readIDB, contigB, posB, strandB. Reads are
mapped separately, sorted by names, then interleaved (rather than mapped in
paired-end mode) to capture the pairs mapping on two different contigs.

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
import sys
import metator.io as mio
import pysam
import subprocess as sp
from metator.log import logger
from os.path import join
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
        "bowtie2 -x {idx} -p {cpus} --quiet --very-sensitive-local {fq} --no-unal"
    ).format(**map_args)

    map_process = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
    sort_process = sp.Popen(
        "samtools sort -n -@ {cpus} -o {bam}".format(**map_args),
        shell=True,
        stdin=map_process.stdout,
    )
    out, err = sort_process.communicate()

    return 0


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

        # Loop while at least one end of one endd is reached. It's possible to
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
    out_file : str
        Path to directory where to write the output file.
    n_cpu : int
        The number of CPUs to use for the alignment.

    Returns
    -------
    list of str:
        List of path of the Files with the Table containing the alignement data
        of the pairs: ReadID, ContigA, Position_startA, Position_endA, StrandA,
        ContigB, Position_startB, Position_endB, StrandB.
    """

    # Check what is the reference. If fasta given build it. If not a bowtie2
    # index or a fasta thraw an error.
    index = mio.check_fasta_index(ref, mode="bowtie2")
    if index is None:
        if mio.check_is_fasta(ref):
            logger.info("Build index from the given fasta.")
            index = join(tmp_dir, "index")
            cmd = "bowtie2-build -q {0} {1}".format(ref, index)
            process = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
            out, err = process.communicate()
        else:
            logger.error(
                "Please give as a reference a bowtie2 index or a fasta."
            )
            sys.exit(1)

    # Iterates on all the fastq files:
    for_fq_list = for_fq_in.split(",")
    rev_fq_list = rev_fq_in.split(",")
    out_file_list = []

    for i in range(len(for_fq_list)):
        for_fq_in = for_fq_list[i]
        rev_fq_in = rev_fq_list[i]
        name = "alignement_" + str(i)

        # Create a temporary file to save the alignment.
        temp_alignment_for = join(out_dir, name + "_for.bam")
        temp_alignment_rev = join(out_dir, name + "_rev.bam")
        filtered_out_for = join(tmp_dir, name + "_for_temp.bed")
        filtered_out_rev = join(tmp_dir, name + "_rev_temp.bed")
        out_file = join(out_dir, name + ".bed")
        out_file_list.append(out_file)

        # Align the forward reads
        logger.info("Alignement of {0}:".format(for_fq_in))
        align(for_fq_in, index, temp_alignment_for, n_cpu)

        # Filters the aligned and non aligned reads
        process_bamfile(temp_alignment_for, min_qual, filtered_out_for)

        # Align the reverse reads
        logger.info("Alignement of {0}:".format(rev_fq_in))
        align(rev_fq_in, index, temp_alignment_rev, n_cpu)

        # Filters the aligned and non aligned reads
        process_bamfile(temp_alignment_rev, min_qual, filtered_out_rev)

        # Merge alignement to create a pairs file
        logger.info("Merging the pairs:")
        merge_alignment(filtered_out_for, filtered_out_rev, out_file)

    return out_file_list


def process_bamfile(alignment, min_qual, filtered_out):
    """Filter alignment BAM files

    Reads all the reads in the input BAM alignment file. Keep reads in the
    output if they are aligned with a good quality saving their uniquely ReadID,
    Contig, Position_start, Position_end, strand to save memory.

    Parameters:
    -----------
    alignment : str
        Path to the input temporary alignment.
    min_qual : int
        Minimum mapping quality required to keep a Hi-C pair.
    filtered_out : str
        Path to the output temporary bed alignement.

    Returns:
    --------
    str:
        Path to the table containing the data of the reads mapping unambiguously
        and with a mapping quality superior to the threshold given. Five
        columns: ReadID, Contig, Position_start, Position_end, strand
    """
    # Check the quality and status of each aligned fragment.
    # Write the ones with good quality in the aligned dataframe.
    # Keep ID of those that do not map unambiguously to be trimmed.

    aligned_reads = 0
    temp_bam = pysam.AlignmentFile(alignment, "rb", check_sq=False)
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
