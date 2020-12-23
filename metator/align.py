#!/usr/bin/env python3
# coding: utf-8
"""Alignement of reads

Align read files onto the assembly and return a 2D BED file with columns:
readIDA, contigA, posA, strandA, readIDB, contigB, posB, strandB. Reads are
mapped separately, sorted by names, then interleaved (rather than mapped in
paired-end mode) to capture the pairs mapping on two different contigs.

If the ligation sites are given, it will make an iteration on the unmapped, or
multi-mapped, cut them at the ligation site and try to align them again to
improve the numbers of reads uniquely mapped. However it's time consumming and
doesn"t improve a lot with short reads (35bp). 

This module contains all these alignment functions:
    - align
    - alignement
    - cutsite_trimming
    - merge_alignement
    - pairs_alignement
    - process_bamfile
    - trim_reads 
"""

import csv
import os
from os.path import join
import pandas as pd
import pysam as ps
import subprocess as sp
import sys
import metator.io as mio
from metator.log import logger


# TODO: Modify the cutsite trimming to try to map both part of the cut read.
# TODO: cutsite trimming in python language
# TODO: adpat trimmed reads
# TODO: Check description of the function


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


def alignment(fq_in, fastq, min_qual, tmp_dir, index, ligation_sites, n_cpu):
    """General function to manage to do the whole alignement with or without
    trimming at the ligation sites. Returns a file with five columns ReadID,
    ContigID, Start, End, Strand.

    Parameters:
    -----------
    fq_in : str
        Path to input fastq file to align. If multiple files are given, list of
        path separated by a comma.
    fastq : str
        "for" or "rev" to precise the fastq file used for the alignement.
    min_qual : int
        Minimum mapping quality required to keep Hi-C pairs.
    tmp_dir : str
        Path where temporary files should be written.
    ref : str
        Path to the index genome.
    ligation_sites : str
        The list of ligations site possible depending on the restriction enzymes
        used separated by a comma. Exemple of the ARIMA kit:
        GATCGATC,GANTGATC,GANTANTC,GATCANTC
    n_cpu : int
        The number of CPUs to use for the alignment.

    Return
    ------
    pandas.core.frame.DataFrame:
        The file will have five columns: ReadID, ContigID, Position_start,
        Position_end, Strand.
    """

    # Create a temporary file to save the alignment.
    temp_alignment_raw = join(tmp_dir, "temp_alignment_raw.bam")
    filtered_out = join(tmp_dir, fastq + "_temp_alignment.bed")

    # Align the reads
    align(fq_in, index, temp_alignment_raw, n_cpu)

    # Filters the aligned and non aligned reads
    unaligned = process_bamfile(temp_alignment_raw, min_qual, filtered_out)

    # Iteration if some ligation sites are given
    if isinstance(ligation_sites, str):

        # Create a temporary fastq with the trimmed reads at the ligation sites.
        fq_trimmed = trimmed_reads(tmp_dir, fq_in, unaligned, ligation_sites)

        # Create a temporary file to save the alignment.
        temp_alignment_trimmed = join(tmp_dir, "temp_alignment_trimmed.bam")

        # Align the trimmed reads
        align(fq_trimmed, index, temp_alignement_trimmed, n_cpu)

        # Filter the aligned reads
        unaligned_trimmed = process_bamfile(
            temp_alignment_trimmed, min_qual, filtered_out
        )

    aligned = pd.DataFrame(csv.reader(open(filtered_out), delimiter="\t"))

    return aligned


# def cutsite_trimming(
#     read,
#     ligation_site
# ):
#     """Cut a read at the ligation site if there is one

#     Parameters
#     ----------
#     read : str
#         Sequence of the read to trimmed.
#     ligation sites: str
#         Sequence of the ligation site to use for the trimming.

#     Return
#     ------
#     str :
#         Trimmed reaads at the ligation site if there is one, else return 0.
#     """
#     return


def merge_alignment(forward_aligned, reverse_aligned, out_file):
    """Merge forward and reverse alignment into one file with pairs which both
    reads are aligned on the genome. The final alignment  dataframe is written
    in the output file.

    Parameters
    ----------
    forward_alignement : pandas.core.frame.DataFrame
        Table containing the data of the forward reads kept after the alignment.
        With five columns: ReadID, Contig, Position_start, Position_end, strand.
    reverse_alignement : pandas.core.frame.DataFrame
        Table containing the data of the forward reads kept after the alignment.
        With five columns: ReadID, Contig, Position_start, Position_end, strand.
    out_file : str
        Path to write the output file.

    Return
    ------
    pandas.core.frame.DataFrame:
        Table conatining the alignement data of the pairs: ReadID, ContigA,
        Position_startA, Position_endA, StrandA, ContigB, Position_startB,
        Position_endB, StrandB
    """

    # Merge the two dataframe on the readID column
    pairs = pd.merge(forward_aligned, reverse_aligned, on=0, how="inner")

    logger.info("{0} pairs aligned.".format(len(pairs)))

    return pairs


def pairs_alignment(
    for_fq_in,
    rev_fq_in,
    min_qual,
    tmp_dir,
    ref,
    ligation_sites,
    out_file,
    n_cpu,
):
    """General function to do the whole alignement of both fastq.
    The Function write at the output file location given as an argument and
    return a 2D bed file of the aligned reads with 9 columns: ReadID, ContigA,
    Position_startA, Position_endA, StrandA, ContigB, Position_startB,
    Position_endB, StrandB

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
    ligation_sites : str
        The list of ligations site possible depending on the restriction enzymes
        used separated by a comma. Exemple of the ARIMA kit:
        GATCGATC,GANTGATC,GANTANTC,GATCANTC
    out_file : str
        Path to write the output file.
    n_cpu : int
        The number of CPUs to use for the alignment.

    Return
    ------
    pandas.core.frame.DataFrame:
        Table conatining the alignement data of the pairs: ReadID, ContigA,
        Position_startA, Position_endA, StrandA, ContigB, Position_startB,
        Position_endB, StrandB
    """

    # Counting reads forward and reverse reads
    total_reads_forward = 0
    with mio.read_compressed(for_fq_in) as inf:
        for _ in inf:
            total_reads_forward += 1
    total_reads_forward /= 4
    total_reads_reverse = 0
    with mio.read_compressed(rev_fq_in) as inf:
        for _ in inf:
            total_reads_reverse += 1
    total_reads_reverse /= 4

    # Safety check: Same numbers of reads in the forward and reverse fastq file.
    if total_reads_forward != total_reads_reverse:
        logger.warning(
            "Different numbers of forward and reverse reads. Please check if \
                your files are not corrupted"
        )

    logger.info(
        "{0} paired-end reads in the library.".format(total_reads_reverse)
    )

    # Throw error if index does not exist
    index = mio.check_fasta_index(ref, mode="bowtie2")
    if index is None:
        logger.error(
            "Reference index is missing, please build the bowtie2 index first."
        )
        sys.exit(1)

    # Align forward and reverse reads
    forward_aligned = alignment(
        for_fq_in, "for", min_qual, tmp_dir, index, ligation_sites, n_cpu
    )

    reverse_aligned = alignment(
        rev_fq_in, "rev", min_qual, tmp_dir, index, ligation_sites, n_cpu
    )

    # Merge alignement to create a pairs file
    pairs = merge_alignment(forward_aligned, reverse_aligned, out_file)

    pairs.to_csv(out_file, sep="\t", index=False, header=False)

    return pairs


def process_bamfile(alignment, min_qual, filtered_out):
    """Filter alignment BAM files

    Reads all the reads in the input BAM alignment file.
    Keep reads in the output if they are aligned with a good
    quality saving their uniquely ReadID, Contig, Position_start, Position_end,
    strand to save memory. Otherwise add their name in a set to stage them in
    order to trim them.

    Parameters
    ----------
    temp_alignment : str
        Path to the input temporary alignment.
    min_qual : int
        Minimum mapping quality required to keep a Hi-C pair.
    filtered_out : str
        Path to the output temporary bed alignement.

    Returns
    -------
    pandas..core.frame.DataFrame:
        Table containing the data of the reads mapping unambiguously and with a
        mapping quality superior to the threshold given. Five columns: ReadID,
        Contig, Position_start, Position_end, strand
    set:
        Contains the names reads that did not align.
    """
    # Check the quality and status of each aligned fragment.
    # Write the ones with good quality in the aligned dataframe.
    # Keep ID of those that do not map unambiguously to be trimmed.

    aligned_reads = 0
    unaligned = set()
    temp_bam = ps.AlignmentFile(alignment, "rb", check_sq=False)
    f = open(filtered_out, "a")
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
            else:
                unaligned.add(r.query_name)
        else:
            unaligned.add(r.query_name)

    f.close()
    temp_bam.close()

    # Display alignement informations
    logger.info("{0} reads aligned.".format(aligned_reads))

    return unaligned


# def trim_reads(
#         tmp_dir,
#         infile,
#         unaligned_set,
#         ligation_sites
#     ):
#     """Trim read ends
#     Trim the fastq sequences in infile to an auxiliary file in the temporary
#     folder using the ligation site sequence provided.
#     Parameters
#     ----------
#     tmp_dir : str
#         Path to the temporary folder.
#     infile : str
#         Path to the fastq file to truncate.
#     unaligned_set : set
#         Contains the names of all reads that did not map unambiguously in
#         previous rounds.
#     ligation_sites : str
#         The list of ligations site possible depending on the restriction
#         enzymes used separated by a comma. Exemple of the ARIMA kit:
#         GATCGATC,GANTGATC,GANTANTC,GATCANTC

#     Returns
#     -------
#     str :
#         Path to the output fastq file containing trimmed reads.
#     """
#     outfile_tmp = "{0}/unmapped.fastq".format(tmp_dir)
#     outfile = "{0}/trimmed.fastq".format(tmp_dir)
#     with ps.FastxFile(infile, "r") as inf, open(outfile_tmp, "w") as outf:
#         for entry in inf:
#             # If the read did not align in previous round or this is the first round
#             if (entry.name in unaligned_set):
#                 entry.sequence = entry.sequence
#                 entry.quality = entry.quality
#                 outf.write(str(entry) + "\n")
#     map_args = {
#         "fq": outfile_tmp,
#         "out": outfile,
#         "cutsite": ligation_sites,
#     }
#     cmd = "cutsite_trimming --fastq {fq} --out {out} --cutsite {cutsite}".format(**map_args)
#     trim_process = sp.Popen(cmd, shell=True)
#     return outfile
