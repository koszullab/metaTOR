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
    - process_bwa_bamfile
"""

import csv
import pysam
import shutil as st
import subprocess as sp
import hicstuff.cutsite as hcc
import hicstuff.io as hio
import hicstuff.iteralign as hci
import metator.io as mio
import metator.network as mtn
from metator.log import logger
from os.path import join
from looseversion import LooseVersion


def align(
    fq_in,
    index,
    aligner,
    bam_out,
    n_cpu,
    tmp_dir,
    iterative=False,
    fq_in_2=None,
):
    """Alignment
    Aligns reads of fq_in with bowtie2. Parameters of bowtie2 are set as
    --very-sensitive-local and only the aligned reads are
    kept in the bam file.
    Multiple files could be given, but only one file will be written as output.

    Parameters:
    -----------
    fq_in : str
        Path to input fastq file to align.
    index : str
        Path to the bowtie2 index genome.
    aligner : str
        Name of the aligner algorithm to use. Either bowtie2 or bwa.
    bam_out : str
        Path where the alignment should be written in BAM format.
    n_cpu : int
        The number of CPUs to use for the alignment.
    tmp_dir:
        Path to directory to store temporary files.
    iterative:
        Wether to use the iterative mapping procedure (truncating reads and
        extending iteratively)
    fq_in_2 : str
        Path to the reverse fastq file for bwa (map both reads at the same time
        but do not take into account their poistion as pairs.)

    Returns:
    --------
    int :
        0
    """

    if iterative:
        iter_tmp_dir = hio.generate_temp_dir(tmp_dir)
        tmp_bam = join(tmp_dir, "tmp.bam")
        hci.iterative_align(
            fq_in,
            tmp_dir=iter_tmp_dir,
            ref=index,
            n_cpu=n_cpu,
            bam_out=tmp_bam,
            min_qual=40,
            aligner=aligner,
            read_len=20,
        )
        st.rmtree(iter_tmp_dir)
        sp.call(
            "samtools sort -n -@ {n_cpu} -o {bam} {tmp}".format(
                n_cpu=n_cpu, tmp=tmp_bam, bam=bam_out
            ),
            shell=True,
        )

    else:
        # Align the reads on the reference genome with the chosen aligner.
        map_args = {
            "cpus": n_cpu,
            "fq": fq_in,
            "idx": index,
            "bam": bam_out,
            "fq_rev": fq_in_2,
        }
        if aligner == "bwa":
            cmd = "bwa mem -5SP -t {cpus} {idx} {fq} {fq_rev}".format(
                **map_args
            )
        elif aligner == "bowtie2":
            cmd = (
                "bowtie2 -x {idx} -p {cpus} --very-sensitive-local {fq} --no-unal"
            ).format(**map_args)

        # Write the outputfile in a temporary bam file.
        map_process = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
        sort_process = sp.Popen(
            "samtools sort -n -@ {cpus} -o {bam}".format(**map_args),
            shell=True,
            stdin=map_process.stdout,
        )
        _out, _err = sort_process.communicate()
        mapping_values = map_process.stderr.read()
        for line in mapping_values.split(b"\n"):
            logger.info(f"{line.decode('utf-8')}")
    return 0


def get_contact_pairs(
    for_in,
    rev_in,
    index,
    assembly,
    aligner,
    aligner_mode,
    min_qual,
    start,
    depth_file,
    enzyme,
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

    Parameters
    ----------
    for_in : str
        Path to input forward fastq or bam file to align. If multiple files are
        given, list of path separated by a comma.
    rev_in : str
        Path to input reverse fastq or bam file to align. If multiple files are
        given, list of path separated by a comma.
    index : str
        Path to the bowtie2 index of the assembly.
    assembly : str
        The initial assembly path acting as the alignment file's reference
        assembly.
    aligner : str
        Either 'bowtie2' or 'bwa' aligner used or to be use to map the reads.
    aligner_mode : str
        Either 'normal', 'iterative' or 'cutsiste'. Mode to align the HiC reads.
    min_qual : int
        Minimum mapping quality required to keep Hi-C pairs.
    start : str
        Either fastq or bam. Starting point for the pipeline.
    depth_file : str or None
        Path to the depth.txt file from jgi_summarize_bam_contig_depths from
        Metabat2 Software.
    enzyme : str or None
        String that contains the names of the enzyme separated by a comma.
    out_dir : str
        Path to directory where to write the output file.
    tmp_dir : str
        Path where temporary files should be written.
    n_cpu : int
        The number of CPUs to use for the alignment.

    Returns
    -------
    list of str:
        List of path of the Files with the table containing the alignement data
        of the pairs: ReadID, ContigA, Position_startA, Position_endA, StrandA,
        ContigB, Position_startB, Position_endB, StrandB.
    dict:
        Dictionnary of the all the contigs from the assembly, the contigs names
        are the keys to the data of the contig available with the following
        keys: "id", "length", "GC", "hit", "coverage". Coverage still at 0 and
        need to be updated later.
    dict:
        Dictionnary for hit information on each contigs.
    """

    # Iterates on all the input files:
    for_list = for_in.split(",")
    rev_list = rev_in.split(",")
    out_file_list = []
    total_aligned_pairs = 0

    # Create the contig data dictionnary and hit from each alignments
    nb_alignment = len(for_list)
    contig_data, hit_data = mtn.create_contig_data(
        assembly, nb_alignment, depth_file, enzyme
    )

    for i in range(len(for_list)):
        for_in = for_list[i]
        try:
            rev_in = rev_list[i]
        except IndexError:
            rev_in = None
        name = "alignment_" + str(i)
        out_file = join(out_dir, "alignment_" + str(i) + ".pairs")

        # Align if necessary
        if start == "fastq":

            iterative = False

            # Digest reads if necessary.
            if aligner_mode == "cutsite":
                digest_for = join(tmp_dir, f"digest_for_{i}.fq.gz")
                digest_rev = join(tmp_dir, f"digest_rev_{i}.fq.gz")
                hcc.cut_ligation_sites(
                    fq_for=for_in,
                    fq_rev=rev_in,
                    digest_for=digest_for,
                    digest_rev=digest_rev,
                    enzyme=enzyme,
                    mode="for_vs_rev",
                    seed_size=20,
                    n_cpu=int(n_cpu),
                )
                for_in, rev_in = digest_for, digest_rev

            elif aligner_mode == "iterative":
                iterative = True

            if iterative or (aligner == "bowtie2"):
                # Create files to save the alignment.
                alignment_for = join(out_dir, name + "_for.bam")
                alignment_rev = join(out_dir, name + "_rev.bam")

                # Align the forward reads
                logger.info(f"Alignment of {for_in}:")
                align(for_in, index, aligner, alignment_for, n_cpu, iterative)

                # Align the reverse reads
                logger.info(f"Alignment of  {rev_in}:")
                align(
                    rev_in,
                    index,
                    aligner,
                    alignment_rev,
                    n_cpu,
                    tmp_dir,
                    iterative,
                )
            elif aligner == "bwa":
                # Create file to save the alignement.
                alignment = join(out_dir, name + ".bam")
                logger.info(f"Alignment of {for_in} and {rev_in}:")
                align(
                    for_in,
                    index,
                    aligner,
                    alignment,
                    n_cpu,
                    tmp_dir,
                    fq_in_2=rev_in,
                )

        elif start == "bam":
            if aligner == "bowtie2":
                logger.info(f"Processing {for_in} and {rev_in}:")
                alignment_for = for_in
                alignment_rev = rev_in
            elif aligner == "bwa":
                alignment = for_in

        else:
            logger.error("Start argument should be either 'fastq' or 'bam'.")
            raise ValueError

        if aligner == "bowtie2":
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
                f"{aligned_reads_for} forward reads aligned and {aligned_reads_rev} reverse reads aligned."
            )

            # Merge alignement to create a pairs file
            logger.info("Merging the pairs:")
            n_pairs = merge_alignment(
                alignment_temp_for, alignment_temp_rev, contig_data, out_file
            )

        # Case where a bam file from bwa is given as input.
        if aligner == "bwa":
            n_pairs = process_bwa_bamfile(
                alignment, min_qual, contig_data, out_file
            )

        logger.info(f"{n_pairs} pairs aligned.\n")
        total_aligned_pairs += n_pairs

    # Sort pairs.
    logger.info(f"Sort and indexed {out_file}")
    out_file = mio.sort_pairs_pairtools(
        out_file,
        threads=n_cpu,
        remove=True,
        force=True
    )
    out_file_list.append(out_file)

    if len(out_file_list) > 1:
        logger.info(f"TOTAL PAIRS MAPPED: {total_aligned_pairs}\n")

    return out_file_list, contig_data, hit_data


def merge_alignment(forward_aligned, reverse_aligned, contig_data, out_file):
    """Merge forward and reverse alignment into one file with only pairs which
    have both reads are aligned on the genome with 9 columns: ReadID, ContigA,
    Position_startA, Position_endA, StrandA, ContigB, Position_startB,
    Position_endB, StrandB.

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
    contig_data : dict
        Dictionary of the all the contigs from the assembly, the contigs names
        are the keys to the data of the contig available with the following
        keys: "id", "length", "GC", "hit", "coverage". Coverage still at 0 and
        need to be updated later.
    out_file : str
        Path to write the output pairs file.

    Returns
    -------
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

        # Write header of the pairs file.
        merged.write("## pairs format v1.0\n")
        merged.write("#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n")
        merged.write("#sorted: readID\n")
        merged.write("#shape: upper triangle\n")
        for contig in contig_data:
            merged.write(
                "#chromsize: {0} {1}\n".format(
                    contig, contig_data[contig]["length"]
                )
            )

        # Loop while at least one end of one end is reached. It's possible to
        # advance like that as the two tsv files are sorted on the id of the
        # reads.
        while n_pairs >= 0:
            # Case of both reads of the pair map.
            if for_read[0] == rev_read[0]:
                # Write read ID
                merged.write(for_read[0] + "\t")
                # Pairs are 1-based so we have to add 1 to 0 based bam position
                for_position = (
                    for_read[1] + "\t" + str(int(for_read[2]) + 1) + "\t"
                )
                rev_position = (
                    rev_read[1] + "\t" + str(int(rev_read[2]) + 1) + "\t"
                )

                # Have upper triangle shape
                if (
                    (
                        for_read[1] == rev_read[1]
                        and int(for_read[2]) <= int(rev_read[2])
                    )
                    or contig_data[for_read[1]]["id"]
                    < contig_data[rev_read[1]]["id"]
                ):

                    merged.write(
                        for_position
                        + rev_position
                        + for_read[3]
                        + "\t"
                        + rev_read[3]
                        + "\n"
                    )
                else:
                    merged.write(
                        rev_position
                        + for_position
                        + rev_read[3]
                        + "\t"
                        + for_read[3]
                        + "\n"
                    )
                n_pairs += 1
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
                names_sorted = sorted(names, key=LooseVersion)
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

    return n_pairs


def process_bamfile(alignment, min_qual, filtered_out):
    """Filter alignment BAM files

    Reads all the reads in the input BAM alignment file. Keep reads in the
    output if they are aligned with a good quality (greater than min quality
    threshold given) saving their only some columns: ReadID, Contig,
    Position_start, Position_end, strand to save memory.

    Parameters
    ----------
    alignment : str
        Path to the input temporary alignment.
    min_qual : int
        Minimum mapping quality required to keep a Hi-C pair.
    filtered_out : str
        Path to the output temporary tsv alignement.

    Returns
    -------
    int:
        Number of reads aligned.
    """

    # Check the quality and status of each aligned fragment.
    aligned_reads = 0
    save = pysam.set_verbosity(0)
    temp_bam = pysam.AlignmentFile(alignment, "rb", check_sq=False)
    pysam.set_verbosity(save)
    with open(filtered_out, "w") as f:
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
                        + "+"
                        + "\n"
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
                        + "-"
                        + "\n"
                    )
                    f.write(read)
    temp_bam.close()

    return aligned_reads


def process_bwa_bamfile(alignment, min_qual, contig_data, out_file):
    """Filter alignment BAM files

    Reads all the reads in the input BAM alignment file. Keep reads in the
    output if they are aligned with a good quality (greater than min quality
    threshold given) saving their only some columns: ReadID, Contig,
    Position_start, Position_end, strand to save memory.

    Parameters
    ----------
    alignment : str
        Path to the input temporary alignment.
    min_qual : int
        Minimum mapping quality required to keep a Hi-C pair.
    contig_data : dict
        Dictionnary of the all the contigs from the assembly, the contigs names
        are the keys to the data of the contig available with the following
        keys: "id", "length", "GC", "hit", "coverage". Coverage still at 0 and
        need to be updated later.
    out_file : str
        Path to the output pairs file.

    Returns
    -------
    int:
        Number of pairs aligned.
    """

    # Read the bam file.
    n_pairs = 0
    save = pysam.set_verbosity(0)
    temp_bam = pysam.AlignmentFile(alignment, "rb", check_sq=False)
    pysam.set_verbosity(save)

    with open(out_file, "w") as merged:

        # Write header of the pairs file.
        merged.write("## pairs format v1.0\n")
        merged.write("#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n")
        merged.write("#sorted: readID\n")
        merged.write("#shape: upper triangle\n")
        for contig in contig_data:
            merged.write(
                "#chromsize: {0} {1}\n".format(
                    contig, contig_data[contig]["length"]
                )
            )

        # Loop until the end of the file. Read the reads by two as the forward
        # and reverse reads should be interleaved.
        while n_pairs >= 0:
            try:
                for_read = next(temp_bam)
                while for_read.is_supplementary:
                    for_read = next(temp_bam)
                rev_read = next(temp_bam)
                while rev_read.is_supplementary:
                    rev_read = next(temp_bam)

                # Check mapping quality
                if (
                    for_read.mapping_quality >= min_qual
                    and rev_read.mapping_quality >= min_qual
                ):

                    # Check flag
                    if not (for_read.is_unmapped or rev_read.is_unmapped):
                        n_pairs += 1

                        # Safety check (forward and reverse are the same reads)
                        if for_read.query_name != rev_read.query_name:
                            logger.error(
                                f"Reads should be paired - {for_read.query_name}\t{rev_read.query_name}"
                            )
                            raise ValueError

                        # Define pairs value.
                        name = for_read.query_name
                        contig1 = for_read.reference_name
                        contig2 = rev_read.reference_name
                        pos1 = for_read.pos + 1
                        pos2 = for_read.pos + 1
                        strand1 = "+"
                        strand2 = "+"
                        if for_read.is_reverse:
                            strand1 = "-"
                        if rev_read.is_reverse:
                            strand2 = "-"

                        # Modify order to have an upper triangle and write
                        # the pair.
                        if (contig1 == contig2 and pos1 <= pos2) or contig_data[
                            contig1
                        ]["id"] < contig_data[contig2]["id"]:
                            merged.write(
                                "\t".join(
                                    [
                                        name,
                                        contig1,
                                        str(pos1),
                                        contig2,
                                        str(pos2),
                                        strand1,
                                        strand2,
                                    ]
                                )
                                + "\n"
                            )
                        else:
                            merged.write(
                                "\t".join(
                                    [
                                        name,
                                        contig2,
                                        str(pos2),
                                        contig1,
                                        str(pos1),
                                        strand2,
                                        strand1,
                                    ]
                                )
                                + "\n"
                            )

            # Exit the loop if no more reads.
            except StopIteration:
                break

    # Close the bam file and return number of pairs
    temp_bam.close()
    return n_pairs
