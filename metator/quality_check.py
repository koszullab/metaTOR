#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Generates quality metrics from the output of MetaTOR.

General utilities functions to extract high  quality contigs which will be used 
to estimate the quality of the metaHiC librairies nad generates some plots to 
illustrate it.

Core functions to assess the quality:
    - extract_hq_contigs
    - extract_pairs
    - hic_quality
    - quality_check
"""

import hicstuff.digest as hcd
import hicstuff.filter as hcf
import hicstuff.io as hio
import metator.figures as mtf
import metator.io as mio
import numpy as np
import pypairix
from Bio import SeqIO
from metator.log import logger
from os.path import join


def extract_hq_contigs(bin_summary, contigs_data):
    """Function to extract the high quality contigs from the metator output.
    These contigs will be the one used to assess the quality of the dataset.

    Parameters
    ----------
    bin_summary : pandas.DataFrame
        Table with the final bin metrics from metator pipeline.
    contigs_data : pandas.DataFrame
        Table with contig information from metator pipeline.

    Returns
    -------
    dict:
        List of the names of the high quality contigs as keys and the final bins
        as values.
    """
    hq_contigs = {}
    # Extract high quality MAGs.
    hq_mags = bin_summary.index[
        np.logical_and(
            (bin_summary["Weighted completeness"] > 0.7),
            (bin_summary["Weighted redundancy"] < 1.15),
        )
    ]
    n_mags = len(hq_mags)
    # Extract contigs bigger than 100kb.
    large_contigs = contigs_data.index[contigs_data["Size"] > 100_000]
    # Build dictionnary of large contigs in HQ MAgs.
    for contig in large_contigs:
        if contigs_data.loc[contig, "Final_bin"] in hq_mags:
            hq_contigs[contig] = contigs_data.loc[contig, "Final_bin"]
    return hq_contigs, n_mags


def extract_pairs(pairs_files, out_file, contigs, contigs_data):
    """Function to extract the pairs fo given contigs from a pairs file and
    write them in a new file.

    Parameters
    ----------
    pairs_files : List of str
        List of path to pairs files. The can have a pypairix index or not.
    out_file : str
        Path where to write the ouput file.
    contigs : List of str
        List of contigs name to extract.
    contigs_data : pandas.DataFrame
        Table with contig information from metator pipeline.
    """
    # Initiation
    n_pairs = 0

    # Write one pair file for all the ones given.
    with open(out_file, "w") as output_pairs:
        # Write the header of the output pairs
        output_pairs.write("## pairs format v1.0\n")
        output_pairs.write(
            "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n"
        )
        for contig in contigs:
            output_pairs.write(
                "#chromsize: {0} {1}\n".format(
                    contig, contigs_data.loc[contig, "Size"],
                )
            )
        for pairs_file in pairs_files:
            pairs_data = mio.get_pairs_data(pairs_file)
            for contig_id1, contig1 in enumerate(contigs):
                # Only need to retrieve the upper triangle.
                for contig2 in contigs[contig_id1:]:
                    pairs_lines = pairs_data.query2D(
                        contig1,
                        0,
                        contigs_data.loc[contig1, "Size"],
                        contig2,
                        0,
                        contigs_data.loc[contig2, "Size"],
                        1,
                    )
                    for pairs_line in pairs_lines:
                        n_pairs += 1
                        output_pairs.write("\t".join(pairs_line[:7]) + "\n")
    logger.info(f"{n_pairs} have been extracted.")
    return n_pairs


def hic_quality(
    contigs,
    n_mags,
    pairs,
    pairs_idx,
    fasta,
    enzyme,
    prefix=None,
    plot=False,
    plot_event=None,
    plot_cam=None,
    threshold=None,
):
    """Function to evaluate the quality of the HiC in the samples. It will
    evaluate the proportion of informative contacts based on the high quality
    contigs metrics.

    Parameters
    ----------
    contigs : dict
        Dictionnary with the list of contigs as keys and the final MAGs as
        values.
    n_mags : int
        Number of MAGs selected as high-quality.
    pairs : str
        Path to the pairs file whith mapping positions of the reads.
    pairs_idx : str
        Path to the file where the pairs with frags indexed will be written.
    fasta : str
        Path to the fasta sequences.
    enzyme : List of str
        List of the enzyme use in the metaHiC experiment.
    prefix : str
        Prefix name of the sample. Default None.
    plot : str
        If True, output some plots.
    plot_event : str
        Path to save the plot with event distribution. Default None.
    plot_cam : str
        Path to save the camembert plot of event distribution. Default None.
    threshold : Tuple of int
        Threshold of religated and loop. By default automatic computation is
        done. Default None.

    Returns
    -------
    List of 6 int:
        Religated loops, weirds, informative, intra_mags, inter_mags counts.
    """
    # Cut the genome on restriction sites to assess quality of the reads.
    restrict_table = {}
    for record in SeqIO.parse(mio.read_compressed(fasta), "fasta"):
        # Get chromosome restriction table
        restrict_table[record.id] = hcd.get_restriction_table(
            record.seq, enzyme, circular=False
        )

    # Add fragment index to pairs (readID, chr1, pos1, chr2,
    # pos2, strand1, strand2, frag1, frag2)
    hcd.attribute_fragments(pairs, pairs_idx, restrict_table)

    if threshold is not None:
        # Thresholds supplied by user beforehand
        uncut_thr = int(threshold[0])
        loop_thr = int(threshold[1])
    else:
        # Threshold defined at runtime
        uncut_thr, loop_thr = hcf.get_thresholds(
            pairs_idx,
            interactive=False,
            plot_events=plot,
            fig_path=plot_event,
            prefix=prefix,
        )
        logger.info(
            "Filtering with thresholds: uncuts={0} loops={1}".format(
                uncut_thr, loop_thr
            )
        )
    # Filter reads and save metrics on informative reads
    n_religated = 0
    n_loops = 0
    n_weirds = 0
    n_informative_intra = 0
    n_informative_inter = 0
    n_intra_contigs = 0
    n_intra_mags = 0
    n_inter_mags = 0

    # open the files for reading and writing
    with open(pairs_idx, "r") as pairs:
        for line in pairs:  # iterate over each line.
            # Skip header lines.
            if line.startswith("#"):
                continue
            # Iterates on pairs.
            p = hcf.process_read_pair(line)
            # Count events.
            if p["chr1"] == p["chr2"]:
                n_intra_contigs += 1
                n_intra_mags += 1
                # Do not report ++ and -- pairs on the same fragment (impossible)
                if p["frag1"] == p["frag2"] and p["strand1"] == p["strand2"]:
                    n_weirds += 1
                elif p["nsites"] <= loop_thr and p["type"] == "-+":
                    n_loops += 1
                elif p["nsites"] <= uncut_thr and p["type"] == "+-":
                    n_religated += 1
                else:
                    n_informative_intra += 1
            else:
                if contigs[p["chr1"]] == contigs[p["chr2"]]:
                    n_intra_mags += 1
                    n_informative_inter += 1
                else:
                    n_inter_mags += 1
    rat_info = (
        100
        * (n_informative_intra + n_informative_inter)
        / (n_intra_mags + n_inter_mags)
    )
    noise_ratio = 100 * n_inter_mags / (n_inter_mags + n_intra_mags)
    if n_mags > 1:
        noise_score = (n_inter_mags / (n_mags * (n_mags - 1) * 0.5)) / (
            n_intra_mags / n_mags + n_inter_mags / (n_mags * (n_mags - 1) * 0.5)
        )
    else:
        noise_score = 0
    logger.info(f"Number of MAGs: {n_mags}")
    logger.info(f"Religated ratio: {100 * n_religated / n_intra_mags:.2f}%.")
    logger.info(f"Loop ratio: {100 * n_loops / n_intra_mags:.2f}%.")
    logger.info(f"Weirds ratio: {100 * n_weirds / n_intra_mags:.2f}%.")
    logger.info(f"Informative contacts estimation: {rat_info:.2f}%.")
    logger.info(
        f"Ratio inter/intra contigs: {n_informative_inter / (n_informative_intra + n_informative_inter):.2f}%."
    )
    logger.info(f"Noise contact ratio: {noise_ratio:.2f}%")
    logger.info(f"Noise score: {noise_score:.2E}")

    if plot:
        mtf.figure_camembert_quality(
            plot_cam,
            prefix,
            n_religated,
            n_loops,
            n_weirds,
            n_informative_intra,
            n_informative_inter,
            n_intra_mags,
            n_inter_mags,
            uncut_thr,
            loop_thr,
        )
    return (
        n_religated,
        n_loops,
        n_weirds,
        n_informative_intra,
        n_informative_inter,
        n_intra_mags,
        n_inter_mags,
    )


def quality_check(
    contig_data_file,
    bin_summary_file,
    fasta_file,
    pairs_files,
    out_dir,
    tmp_dir,
    prefix,
    plot,
    enzyme,
    threshold,
):
    """Main function to compute the quality of the metaHiC library and to
    display some metrics about it.

    Parameters
    ----------
    contig_data_file : str
    bin_summary_file : str
    fasta_file : str
    pairs_files : List of str
    out_dir : str
    tmp_dir : str
    prefix : str
    plot : bool
    enzyme : List of str
    threshold : Tuple of int
    """
    # Import contigs data.
    contigs_data = mio.read_contig_data(contig_data_file)
    # Import bin summary.
    bin_summary = mio.read_bin_summary(bin_summary_file)

    # Define temporary and output files.
    if plot:
        plot_event = join(out_dir, f"{prefix}_event.pdf")
        plot_cam = join(out_dir, f"{prefix}_camembert_plot.pdf")
    else:
        plot_event = None
        plot_cam = None

    pairs = join(tmp_dir, f"{prefix}.pairs")
    pairs_idx = join(tmp_dir, f"{prefix}_idx.pairs")

    # Extract high quality contigs.
    hq_contigs, n_mags = extract_hq_contigs(bin_summary, contigs_data)
    # Extract pairs
    n_pairs = extract_pairs(
        pairs_files, pairs, list(hq_contigs.keys()), contigs_data
    )

    # Estimate HiC quality.
    if n_pairs > 0:
        _ = hic_quality(
            hq_contigs,
            n_mags,
            pairs,
            pairs_idx,
            fasta_file,
            enzyme,
            prefix=prefix,
            plot=plot,
            plot_event=plot_event,
            plot_cam=plot_cam,
            threshold=threshold,
        )
    else:
        logger.info("No MAGs have been found.")
    # Plot normalised contact map.
    # plot_contact_map()
    return 0
