#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Generates quality metrics from the output of MetaTOR.

General utilities functions to extract high  quality contigs which will be used
to estimate the quality of the metaHiC librairies nad generates some plots to
illustrate it.

From a .metacool contact map instead of the pairs files, the noise of the
library is measured, and only it: the contacts within contigs, within MAGs
and between MAGs, the noise ratio and the noise score, which are those of
the pairs that the map holds. The events of the pairs within a contig need
the positions and strands of their reads, which a map does not hold.

Core functions to assess the quality:
    - extract_hq_contigs
    - extract_pairs, extract_map_pairs
    - hic_quality, map_quality
    - quality_check
"""

import hicstuff.digest as hcd
import hicstuff.filter as hcf
import metator.figures as mtf
import metator.io as mio
import metator.metacool as mtmc
import numpy as np
from Bio import SeqIO
from metator.log import logger
from os.path import join

# What a map does not give, which the log says when the quality is measured
# from one.
NOT_FROM_A_MAP = (
    "The religated, loop and weirds ratios, the informative contacts estimation and the ratio inter/intra contigs are "
    "not measured from a contact map: they need the positions and strands of the reads, which only pairs files hold."
)


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
        output_pairs.write("#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n")
        for contig in contigs:
            output_pairs.write(
                "#chromsize: {0} {1}\n".format(
                    contig,
                    contigs_data.loc[contig, "Size"],
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
    logger.info(f"{n_pairs} pairs have been extracted.")
    return n_pairs


def extract_map_pairs(metacool, contigs, contigs_data, contig_data_file):
    """The pairs between given contigs in a .metacool contact map: the
    counterpart of extract_pairs, without the reads. The map holds the
    number of pairs between each two bins: those between the bins of two
    contigs, or within the bins of one, are summed (see
    metator.metacool.contig_contacts). They are the pairs that extract_pairs
    extracts from the pairs files that the map was made of: each pair whose
    two reads lie on the contigs, within their length (a read at position 0
    counts on the first base).

    Parameters
    ----------
    metacool : str
        Path to the .metacool file, or URI of an entry of a .mmetacool
        collection: ``file.mmetacool::entries/<name>``.
    contigs : List of str
        List of contigs name to extract.
    contigs_data : pandas.DataFrame
        Table with contig information from metator pipeline: the lengths of
        the contigs (Size), which must be theirs in the map.
    contig_data_file : str
        The file of contigs_data, for the messages.

    Returns
    -------
    pandas.DataFrame
        One row per pair of contigs in contact: contig1, contig2 (the same
        contig for the pairs within a contig) and count, their pairs.

    Raises
    ------
    metator.metacool.MetacoolError
        If the map cannot be read, stores both triangles of the matrix, or
        gives one of the contigs another length than contigs_data. The
        contigs that it lacks have no pairs, with a warning.
    """
    sizes = {name: int(contigs_data.loc[name, "Size"]) for name in contigs}
    contacts = mtmc.contig_contacts(metacool, sizes, contig_data_file, "the contigs longer than 100 kb of the high-quality MAGs")
    logger.info(f"{int(contacts['count'].sum())} pairs have been extracted from {metacool}.")
    return contacts


def parse_fasta(fasta):
    """Iterate over the records of a FASTA file, compressed or not.

    Parameters
    ----------
    fasta : str
        Path to the fasta sequences.

    Returns
    -------
    Iterator of Bio.SeqRecord.SeqRecord:
        The records of the file.

    Raises
    ------
    ValueError:
        If the file does not start with a '>' header line. Biopython >= 1.87
        raises it on any text before the first record, blank or '#' lines
        included; it is logged here as an error on the input file.
    UnicodeDecodeError:
        If the file holds a byte that cannot be decoded. It is a ValueError
        too, but says nothing about the header, so it is passed on as is.
    """
    try:
        yield from SeqIO.parse(mio.read_compressed(fasta), "fasta")
    except UnicodeDecodeError:
        raise
    except ValueError as err:
        logger.error(
            f"{fasta} is not a FASTA file: it must start with a '>' header line, with no blank or '#' comment line before it."
        )
        raise ValueError(f"Not a FASTA file: {fasta}") from err


def _percent(part, whole):
    """part as a percentage of whole, 0 when whole is 0."""
    return 100 * part / whole if whole else 0


def noise_metrics(n_mags, n_intra_mags, n_inter_mags):
    """The noise of a library: the noise contact ratio, the pairs between
    two MAGs in percent of the pairs within or between the MAGs, and the
    noise score, x / (y + x), with x the pairs between MAGs per pair of
    MAGs and y the pairs within MAGs per MAG. Both are 0 without a pair,
    and the score is 0 with fewer than two MAGs.

    Parameters
    ----------
    n_mags : int
        Number of MAGs selected as high-quality.
    n_intra_mags : int
        The pairs within the MAGs: within a contig, or between two contigs
        of the same MAG.
    n_inter_mags : int
        The pairs between two MAGs.

    Returns
    -------
    tuple of float
        The noise contact ratio, in percent, and the noise score.
    """
    if n_intra_mags + n_inter_mags == 0:
        return 0, 0
    noise_ratio = 100 * n_inter_mags / (n_inter_mags + n_intra_mags)
    if n_mags < 2:
        return noise_ratio, 0
    pairs_of_mags = n_mags * (n_mags - 1) * 0.5
    return noise_ratio, (n_inter_mags / pairs_of_mags) / (n_intra_mags / n_mags + n_inter_mags / pairs_of_mags)


def _log_contacts(n_mags, n_intra_contigs, n_intra_mags, n_inter_mags):
    """Log the MAGs and their pairs: those within contigs, within MAGs (a
    contig, or two contigs of a MAG) and between MAGs."""
    logger.info(f"Number of MAGs: {n_mags}")
    logger.info(f"Contacts within contigs: {n_intra_contigs}.")
    logger.info(f"Contacts within MAGs: {n_intra_mags}.")
    logger.info(f"Contacts between MAGs: {n_inter_mags}.")


def _log_noise(n_mags, n_intra_mags, n_inter_mags):
    """Log the noise contact ratio and the noise score (see noise_metrics)."""
    noise_ratio, noise_score = noise_metrics(n_mags, n_intra_mags, n_inter_mags)
    logger.info(f"Noise contact ratio: {noise_ratio:.2f}%")
    logger.info(f"Noise score: {noise_score:.2E}")


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
    for record in parse_fasta(fasta):
        # Get chromosome restriction table
        restrict_table[record.id] = hcd.get_restriction_table(record.seq, enzyme, circular=False)

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
    logger.info("Filtering with thresholds: uncuts={0} loops={1}".format(uncut_thr, loop_thr))
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
    if n_intra_mags + n_inter_mags > 0:
        rat_info = 100 * (n_informative_intra + n_informative_inter) / (n_intra_mags + n_inter_mags)
    else:
        logger.warning("No pairs have been extracted. All scores set to 0.")
        rat_info = 0
    # No pair at all (or none within the MAGs, or none informative): the
    # scores and ratios that would divide by 0 are 0.
    n_informative = n_informative_intra + n_informative_inter
    _log_contacts(n_mags, n_intra_contigs, n_intra_mags, n_inter_mags)
    logger.info(f"Religated ratio: {_percent(n_religated, n_intra_mags):.2f}%.")
    logger.info(f"Loop ratio: {_percent(n_loops, n_intra_mags):.2f}%.")
    logger.info(f"Weirds ratio: {_percent(n_weirds, n_intra_mags):.2f}%.")
    logger.info(f"Informative contacts estimation: {rat_info:.2f}%.")
    logger.info(f"Ratio inter/intra contigs: {_percent(n_informative_inter, n_informative):.2f}%.")
    _log_noise(n_mags, n_intra_mags, n_inter_mags)

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


def map_quality(contigs, n_mags, contacts):
    """The noise of the library, from the pairs that a .metacool contact
    map holds between the contigs of the high-quality MAGs (see
    extract_map_pairs): the pairs within contigs, within MAGs (a contig, or
    two contigs of a MAG) and between MAGs, the noise contact ratio and the
    noise score, logged as hic_quality logs them from the pairs files, which
    give the same numbers. The events of the pairs within a contig (uncuts,
    loops, weirds and informative pairs) need the positions and strands of
    their reads, which a map does not hold: neither they nor the ratios
    made of them are measured, which the log says.

    Parameters
    ----------
    contigs : dict
        Dictionnary with the list of contigs as keys and the final MAGs as
        values.
    n_mags : int
        Number of MAGs selected as high-quality.
    contacts : pandas.DataFrame
        The pairs between the contigs, by pair of contigs: contig1, contig2
        and count (see extract_map_pairs).

    Returns
    -------
    tuple of 3 int:
        The pairs within contigs, within MAGs and between MAGs.
    """
    count = contacts["count"].to_numpy()
    same_contig = (contacts["contig1"] == contacts["contig2"]).to_numpy()
    same_mag = (contacts["contig1"].map(contigs) == contacts["contig2"].map(contigs)).to_numpy()
    n_intra_contigs = int(count[same_contig].sum())
    n_intra_mags = int(count[same_mag].sum())
    n_inter_mags = int(count[~same_mag].sum())
    _log_contacts(n_mags, n_intra_contigs, n_intra_mags, n_inter_mags)
    _log_noise(n_mags, n_intra_mags, n_inter_mags)
    logger.info(NOT_FROM_A_MAP)
    return n_intra_contigs, n_intra_mags, n_inter_mags


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
    metacool=None,
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
    metacool : str, optional
        A .metacool contact map of the run (or an entry of a .mmetacool
        collection), read instead of the pairs files: the noise of the
        library only (see map_quality), which is that of the pairs the map
        holds. fasta_file, enzyme, threshold and plot are then not used.

    Raises
    ------
    metator.metacool.MetacoolError
        If the map cannot be read, stores both triangles of the matrix, or
        gives a contig of the high-quality MAGs another length than the
        contig data.
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
    # Extract pairs, from the pairs files or from the map.
    contacts = None
    if metacool is None:
        n_pairs = extract_pairs(pairs_files, pairs, list(hq_contigs.keys()), contigs_data)
    else:
        contacts = extract_map_pairs(metacool, list(hq_contigs.keys()), contigs_data, contig_data_file)
        n_pairs = int(contacts["count"].sum())

    # Estimate HiC quality: only the noise from a map.
    if n_pairs > 0 and contacts is not None:
        map_quality(hq_contigs, n_mags, contacts)
    elif n_pairs > 0:
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
    elif n_mags == 0:
        logger.info("No high-quality MAG (weighted completeness above 70%, weighted redundancy under 1.15): no quality check.")
    elif not hq_contigs:
        logger.info(f"The {n_mags} high-quality MAGs have no contig longer than 100 kb: no quality check.")
    else:
        logger.info("No pair between the contigs longer than 100 kb of the high-quality MAGs: no quality check.")
    # Plot normalised contact map.
    # plot_contact_map()
    return 0
