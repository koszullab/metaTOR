#!/usr/bin/env python3
# coding: utf-8

"""Module with the mges bacterial host detection functions. 

It detects mges bacterial host MAGs from the output binning of MetaTOR and the
metaHiC network. The module extract the subnetwork and then identify the bins 
which interaction with the mge to rank potential bacterial host.

A core class Subnetwork is build to handle the mge contigs interaction.

Main function to call the class and build the output is host_detection.
"""


import networkx
import numpy as np
import pandas as pd
from metator.log import logger
from typing import List


class Subnetwork:
    """Class to handle subnetwork information of one bin."""

    def __init__(self, subnetwork, threshold):
        """Initialize nodes list and their weights (number of hits with the
        associated contigs normalized by metaTOR).

        Parameters:
        -----------
        subnetwork : networkx.classes.reportviews.EdgeDataView
            Subnetwork with all edges of a given contig.
        threshold : float
            Threshold to consider an association.
        """
        self.nodes = []
        self.weights = []
        for edge in subnetwork:
            self.id = edge[0]
            self.nodes.append(edge[1])
            self.weights.append(edge[2])
        self.len = len(self.weights)
        self.scored = False
        self.threshold = threshold

    def setScore(self):
        """Set scores for each associated contigs. The score is just a ratio of
        the sum of the weigths.
        """
        self.score = [x / sum(self.weights) for x in self.weights]

    def getMaxScore(self):
        """Return highest score and the associated contig.

        Returns:
        --------
        str:
            Contig name with the highest score (maximum connectivity with the
            given contig).
        float:
            Highest score value (contig score).
        """
        return self.nodes[self.score.index(max(self.score))], max(self.score)

    def setBinScore(self, contig_data):
        """Set scores for each associated bins.

        Parameters:
        -----------
        contig_data : pandas.DataFrame:
            Table with the contig information given with more column with the given
            anvio binning and the mge annotation.
        """
        self.bins = dict()
        for i in range(self.len):
            node = self.nodes[i] - 1
            score = self.score[i]
            # Uncomment to add virus as bins.
            # if contig_data[node]["virus"]:
            #     try:
            #         self.bins[contig_data[node]["name"]]["score"] += score
            #     except KeyError:
            #         self.bins[contig_data[node]["name"]] = {"score" : score}
            if (
                contig_data.loc[node, "Binned"]
                and not contig_data.loc[node, "MGE"]
            ):
                try:
                    self.bins[contig_data.loc[node, "Final_bin"]][
                        "score"
                    ] += score
                except KeyError:
                    self.bins[contig_data.loc[node, "Final_bin"]] = {
                        "score": score
                    }
        self.scored = True

    def getMaxBinScore(self, contig_data=None):
        """Return highest score and the associated bin.

        Returns:
        --------
        str:
            Bin name with the highest score (maximum connectivity with the given
            contig).
        float:
            Highest score value (bin score).
        """
        if not self.scored:
            self.setBinScore(contig_data)
        max_score = 0
        max_bin = "-"
        for bin_name in self.bins:
            score = self.bins[bin_name]["score"]
            if score > max_score:
                max_score = score
                max_bin = bin_name
        return max_bin, max_score

    def getBinScore(self):
        """Return the count of connected bin based on a threshold.

        Returns:
        --------
        int:
            Count of connected bins.
        """
        c = 0
        for bin_name in self.bins:
            score = self.bins[bin_name]["score"]
            if score >= self.threshold:
                c += 1
        return c

    def getScoreList(self):
        """Return the list of connected bin with respective scores.

        Returns:
        --------
        str:
            List of bins with scores separated by ';' (bin_name|score).
        """
        score_list = []
        if len(self.bins) > 1:
            for bin_name in self.bins:
                score_list.append(
                    f'{bin_name}|{self.bins[bin_name]["score"]:.1e}'
                )
            return ";".join(score_list)
        else:
            return "NA"


def associate_bin(
    bin_contigs: dict,
    network: "networkx.classes.graph.Graph",
    contig_data: "pandas.DataFrame",
    threshold: float,
) -> dict:
    """Function to associate one bin to one MAG.

    Parameters
    ----------
    bin_contigs : dict
        Dictionary with the name of the contigs.
    network : networkx.classes.graph.Graph
        MetaTOR network of the HiC data.
    contig_data : dict
        Dictionary with the contig name as keys and with the values of the
        contig id the associated bin name, and either if the contig is binned
        and if it's a mge contig. The name of the keys are "id", "bin",
        "binned", "mge".
    threshold : float
        Threshold to consider an association.

    Return
    ------
    dict :
        Updated dictionary with the associated MAG.
    """
    # Transform from contig to id
    contigs_id = [
        list(contig_data.index[contig_data.Name == name])[0]
        for name in bin_contigs["Contigs"]
    ]
    network_id = [i + 1 for i in contigs_id]
    # Do not compute subnetwork with no HiC contacts.
    if (
        np.sum([contig_data.loc[contig_id, "Hit"] for contig_id in contigs_id])
        > 0
    ):
        # Manage the case of the self interacting contigs.
        try:
            subnetwork = network.edges(network_id, data="weight")
            # Remove edges inside the bin. Edge 0 is always the asked bins.
            subnetwork = [
                edge for edge in subnetwork if (edge[1] not in network_id)
            ]
            if len(subnetwork) > 0:
                sub = Subnetwork(subnetwork, threshold)
                sub.setScore()
                sub.setBinScore(contig_data)
                bin_name, score = sub.getMaxBinScore()
                count = sub.getBinScore()
                if count == 0:
                    bin_contigs["Bin"] = "None"
                    bin_contigs["AssociationScore"] = score
                    bin_contigs["BinList"] = sub.getScoreList()
                elif count == 1:
                    bin_contigs["Bin"] = bin_name
                    bin_contigs["AssociationScore"] = score
                    bin_contigs["BinList"] = sub.getScoreList()
                else:
                    bin_contigs["Bin"] = "Multiple"
                    bin_contigs["AssociationScore"] = score
                    bin_contigs["BinList"] = sub.getScoreList()
            else:
                bin_contigs["Bin"] = "None"
                bin_contigs["AssociationScore"] = np.nan
                bin_contigs["BinList"] = "NA"
        # No contacts with others contigs.
        except networkx.exception.NetworkXError:
            bin_contigs["Bin"] = "None"
            bin_contigs["AssociationScore"] = np.nan
            bin_contigs["BinList"] = "NA"
    else:  # No HiC coverage.
        bin_contigs["Bin"] = "None"
        bin_contigs["AssociationScore"] = np.nan
        bin_contigs["BinList"] = "NA"
    return bin_contigs


def host_detection(
    network, contig_data, mges_list, mges_list_id, outfile, threshold
):
    """Main function to detect the host.

    Parameters:
    -----------
    network : networkx.classes.graph.Graph
        MetaTOR network of the HiC data.
    contig_data : pandas.DataFrame
        Table with the contig name as keys and with the values of the
        contig id the associated bin name, and either if the contig is binned
        and if it's a mge contig. The name of the keys are "id", "bin",
        "binned", "mge".
    mges_list : list
        List of the mges contigs names.
    mges_list_id : list
        List of the mges contigs IDs.
    outfile : str
        Path to write the output table.
    threshold : float
        Threshold to consider an association.

    Returns:
    --------
    dict:
        Dictionary with the mge contig as keys and with the associated host
        as value.
    """
    # Compute the score with the subnetwork and return bins in each categories
    # and build the associated table.
    mge_data = pd.DataFrame(contig_data.loc[mges_list_id, :])
    mge_data["Host"] = "ND"
    A, B, C = 0, 0, 0
    for contig_id in mges_list_id:
        network_id = contig_id + 1  # 1-based in network
        contig = contig_data.loc[contig_id, "Name"]
        # Do not compute subnetwork with no HiC contacts.
        if contig_data.loc[contig_id, "Hit"] > 0:
            # Manage the case of the self interacting contigs.
            try:
                subnetwork = network.edges(network_id, data="weight")
                sub = Subnetwork(subnetwork, threshold)
                sub.setScore()
                sub.setBinScore(contig_data)
                bin_name, score = sub.getMaxBinScore()
                count = sub.getBinScore()
                if count == 0:
                    C += 1
                    mge_data.loc[contig_id, "Host"] = "None"
                elif count == 1:
                    A += 1
                    mge_data.loc[contig_id, "Host"] = bin_name
                else:
                    B += 1
                    mge_data.loc[contig_id, "Host"] = "Multiple"
            # No contacts with others contigs.
            except networkx.exception.NetworkXError:
                C += 1
                mge_data.loc[contig_id, "Host"] = "None"
        else:  # No HiC coverage.
            C += 1
            mge_data.loc[contig_id, "Host"] = "None"

    logger.info("{0} mges associated with one bin.".format(A))
    logger.info("{0} mges associated with more than one bin.".format(B))
    logger.info("{0} mges associated with no bin.".format(C))

    # Save the host table.
    mge_data.to_csv(outfile, sep="\t")
    return mge_data
