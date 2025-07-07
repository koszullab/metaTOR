#!/usr/bin/env python3
# coding: utf-8
"""
This module is used to detect and associate an MGE with its bacterial MAG host.

It defines two main classes, `MAG` and `MgeMag`, along with their methods.
Starting from a `MAG` or `MgeMag` object, the module retrieves all contigs belonging to a given MAG or MGE-MAG based on binning results.

The module allows evaluation of background noise by calculating intra- and inter-MAG interaction signals, as well as interactions between MGE-MAGs and MAGs. This helps identify the host MAG for a given MGE-MAG.

Main functions and methods include:
    - Class methods:
        - `__init__()`: Initializes a MAG or MGE-MAG object.
        - `add_contig()`: Adds a contig to a MAG or MGE-MAG.
        - `has_contig()`: Checks if a MAG or MGE-MAG contains a specific contig.
        - `list_all_mags()` and `list_all_mge_mags()`: List all existing MAG and MGE-MAG objects.
    - `create_mags()` and `create_mge_mags()`: Instantiate MAG and MGE-MAG objects from input data.
    - `evaluate_experience_noise()`: Calculates intra-MAG and inter-MAG interaction signals, writes them to a file, and plots a normalized log10 boxplot.
       Also compares the distributions of intra- and inter-MAG interactions on a violin plot by quality of mags.
    - `compute_mge_mag_interactions()`: Compute MGE-MAG to MAG interaction signals with signal filtering by percentage and minimum contact threshold.
"""

# import checkv
# import metator.figures as mtf
# import metator.io as mio
from collections import Counter
from typing import Literal
import networkx as nx
import matplotlib.pyplot as plt
from metator.log import logger
import pandas as pd
import numpy as np
import seaborn as sns
import math
from typing import Union

class Bin:
    """
        Classe de base représentant un ensemble de contigs appartenant à un même groupe (ex: MAG ou MGE-MAG).
    """

    def __init__(self, name):
        self.name: str = name
        self.info: dict = {}
        self.contigs: dict = {}
        self.completeness: float = None
        self.contamination: float = None
        self.quality: str = None
        self.intra_signal: float = None
        self.inter_signals: dict = {}

    def add_contig(self, contig_name, contig) -> None:
        self.contigs[contig_name] = contig

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{self.name}', {len(self.contigs)} contigs)"


class Mag(Bin):
    """
        Représente un MAG (Metagenome-Assembled Genome). Hérite de Bin.
    """

    def __init__(self, name):
        super().__init__(name)
        self.mges: dict = {}

    def add_mge(self, mge_name, mge_mag) -> None:
        self.mges[mge_name] = mge_mag

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{self.name}', {len(self.contigs)} contigs, {len(self.mges)} MGEs)"


class MgeMag(Bin):
    """
        Représente un MGE-MAG (intégrant des éléments génétiques mobiles). Hérite de Bin.
    """

    def __init__(self, name):
        super().__init__(name)
        self.hosts: dict = {}

    def add_host(self, mag_name, mag) -> None:
        self.hosts[mag_name] = mag

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{self.name}', {len(self.contigs)} contigs, {len(self.hosts)} hosts)"


def update_final_bin_with_mge(contig_data: pd.DataFrame, mges_bin_summary: pd.DataFrame) -> pd.DataFrame :
    """
    Updates the 'Final_bin' column of df1 by replacing 'ND' values or metator_ MAGs with the corresponding MGEs (from df2), if available.

    Args:

        contig_data (pd.DataFrame): DataFrame containing at least the columns ['Name', 'Final_bin']

        mges_bin_summary (pd.DataFrame): DataFrame containing the columns ['Contigs', 'BinName']

    Returns:

        pd.DataFrame: An updated copy of df1 with the 'Final_bin' column modified when applicable
    """

    contig_data = contig_data.copy()
    mges_bin = mges_bin_summary.copy()

    # Construire dictionnaire contig → MGE
    mges_bin = mges_bin.assign(Contigs=mges_bin["Contigs"].str.split(",")).explode("Contigs")
    mges_bin["Contigs"] = mges_bin["Contigs"].str.strip()
    contig_to_mge = dict(zip(mges_bin["Contigs"], mges_bin["BinName"]))

    def update(row):
        name = row.get("Name")
        final_bin = row.get("Final_bin")
        mge = contig_to_mge.get(name)

        if pd.isna(final_bin):
            final_bin = "ND"  # pour que les NaN soient aussi traités comme "ND"

        if final_bin == "ND" and mge:
            return mge
        elif isinstance(final_bin, str) and final_bin.startswith("metator_") and mge:
            return mge
        else:
            return final_bin

    contig_data["Final_bin"] = contig_data.apply(update, axis=1)
    return contig_data



def create_bins(contig_data: pd.DataFrame, mges_bin_summary: pd.DataFrame) -> tuple[dict, dict]:
    """
    Creates MGE-MAGs from contig data.

    Args:
        contig_data (pd.DataFrame): DataFrame containing information about contigs.

    Returns:
        dict: A dictionary mapping MGE-MAG names to their corresponding MgeMag objects.
    """

    mags = {}
    mge_mags = {}
    # Loading data and updating the contig_data_final  to a full version with MAGs and MGEs
    logger.info("Loading data and generating the contig_data")
    contig_data = update_final_bin_with_mge(contig_data, mges_bin_summary)

    for _, row in contig_data.iterrows():
        contig_id, contig_name, bin_name = row["ID"], row["Name"], row["Final_bin"]

        if bin_name.startswith("MetaTOR_MGE_"):  # Identify MGE-MAGs
            if bin_name not in mge_mags:
                mge_mags[bin_name] = MgeMag(bin_name)  # create a new MGE-MAG
            mge_mags[bin_name].add_contig(contig_name, contig_id)

        elif bin_name.startswith("metator_"):
            if bin_name not in mags:
                mags[bin_name] = Mag(bin_name)  # create a new MGE-MAG
            mags[bin_name].add_contig(contig_name, contig_id)
    logger.info(f"{len(mags)} MAGs and {len(mge_mags)} MGE-MAGs have been parsed.")

    return mags, mge_mags


def map_contigs_to_mags(mags: dict) -> dict:
    """
    Creates a mapping from each contig to its corresponding MAG name.

    Iterates through a dictionary of MAG objects and associates each contig with the name of the MAG it belongs to.

    Args:
    - mags (dict): Dictionary in the form {mag_name: mag_object}, where each mag_object has a `contigs` attribute
      that returns a dictionary of contigs.

    Returns:
    - contig_to_mag (dict): Dictionary mapping each contig ID to its corresponding MAG name.
    """
    contig_to_mag = {contig: mag_name for mag_name, mag_obj in mags.items() for contig in mag_obj.contigs.values()}
    logger.info(f"{len(contig_to_mag)} contigs mapped to {len(mags)} MAGs.")

    return contig_to_mag


def build_contig_graph(network_data: pd.DataFrame) -> nx.Graph:
    """
        Builds an undirected interaction graph from a DataFrame containing contig pairs and interaction signals.

        Each row of the DataFrame should contain two contig IDs and a signal value representing the interaction strength.
        The function constructs a NetworkX graph where each contig is a node, and each edge is weighted by the signal.

        Args:
        - network_data (pd.DataFrame): DataFrame with columns ['contig1', 'contig2', 'signal'].

        Returns:
        - G (nx.Graph): An undirected graph where edges are weighted by the interaction signal.
    """
    logger.info("Creating contigs contacts graph...")

    G = nx.Graph()
    edges = [[int(u), int(v), w] for u, v, w in zip(network_data["contig1"], network_data["contig2"], network_data["signal"])]
    G.add_weighted_edges_from(edges)
    logger.info(f"Contigs graph has {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")
    return G


def estimate_noise(mags, network_data) -> pd.DataFrame:
    """
        Evaluates inter- and intra-MAG interaction signals from a contig interaction graph and
        MAG quality data, in order to detect potential noise or assembly artifacts.

        Parameters
        ----------
        mags : dict
            Dictionary where keys are MAG names and values are objects representing MAGs. Each object must have
            a `contigs` attribute that returns a dictionary of contigs.

        network_data : pd.DataFrame
            DataFrame containing interaction data with columns ['contig1', 'contig2', 'signal'].

        Returns
        -------
        pd.DataFrame
            DataFrame containing the average intra- and inter-MAG interaction signals for each MAG.

        Notes
        -----
        The function:
        - Builds a graph of contig interactions.
        - Maps each contig to its corresponding MAG.
        - Calculates average intra- and inter-MAG interaction signals.
        - Groups MAG pairs by quality (e.g., HQ-HQ, HQ-MQ, HQ-contaminated).
        - Applies log transformation and normalization to the interaction signals for downstream analysis.
    """
    logger.info("Computing the background noise of the experience.")

    # Get contigs <-> MAG mapping
    contig_to_mag = map_contigs_to_mags(mags)
    contig_to_mag_df = pd.DataFrame({"contig": list(contig_to_mag.keys()), "MAG": list(contig_to_mag.values())})

    # Merge network data with contig to MAG mapping
    mags_df = (
        network_data.merge(contig_to_mag_df, left_on="contig1", right_on="contig", how="left")
        .rename(columns={"MAG": "MAG1"})[["contig1", "contig2", "signal", "MAG1"]]
        .merge(contig_to_mag_df, left_on="contig2", right_on="contig", how="left")
        .rename(columns={"MAG": "MAG2"})[["MAG1", "MAG2", "contig1", "contig2", "signal"]]
    )

    # Define intra/inter contacts
    mags_df["type"] = np.where(mags_df["MAG1"] == mags_df["MAG2"], "intra", "inter")

    # Remove rows with NaN values
    mags_df = mags_df.dropna(subset=["MAG1", "MAG2"])

    ###################### INTRA-MAG INTERACTIONS COMPUTATION #################
    for _, mag in mags.items():
        intra_signals = mags_df.query("(MAG1 == @mag.name or MAG2 == @mag.name) and type == 'intra'")["signal"]
        mag.intra_signal = np.mean(intra_signals)
        if np.isnan(mag.intra_signal):
            mag.intra_signal = 0.0

    ########### INTER-MAG INTERACTIONS COMPUTATION #################
    def get_other_mag(row, target_mag):
        if row["MAG1"] == target_mag:
            return row["MAG2"]
        else:
            return row["MAG1"]

    for _, mag in mags.items():
        subset = ((mags_df["MAG1"] == mag.name) | (mags_df["MAG2"] == mag.name)) & (mags_df["type"] == "inter")
        if subset.sum() == 0:
            mag.inter_signals = 0
            continue
        else:
            df_inters = mags_df.copy()[subset]
            df_inters.loc[:, "mag0"] = mag.name
            df_inters.loc[:, "mag"] = df_inters.apply(get_other_mag, axis=1, target_mag=mag.name)
            scores = list(
                (
                    df_inters.groupby(["mag0", "mag"])["signal"]
                    .mean()
                    .reset_index()
                    .rename(columns={"mag0": "MAG1", "mag": "MAG2"})
                    .sort_values("MAG2")
                )["signal"]
            )
            while len(scores) < (len(mags) - 1):
                scores.append(0)
            mag.inter_signals = scores

    return mags


def plot_mag_noise(mags, image_file="background_noise.png") -> None:
    """
        Plots the background noise of MAG interactions.

        Args:
            mags (dict): Dictionary of MAG objects.
            image_file (str): Path to save the plot image.

        Returns:
            None
    """
    logger.info("Plotting background noise...")

    intra_signals = pd.DataFrame(
        {
            "MAG": [mag.name for mag in mags.values()],
            "type": ["intra"] * len(mags),
            "quality": [mag.quality for mag in mags.values()],
            "signal": [mag.intra_signal for mag in mags.values()],
        }
    )
    # Ajoute ce typage ici :
    intra_signals["signal"] = intra_signals["signal"].astype(float)
    intra_signals["log_signal"] = np.log10(intra_signals["signal"] + 1e-10)
    min_val = intra_signals["log_signal"].min()
    max_val = intra_signals["log_signal"].max()
    intra_signals["norm_log"] = (intra_signals["log_signal"] - min_val) / (max_val - min_val)

    inter_signals = pd.DataFrame(
        {
            "MAG": [mag.name for mag in mags.values()],
            "type": ["inter"] * len(mags),
            "quality": [mag.quality for mag in mags.values()],
            "signal": [mag.inter_signals for mag in mags.values()],
        }
    ).explode("signal")
    inter_signals["signal"] = inter_signals["signal"].astype(float)
    inter_signals["log_signal"] = np.log10(inter_signals["signal"] + 1e-10)
    min_val = inter_signals["log_signal"].min()
    max_val = inter_signals["log_signal"].max()
    inter_signals["norm_log"] = (inter_signals["log_signal"] - min_val) / (max_val - min_val)

    # Create a DataFrame for plotting
    data = pd.concat([intra_signals, inter_signals], ignore_index=True)
   
    
    # Plotting a split violin and boxplots plot; each violin show distribution of intra (left) and inter (right) signals, per quality
    #plt.figure(figsize=(12, 6))
    
    fig, axes = plt.subplots(2, 1, figsize=(16, 16))
    #plt.title("Distribution of Intra- and Inter-MAG Interaction Signals")
    # Boxplot
    sns.boxplot(data=data, x="quality", y="norm_log", hue="type", ax=axes[0])
    axes[0].set_title("Grouped Boxplots -  Normalised Log Signal")
    axes[0].legend(title="Interaction Type")

    sns.violinplot(data=data, x="quality", y="norm_log", hue="type", inner="quartile", ax=axes[1])
    axes[1].set_title("Violin Plot - Normalised Log Signal")
    #plt.legend(title="Interaction Type")
    axes[1].legend(title="Interaction Type")

    for ax in axes:
        ax.set_xlabel("MAGs Quality")
        ax.set_ylabel("Normalised interaction Signal (log scale)")
    
    
    plt.tight_layout()
    plt.savefig(image_file)
    plt.close(fig)

def classify_mags(mags, bin_summary) -> dict:
    """
        Classifies MAGs based on their completeness and contamination values.

        Args:
            mags (dict): Dictionary of MAG objects.
            bin_summary (pd.DataFrame): DataFrame containing MAG quality metrics.

        Returns:
            dict: Updated dictionary of MAG objects with quality classifications.
    """
    logger.info("Classifying MAGs based on their completeness and contamination values...")

    def classify_mag( mag,) -> Union[str, None]:
        
        c, r = mag.completeness, mag.contamination
        if pd.isna(c) or pd.isna(r):
            return "Unknown"
        elif r > 1.1:
            return "Contaminated"
        elif c > 0.9 and r < 1.05:
            return "Complete"
        elif c > 0.9 and 1.05 <= r < 1.1:
            return "HQ"
        elif c > 0.7 and r <= 1.1:
            return "MQ"
        elif c > 0.5 and r <= 1.1:
            return "LQ"
        elif c <= 0.5 and r <= 1.1:
            return "PQ"

    for mag in mags.values():
        if mag.name in bin_summary["MAG"].values:
            mag_data = bin_summary[bin_summary["MAG"] == mag.name].iloc[0]
            mag.info = dict(mag_data)
            mag.completeness = mag_data["Completeness"]
            mag.contamination = mag_data["Redundancy"]
            mag.quality = classify_mag(mag)
        else:
            logger.warning(f"Warning: {mag.name} not found in bin summary.")

    quals = [x.quality for x in mags.values()]
    cnts = Counter(quals)
    for value, count in cnts.items():
        print(f"{value} MAGs: {count}")

    return mags


def compute_mge_mag_interactions(  
    network_data,
    mags,
    mge_mags,
    interaction_threshold=10.0,
    min_interacting_contigs=5,
    output_file="mge_mag_interactions.tsv",
    image_file="mge_mag_histogram.png",
):
    """
        Computes interactions between contigs from mgeMAGs and those from MAGs,
        and filters them based on a percentage threshold.
        Generates two histograms: one showing the total interaction signal (log scale)
        and another showing the number of MAG contigs interacting with an mgeMAG.

        Args:
        - network_data (pd.DataFrame): DataFrame with columns ['contig1', 'contig2', 'signal'].
        - mge_mags (dict): Dictionary {mgeMAG_name: MAG_object} representing the mgeMAGs.
        - mags (dict): Dictionary {MAG_name: MAG_object} representing the MAGs.
        - output_file (str): Name of the output text file.
        - image_file (str): Name of the output image file for the histograms.
        - interaction_threshold (float): Threshold (%) to validate an interaction.
        - min_interacting_contigs (int): Minimum number of MAG contigs interacting with the mgeMAG.

        Returns:
        - image_file (str): Name of the saved image file.
    """

    G = build_contig_graph(network_data)

    contig_to_mag = map_contigs_to_mags(mags)

    interaction_results = []

    for mgemag in mge_mags.values():
        interaction_signals = {}
        interacting_contigs_by_mag = {}
        for contig1 in mgemag.contigs.values():
            if contig1 in G:
                for contig2 in G.neighbors(contig1):
                    if contig2 in contig_to_mag:
                        mag_name = contig_to_mag[contig2]
                        signal = G[contig1][contig2]["weight"]
                        interaction_signals[mag_name] = interaction_signals.get(mag_name, 0) + signal
                        interacting_contigs_by_mag.setdefault(mag_name, set()).add(contig2)

        total_signal_mgemag = sum(interaction_signals.values())

        for mag_name, signal_sum in interaction_signals.items():
            percent_signal = (signal_sum / total_signal_mgemag) * 100 if total_signal_mgemag > 0 else 0
            if percent_signal >= interaction_threshold:
                interacting_contig_count = len(interacting_contigs_by_mag[mag_name])
                if interacting_contig_count >= min_interacting_contigs:
                    interaction_results.append((mgemag.name, mag_name, signal_sum, percent_signal, interacting_contig_count))

    with open(output_file, "w") as f:
        f.write("MGE\tMAG\tTotal Signal\tSignal rate\tInteracting Contigs\n")
        for res in interaction_results:
            mge_name, mag_name, signal_sum, percent_signal, contig_count = res
            f.write(f"{mge_name:<24}\t{mag_name:<24}\t{signal_sum:.6f}\t{percent_signal:.2f}\t{contig_count}\n")

    filtered_signals = [x[2] for x in interaction_results if x[2] > 0]
    percentage_counts = [x[3] for x in interaction_results]

    if filtered_signals and percentage_counts:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Histogramm of total interaction signal  (log scale)
        axes[0].hist(
            filtered_signals,
            bins=np.logspace(np.log10(min(filtered_signals)), np.log10(max(filtered_signals)), 20),
            edgecolor="black",
        )
        axes[0].set_xscale("log")
        axes[0].set_xlabel("MGE-MAG total interaction Signal  (log scale)")
        axes[0].set_ylabel("Frequences")
        axes[0].set_title(f"Signal Distribution with threshold = (\u2265 {interaction_threshold}%)")

        # Histogramm of interaction rate
        bins = range(1, math.ceil(max(percentage_counts)) + 1 + 1)  # +2 pour être équivalent à +2 dans ta ligne
        axes[1].hist(percentage_counts, bins=bins, edgecolor="black", align="left")
        axes[1].set_xlabel("Association Rate")
        axes[1].set_ylabel("Frequence")
        axes[1].set_title("Distribution of association rates of MgeMags with their MAGs.")

        plt.tight_layout()
        fig.savefig(image_file)
        plt.close(fig)

    for res in interaction_results:
        mge_name, mag_name, _, _, _ = res
        mge_mags[mge_name].add_host(mag_name=mag_name, mag=mags[mag_name])
        mags[mag_name].add_mge(mge_name=mge_name, mge_mag=mge_mags[mge_name])

    print(f"MAGs with: 0  MGEs: {len([x for x in mags.values() if len(x.mges) == 0])}")
    print(f"           1  MGEs: {len([x for x in mags.values() if len(x.mges) == 1])}")
    print(f"           2  MGEs: {len([x for x in mags.values() if len(x.mges) == 2])}")
    print(f"           3+ MGEs: {len([x for x in mags.values() if len(x.mges) >= 3])}")
    print("---")
    print(f"MGEs with: 0  hosts: {len([x for x in mge_mags.values() if len(x.hosts) == 0])}")
    print(f"           1  hosts: {len([x for x in mge_mags.values() if len(x.hosts) == 1])}")
    print(f"           2  hosts: {len([x for x in mge_mags.values() if len(x.hosts) == 2])}")
    print(f"           3  hosts: {len([x for x in mge_mags.values() if len(x.hosts) == 3])}")
    print(f"           4  hosts: {len([x for x in mge_mags.values() if len(x.hosts) == 4])}")
    print(f"           5+ hosts: {len([x for x in mge_mags.values() if len(x.hosts) >= 5])}")
    
    return mags, mge_mags


def annotate_hosts(                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
    contig_data: pd.DataFrame,
    network_data: pd.DataFrame,
    bin_summary: pd.DataFrame,
    mges_bin_summary: pd.DataFrame,
    interaction_threshold: int,
    min_interacting_contigs: int
) -> None:

    
    # Instantiating objects.
    logger.info("Instantiating objects...")
    mags, mge_mags = create_bins(contig_data, mges_bin_summary)
    mags = classify_mags(mags, bin_summary)

    # Evaluate background noise
    logger.info("Evaluating the background noise of the experience...")
    mags = estimate_noise(mags, network_data)

    # Plot background noise
    #logger.info("Plotting background noise...")
    plot_mag_noise(mags, image_file="background_noise.png")

    # Associate each MGE to their most likely MAG(s)
    logger.info("Computing MGE and MAG interactions and associate each MGE to potential host(s)...")
    mags, mge_mags = compute_mge_mag_interactions(
        network_data,
        mags,
        mge_mags,
        interaction_threshold=interaction_threshold,
        min_interacting_contigs=min_interacting_contigs,
    )
    logger.info("Host association done !!!")

    return mags, mge_mags
