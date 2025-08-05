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
import os
from scipy.sparse import coo_matrix
from collections import Counter


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
        Represent a MAG (Metagenome-Assembled Genome). inherit from Bin.
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
        Represent a MGE-MAG (mobile genetic element). inherited from Bin.
    """

    def __init__(self, name):
        super().__init__(name)
        self.hosts: dict = {}

    def add_host(self, mag_name, mag) -> None:
        self.hosts[mag_name] = mag

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{self.name}', {len(self.contigs)} contigs, {len(self.hosts)} hosts)"


def update_final_bin_with_mge_legacy(contig_data: pd.DataFrame, mges_bin_summary: pd.DataFrame) -> pd.DataFrame :
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

def update_final_bin_with_mge_lagacy2(
    contig_data: pd.DataFrame,
    mges_bin_summary: pd.DataFrame
) -> tuple[pd.DataFrame, list[str]]:
    """
    Updates the 'Final_bin' column of contig_data1 by replacing 'ND' values with
    the corresponding MGEs (from mges_bin_summary), if available.
    Also returns the list of contigs already assigned to a MAG and found in a MGE.

    Args:
        contig_data1 (pd.DataFrame): DataFrame containing at least ['Name', 'Final_bin']
        mges_bin_summary (pd.DataFrame): DataFrame with ['Contigs', 'BinName']

    Returns:
        tuple:
            - pd.DataFrame: Updated contig_data
            - list[str]: List of contigs assigned both to a MAG and present in MGE binning
    """

    contig_data = contig_data.copy()
    mges_bin = mges_bin_summary.copy()

    # Construire dictionnaire contig → MGE
    mges_bin = mges_bin.assign(Contigs=mges_bin["Contigs"].str.split(",")).explode("Contigs")
    mges_bin["Contigs"] = mges_bin["Contigs"].str.strip()
    contig_to_mge = dict(zip(mges_bin["Contigs"], mges_bin["BinName"]))

    # Liste pour suivre les contigs ambigus
    ambiguous_contigs = []

    def update(row):
        name = row.get("Name")
        final_bin = row.get("Final_bin")
        mge = contig_to_mge.get(name)

        if pd.isna(final_bin) or final_bin == "ND":
            return mge if mge else final_bin
        elif isinstance(final_bin, str) and final_bin.startswith("metator_") and mge:
            ambiguous_contigs.append(name)
            return final_bin  # on ne modifie pas, mais on signale
        else:
            return final_bin

    contig_data["Final_bin"] = contig_data.apply(update, axis=1)

   
    return contig_data


def update_final_bin_with_mge(contig_data: pd.DataFrame, mges_bin_summary: pd.DataFrame) -> pd.DataFrame:
    """
    Updates the 'Final_bin' column of contig_data by replacing 'ND' values or MAGs containing 'metator_'
    with the corresponding MGEs (from mges_bin_summary), if available.

    Args:
        contig_data (pd.DataFrame): DataFrame containing at least the columns ['Name', 'Final_bin']
        mges_bin_summary (pd.DataFrame): DataFrame containing the columns ['Contigs', 'BinName']

    Returns:
        pd.DataFrame: An updated copy of contig_data with the 'Final_bin' column modified when applicable
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
            final_bin = "ND"  # pour traiter aussi les NaN

        if final_bin == "ND" and mge:
            return mge
        elif isinstance(final_bin, str) and "metator_" in final_bin and mge:
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

        #if bin_name.startswith("MetaTOR_MGE_"):  
        if isinstance(bin_name, str) and "MetaTOR_MGE_" in bin_name: # Identify MGE-MAGs
            if bin_name not in mge_mags:
                mge_mags[bin_name] = MgeMag(bin_name)  # create a new MGE-MAG
            mge_mags[bin_name].add_contig(contig_name, contig_id)

        elif isinstance(bin_name, str) and "metator_" in bin_name:
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


def estimate_noise_legacy(mags, network_data) -> pd.DataFrame:
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

def estimate_noise(mags, network_data):
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
        ------
            Tuple (mags, background_noise) where:
            - mags (dict): MAGs enriched with intra- and inter-signal information
            - mean_intra (float): Mean of intra-MAG interaction signals
            - mean_inter (float): Mean of inter-MAG interaction signals
            - sparse_matrix (coo_matrix): Sparse full (non-symmetric, non-triangular) matrix of interactions
            - mag_names (list): List of MAG names in order of matrix indices
    """

    logger.info("Computing the background noise of the experience.")

    contig_to_mag = map_contigs_to_mags(mags)
    contig_to_mag_df = pd.DataFrame({"contig": list(contig_to_mag.keys()), "MAG": list(contig_to_mag.values())})

    mags_df = (
        network_data.merge(contig_to_mag_df, left_on="contig1", right_on="contig", how="left")
        .rename(columns={"MAG": "MAG1"})[["contig1", "contig2", "signal", "MAG1"]]
        .merge(contig_to_mag_df, left_on="contig2", right_on="contig", how="left")
        .rename(columns={"MAG": "MAG2"})[["MAG1", "MAG2", "contig1", "contig2", "signal"]]
    )

    mags_df["type"] = np.where(mags_df["MAG1"] == mags_df["MAG2"], "intra", "inter")
    mags_df = mags_df.dropna(subset=["MAG1", "MAG2"])

    for _, mag in mags.items():
        intra_signals = mags_df.query("(MAG1 == @mag.name and MAG2 == @mag.name) and type == 'intra'")["signal"]
        mag.intra_signal = round(np.sum(intra_signals), 8)
        if np.isnan(mag.intra_signal):
            mag.intra_signal = 0.0

    def get_other_mag(row, target_mag):
        return row["MAG2"] if row["MAG1"] == target_mag else row["MAG1"]

    for _, mag in mags.items():
        subset = ((mags_df["MAG1"] == mag.name) | (mags_df["MAG2"] == mag.name)) & (mags_df["type"] == "inter")
        if subset.sum() == 0:
            mag.inter_signals = [0.0] * (len(mags) - 1)
            mag.inter_signal_sum = 0.0
            continue
        df_inters = mags_df.loc[subset].copy()
        df_inters["mag0"] = mag.name
        df_inters["mag"] = df_inters.apply(get_other_mag, axis=1, target_mag=mag.name)
        scores = list(
            df_inters.groupby(["mag0", "mag"])["signal"]
            .sum()
            .reset_index()
            .rename(columns={"mag0": "MAG1", "mag": "MAG2"})
            .sort_values("MAG2")["signal"]
        )
        while len(scores) < (len(mags) - 1):
            scores.append(0)
        mag.inter_signals = scores

    # --- Build full sparse matrix (non-symmetric, non-triangular) ---
    mag_names = sorted(map(str, mags.keys()))

    mag_to_idx = {name: i for i, name in enumerate(mag_names)}

    data, row, col = [], [], []

    # Diagonal (intra)
    for mag_name in mag_names:
        idx = mag_to_idx[mag_name]
        signal = getattr(mags[mag_name], "intra_signal", 0.0)
        data.append(signal)
        row.append(idx)
        col.append(idx)

    # All inter-signals as-is (non-symmetric)
    for i, mag_name in enumerate(mag_names):
        if not hasattr(mags[mag_name], "inter_signals"):
            continue
        other_names = [name for name in mag_names if name != mag_name]
        for j, other_name in enumerate(other_names):
            source_idx = mag_to_idx[mag_name]
            target_idx = mag_to_idx[other_name]
            signal = mags[mag_name].inter_signals[j]
            data.append(signal)
            row.append(source_idx)
            col.append(target_idx)

    sparse_matrix = coo_matrix((data, (row, col)), shape=(len(mag_names), len(mag_names)))

    diag_mask = np.array(row) == np.array(col)
    intra_values = np.array(data)[diag_mask]
    mean_intra = np.mean(intra_values) if intra_values.size > 0 else 0.0

    off_diag_mask = np.array(row) != np.array(col)
    inter_values = np.array(data)[off_diag_mask]
    mean_inter = np.mean(inter_values) if inter_values.size > 0 else 0.0

    logger.info(f"Moyenne intra-MAG : {mean_intra:.8f}")
    logger.info(f"Moyenne inter-MAG : {mean_inter:.8f}")

    return mags, mean_inter, sparse_matrix, mag_names

def plot_mag_noise_legacy(mags, image_file: str = "background_noise.png") -> None:

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

def plot_mag_noise(sparse_matrix, mag_names, prefix="mag_interactions", matrix_file="inter_mag_matrix.tsv"):
    """
    Visualizes the MAG interaction sparse matrix:
    - Boxplot of intra- vs inter-MAG signals (log10 transformed)
    - Heatmap of log-transformed matrix (dense version)
    - Save matrix as TSVmean_inter

    Parameters
    ----------
    sparse_matrix : scipy.sparse.coo_matrix
        Sparse interaction matrix (non-symmetric, non-triangular).
    mag_names : list of str
        List of MAG names corresponding to matrix rows/columns.
    prefix : str
        Prefix for output image filenames.
    matrix_file : str
        Path to save dense matrix as TSV.
    """
    logger.info("Plotting background noise from sparse matrix...")

    epsilon = 1e-6
    data = sparse_matrix.data
    row = sparse_matrix.row
    col = sparse_matrix.col

    # Identify intra/inter based on diagonal
    intra_idx = row == col
    inter_idx = row != col

    intra_values = data[intra_idx]
    inter_values = data[inter_idx]
    mean_intra = np.mean(intra_values) if intra_values.size > 0 else 0.0
    mean_inter = np.mean(inter_values) if inter_values.size > 0 else 0.0
    # Log10 transform
    intra_log = np.log10(intra_values + epsilon)
    inter_log = np.log10(inter_values + epsilon)

    # --- BOXPLOT ---
    df_box = pd.DataFrame({
        "Type": ["Intra-MAG"] * len(intra_log) + ["Inter-MAG"] * len(inter_log),
        "log10(Signal)": np.concatenate([intra_log, inter_log])
    })
    

    plt.figure(figsize=(12, 9))
    sns.boxplot(data=df_box, x="Type", y="log10(Signal)", hue="Type", palette="Set2", legend=False)
    plt.title("Interaction signal Distribution", fontsize=15)
    plt.text(0, 1.02, f"Intra-mag mean signal = {mean_intra:.8f}", ha='left', va='bottom', transform=plt.gca().transAxes, fontsize=12)
    plt.text(1, 1.02, f"Inter-mag mean signal = {mean_inter:.8f}", ha='right', va='bottom', transform=plt.gca().transAxes, fontsize=12)

    plt.xlabel("Interaction Type")
    plt.ylabel("log10(Signal)")
    plt.tight_layout()
    os.makedirs(os.path.dirname(prefix), exist_ok=True)
    plt.savefig(f"{prefix}_boxplot.png")
    plt.close()
    logger.info(f"Boxplot saved in : {prefix}_boxplot.png")

    # # --- HEATMAP ---
    # dense_matrix = sparse_matrix.toarray()
    # matrix_log = np.log10(dense_matrix + epsilon) # matrix_log = np.log10(dense_matrix + epsilon)
    # df_matrix = pd.DataFrame(matrix_log, index=mag_names, columns=mag_names)
    # plt.figure(figsize=(30, 30))
    # sns.heatmap(df_matrix, cmap="viridis", square=True, cbar_kws={"label": "log10(Signal)"})
    # plt.title("MAG-MAG log10 interactions Heatmap", fontsize=15)
    # plt.text(0, 1.02, f"Intra-mag mean signal = {mean_intra:.8f}", ha='left', va='bottom', transform=plt.gca().transAxes, fontsize=13)
    # plt.text(1, 1.02, f"Inter-mag mean signal = {mean_inter:.8f}", ha='right', va='bottom', transform=plt.gca().transAxes, fontsize=13)
    # plt.xlabel("MAG")
    # plt.ylabel("MAG")
    # plt.tight_layout()
    # plt.savefig(f"{prefix}_heatmap.png")
    # plt.close()
    # logger.info(f"Heatmap saved in : {prefix}_heatmap.png")

    
    # df_export = pd.DataFrame(dense_matrix, index=mag_names, columns=mag_names)
    # with open(matrix_file, "w") as f:
    #     f.write(f"# intra-MAG mean: {mean_intra:.4f}\tinter-MAG mean: {mean_inter:.4f}\n")
    #     df_export.to_csv(f, sep="\t")
    # logger.info(f"Interaction matrix saved in: {matrix_file}")

def classify_mags_legacy(mags, bin_summary) -> dict:
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
        if pd.isna(c) or pd.isna(r) or c < 0.1:
            return "Unknown"
        if r > 1.1:
            return "Contaminated"
        elif c > 0.9 and r < 1.05:
            return "Complete"
        elif c > 0.9 and 1.05 <= r < 1.10:
            return "HQ"
        elif 0.7 < c <= 0.9 and r <= 1.10:
            return "MQ"
        if 0.5 < c <= 0.7 and r <= 1.10:
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

def classify_mags(mags, bin_summary, plot_prefix: str = None) -> dict:
    """
    Classifies MAGs based on their completeness and contamination values.

    Args:
        mags (dict): Dictionary of MAG objects.
        bin_summary (pd.DataFrame): DataFrame containing MAG quality metrics.
        plot_prefix (str, optional): Path prefix for saving output plots.

    Returns:
        dict: Updated dictionary of MAG objects with quality classifications.
    """
    logger.info("Classifying MAGs based on their completeness and contamination values...")

    def classify_mag(mag) -> Union[str, None]:
        c, r = mag.completeness, mag.contamination
        if pd.isna(c) or pd.isna(r) or c < 0.1:
            return "Unknown"
        if r > 1.1:
            return "Contaminated"
        elif c > 0.9 and r < 1.05:
            return "Complete"
        elif c >= 0.9 and 1.05 <= r < 1.10:
            return "HQ"
        elif 0.7 < c <= 0.9 and r <= 1.10:
            return "MQ"
        if 0.5 < c <= 0.7 and r <= 1.10:
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

    # Statistiques
    quals = [x.quality for x in mags.values()]
    cnts = Counter(quals)
    for value, count in cnts.items():
        print(f"{value} MAGs: {count}")

    # Plot
    if plot_prefix:
        df_qual = pd.DataFrame(cnts.items(), columns=["Quality", "Count"]).sort_values("Count", ascending=False)

        plt.figure(figsize=(10, 6))
        ax = sns.barplot(x="Quality", y="Count", data=df_qual)

        # Ajouter les valeurs exactes au-dessus des barres
        for p in ax.patches:
            height = p.get_height()
            ax.text(
                p.get_x() + p.get_width() / 2,
                height,
                f"{int(height)}",
                ha="center",
                va="bottom",
                fontsize=11
            )

        plt.title("Distribution of MAG quality categories")
        plt.xlabel("MAG quality")
        plt.ylabel("Number of MAGs")
        plt.tight_layout()

        # Créer le dossier si nécessaire
        os.makedirs(os.path.dirname(plot_prefix), exist_ok=True)
        plt.savefig(f"{plot_prefix}_mag_quality_distribution.png")
        plt.close()
    # Nouveau plot : MAGs qui passent les filtres (c > 0.5 et r < 1.1)
    if plot_prefix:
        passes = 0
        fails = 0
        for mag in mags.values():
            if hasattr(mag, "completeness") and hasattr(mag, "contamination"):
                if mag.completeness > 0.5 and mag.contamination < 1.1:
                    passes += 1
                else:
                    fails += 1

        df_filter = pd.DataFrame({
            "Filter": ["Valids", "Others"],
            "Count": [passes, fails]
        })

        plt.figure(figsize=(10, 7))
        ax = sns.barplot(
                    x="Filter",
                    y="Count",
                    hue="Filter",
                    data=df_filter,
                    palette={"Valids": "green", "Others": "red"},
                    legend=False
                )


        # Ajouter le nombre sur chaque barre
        for p in ax.patches:
            height = p.get_height()
            ax.text(
                p.get_x() + p.get_width() / 2,
                height,
                f"{int(height)}",
                ha="center",
                va="bottom",
                fontsize=11
            )

        plt.title("MAGs with completeness > 50% and contamination < 1.1 Vs Others")
        plt.xlabel("Status")
        plt.ylabel("Number of MAGs")
        plt.tight_layout()
        plt.savefig(f"{plot_prefix}_HQ_mag_vs_others.png")
        plt.close()
    return mags

def compute_mge_mag_interactions(  
    mean_inter,
    network_data,
    mags,
    mge_mags,
    interaction_threshold=10.0,
    min_interacting_contigs=5,
    output_file: str = "mge_mag_interactions.tsv",
    plot_prefix: str ="mge_mag_interactions",
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
    all_signals = []  
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
            all_signals.append(signal_sum)
            interacting_contig_count = len(interacting_contigs_by_mag[mag_name])
            if ((percent_signal >= interaction_threshold) and (interacting_contig_count >= min_interacting_contigs) and (signal_sum > mean_inter)):
                interaction_results.append((mgemag.name, mag_name, signal_sum, percent_signal, interacting_contig_count))

    with open(output_file, "w") as f:
        f.write("MGE\tMAG\tTotal Signal\tSignal rate\tInteracting Contigs\n")
        for res in interaction_results:
            mge_name, mag_name, signal_sum, percent_signal, contig_count = res
            f.write(f"{mge_name:<24}\t{mag_name:<24}\t{signal_sum:.6f}\t{percent_signal:.2f}\t{contig_count}\n")


    df = pd.DataFrame(interaction_results, columns=["MGE", "MAG", "Signal", "Percentage", "Contigs"])
    df["MAG"]=df["MAG"].astype(str)
    df["MGE"]=df["MGE"].astype(str)
    if not df.empty:
        # Boxplot: Raw vs Filtered Signals
        df_all = pd.DataFrame({
            "Signal": all_signals + list(df["Signal"]),
            "Type": ["All"] * len(all_signals) + ["Filtered"] * len(df["Signal"])
        })
        plt.figure(figsize=(10, 7))
        sns.boxplot(data=df_all, x="Type", y="Signal", showfliers=False)
        plt.yscale("log")
        plt.title("Log10-scaled distribution of MGE-MAG interaction signals")
        plt.ylabel("Signal")
        plt.tight_layout()
        plt.savefig(f"{plot_prefix}_boxplot_signal_filtered.png")
        plt.close()

        # Combined Histogram: MAGs per MGE & MGEs per MAG
        
        mag_count = df["MGE"].value_counts()
        mge_count = df["MAG"].value_counts()

        fig, axs = plt.subplots(1, 2, figsize=(16, 8))

        sns.countplot(x=mag_count.values, ax=axs[0])
        axs[0].set_title("Number of MAGs linked to each MGE")
        axs[0].set_xlabel("Number of MAGs")
        axs[0].set_ylabel("Number of MGEs")

        sns.countplot(x=mge_count.values, ax=axs[1])
        axs[1].set_title("Number of MGEs linked to each MAG")
        axs[1].set_xlabel("Number of MGEs")
        axs[1].set_ylabel("Number of MAGs")

        plt.tight_layout()
        plt.savefig(f"{plot_prefix}_counts_summary.png")
        plt.close()

        # Optional Heatmap if matrix is not too sparse or large
        if len(df["MGE"].unique()) <= 100 and len(df["MAG"].unique()) <= 100:
            heatmap_matrix = df.pivot(index="MGE", columns="MAG", values="Signal").fillna(0)
            plt.figure(figsize=(max(16, 1 * heatmap_matrix.shape[1]), max(16, 1 * heatmap_matrix.shape[0])))
            sns.heatmap(np.log10(heatmap_matrix + 1e-6), cmap="viridis", cbar_kws={"label": "log10(Signal)"})
            plt.title("Heatmap of MGE-MAG association signals (log10 scale)")
            plt.xlabel("MAG")
            plt.ylabel("MGE")
            plt.tight_layout()
            plt.savefig(f"{plot_prefix}_heatmap.png")
            plt.close()

            # --- HEATMAP automatique si taille raisonnable ---
        if df.shape[0] > 0 and df["MGE"].nunique() <= 100 and df["MAG"].nunique() <= 100:
            pivot = df.pivot(index="MGE", columns="MAG", values="Signal").fillna(0)
            
            # Log10 transformation
            pivot_log = np.log10(pivot + 1e-6)

            # Tri par somme décroissante
            mge_order = pivot_log.sum(axis=1).sort_values(ascending=False).index
            mag_order = pivot_log.sum(axis=0).sort_values(ascending=False).index
            pivot_log_sorted = pivot_log.loc[mge_order, mag_order]

            plt.figure(figsize=(min(1 * len(mag_order), 18), min(1 * len(mge_order), 18)))
            sns.heatmap(pivot_log_sorted, cmap="viridis", cbar_kws={"label": "log10(Signal)"})
            plt.title("Heatmap log10 des signaux MGE–MAG (triés)")
            plt.xlabel("MAG")
            plt.ylabel("MGE")
            plt.tight_layout()
            plt.savefig(f"{plot_prefix}_heatmap.png")
            plt.close()

            # --- Heatmap avec clustering (facultatif mais informatif) ---
            sns.clustermap(
                pivot_log,
                cmap="viridis",
                figsize=(min(1 * pivot_log.shape[1], 25), min(1 * pivot_log.shape[0], 25)),
                cbar_kws={"label": "log10(Signal)"},
                metric="euclidean",
                method="average"
            )
            plt.savefig(f"{plot_prefix}_heatmap_clustermap.png")
            plt.close()

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
    min_interacting_contigs: int,
    outdir
) -> None:

    
    # Instantiating objects.
    logger.info("Instantiating objects...")
    mags, mge_mags = create_bins(contig_data, mges_bin_summary)
    mags = classify_mags(mags, bin_summary, plot_prefix =os.path.join(outdir, "plot"))

    # Evaluate background noise
    logger.info("Evaluating the background noise of the experience...")
    mags, mean_inter, interactions_matrix, mag_names = estimate_noise(mags, network_data)

    # Plot background noise
    #logger.info("Plotting background noise...")
    plot_mag_noise(interactions_matrix, mag_names, prefix=os.path.join(outdir, "mag_interactions"), matrix_file=os.path.join(outdir, "inter_mag_matrix.tsv") )

    # Associate each MGE to their most likely MAG(s)
    logger.info("Computing MGE and MAG interactions and associate each MGE to potential host(s)...")
    mags, mge_mags = compute_mge_mag_interactions(
        mean_inter,
        network_data,
        mags,
        mge_mags,
        interaction_threshold=interaction_threshold,
        min_interacting_contigs=min_interacting_contigs,
        output_file=os.path.join(outdir, "mge_mag_interactions.tsv"),
        plot_prefix=os.path.join(outdir, "mge_mag_interactions")
        
    )
    logger.info("Host association done !!!")

    return mags, mge_mags