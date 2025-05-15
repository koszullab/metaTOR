#!/usr/bin/env python3
# coding: utf-8
"""
This module is used to detect and associate an MGE with its bacterial MAG host.

It defines two main classes, `MAG` and `MGEMAG`, along with their methods.
Starting from a `MAG` or `MGEMAG` object, the module retrieves all contigs belonging to a given MAG or MGE-MAG based on binning results.

The module allows evaluation of background noise by calculating intra- and inter-MAG interaction signals, as well as interactions between MGE-MAGs and MAGs. This helps identify the host MAG for a given MGE-MAG.

Main functions and methods include:
    - Class methods:
        - `__init__()`: Initializes a MAG or MGE-MAG object.
        - `add_contig()`: Adds a contig to a MAG or MGE-MAG.
        - `get_contigs()`: Retrieves contigs from a MAG or MGE-MAG.
        - `add_intra_signal()`: Adds an intra-MAG interaction signal.
        - `add_inter_signal()`: Adds an inter-MAG interaction signal.
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
import networkx as nx  # type: ignore
import matplotlib.pyplot as plt  # type: ignore

# import pypairix
# from metator.log import logger
from os.path import join
import pandas as pd  # type: ignore
import numpy as np  # type: ignore
import seaborn as sns  # type: ignore
from itertools import combinations
import math

import logging

# Configuration du logger
logging.basicConfig(
    level=logging.INFO,  # Ou DEBUG pour plus de détails
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.FileHandler("mags_mge_analysis.log"), logging.StreamHandler()],
)

logger = logging.getLogger(__name__)


#################### MAG CLASS AND FUNCTIONS ########################################################################
#####################################################################################################################
###################################################################################################################
class Bin:
    """
    Classe de base représentant un ensemble de contigs appartenant à un même groupe (ex: MAG ou MGE-MAG).
    """

    def __init__(self, name):
        self.name = name
        self.contigs = {}
        self.intra_signal = []
        self.inter_signal = {}

    def add_contig(self, contig_name, contig):
        self.contigs[contig_name] = contig

    def get_contigs(self):
        return self.contigs

    def add_intra_signal(self, signal):
        self.intra_signal.append(signal)

    def get_intra_signal(self):
        return self.intra_signal

    def add_inter_signal(self, partner, signal):
        self.inter_signal.setdefault(partner, []).append(signal)

    def get_inter_signal(self):
        return self.inter_signal

    def __repr__(self):
        return f"{self.__class__.__name__}('{self.name}', {len(self.contigs)} contigs)"


class MAG(Bin):
    """Représente un MAG (Metagenome-Assembled Genome). Hérite de Bin."""

    all_mags = {}

    def __init__(self, name):
        super().__init__(name)
        MAG.all_mags[name] = self

    @classmethod
    def list_all_mags(cls):
        return [(mag_name, mag_obj.contigs) for mag_name, mag_obj in cls.all_mags.items()]


class MGEMAG(Bin):
    """Représente un MGE-MAG (intégrant des éléments génétiques mobiles). Hérite de Bin."""

    all_mge_mags = {}

    def __init__(self, name):
        super().__init__(name)
        MGEMAG.all_mge_mags[name] = self

    @classmethod
    def list_all_mge_mags(cls):
        return [(mge_name, mge_obj.contigs) for mge_name, mge_obj in cls.all_mge_mags.items()]


def create_bin(contig_data):
    """
    Creates MGE-MAGs from contig data.

    Args:
        contig_data (pd.DataFrame): DataFrame containing information about contigs.

    Returns:
        dict: A dictionary mapping MGE-MAG names to their corresponding MGEMAG objects.
    """

    mags = {}
    mge_mags = {}
    logger.info("Intantiating bins (MAGs and MGE-MAGs)...")

    for _, row in contig_data.iterrows():
        contig_id, contig_name, bin_name = row["ID"], row["Name"], row["Final_bin"]

        if bin_name.startswith("MetaTOR_MGE_"):  # Identify MGE-MAGs
            if bin_name not in mge_mags:
                mge_mags[bin_name] = MGEMAG(bin_name)  # create a new MGE-MAG
            mge_mags[bin_name].add_contig(contig_name, contig_id)

        elif bin_name.startswith("metator_"):
            if bin_name not in mags:
                mags[bin_name] = MAG(bin_name)  # create a new MGE-MAG
            mags[bin_name].add_contig(contig_name, contig_id)
    logger.info(f"{len(mags)} MAGs and {len(mge_mags)} MGE-MAGs have been created.")

    return mags, mge_mags


###########################  OTHER FUNCTIONS ########################################################
#####################################################################################################
#####################################################################################################


def map_contigs_to_mags(mags: dict) -> dict:
    """
    Creates a mapping from each contig to its corresponding MAG name.

    Iterates through a dictionary of MAG objects and associates each contig with the name of the MAG it belongs to.

    Args:
    - mags (dict): Dictionary in the form {mag_name: mag_object}, where each mag_object has a get_contigs() method
      that returns a dictionary of contigs.

    Returns:
    - contig_to_mag (dict): Dictionary mapping each contig ID to its corresponding MAG name.
    """
    contig_to_mag = {contig: mag_name for mag_name, mag_obj in mags.items() for contig in mag_obj.get_contigs().values()}
    logger.info(f"{len(contig_to_mag)} contigs mapped to their MAGs.")

    return contig_to_mag


def build_interaction_graph(network_data: pd.DataFrame) -> nx.Graph:
    """
    Builds an undirected interaction graph from a DataFrame containing contig pairs and interaction signals.

    Each row of the DataFrame should contain two contigs and a signal value representing the interaction strength.
    The function constructs a NetworkX graph where each contig is a node, and each edge is weighted by the signal.

    Args:
    - network_data (pd.DataFrame): DataFrame with columns ['contig1', 'contig2', 'signal'].

    Returns:
    - G (nx.Graph): An undirected graph where edges are weighted by the interaction signal.
    """
    logger.info("Creating contigs contacts graph...")

    G = nx.Graph()
    for _, row in network_data.iterrows():
        contig1, contig2, signal = int(row["contig1"]), int(row["contig2"]), row["signal"]
        G.add_edge(contig1, contig2, weight=signal)
    logger.info(f"Graph created with {G.number_of_nodes()} nodes et {G.number_of_edges()} edges.")

    return G


def evaluate_experience_noise(
    network_data,
    mags,
    completeness_data,
    inter_output_file="inter_mag_interactions.txt",
    intra_output_file="intra_mag_interactions.txt",
    image_file="mags_interactions_by_quality.png",
):
    """
    Evaluates inter- and intra-MAG interaction signals from a contig interaction graph and
    MAG quality data, in order to detect potential noise or assembly artifacts.

    Parameters
    ----------
    network_data : pandas.DataFrame
        DataFrame containing the contig interaction graph, with columns "contig1", "contig2", and "signal".

    mags : dict
        Dictionary where keys are MAG names and values are objects representing MAGs. Each object must have
        a `get_contigs()` method that returns a dictionary of contigs.

    completeness_data : pandas.DataFrame
        DataFrame containing completeness and contamination data for each MAG,
        with at least the columns "Mag" (MAG name), "Completeness", and "Redundancy".

    inter_output_file : str, optional
        Name of the output text file for average inter-MAG interaction signals. Default is "inter_mag_interactions.txt".

    intra_output_file : str, optional
        Name of the output text file for average intra-MAG interaction signals. Default is "intra_mag_interactions.txt".

    image_file : str, optional
        Name of the image file for plotting MAG interaction signals by quality.
        This parameter is declared but not currently used in the function. Default is "mags_interactions_by_quality.png".

    Side Effects
    ------------
    - Writes inter-MAG average signals to a text file.
    - Writes intra-MAG average signals to a text file.
    - Adds inter-MAG interaction signal information to each MAG object using its `add_inter_signal()` method, if available.

    Returns
    -------
    None

    Notes
    -----
    The function:
    - Builds a graph of contig interactions.
    - Maps each contig to its corresponding MAG.
    - Calculates average intra- and inter-MAG interaction signals.
    - Groups MAG pairs by quality (e.g., HQ-HQ, HQ-MQ, HQ-contaminated).
    - Applies log transformation and normalization to the interaction signals for downstream analysis.
    """
    logger.info("Calculating the backgroung noise of the experience.")

    # Construction of the contigs graph
    logger.info("Construction of the interaction's graph...")

    G = build_interaction_graph(network_data)

    contig_to_mag = map_contigs_to_mags(mags)

    ########### INTER-MAG INTERACTIONS COMPUTATION #################
    # Caculate inter-MAGs signals
    inter_mag_signals = {}
    for contig1, contig2, data in G.edges(data=True):
        signal = data["weight"]
        mag1 = contig_to_mag.get(contig1)
        mag2 = contig_to_mag.get(contig2)
        if mag1 and mag2 and mag1 != mag2:
            key = tuple(sorted([mag1, mag2]))
            inter_mag_signals.setdefault(key, []).append(signal)
            mags[mag1].add_inter_signal(mag2, signal)
            mags[mag2].add_inter_signal(mag1, signal)
    logger.info(f"{len(inter_mag_signals)} paires of MAGs have inter-MAG interactions.")

    # Write in the text file
    logger.info(f"Writting inter-MAG signals results in {inter_output_file}")
    with open(inter_output_file, "w") as f:
        f.write("MAG1\tMAG2\tinter-MAG mean signal\n")
        for (mag1, mag2), signals in inter_mag_signals.items():
            mean_signal = sum(signals) / len(signals) if signals else 0

            f.write(f"{mag1:<24}\t{mag2:<24}\t{mean_signal:.6f}\n")

    # Construction of a  dataframe for analyses
    all_mag_pairs = list(combinations(mags.keys(), 2))
    rows = []

    for mag1, mag2 in all_mag_pairs:
        key = tuple(sorted([mag1, mag2]))
        signals = inter_mag_signals.get(key, [])
        mean_signal = sum(signals) / len(signals) if signals else 0
        comp1 = completeness_data[completeness_data["Unnamed: 0"] == mag1]
        comp2 = completeness_data[completeness_data["Unnamed: 0"] == mag2]

        c1 = comp1["Completeness"].values[0] if not comp1.empty else np.nan
        c2 = comp2["Completeness"].values[0] if not comp2.empty else np.nan
        r1 = comp1["Redundancy"].values[0] if not comp1.empty else np.nan
        r2 = comp2["Redundancy"].values[0] if not comp2.empty else np.nan

        rows.append(
            {
                "MAG1": mag1,
                "MAG2": mag2,
                "Total_Signal": mean_signal,
                "Completeness1": c1,
                "Completeness2": c2,
                "Contamination1": r1,
                "Contamination2": r2,
            }
        )

    inter_df = pd.DataFrame(rows)

    # Calculation of the mean of  Completeness/Contamination of two MAGs
    inter_df["Completeness"] = inter_df[["Completeness1", "Completeness2"]].mean(axis=1)
    inter_df["Contamination"] = inter_df[["Contamination1", "Contamination2"]].mean(axis=1)

    def classify_mag(row):
        c, r = row["Completeness"], row["Contamination"]
        if pd.isna(c) or pd.isna(r):
            return "Unknown"
        elif r > 1.1:
            return "Contaminated"
        elif c > 0.9 and r <= 1.05:
            return "Complete"
        elif c > 0.9 and r < 1.1 and r > 1.05:
            return "HQ"
        elif c > 0.7 and r < 1.1:
            return "MQ"
        elif c > 0.5 and r < 1.1:
            return "LQ"
        elif c < 0.5 and r < 1.1:
            return "PQ"

    inter_df["MAG_Quality"] = inter_df.apply(classify_mag, axis=1)

    print(inter_df.head())
    # Log Transformation and normalisation
    inter_df["Log_Signal"] = np.log10(inter_df["Total_Signal"].replace(0, np.nan))
    inter_df["Log_Signal"] = inter_df["Log_Signal"].fillna(-10)
    inter_df["Normalized_Log"] = (inter_df["Log_Signal"] - inter_df["Log_Signal"].min()) / (
        inter_df["Log_Signal"].max() - inter_df["Log_Signal"].min()
    )

    ######################INTRA-MAG INTERACTIONS COMPUTATION#################
    intra_mag_means = []
    mag_names = []
    completeness_list = []
    contamination_list = []
    logger.info(f"Writting intra-MAG signals results in {intra_output_file}")
    with open(intra_output_file, "w") as f:
        f.write("MAG\tIntra-MAG mean interaction\n")
        for mag_name, mag in mags.items():
            intra_signals = []
            contigs = set(mag.get_contigs().values())
            for contig1, contig2 in combinations(contigs, 2):
                if G.has_edge(contig1, contig2):
                    intra_signals.append(G[contig1][contig2]["weight"])
            mean_signal = sum(intra_signals) / len(intra_signals) if intra_signals else 0

            intra_mag_means.append(mean_signal)
            mag_names.append(mag_name)
            f.write(f"{mag_name:<24}\t{mean_signal:.6f}\n")
            mag.intra_signals = intra_signals

            comp_info = completeness_data[completeness_data["Unnamed: 0"] == mag_name]
            if not comp_info.empty:
                completeness_list.append(comp_info["Completeness"].values[0])
                contamination_list.append(comp_info["Redundancy"].values[0])
            else:
                completeness_list.append(np.nan)
                contamination_list.append(np.nan)

    intra_df = pd.DataFrame(
        {
            "MAG": mag_names,
            "Intra_MAG_signal": intra_mag_means,
            "Completeness": completeness_list,
            "Contamination": contamination_list,
        }
    )

    intra_df["MAG_Quality"] = intra_df.apply(classify_mag, axis=1)
    # Log10 normalisation
    intra_df["Log_Mean"] = np.log10(intra_df["Intra_MAG_signal"] + 1e-10)
    min_val = intra_df["Log_Mean"].min()
    max_val = intra_df["Log_Mean"].max()
    intra_df["Normalized_Log"] = (intra_df["Log_Mean"] - min_val) / (max_val - min_val)
    # Adding a  "Type" column to the dataframe
    intra_df["Type"] = "Intra"
    inter_df["Type"] = "Inter"
    # Harmonisation of the columns
    inter_df_plot = inter_df[["MAG_Quality", "Normalized_Log", "Type"]]
    intra_df_plot = intra_df[["MAG_Quality", "Normalized_Log", "Type"]]
    # Concatenate both intra and inter dataframes
    combined_df = pd.concat([intra_df_plot, inter_df_plot], ignore_index=True)

    # Plotting the  boxplots and the violin plot
    fig, axes = plt.subplots(2, 1, figsize=(16, 16))
    # Boxplot
    sns.boxplot(data=combined_df, x="MAG_Quality", y="Normalized_Log", hue="Type", ax=axes[0])
    axes[0].set_title("Grouped Boxplots -  Normalised Signal")

    # Violin plot
    sns.violinplot(
        data=combined_df, x="MAG_Quality", y="Normalized_Log", hue="Type", split=True, inner="quart", palette="Set2", ax=axes[1]
    )
    axes[1].set_title("Violin Plot - Normalised Signal")

    # Ajustements finaux
    for ax in axes:

        ax.set_xlabel("MAGs Quality")
        ax.set_ylabel("log10 normalised Signal")

    plt.legend(title="Interaction Type")
    plt.tight_layout()
    plt.savefig(image_file)
    plt.close(fig)
    return image_file


def compute_mge_mag_interactions(
    network_data,
    mge_mags,
    mags,
    output_file="mge_mag_interactions.txt",
    image_file="mge_mag_histogram.png",
    interaction_threshold=10.0,
    min_interacting_contigs=5,
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

    G = build_interaction_graph(network_data)

    contig_to_mag = map_contigs_to_mags(mags)

    interaction_results = []

    for mgemag_name, mgemag_obj in mge_mags.items():
        mgemag_contigs = set(mgemag_obj.get_contigs().values())
        interaction_signals = {}
        interacting_contigs_by_mag = {}

        for contig1 in mgemag_contigs:
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
                    interaction_results.append((mgemag_name, mag_name, signal_sum, percent_signal, interacting_contig_count))

    with open(output_file, "w") as f:
        f.write("mgeMAG\tMAG\tTotal Signal\t% of signal\tInteracting Contigs\n")
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
        axes[0].set_xlabel("mge-MAG total interaction Signal  (log scale)")
        axes[0].set_ylabel("Frequences")
        axes[0].set_title(f"Signaux Distribution (\u2265 {interaction_threshold}%)")

        # Histogramm of interaction percentage
        bins = range(1, math.ceil(max(percentage_counts)) + 1 + 1)  # +2 pour être équivalent à +2 dans ta ligne
        axes[1].hist(percentage_counts, bins=bins, edgecolor="black", align="left")
        axes[1].set_xlabel("Taux d'association")
        axes[1].set_ylabel("Fréquence")
        axes[1].set_title("Distribution des taux d'association des MGEMAGs à leur MAGs")

        plt.tight_layout()
        fig.savefig(image_file)
        plt.close(fig)

    return image_file


def main(
    contig_data_file, network_data_file, completeness_data_file, interaction_threshold, min_interacting_contigs
):  # contig_data_file = contig_data_final_reformarted, network_data_file = network_{0...}.txt, completeness_data_file = "bin_summary.txt

    # Loading data and instantiating objects.
    print("Loading data and instantiating objects...")
    contig_data = pd.read_csv(contig_data_file, sep="\t")
    completeness_data = pd.read_csv(completeness_data_file, sep="\t")
    network_data = pd.read_csv(network_data_file, sep="\t", names=["contig1", "contig2", "signal"])
    mags, mge_mags = create_bin(contig_data)

    # Evaluation of intra- and inter-MAG interactions.
    print("Evaluating of the background noise of the experience...")
    evaluate_experience_noise(
        network_data,
        mags,
        completeness_data,
        inter_output_file="inter_mag_interactions.txt",
        intra_output_file="intra_mag_interactions.txt",
        image_file="mags_interactions_by_quality.png",
    )
    print("Evaluation terminated!!!")
    print("Computing mgeMAG and MAG interactions and association of the mge to it host...")
    compute_mge_mag_interactions(network_data, mge_mags, mags, interaction_threshold=10.0, min_interacting_contigs=5)
    print("Association terminated!!!")


# Call to the main function
if __name__ == "__main__":
    # data loading
    network_data_file = "../host_test/network_7.txt"
    contig_data_file = "../host_test/contig_data_final_reformated.txt"
    completeness_data_file = "../host_test/bin_summary.txt"
    interaction_threshold = 10
    min_interacting_contigs = 5
    # calling main
    main(contig_data_file, network_data_file, completeness_data_file, interaction_threshold, min_interacting_contigs)
