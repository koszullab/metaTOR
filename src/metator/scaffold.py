#!/usr/bin/env python3
# coding: utf-8

"""Scaffolds bins using HiC contact maps.

General fonctions and classes to reorganize and scaffold genomes.

Function in this module:
    - get_scaffolds
    - parallel_scaffold
    - write_info_frags
    - write_scaffolds_fasta

Class to handle the scaffolding:
    Scaffolder
"""


from dataclasses import dataclass, field
import metator.io as mio
import networkx as nx
import numpy as np
import pypairix
import pyfastx
from metator.regions import Bin, Scaffold
from metator.log import logger
from typing import List, Optional
from os.path import join
import warnings
import sys
from io import StringIO


def get_scaffolds(
    bin_name: str,
    fasta_bin: str,
    pair_files: List[str],
    out_fasta: str,
    out_info: str,
    threshold: float = 0.1,
    threads: int = 1,
    window_size: int = 5000,
    junctions: Optional[str] = None,
):
    """
    Function to scaffold the a bin.

    Parameters
    ----------
    bin_name : str
        Name of the bin
    fasta_bin : str
        Path to the fasta sequences of the bin contigs.
    pair_files : str
        List of path to the pair file to pair file. If not sorted and indexed it
        will be sorted and indexed.
    out_fasta : str
        Path to write the fasta output.
    out_info : str
        Path to write the info frags file.
    threshold : float
        Threshold score to consider an association between two contigs.
        [Default 0.1]
    threads : int
        Number of threads to allocate. [Default: 1]
    window_size : int
        Size of the window in base pair to use as edge window to scaffold.
        [Default: 5000]
    junctions : Optional[str]
        Sequence to put between two contigs. [Default: None]
    """
    fasta = pyfastx.Fasta(fasta_bin)
    bin = Bin(bin_name, fasta)
    pairs_data = []
    for pair_file in pair_files:
        pairs_data.append(mio.get_pairs_data(pair_file, threads))
    scaffolder = Scaffolder(bin, pairs_data, window_size, threshold,)
    scaffolder.get_matrix()
    scaffolder.get_normalize_values()
    scaffolder.normalize_matrix()
    scaffolder.build_graph()
    scaffolds = scaffolder.solve_graph()
    if junctions is None:
        junctions = ""
    write_info_frags(scaffolds, out_info)
    write_scaffolds_fasta(scaffolds, out_fasta, fasta, junctions)


def parallel_scaffold(
    bin_name, final_fasta_dir, alignment_files, scaffold_fasta_dir, junctions,
):
    """
    Task function for parallel worker to scaffold the a bin.

    Parameters
    ----------
    bin_name : str
        Name of the bin
    final_fasta_dir : str
        Folder path to the fasta sequences of the bin contigs.
    alignment_files : str
        List of path to the pair file to pair file. If not sorted and indexed it
        will be sorted and indexed.
    scaffold_fasta_dir : str
        Folder path to write the fasta and table output.
    out_info : str
        Path to write the info frags file.
    junctions : Optional[str]
        Sequence to put between two contigs. [Default: None]
    """
    logger.info(f"Scaffolding {bin_name}...")
    get_scaffolds(
        bin_name,
        join(final_fasta_dir, f"{bin_name}.fa"),
        alignment_files,
        join(scaffold_fasta_dir, f"{bin_name}_scaffolded.fa"),
        join(scaffold_fasta_dir, f"{bin_name}_info_frags.txt"),
        threshold=0.05,
        threads=1,
        window_size=5000,
        junctions=junctions,
    )


def write_info_frags(scaffolds: List[Scaffold], outfile: str):
    """Write the info of the scaffolded fragments.

    Parameters
    ----------
    scaffolds : List[Scaffold]
        List of the scaffold to write.
    outfile : str
        Path where to write the output table.
    """
    # Write th eoutput one entry for each scaffold. Then one entry whit each
    # contig with their name, positions and orientations.
    with open(outfile, "w") as out:
        for scaffold in scaffolds:
            out.write(f"{scaffold.get_fasta_entry()}\n")
            for contig in scaffold.contigs:
                orientation = -1 if contig.is_reverse else 1
                out.write(f"{contig.name}\t{0}\t{len(contig)}\t{orientation}\n")


def write_scaffolds_fasta(
    scaffolds: List[Scaffold],
    outfile: str,
    fasta: pyfastx.Fasta,
    junctions: Optional[str] = None,
):
    """Write the scaffolds in a fasta file.

    Parameters
    ----------
    scaffolds : List[Scaffold]
        List of the scaffold to write.
    outfile : str
        Path where to write the output fasta.
    fasta : pyfastx.Fasta
        Object containing the initial fasta.
    junctions : Optional[str]
        Sequence to put between two contigs. [Default: None]
    """
    with open(outfile, "w") as out:
        for scaffold in scaffolds:
            out.write(f"{scaffold.get_fasta_entry()}\n")
            is_junction = False
            for seq in scaffold.get_seq(fasta):
                if is_junction:
                    out.write(f"{junctions}{seq}\n")
                else:
                    out.write(f"{seq}\n")
                is_junction = True


@dataclass()
class Scaffolder:
    """
    Provides a class to scaffold contigs from a metator bin. The scaffolding is
    based on the physicla contacts between

    Attributes
    ----------
    bin : Bin
        Bin object from metator.regions.
    pairs_data : List[pypairix.open]
        pypairix object wth pair information.
    window_size : int
        Size in base pair of the window to take to bridge contigs together.
        [Default: 5000]
    threshold : float
        Threshold to consider a bridge. [Default: 0.1]
    default_norm : float
        Default value to use if contigs too small to normalize. [Default: None]
    """

    bin: Bin
    pairs_data: List[pypairix.open]
    window_size: int = field(default=5000)
    threshold: float = field(default=0.1)
    default_norm: Optional[float] = field(default=None)

    def __post_init__(self):
        # Set the length of the bin (contigs number).
        self.set_length()
        # Set length of the genome in base pair.
        self.size = len(self.bin)
        # Set the contigs edges bed positions.
        self.set_contigs_edges()

    def set_length(self):
        """Set number of breakpoints."""
        self.length = len(self.bin.contigs)

    def set_contigs_edges(self):
        """Set positions in bed format of each breakpoints."""
        self.positions = np.empty(2 * self.length, dtype="S256")
        for i, contig in enumerate(self.bin.contigs):
            # Start/end are defined for small fragments
            end = (
                contig.start + self.window_size
                if len(contig) > self.window_size
                else len(contig)
            )
            self.positions[2 * i] = f"{contig.name}:{0}-{end}"
            start = (
                contig.end - self.window_size
                if len(contig) > self.window_size
                else 0
            )
            self.positions[2 * i + 1] = f"{contig.name}:{start}-{contig.end}"

    def get_matrix(self):
        """Get the contact count between breakpoints."""
        # Compute the score between breakpoints.
        self.score_matrice = np.zeros((2 * self.length, 2 * self.length))
        for i in range(self.length * 2):
            for j in range(self.length * 2):
                if i != j:
                    for pairs_data_object in self.pairs_data:
                        pairs = pairs_data_object.querys2D(
                            f"{(self.positions[i]).decode()}|{(self.positions[j]).decode()}",
                            1,
                        )
                        for _pair in pairs:
                            self.score_matrice[i, j] += 1
        self.score_matrice = np.triu(
            self.score_matrice + self.score_matrice.T, k=1
        )
        return self.score_matrice

    def get_normalize_values(self):
        """Get the values to normalize"""
        # Compute the score to normalize values.
        self.normalize = np.zeros(2 * self.length)
        if self.default_norm is None:
            if not hasattr(self, "score_matrice"):
                self.get_matrix()
            self.default_norm = np.max(self.score_matrice) * np.sqrt(2)
        for i, contig in enumerate(self.bin.contigs):
            # Compute the score at the same distance in intra to have a local
            # approximation of the signal to use to normalize.
            if len(contig) > 2 * self.window_size:
                # Left Breakpoint.
                start, end = self.window_size, 2 * self.window_size
                query_left_bp1 = (self.positions[2 * i]).decode()
                query_right_bp1 = f"{contig.name}:{start}-{end}"
                # Right Breakpoint.
                start, end = (
                    contig.end - 2 * self.window_size,
                    contig.end - self.window_size,
                )
                query_left_bp2 = f"{contig.name}:{start}-{end}"
                query_right_bp2 = (self.positions[2 * i + 1]).decode()

                # Extract pairs.
                for pairs_data_object in self.pairs_data:
                    pairs = pairs_data_object.querys2D(
                        f"{query_left_bp1}|{query_right_bp1}", 1
                    )
                    for _pair in pairs:
                        self.normalize[2 * i] += 1
                    pairs = pairs_data_object.querys2D(
                        f"{query_left_bp2}|{query_right_bp2}", 1
                    )
                    for _pair in pairs:
                        self.normalize[2 * i + 1] += 1

            # Defined to default if contigs are too small to have equivalent
            # intra signal.
            else:
                self.normalize[2 * i] = self.default_norm
                self.normalize[2 * i + 1] = self.default_norm
        return self.normalize

    def normalize_matrix(self):
        """Normalize the matrix."""
        # Test whether matrix and normlized values have been computed.
        if not hasattr(self, "normalize"):
            self.get_normalize_values()
        if not hasattr(self, "score_matrice"):
            self.get_matrix()

        # Normalized the matrix.
        self.score_matrice_norm = np.zeros((2 * self.length, 2 * self.length))
        for i in range(self.length * 2):
            for j in range(i + 1, self.length * 2):
                if (j - i > 2) or (i % 2 == 1):
                    self.score_matrice_norm[i, j] = self.score_matrice[
                        i, j
                    ] / np.sqrt(self.normalize[i] * self.normalize[j])
        self.score_matrice_norm = (
            self.score_matrice_norm + self.score_matrice_norm.T
        )
        return self.score_matrice_norm

    def build_graph(self):
        """Transform score matrix in graph. Each represents a breakpoint and can
        have only two edges."""
        # Test whether the matrix have been normalized or already resolved.
        try:
            if np.max(self.score_matrice_norm) < self.threshold:
                self.normalize_matrix()
        except AttributeError:
            self.normalize_matrix()

        self.graph = nx.Graph()
        # Add edges between the breakpoints of the same contigs (the two edges).
        for i in range(self.length):
            self.graph.add_edge(2 * i, 2 * i + 1)
            self.score_matrice_norm[2 * i, 2 * i + 1] = 0

        while np.max(self.score_matrice_norm) > self.threshold:
            # If one breakpoint have only one possible association. Take it
            # first even if there are better association.
            if (
                np.sum(
                    (
                        np.sum(self.score_matrice_norm > self.threshold, axis=1)
                        == 1
                    )
                )
                > 0
            ):
                if (
                    np.max(
                        self.score_matrice_norm[
                            (
                                np.sum(
                                    self.score_matrice_norm > self.threshold,
                                    axis=1,
                                )
                                == 1
                            ),
                            :,
                        ]
                    )
                    > 5 * self.threshold
                ):
                    row1, col1 = np.where(
                        self.score_matrice_norm
                        == np.max(
                            self.score_matrice_norm[
                                (
                                    np.sum(
                                        self.score_matrice_norm
                                        > self.threshold,
                                        axis=1,
                                    )
                                    == 1
                                ),
                                :,
                            ]
                        )
                    )
                else:
                    row1, col1 = np.where(
                        self.score_matrice_norm
                        == np.max(self.score_matrice_norm)
                    )
            # Else take the highest association score
            else:
                row1, col1 = np.where(
                    self.score_matrice_norm == np.max(self.score_matrice_norm)
                )
            # Always take first indice first if multiple positions found.
            row1 = row1[0]
            col1 = col1[0]
            # Add the edge in the graph.
            self.graph.add_edge(row1, col1)
            # Remove the  other scores as an association have been found.
            self.score_matrice_norm[row1, :] = 0
            self.score_matrice_norm[col1, :] = 0
            self.score_matrice_norm[:, row1] = 0
            self.score_matrice_norm[:, col1] = 0

    def solve_graph(self):
        """Solve the connected components in the graph into scaffolds."""
        # Test whether the scaffolded graph have been computed.
        if not hasattr(self, "graph"):
            self.build_graph()

        # Iterates on each scaffold as a connected component.
        n = 0
        for components in nx.connected_components(self.graph):
            n += 1
            scaf_name = f"metator_scaffold_{n:04d}"
            subgraph = nx.subgraph(self.graph, components)
            # Case of an isolated contigs add it alone as a scaffold.
            # We choose to keep the forward strand of the contig.
            if len(subgraph.nodes) == 2:
                contig_id = [node for node in subgraph.nodes][0] // 2
                self.bin.contigs[contig_id].is_reverse = False
                scaffold = Scaffold(
                    scaf_name, [self.bin.contigs[contig_id]], False
                )
                self.bin.scaffolds.append(scaffold)
            # Look at dangling residues to search for edges of the scaffold.
            else:
                ordered_nodes = np.zeros(len(subgraph.nodes))
                dangling_residues = [
                    node
                    for node in subgraph.nodes
                    if subgraph.degree(node) == 1
                ]
                # Either two dangling residues (linear)
                if len(dangling_residues) == 2:
                    circular = False
                    # Choose arbitrarly the start and end.
                    start = dangling_residues[0]
                    end = dangling_residues[1]
                    ordered_nodes[0] = start
                    k = 0
                # Or none (circular)
                else:
                    # Take first contig as a start as it's circular.
                    circular = True
                    start = list(subgraph.nodes)[0]
                    end = start + 1 if start % 2 == 0 else start - 1
                    ordered_nodes[0], ordered_nodes[1] = end, start
                    k = 1
                prev = end
                # Build a list of the breakpoint ID.
                while start != end:
                    k += 1
                    for i, j in subgraph.edges(start):
                        if (i != start) & (i != prev):
                            prev = start
                            start = i
                            # Avoid to add the last one in circular case.
                            if (start != end) or not circular:
                                ordered_nodes[k] = start
                            break
                        elif (j != start) & (j != prev):
                            prev = start
                            start = j
                            # Avoid to add the last one in circular case.
                            if (start != end) or not circular:
                                ordered_nodes[k] = start
                            break
                # Iterates on the ordered nodes to build the scaffold
                len_scaffold = len(ordered_nodes) // 2
                contigs = []
                for i in range(len_scaffold):
                    bp_id_1 = ordered_nodes[2 * i]
                    bp_id_2 = ordered_nodes[2 * i + 1]
                    if bp_id_1 > bp_id_2:
                        reverse = True
                    else:
                        reverse = False
                    contig_id = int(bp_id_1 // 2)
                    self.bin.contigs[contig_id].is_reverse = reverse
                    contigs.append(self.bin.contigs[contig_id])
                scaffold = Scaffold(scaf_name, contigs, circular)
                self.bin.scaffolds.append(scaffold)
        return self.bin.scaffolds
