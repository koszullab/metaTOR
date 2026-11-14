# Scaffolding tests on hand-made breakpoint graphs: how the connected
# components of the graph are read into scaffolds (no Hi-C pairs needed).

import igraph as ig
import numpy as np
import pyfastx

import metator.scaffold as mts
from metator.regions import Bin


def _scaffolder(tmp_path, n_contigs):
    """A Scaffolder of ``n_contigs`` contigs c0, c1...: contig i has the
    breakpoints 2i (start) and 2i + 1 (end)."""
    fasta = tmp_path / "bin.fa"
    fasta.write_text("".join(f">c{i}\n{'ACGT' * 5}\n" for i in range(n_contigs)))
    return mts.Scaffolder(Bin("bin", pyfastx.Fasta(str(fasta))), [])


def test_build_graph(tmp_path):
    """Each breakpoint is joined to the other end of its contig and to at
    most one breakpoint of another contig."""
    scaffolder = _scaffolder(tmp_path, 4)
    scores = np.zeros((8, 8))
    scores[1, 2], scores[3, 4], scores[5, 6] = 2.0, 1.5, 1.2  # c0 - c1 - c2 - c3
    scores[1, 4] = 0.5  # a weaker contact, between taken breakpoints
    scaffolder.score_matrice_norm = scores + scores.T
    scaffolder.build_graph()
    graph = scaffolder.graph
    assert graph.vcount() == 8
    assert sorted(tuple(sorted(edge)) for edge in graph.get_edgelist()) == [
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7)
    ]  # fmt: skip
    (scaffold,) = scaffolder.solve_graph()
    assert not scaffold.circular
    assert [repr(c) for c in scaffold.contigs] == ["c0:0-20:+", "c1:0-20:+", "c2:0-20:+", "c3:0-20:+"]


def _scaffolder_on_graph(tmp_path, n_contigs, joins):
    """A Scaffolder (see _scaffolder) whose breakpoint graph links the two ends
    of each contig and the breakpoints of ``joins``."""
    scaffolder = _scaffolder(tmp_path, n_contigs)
    edges = [(2 * i, 2 * i + 1) for i in range(n_contigs)] + joins
    scaffolder.graph = ig.Graph(n=2 * n_contigs, edges=edges)
    return scaffolder


def test_solve_graph_order_and_orientation(tmp_path):
    """Scaffolds are numbered in the order of their first contig in the bin; a
    linear scaffold is read from its end with the lowest breakpoint, a circular
    one from the start of its first contig."""
    joins = [
        (9, 5),  # end of c4 - end of c2
        (1, 6),  # end of c0 - start of c3
        (7, 10),  # end of c3 - start of c5
        (11, 0),  # end of c5 - start of c0
    ]
    scaffolds = _scaffolder_on_graph(tmp_path, 6, joins).solve_graph()
    assert [(s.name, s.circular) for s in scaffolds] == [
        ("metator_scaffold_0001", True),
        ("metator_scaffold_0002", False),
        ("metator_scaffold_0003", False),
    ]
    assert [repr(c) for c in scaffolds[0].contigs] == ["c0:0-20:-", "c5:0-20:-", "c3:0-20:-"]
    assert [repr(c) for c in scaffolds[1].contigs] == ["c1:0-20:+"]
    assert [repr(c) for c in scaffolds[2].contigs] == ["c2:0-20:+", "c4:0-20:-"]


def test_solve_graph_start_of_a_small_scaffold(tmp_path):
    """The start of a scaffold holding few of the bin's contigs follows the
    same rule. metaTOR 1.3.10 iterated their breakpoints as a set: these two
    scaffolds started from c8."""
    # Linear: end of c8 - end of c5.
    scaffolds = _scaffolder_on_graph(tmp_path, 20, [(17, 11)]).solve_graph()
    (linear,) = [s for s in scaffolds if len(s.contigs) == 2]
    assert [repr(c) for c in linear.contigs] == ["c5:0-20:+", "c8:0-20:-"]
    # Circular: the ends of c5 and c8 joined pairwise.
    scaffolds = _scaffolder_on_graph(tmp_path, 20, [(17, 11), (16, 10)]).solve_graph()
    (circular,) = [s for s in scaffolds if s.circular]
    assert [repr(c) for c in circular.contigs] == ["c5:0-20:-", "c8:0-20:+"]
