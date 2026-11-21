"""Seeded, deterministic mock metagenome for metaTOR tests."""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
_BASES = np.frombuffer(b"ACGT", dtype="S1")


@dataclass(frozen=True)
class GenomeSpec:
    """A genome of the mock community."""

    name: str
    gc: float
    abundance: float
    length: int
    circular: bool = True


@dataclass(frozen=True)
class MgeSpec:
    """A mobile genetic element: a circular molecule bound to a host genome."""

    name: str
    host: str
    gc: float
    abundance: float
    contig_lengths: Tuple[int, ...]


@dataclass
class Contig:
    """One assembled contig: a strand of a slice of a molecule."""

    name: str
    molecule: str
    start: int
    end: int
    reverse: bool
    seq: str

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class Molecule:
    """A genome or an MGE, with its contigs in molecule order."""

    name: str
    kind: str
    gc: float
    abundance: float
    circular: bool
    seq: str
    host: Optional[str] = None
    contigs: List[Contig] = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.seq)


# Five genomes from 30% to 70% GC, of distinct abundances, ~560 kb in total.
DEFAULT_GENOMES = (
    GenomeSpec("genome_1", gc=0.30, abundance=1.0, length=120_000),
    GenomeSpec("genome_2", gc=0.40, abundance=2.0, length=100_000),
    GenomeSpec("genome_3", gc=0.50, abundance=4.0, length=140_000),
    GenomeSpec("genome_4", gc=0.60, abundance=6.0, length=90_000),
    GenomeSpec("genome_5", gc=0.70, abundance=3.0, length=110_000),
)

# A phage assembled in two contigs and a plasmid in one, with different hosts.
DEFAULT_MGES = (
    MgeSpec("phage_1", host="genome_2", gc=0.45, abundance=6.0, contig_lengths=(7_000, 5_000)),
    MgeSpec("plasmid_1", host="genome_4", gc=0.55, abundance=10.0, contig_lengths=(6_000,)),
)


def random_sequence(rng: np.random.Generator, length: int, gc: float) -> str:
    """Random DNA of the given GC fraction."""
    weights = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return _BASES[rng.choice(4, size=length, p=weights)].tobytes().decode()


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def gc_percent(seq: str) -> float:
    return 100 * (seq.count("G") + seq.count("C")) / len(seq)


def cut_lengths(
    rng: np.random.Generator,
    total: int,
    median: float,
    sigma: float,
    bounds: Tuple[int, int],
) -> List[int]:
    """Cut ``total`` bp into log-normal lengths clipped to ``bounds``."""
    low, high = bounds
    lengths: List[int] = []
    remaining = total
    while remaining > 0:
        length = int(np.clip(rng.lognormal(math.log(median), sigma), low, high))
        if remaining - length < low:
            length = remaining
        lengths.append(length)
        remaining -= length
    return lengths


class Metagenome:
    """A simulated metaHiC experiment with its ground truth."""

    def __init__(self, molecules: List[Molecule], contigs: List[Contig], pairs: pd.DataFrame):
        self.molecules = molecules
        self.contigs = contigs
        self.pairs = pairs

    # -- Truth -------------------------------------------------------------

    @property
    def genomes(self) -> List[Molecule]:
        return [m for m in self.molecules if m.kind == "genome"]

    @property
    def mges(self) -> List[Molecule]:
        return [m for m in self.molecules if m.kind == "mge"]

    def molecule(self, name: str) -> Molecule:
        return next(m for m in self.molecules if m.name == name)

    @property
    def truth(self) -> pd.DataFrame:
        """One row per contig, in assembly order."""
        rows = []
        for molecule in self.molecules:
            for order, contig in enumerate(molecule.contigs):
                rows.append(
                    {
                        "Name": contig.name,
                        "Length": contig.length,
                        "Genome": molecule.name,
                        "Kind": molecule.kind,
                        "Host": molecule.host or molecule.name,
                        "Start": contig.start,
                        "End": contig.end,
                        "Reverse": contig.reverse,
                        "Order": order,
                        "GC": gc_percent(contig.seq),
                        "Abundance": molecule.abundance,
                    }
                )
        table = pd.DataFrame(rows).set_index("Name", drop=False)
        table.index.name = None
        return table.loc[[c.name for c in self.contigs]]

    def genome_of(self) -> Dict[str, str]:
        """Contig name -> source molecule name."""
        return {c.name: c.molecule for c in self.contigs}

    # -- Files -------------------------------------------------------------

    def write_fasta(self, path, names: Optional[Sequence[str]] = None, width: int = 80) -> Path:
        """Write the assembly, or only the contigs ``names`` in that order."""
        by_name = {c.name: c for c in self.contigs}
        selected = [by_name[n] for n in names] if names is not None else self.contigs
        with open(path, "w") as fasta:
            for contig in selected:
                fasta.write(f">{contig.name}\n")
                for i in range(0, contig.length, width):
                    fasta.write(contig.seq[i : i + width] + "\n")
        return Path(path)

    def write_depth(self, path, rng_seed: int = 1) -> Path:
        """Write a jgi_summarize_bam_contig_depths table (assembly order)."""
        rng = np.random.default_rng(rng_seed)
        abundance = {m.name: m.abundance for m in self.molecules}
        with open(path, "w") as depth:
            depth.write("contigName\tcontigLen\ttotalAvgDepth\tsynthetic.bam\tsynthetic.bam-var\n")
            for contig in self.contigs:
                value = 10 * abundance[contig.molecule] * (1 + 0.05 * rng.standard_normal())
                depth.write(f"{contig.name}\t{contig.length}\t{value:.4f}\t{value:.4f}\t{value / 2:.4f}\n")
        return Path(path)

    def write_pairs(self, path) -> Path:
        """Write the pairs file with the header that metaTOR 1.x wrote."""
        with open(path, "w") as out:
            out.write("## pairs format v1.0\n")
            out.write("#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n")
            out.write("#sorted: readID\n")
            out.write("#shape: upper triangle\n")
            for contig in self.contigs:
                out.write(f"#chromsize: {contig.name} {contig.length}\n")
            self.pairs.to_csv(out, sep="\t", header=False, index=False)
        return Path(path)

    def write_truth(self, path) -> Path:
        self.truth.to_csv(path, sep="\t", index=False)
        return Path(path)

    def write_all(self, directory) -> Dict[str, Path]:
        """Write assembly, depth, pairs and truth in ``directory``."""
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        return {
            "assembly": self.write_fasta(directory / "assembly.fa"),
            "depth": self.write_depth(directory / "depth.txt"),
            "pairs": self.write_pairs(directory / "synthetic.pairs"),
            "truth": self.write_truth(directory / "truth.tsv"),
        }


def make_metagenome(
    seed: int = 0,
    genomes: Sequence[GenomeSpec] = DEFAULT_GENOMES,
    mges: Sequence[MgeSpec] = DEFAULT_MGES,
    n_pairs: int = 60_000,
    contig_median: float = 10_000,
    contig_sigma: float = 0.7,
    contig_bounds: Tuple[int, int] = (2_000, 40_000),
    p_noise: float = 0.02,
    mge_host_share: float = 0.4,
    min_distance: int = 100,
    decay: float = 1.0,
    intra_uniform: float = 0.3,
    shuffle: bool = True,
    strands: str = "random",
    prefix: str = "contig",
) -> Metagenome:
    """Simulate a metaHiC experiment on a mock community."""
    rng = np.random.default_rng(seed)

    molecules: List[Molecule] = []
    for spec in genomes:
        seq = random_sequence(rng, spec.length, spec.gc)
        molecule = Molecule(spec.name, "genome", spec.gc, spec.abundance, spec.circular, seq)
        lengths = cut_lengths(rng, spec.length, contig_median, contig_sigma, contig_bounds)
        _add_contigs(molecule, lengths, rng, strands)
        molecules.append(molecule)
    for spec in mges:
        seq = random_sequence(rng, sum(spec.contig_lengths), spec.gc)
        molecule = Molecule(spec.name, "mge", spec.gc, spec.abundance, True, seq, host=spec.host)
        _add_contigs(molecule, list(spec.contig_lengths), rng, strands)
        molecules.append(molecule)

    # Name the contigs in assembly order.
    contigs = [c for m in molecules for c in m.contigs]
    order = rng.permutation(len(contigs)) if shuffle else np.arange(len(contigs))
    contigs = [contigs[i] for i in order]
    for rank, contig in enumerate(contigs, start=1):
        contig.name = f"{prefix}_{rank:03d}"

    pairs = _simulate_pairs(molecules, contigs, rng, n_pairs, p_noise, mge_host_share, min_distance, decay, intra_uniform)
    return Metagenome(molecules, contigs, pairs)


def _add_contigs(molecule: Molecule, lengths: List[int], rng: np.random.Generator, strands: str) -> None:
    if strands not in ("random", "forward", "alternate"):
        raise ValueError(f"Unknown strands mode: {strands}")
    start = 0
    for rank, length in enumerate(lengths):
        draw = rng.random()
        reverse = bool({"random": draw < 0.5, "forward": False, "alternate": rank % 2 == 1}[strands])
        seq = molecule.seq[start : start + length]
        if reverse:
            seq = reverse_complement(seq)
        molecule.contigs.append(Contig("", molecule.name, start, start + length, reverse, seq))
        start += length


def _power_law(rng, low, high, size, exponent):
    """Distances in [low, high] with a density proportional to s^-exponent."""
    u = rng.uniform(0, 1, size)
    if exponent == 1:
        return np.exp(math.log(low) + u * (math.log(high) - math.log(low))).astype(int)
    a, b = low ** (1 - exponent), high ** (1 - exponent)
    return ((a + u * (b - a)) ** (1 / (1 - exponent))).astype(int)


def _simulate_pairs(molecules, contigs, rng, n_pairs, p_noise, mge_host_share, min_distance, decay, intra_uniform):
    rank = {id(c): i for i, c in enumerate(contigs)}
    names = np.array([c.name for c in contigs])
    weights = np.array([m.abundance * m.length for m in molecules], dtype=float)
    weights /= weights.sum()
    index = {m.name: i for i, m in enumerate(molecules)}
    lengths = np.array([m.length for m in molecules])

    # Contig boundaries on each molecule, to map a position to a contig.
    starts = [np.array([c.start for c in m.contigs]) for m in molecules]
    ranks = [np.array([rank[id(c)] for c in m.contigs]) for m in molecules]
    reverse = np.array([c.reverse for c in contigs])
    offsets = np.array([c.start for c in contigs])
    sizes = np.array([c.length for c in contigs])

    def to_contig(mol_ids, positions):
        """Molecule positions (0-based) -> contig ranks, 1-based positions."""
        out_rank = np.empty(len(positions), dtype=int)
        for mol in np.unique(mol_ids):
            sel = mol_ids == mol
            k = np.searchsorted(starts[mol], positions[sel], side="right") - 1
            out_rank[sel] = ranks[mol][k]
        local = positions - offsets[out_rank]
        local = np.where(reverse[out_rank], sizes[out_rank] - 1 - local, local)
        return out_rank, local + 1

    mol1, x1, mol2, x2 = [], [], [], []
    counts = rng.multinomial(n_pairs, weights)
    for mol, n in enumerate(counts):
        molecule = molecules[mol]
        length = molecule.length
        pos1 = rng.integers(0, length, n)
        partner = np.full(n, mol)
        pos2 = np.empty(n, dtype=int)
        draw = rng.random(n)
        noise = draw < p_noise
        to_host = np.zeros(n, dtype=bool)
        if molecule.host is not None:
            to_host = (draw >= p_noise) & (draw < p_noise + mge_host_share)
        cis = ~(noise | to_host)

        # Noise: a random molecule, weighted like the first ends.
        partner[noise] = rng.choice(len(molecules), size=noise.sum(), p=weights)
        pos2[noise] = rng.integers(0, lengths[partner[noise]])

        # MGE <-> host contacts land anywhere on the host.
        if to_host.any():
            host = index[molecule.host]
            partner[to_host] = host
            pos2[to_host] = rng.integers(0, lengths[host], to_host.sum())

        # Intra-molecule contacts: the contact probability of two loci decays
        # as s^-decay. On a circle every distance up to L/2 has L pairs of
        # loci, so the distance of a pair follows that power law; on a line
        # only L - s pairs of loci are s apart, hence the rejection, and the
        # first end is then placed where the second one fits.
        n_cis = int(cis.sum())
        sign = rng.choice([-1, 1], size=n_cis)
        if molecule.circular:
            distance = _power_law(rng, min_distance, length // 2, n_cis, decay)
            pos2[cis] = (pos1[cis] + sign * distance) % length
        else:
            distance = _power_law(rng, min_distance, length - 1, n_cis, decay)
            rejected = rng.random(n_cis) >= (length - distance) / length
            while rejected.any():
                distance[rejected] = _power_law(rng, min_distance, length - 1, rejected.sum(), decay)
                rejected[rejected] = rng.random(rejected.sum()) >= (length - distance[rejected]) / length
            left = rng.integers(0, length - distance)
            pos1[cis] = np.where(sign > 0, left, left + distance)
            pos2[cis] = np.where(sign > 0, left + distance, left)
        # Long-range plateau: some pairs join two random loci of the molecule,
        # whatever their distance (compaction of the nucleoid).
        plateau = np.flatnonzero(cis)[rng.random(n_cis) < intra_uniform]
        pos2[plateau] = rng.integers(0, length, len(plateau))

        mol1.append(np.full(n, mol))
        x1.append(pos1)
        mol2.append(partner)
        x2.append(pos2)

    mol1, x1 = np.concatenate(mol1), np.concatenate(x1)
    mol2, x2 = np.concatenate(mol2), np.concatenate(x2)
    shuffle = rng.permutation(len(x1))
    rank1, pos1 = to_contig(mol1[shuffle], x1[shuffle])
    rank2, pos2 = to_contig(mol2[shuffle], x2[shuffle])
    strands = np.array(["+", "-"])
    strand1 = strands[rng.integers(0, 2, len(pos1))]
    strand2 = strands[rng.integers(0, 2, len(pos1))]

    # Upper-triangle order, as in the pairs that metaTOR 1.x wrote.
    swap = (rank1 > rank2) | ((rank1 == rank2) & (pos1 > pos2))
    rank1, rank2 = np.where(swap, rank2, rank1), np.where(swap, rank1, rank2)
    pos1, pos2 = np.where(swap, pos2, pos1), np.where(swap, pos1, pos2)
    strand1, strand2 = np.where(swap, strand2, strand1), np.where(swap, strand1, strand2)

    return pd.DataFrame(
        {
            "readID": [f"read{i:07d}" for i in range(1, len(pos1) + 1)],
            "chr1": names[rank1],
            "pos1": pos1,
            "chr2": names[rank2],
            "pos2": pos2,
            "strand1": strand1,
            "strand2": strand2,
        }
    )


# -- Stub miComplete --------------------------------------------------------

MICOMPLETE_COLUMNS = [
    "Length",
    "GC-content",
    "Present Markers",
    "Completeness",
    "Redundancy",
    "Weighted completeness",
    "Weighted redundancy",
    "Contigs",
    "N50",
    "L50",
    "N90",
    "L90",
    "CDs",
    "Markers",
]


def read_fasta(path) -> Dict[str, str]:
    """Minimal FASTA reader: name -> sequence."""
    sequences: Dict[str, List[str]] = {}
    name = None
    with open(path) as fasta:
        for line in fasta:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                sequences[name] = []
            elif name is not None:
                sequences[name].append(line)
    return {name: "".join(parts) for name, parts in sequences.items()}


def bin_quality(contig_names: Sequence[str], metagenome: Metagenome) -> Tuple[float, float]:
    """Completeness and redundancy of a bin, as an oracle miComplete would see
    them: the completeness is the share of the best represented genome present
    in the bin, and every other genome present adds its share to the redundancy
    (two complete genomes: redundancy 2)."""
    truth = metagenome.truth
    genome_length = {g.name: g.length for g in metagenome.genomes}
    present: Dict[str, int] = {}
    for name in contig_names:
        row = truth.loc[name]
        if row["Kind"] == "genome":
            present[row["Genome"]] = present.get(row["Genome"], 0) + int(row["Length"])
    if not present:
        return 0.0, 1.0
    shares = sorted((length / genome_length[g] for g, length in present.items()), reverse=True)
    return shares[0], 1.0 + sum(shares[1:])


def _nx(lengths: List[int], fraction: float) -> Tuple[int, int]:
    """Nx and Lx of contig lengths."""
    lengths = sorted(lengths, reverse=True)
    target, total = fraction * sum(lengths), 0
    for rank, length in enumerate(lengths, start=1):
        total += length
        if total >= target:
            return length, rank
    return 0, 0


def write_micomplete_summary(fasta_dir, outfile, metagenome: Metagenome) -> pd.DataFrame:
    """Write the merged miComplete table that ``validation.merge_micomplete``
    writes, for the bins (fasta files) of ``fasta_dir``, with the oracle
    qualities of :func:`bin_quality`."""
    rows = {}
    for path in sorted(os.listdir(fasta_dir)):
        if not path.endswith((".fa", ".fasta", ".fna")):
            continue
        sequences = read_fasta(os.path.join(fasta_dir, path))
        completeness, redundancy = bin_quality(list(sequences), metagenome)
        lengths = [len(s) for s in sequences.values()]
        n50, l50 = _nx(lengths, 0.5)
        n90, l90 = _nx(lengths, 0.9)
        total = sum(lengths)
        gc = sum(gc_percent(s) * len(s) for s in sequences.values()) / total
        rows[os.path.splitext(path)[0]] = [
            float(total),
            round(gc, 2),
            float(round(105 * completeness)),
            round(completeness, 4),
            round(redundancy, 4),
            round(completeness, 4),
            round(redundancy, 4),
            float(len(lengths)),
            float(n50),
            float(l50),
            float(n90),
            float(l90),
            float(total // 1000),
            "Bacteria",
        ]
    with open(outfile, "w") as out:
        out.write("## miComplete\n")
        out.write("## v1.1.1\n")
        out.write("## Weights: Bact105 and Arch131\n")
        out.write("Name\t" + "\t".join(MICOMPLETE_COLUMNS) + "\n")
        for name, values in rows.items():
            out.write(name + "\t" + "\t".join(map(str, values)) + "\n")
    return pd.DataFrame.from_dict(rows, orient="index", columns=MICOMPLETE_COLUMNS)


def fake_micomplete(metagenome: Metagenome, calls: Optional[list] = None) -> Callable:
    """A drop-in replacement of ``metator.validation.micomplete_quality``."""

    def micomplete_quality(fasta_dir, outfile, threads):
        if calls is not None:
            calls.append((str(fasta_dir), str(outfile)))
        write_micomplete_summary(fasta_dir, outfile, metagenome)

    return micomplete_quality


def first_cluster_contaminated(metagenome: Metagenome, kept_whole: bool = False) -> Callable:
    """A stub miComplete (see :func:`fake_micomplete`) that finds the first
    overlapping cluster contaminated (a weighted redundancy of 1.5): the
    validation partitions it again, which makes the round recursive_1, whose
    clusters replace it; with kept_whole, they are too incomplete for that (a
    weighted completeness of 0.1), and the cluster is kept whole."""

    def micomplete_quality(fasta_dir, outfile, threads):
        write_micomplete_summary(fasta_dir, outfile, metagenome)
        # Three comment lines, the header, then a row per cluster.
        with open(outfile) as summary:
            lines = [line.split("\t") for line in summary.read().splitlines()]
        if os.path.basename(os.path.normpath(fasta_dir)) == "overlapping_clusters":
            lines[4][lines[3].index("Weighted redundancy")] = "1.5"
        elif kept_whole:
            for line in lines[4:]:
                line[lines[3].index("Weighted completeness")] = "0.1"
        with open(outfile, "w") as summary:
            summary.write("".join("\t".join(line) + "\n" for line in lines))

    return micomplete_quality


def contaminated_twice(metagenome: Metagenome) -> Callable:
    """A stub miComplete (see :func:`fake_micomplete`) whose validation goes two
    recursive steps deep: the first overlapping cluster is contaminated (a
    weighted completeness of 1, a weighted redundancy of 1.5), its parts
    complete enough to replace it, the first of which is contaminated again at
    step 1, and partitioned again at step 2 (the rounds recursive_1 and
    recursive_2)."""

    def micomplete_quality(fasta_dir, outfile, threads):
        write_micomplete_summary(fasta_dir, outfile, metagenome)
        step = os.path.basename(os.path.normpath(fasta_dir))
        with open(outfile) as summary:
            lines = [line.split("\t") for line in summary.read().splitlines()]
        completeness = lines[3].index("Weighted completeness")
        redundancy = lines[3].index("Weighted redundancy")
        for row, line in enumerate(lines[4:]):
            if step == "overlapping_clusters":
                if row == 0:
                    line[completeness], line[redundancy] = "1.0", "1.5"
            else:
                contaminated = step == "step_1" and row == 0
                line[completeness], line[redundancy] = "0.9", ("1.5" if contaminated else "1.0")
        with open(outfile, "w") as summary:
            summary.write("".join("\t".join(line) + "\n" for line in lines))

    return micomplete_quality


# -- Comparing a binning with the truth -------------------------------------


def _codes(labels: Sequence) -> List[int]:
    index: Dict = {}
    return [index.setdefault(label, len(index)) for label in labels]


def adjusted_rand(truth: Sequence, found: Sequence) -> float:
    """Adjusted Rand index of two labelings of the same contigs."""
    import igraph

    return igraph.compare_communities(_codes(truth), _codes(found), method="adjusted_rand")


def purity(truth: Sequence, found: Sequence) -> float:
    """Share of the contigs that belong to the majority genome of their bin."""
    table = pd.DataFrame({"truth": list(truth), "found": list(found)})
    majority = table.groupby("found")["truth"].agg(lambda genomes: genomes.value_counts().iloc[0])
    return majority.sum() / len(table)


# -- metaTOR 1.x ------------------------------------------------------------


def as_1x(path, out) -> str:
    """A copy, at out, of the contig data file path, its columns named as metaTOR
    1.x named them, after bins (metator.io.CONTIG_DATA_1X_COLUMNS): the table
    of a 1.x run."""
    import metator.io as mio

    names = {new: old for old, new in mio.CONTIG_DATA_1X_COLUMNS.items()}
    header, *rows = Path(path).read_text().split("\n")
    renamed = "\t".join(names.get(name, name) for name in header.split("\t"))
    assert renamed != header, f"{path} has no column that metaTOR 1.x named otherwise"
    Path(out).write_text("\n".join([renamed, *rows]))
    return str(out)
