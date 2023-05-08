"""Structures and methods to represent SV breakpoints and contig. This submodule
defines a BreakPoint class representing connections between positions on 
separate contig. Each breakpoint consists of two positions. Each of position is 
located on a contigment. Positions have a sign representing what end of a contig
they're on.

Classes to handle genomic regions:
    - Position
    - Contig
    - Breakpoint
    - Bin 
    - Scaffold
"""

from __future__ import annotations
import copy
from dataclasses import dataclass, field
import numpy as np
import pyfastx
from typing import List, Iterator, Optional, Tuple


@dataclass(order=True)
class Position:
    """A single position in the assembly, defined by a bin (Optional), a
    scaffold (Optional) a contig, genomic position and optionally a sign.

    Attributes
    ----------
    contig : str
        The contig name.
    coord : int
        The 0-based coordinate on the contig.
    sign : Optional[bool]
        Whether the position is on the 5' (-) or 3' (+) side.
    """

    contig: str
    coord: int
    #  None means unknown, True means 3'
    sign: Optional[bool] = field(default=None, compare=False)

    def __repr__(self) -> str:
        sign_symbol = {None: "", True: "+", False: "-"}
        return f"{self.contig}:{self.coord}:{sign_symbol[self.sign]}"

    def has_sign(self) -> bool:
        """Whether sign information is available."""
        return not self.sign is None


@dataclass
class Contig:
    """A region representing an assembled DNA sequence. Coordinates are 0-based
    and right-open.

    Attributes
    ----------
    name : str
        Name of the contig sequence.
    start : int
        Coordinate where the contig starts in the bin (smallest).
    end : int
        Coordinate where the contig ends (largest).
    is_reverse : bool
    bin : Optional[str]
        Bin in which the contig is binned. None if unbinned.
    scaffold : Optional[str]
        Scaffold in which the contig is. None if unbinned.
    """

    name: str
    start: int
    end: int
    is_reverse: bool = field(default=False)
    # None means unbinned.
    bin: Optional[str] = field(default=None, compare=False)
    # None means unscaffolded.
    scaffold: Optional[str] = field(default=None, compare=False)

    def __post_init__(self,):
        if self.end < self.start:
            raise ValueError("end cannot be smaller than start.")
        if self.start < 0:
            raise ValueError("Coordinates cannot be negative.")

    def __len__(self) -> int:
        return self.end - self.start

    def __repr__(self) -> str:
        sign = "-" if self.is_reverse else "+"
        return f"{self.name}:{self.start}-{self.end}:{sign}"

    def __hash__(self) -> int:
        return hash(str(self))

    def middle(self) -> int:
        return (self.start + self.end) // 2

    def intersect(self, other: Contig) -> int:
        """Return the length of intersection with another contigment. If there is
        no intersection, returns 0."""
        same_contig = self.name == other.name
        separate = (self.start > other.end) or (self.end < other.start)
        if same_contig and not separate:
            overlap_start = max(self.start, other.start)
            overlap_end = min(self.end, other.end)
            return overlap_end - overlap_start
        return 0

    def split(self, rel_pos: int) -> Tuple[Contig, Contig]:
        """Split the contigment at a positiion relative from the start."""
        if rel_pos < 0:
            raise ValueError("Cannot split on negative coord")
        if rel_pos > (self.end - self.start):
            raise ValueError(
                f"Cannot split {self} at position {rel_pos}: Beyond contigment end."
            )
        left_contig, right_contig = copy.copy(self), copy.copy(self)
        if not self.is_reverse:
            left_contig.end = self.start + rel_pos
            right_contig.start = self.start + rel_pos
            return (left_contig, right_contig)
        else:
            left_contig.start = self.end - rel_pos
            right_contig.end = self.end - rel_pos
            return (left_contig, right_contig)

    def flip(self):
        """Change contigment's sign."""
        self.is_reverse = not self.is_reverse


class BreakPoint:
    """Defines a breakpoint in the genome.
    A breakpoint associates 2 different genomic positions from independent
    contigs. pos1 and contig1 represent the left side of the breakpoint.

    Attributes
    ----------
    pos1, pos2:
        The genomic positions connected by the breakpoint.
    contig1, contig2:
        Segments of DNA connected by the breakpoint
    """

    def __init__(
        self, pos1: Position, pos2: Position, type: str = "UNK",
    ):

        # Properties are use to set/get contigments instead
        self._contig1 = None
        self._contig2 = None
        self.pos1, self.pos2 = pos1, pos2
        self.type = type

    def __repr__(self) -> str:
        p1 = self.pos1
        p2 = self.pos2
        if self.has_contigs():
            p1 = f"[{self.contig1}]{p1}"
            p2 = f"{p2}[{self.contig2}]"
        return f"{p1}|{p2} {self.type}"

    @property
    def signs(self) -> Optional[Tuple[str, str]]:
        """Return signs if available."""
        if self.has_signs():
            sign1 = "+" if self.pos1.sign else "-"
            sign2 = "+" if self.pos2.sign else "-"
            return sign1, sign2

    @staticmethod
    def _check_pos_at_contig_end(pos: Position, contig: Contig) -> bool:
        """Check whether a position is at the correct end of a contig."""
        right_contig = pos.contig == contig.name
        if pos.has_sign():
            if pos.sign:
                if contig.is_reverse:
                    right_pos = pos.coord == contig.start
                else:
                    right_pos = pos.coord == contig.end
            else:
                if contig.is_reverse:
                    right_pos = pos.coord == contig.end
                else:
                    right_pos = pos.coord == contig.start
        else:
            right_pos = pos.coord in [contig.start, contig.end]
        return right_contig & right_pos

    @property
    def contig1(self):
        """Get the left contig of the breakpoint"""
        return self._contig1

    @contig1.setter
    def contig1(self, contig: Contig):
        if self._check_pos_at_contig_end(self.pos1, contig):
            self._contig1 = contig
        else:
            raise ValueError(
                f"Attempted to set contig1 to {contig} which does not match "
                f"breakpoint pos1 of {self.pos1}."
            )

    @property
    def contig2(self) -> Contig:
        """Get the right contig of the breakpoint"""
        return self._contig2

    @contig2.setter
    def contig2(self, contig: Contig):
        if self._check_pos_at_contig_end(self.pos2, contig):
            self._contig2 = contig
        else:
            raise ValueError(
                f"Attempted to set contig2 to {contig}, which does not match "
                f"breakpoint pos2 of {self.pos2}."
            )

    def has_signs(self) -> bool:
        """Whether sign information is available"""
        return self.pos1.has_sign() & self.pos2.has_sign()

    def has_contigs(self) -> bool:
        """Whether contigments information is available."""
        return (self.contig1 is not None) & (self.contig2 is not None)

    def can_connect(self, other: BreakPoint, min_intersect: int = 1000) -> bool:
        """Whether two breakpoints could be connected (i.e. whether they could
        share a contigment). This only checks for one-way connection: the second
        contigment of self should overlap the first contigment of other."""
        both_contigs = self.has_contigs() and other.has_contigs()
        both_signs = self.has_signs() and other.has_signs()
        compat_contigs = self.pos2.contig == other.pos1.contig
        if not (both_contigs and both_signs):
            raise ValueError(
                "Cannot connect breakpoints without contigment or sign information."
            )
        compat_coords = False
        if self.pos2.sign != other.pos1.sign:
            if self.intersect(other) >= min_intersect:
                compat_coords = True

        return compat_contigs & compat_coords

    def intersect(self, other: BreakPoint) -> int:
        """Amount of overlap between two breakpoint's contigments"""
        if not (self.has_contigs() and other.has_contigs()):
            raise ValueError(
                "Cannot compute overlap if without contigment information."
            )
        # TODO: Could use signs to skip some comparisons
        intersect = 0
        for source in [self.contig1, self.contig2]:
            for target in [other.contig1, other.contig2]:
                intersect = max(intersect, source.intersect(target))
        return intersect


@dataclass
class Bin:
    """Representation of a bin as a collection of contigs. Each contig
    represents a (0-based, right-open) region of the binned genome.

    Attributes
    ----------

    bin_name : str
        Name of the bin.
    fasta : pyfastx.Fasta
        Fasta sequences of the bin.
    """

    bin_name: str
    fasta: pyfastx.Fasta

    def __post_init__(self):
        # list rather than np.array, due to better insert/append performance.
        self.contigs = []
        for seq in self.fasta:
            self.contigs.append(
                Contig(seq.name, 0, len(seq), bin=self.bin_name)
            )

        # Compute breakpoint positions.
        self.breakpoints = np.cumsum(
            [0] + [len(contig) for contig in self.contigs]
        )
        self.scaffolds = []

    def __len__(self) -> int:
        """Returns the total bin length."""
        return self.breakpoints[-1]

    def get_contigs_counts(self) -> int:
        """Returns the number of bins."""
        return len(self.contigs)

    @property
    def boundaries(self) -> List[int]:
        """Get array of contig boundaries, from the start to the end of the
        bin."""
        # Memorize whether contigs have changed to avoid recomputing the same
        # values.
        contig_hash = hash(tuple(self.contigs))
        try:
            if self._contig_hash == contig_hash:
                changed = False
            else:
                changed = True
        # On first access, required attrs are generated
        except AttributeError:
            self._contig_hash = contig_hash
            changed = True
        if changed:
            self._bds = np.cumsum(
                [0] + [len(contig) for contig in self.contigs]
            )
        return self._bds

    def get_contig_bounds(self, coord: int) -> Tuple[int, Tuple[int, int]]:
        """Returns the index and boundaries of the contig in which input
        coordinate falls. Return format is (id, (start, end))."""
        bounds = self.boundaries
        if coord >= bounds[-1]:
            raise ValueError(
                f"Coordinate out of bounds: {self.bin_name}:{coord}"
            )
        contig_id = max(0, np.searchsorted(bounds, coord, side="right") - 1)
        return (contig_id, (bounds[contig_id], bounds[contig_id + 1]))

    def get_contigs_between(self, start: int, end: int) -> List[Contig]:
        """Returns a list of the contigs between 2 positions in the bin."""
        result = []
        bounds = self.boundaries
        if start < 0 or end > bounds[-1]:
            raise ValueError(
                f"Coordinate out of bounds: {self.bin_name}:{start}-{end}"
            )

        (
            contig_id_left,
            (start_bound_left, end_bound_left,),
        ) = self.get_contig_bounds(start)
        (
            contig_id_right,
            (start_bound_right, end_bound_right,),
        ) = self.get_contig_bounds(end)

        # If the start and end falls in the same contig
        if contig_id_left == contig_id_right:
            return [
                Contig(
                    self.bin_name,
                    self.contigs[contig_id_left].start
                    + (start - start_bound_left),
                    self.contigs[contig_id_left].start
                    + (end - start_bound_left),
                    self.contigs[contig_id_left].is_reverse,
                )
            ]
        # Else
        # Start falls in a contig, trim it
        print("start index: ", start)
        print("start bound: ", self.get_contig_bounds(start))
        left_contig = self.contigs[contig_id_left]
        print("left_contig = ", left_contig)

        left_contig_trimmed = Contig(
            self.bin_name,
            left_contig.start + (start - start_bound_left),
            left_contig.end,
            left_contig.is_reverse,
        )
        print("left contig trimmed = ", left_contig_trimmed)
        result.append(left_contig_trimmed)

        # Ends falls in a contig, trim it
        right_contig = self.contigs[contig_id_right]
        print("right_contig = ", right_contig)

        right_contig_trimmed = Contig(
            self.bin_name,
            right_contig.start,
            right_contig.start + end - start_bound_right,
            right_contig.is_reverse,
        )
        print("right contig trimmed = ", right_contig_trimmed)
        # Add all contigs between the first and the last
        for contig_i in range(contig_id_left + 1, contig_id_right):
            contig = Contig(
                self.bin_name,
                self.contigs[contig_i].start,
                self.contigs[contig_i].end,
                self.contigs[contig_i].is_reverse,
            )
            result.append(contig)

        # Finally add last contig
        result.append(right_contig_trimmed)

        return result

    def insert_contigs(self, position: int, contigs: List[Contig]):
        insert_pos = position
        for i in range(len(contigs)):
            self.insert(insert_pos, contigs[i])
            insert_pos += len(contigs[i])

    def clean_contigs(self):
        """Purge 0-length contigs."""
        self.contigs = [contig for contig in self.contigs if len(contig)]

    def insert(self, position: int, contig_ins: Contig):
        """Updates contigs by inserting a sequence in the bin."""
        bounds = self.boundaries
        # Append after the end of bin.
        if position == len(self):
            contig_id = len(bounds)
        else:
            contig_id, (contig_start, _) = self.get_contig_bounds(position)
        if position in bounds:
            # Insertion right between two contigs, add a contig.
            self.contigs.insert(contig_id, contig_ins)
        else:
            # Insertion inside a contig, split it and add contig in between.
            contig_l, contig_r = self.contigs.pop(contig_id).split(
                position - contig_start
            )
            for contig in [contig_r, contig_ins, contig_l]:
                self.contigs.insert(contig_id, contig)
        bp = BreakPoint(
            Position(self.bin_name, contig_ins.start,),
            Position(self.bin_name, contig_ins.end),
            "INS",
        )
        self.breakpoints.append(bp)

    def invert(self, start: int, end: int):
        """Updates contigs by inverting a portion of the bin.
        The interval is 0-based and right open [start;end[."""
        s_contig_id, (s_contig_start, _) = self.get_contig_bounds(start)
        e_contig_id, (e_contig_start, _) = self.get_contig_bounds(end - 1)
        s_start_dist = start - s_contig_start
        e_start_dist = end - e_contig_start

        # Inversion inside a single contig.: Split it in 3 and invert middle.
        if s_contig_id == e_contig_id:
            inv_size = end - start
            contig_l, contig_mr = self.contigs.pop(s_contig_id).split(
                s_start_dist
            )
            contig_m, contig_r = contig_mr.split(inv_size)
            contig_m.flip()
            for contig in [contig_r, contig_m, contig_l]:
                self.contigs.insert(s_contig_id, contig)
        else:
            # Split contig where inversion starts, we'll flip the right part.
            start_l, start_r = self.contigs.pop(s_contig_id).split(s_start_dist)
            for contig in [start_r, start_l]:
                self.contigs.insert(s_contig_id, contig)
            s_contig_id += 1
            e_contig_id += 1
            # Split contig where inversion ends we'll flip the left part.
            end_l, end_r = self.contigs.pop(e_contig_id).split(e_start_dist)
            for contig in [end_r, end_l]:
                self.contigs.insert(e_contig_id, contig)
            e_contig_id += 1
            # If contigs are contained in the inversion, invert and flip them.
            for contig_id in range(s_contig_id, e_contig_id):
                self.contigs[contig_id].flip()
            self.contigs[s_contig_id:e_contig_id] = self.contigs[
                e_contig_id - 1 : s_contig_id - 1 : -1
            ]
        self.clean_contigs()
        bp = BreakPoint(
            Position(self.bin_name, start), Position(self.bin_name, end), "INV"
        )
        self.breakpoints.append(bp)

    def delete(self, start: int, end: int):
        """Updates contigs by deleting a portion of the chromosome.
        The interval is 0-based and right open [start;end[."""
        s_contig_id, (s_contig_start, _) = self.get_contig_bounds(start)
        e_contig_id, (_, e_contig_end) = self.get_contig_bounds(end - 1)
        del_size = end - start
        start_dist = start - s_contig_start
        end_dist = e_contig_end - end
        # Deletion contained in a single contig: split it and trim right part
        if e_contig_id == s_contig_id:
            start_l, start_r = self.contigs.pop(s_contig_id).split(start_dist)
            start_r.start += del_size
            for contig in [start_r, start_l]:
                self.contigs.insert(s_contig_id, contig)
        # Deletion spans multiple contigs
        else:
            # Deletion starts in contig, end gets trimmed
            self.contigs[s_contig_id].end = (
                self.contigs[s_contig_id].start + start_dist
            )

            # Contigs contained in deletion disappear
            for contig_id in range(s_contig_id + 1, e_contig_id):
                curr_start = self.contigs[contig_id].start
                self.contigs[contig_id].end = curr_start

            from copy import copy

            ori_end = copy(self.contigs[e_contig_id])
            # Deletion ends in contig, trim left side
            self.contigs[e_contig_id].start = (
                self.contigs[e_contig_id].end - end_dist
            )

        bp = BreakPoint(
            Position(self.bin_name, start), Position(self.bin_name, end), "DEL"
        )
        self.breakpoints.append(bp)

    def get_seq(self, fasta: pyfastx.Fasta) -> Iterator[str]:
        """Retrieve the bin sequence, as a generator yielding the sequence by
        contig."""
        self.clean_contigs()
        for contig in self.contigs:
            strand = "-" if contig.is_reverse else "+"
            # Note: fasta.fetch is 1-based...
            yield fasta.fetch(
                contig.name,
                (int(contig.start + 1), (contig.end)),
                strand=strand,
            )

    def get_breakpoints(self) -> Iterator[BreakPoint]:
        """Retrieve a generator yielding breakpoints in the chromosome."""
        self.clean_contigs()
        for contig_id in range(0, len(self.contigs) - 1):
            contig1, contig2 = self.contigs[contig_id : contig_id + 2]
            p1 = Position(contig1.bin, contig1.end, not contig1.is_reverse)
            p2 = Position(contig2.bin, contig2.start, contig2.is_reverse)
            breakpoint = BreakPoint(p1, p2)
            breakpoint.contig1 = contig1
            breakpoint.contig2 = contig2
            yield breakpoint


@dataclass
class Scaffold:
    """A collection of contigs scaffolded together. The contigs are ordered and
    oriented according to the genome.

    Attributes
    ----------
    name : str
        Name of the scaffold.
    contigs : List of Contig
        List of the contigs which are oredered and orientated.
    circular : bool
        Either the scaffold is circular or not. [Default: False]
    """

    name: str
    contigs: List[Contig]
    circular: bool = field(default=False)

    def __len__(self) -> int:
        """Length defined as the length of all contigs"""
        return np.abs(np.sum([len(contig) for contig in self.contigs]))

    def clean_contigs(self):
        """Purge 0-length contigs."""
        self.contigs = [contig for contig in self.contigs if len(contig)]

    def get_seq(self, fasta: pyfastx.Fasta) -> Iterator[str]:
        """Retrieve the scaffold sequence, as a generator yielding the sequence
        by contig."""
        self.clean_contigs()
        for contig in self.contigs:
            strand = "-" if contig.is_reverse else "+"
            # Note: fasta.fetch is 1-based...
            yield fasta.fetch(
                contig.name,
                (int(contig.start + 1), (contig.end)),
                strand=strand,
            )

    def get_fasta_entry(self) -> str:
        """Generate the fasta entry for the scaffold sequence."""
        circular = "circular" if self.circular else "linear"
        length = self.__len__()
        return f">{self.name} {circular} length: {length}"
