"""The files of a run beside its map (see metator.outputs): the FASTA files of
its genomes, as their contigs and stitched, and the cool files of its MAGs;
and what cannot write them (metator.cli.write_map_outputs)."""

import cooler
import metacooler
import pandas as pd
import pytest

import metator.cli as mtcli
import metator.metacool as mtmc
import metator.outputs as mto
from tests import synthetic

# The contigs of the map of these tests, and their sequences: c1 and c3 in the
# genome mag_1, c2 and c4 in none.
SEQUENCES = {"c1": "A" * 100, "c2": "C" * 50, "c3": "G" * 80, "c4": "T" * 80}


def _map(directory):
    """A .metacool map of the contigs of SEQUENCES, one bin per contig, c1 and c3
    in the genome mag_1, and the FASTA file of their assembly, in directory."""
    pairs = directory / "a.pairs"
    pairs.write_text("## pairs format v1.0\nr1\tc1\t10\tc3\t20\t+\t-\nr2\tc1\t10\tc1\t90\t+\t-\nr3\tc2\t5\tc3\t60\t+\t-\n")
    bins = pd.DataFrame({"chrom": list(SEQUENCES), "start": 0, "end": [len(sequence) for sequence in SEQUENCES.values()]})
    metacool = directory / "map.metacool"
    assignment = {"c1": "mag_1", "c2": "ND", "c3": "mag_1", "c4": "ND"}
    metacooler.from_pairs(str(metacool), pairs=str(pairs), bins=bins, contig_genome_map=assignment, ps=False)
    assembly = directory / "assembly.fa"
    assembly.write_text("".join(f">{name} a contig\n{sequence}\n" for name, sequence in SEQUENCES.items()))
    return metacool, assembly


def _logged(logged):
    """What a test logged, but the bookkeeping (DEBUG)."""
    return [record for record in logged() if record[0] != "DEBUG"]


def test_write_fasta_files(tmp_path):
    """Each genome, side by side: its contigs, the records of the assembly as
    they are, in the order of the map, and stitched, one record named after it,
    its contigs the longest first, joined by 200 N."""
    metacool, assembly = _map(tmp_path)
    out = tmp_path / "mags"
    assert mto.write_fasta_files(str(metacool), str(assembly), str(out), genomes=["mag_1"]) == ["mag_1"]
    assert sorted(path.name for path in out.iterdir()) == ["mag_1.fa", "mag_1.stitched.fa"]
    assert (out / "mag_1.fa").read_text() == f">c1 a contig\n{SEQUENCES['c1']}\n>c3 a contig\n{SEQUENCES['c3']}\n"
    assert synthetic.read_fasta(out / "mag_1.stitched.fa") == {"mag_1": SEQUENCES["c1"] + "N" * 200 + SEQUENCES["c3"]}


def test_write_fasta_files_as_contigs_only(tmp_path):
    """Without stitched, the FASTA file of each group as its contigs only, as the
    validation writes the final clusters (fasta/mags_with-mges): in the order
    given."""
    metacool, assembly = _map(tmp_path)
    out = tmp_path / "final"
    written = mto.write_fasta_files(str(metacool), str(assembly), str(out), groups={"f1": ["c3", "c2"]}, stitched=False)
    assert written == ["f1"] and [path.name for path in out.iterdir()] == ["f1.fa"]
    assert (out / "f1.fa").read_text() == f">c3 a contig\n{SEQUENCES['c3']}\n>c2 a contig\n{SEQUENCES['c2']}\n"


def test_round_and_mge_groups(tmp_path):
    """The clusters of a round of a map, by name, their contigs in the order of
    the map; those of an MGE round, the longest first, ties by name, as
    metator mge writes them; None for a round that the map lacks."""
    metacool, _ = _map(tmp_path)
    rounds = [("final", "mag_1", ["c1", "c3"]), ("mge", "MetaTOR_MGE_00001", ["c4", "c2", "c3"])]
    metacooler.write_clusters(
        str(metacool),
        [{"name": name, "tool": "metaTOR", "step": name} for name, _, _ in rounds],
        pd.DataFrame({"round": [name for name, _, _ in rounds], "name": [cluster for _, cluster, _ in rounds]}),
        pd.DataFrame(
            [(name, cluster, contig) for name, cluster, contigs in rounds for contig in contigs],
            columns=["round", "cluster", "contig"],
        ),
    )
    mc = mtmc.open_map(str(metacool))
    order = [str(name) for name in mc.contigs["name"]]
    assert mto.round_groups(mc, "final") == {"mag_1": sorted(["c1", "c3"], key=order.index)}
    assert mto.mge_groups(mc, "mge") == {"MetaTOR_MGE_00001": ["c3", "c4", "c2"]}
    assert mto.round_groups(mc, "mge_random") is None and mto.mge_groups(mc, "mge_random") is None


def test_write_fasta_files_of_nothing(tmp_path):
    """No genome, or no group: nothing is written, not even the directory (for
    no genome named, metacooler would write every genome of the map)."""
    metacool, assembly = _map(tmp_path)
    assert mto.write_fasta_files(str(metacool), str(assembly), str(tmp_path / "genomes"), genomes=[]) == []
    assert mto.write_fasta_files(str(metacool), str(assembly), str(tmp_path / "groups"), groups={}) == []
    assert not (tmp_path / "genomes").exists() and not (tmp_path / "groups").exists()


def test_cool_files_that_cooler_balance_cannot_balance(tmp_path, monkeypatch, logged):
    """A cool file whose every bin cooler balance masks (here, too few bins) has
    weights all NaN; one that cooler balance fails on keeps its raw contacts,
    without weights. Both are kept, with a warning."""
    metacool, _ = _map(tmp_path)
    path = tmp_path / "masked" / "mag_1.cool"
    assert mto.write_cool_files(str(metacool), str(path.parent), ["mag_1"]) == [str(path)]
    clr = cooler.Cooler(str(path))
    assert clr.chromnames == ["c1", "c3"] and clr.bins()[:]["weight"].isna().all()
    masked = (
        f"Every bin of {path} is masked by its balancing (too few bins, or contacts, for cooler balance): its weights "
        "are all NaN, and its raw contacts are kept."
    )
    assert _logged(logged) == [("WARNING", masked)]

    def fails(*args, **kwargs):
        raise ValueError("cannot")

    monkeypatch.setattr("cooler.balance_cooler", fails)
    path = tmp_path / "unbalanced" / "mag_1.cool"
    assert mto.write_cool_files(str(metacool), str(path.parent), ["mag_1"]) == [str(path)]
    clr = cooler.Cooler(str(path))
    assert clr.chromnames == ["c1", "c3"] and "weight" not in clr.bins()[:].columns
    unbalanced = f"{path} is not balanced, and keeps its raw contacts only: ValueError: cannot"
    assert _logged(logged) == [("WARNING", masked), ("WARNING", unbalanced)]


def test_write_map_outputs_refuses_what_it_cannot_write(tmp_path, logged):
    """What the outputs of a run that write_map_outputs writes cannot be written
    from (here, the FASTA files of a MAG, from an assembly that lacks one of its
    contigs) is raised, naming the reason; the caller (metator export, the
    step export of pipeline) logs it (see run_export)."""
    metacool, _ = _map(tmp_path)
    assembly = tmp_path / "other.fa"
    assembly.write_text(f">c1\n{SEQUENCES['c1']}\n>c2\n{SEQUENCES['c2']}\n")
    with pytest.raises(mtmc.MetacoolError, match="c3"):
        mtcli.write_map_outputs(str(metacool), str(assembly), str(tmp_path / "out"), str(tmp_path))
    assert not [record for record in _logged(logged) if record[0] == "ERROR"]
