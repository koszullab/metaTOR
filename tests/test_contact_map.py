# Tests for the contact_map module and the contactmap command.

import shutil
from pathlib import Path

import metator.contact_map as mtcm
import pandas as pd
import pytest
from tests import synthetic

CONTIG_DATA = "tests_data/contig_data_final.txt"


@pytest.mark.parametrize(
    "metator_object, name, column",
    [
        ("core_bin", "1", "Core_bin_ID"),
        ("core_bin", 1, "Core_bin_ID"),
        ("overlapping_bin", "1", "Overlapping_bin_ID"),
        ("overlapping_bin", "metator_00001_00000", "Overlapping_bin_ID"),
    ],
)
def test_bin_ids_match_zero_padded_ids(metator_object, name, column):
    """Bin ids are written zero-padded ("00001"): "1" must select that bin."""
    data = pd.read_csv(CONTIG_DATA, sep="\t")
    expected = data.loc[pd.to_numeric(data[column], errors="coerce") == 1, "Name"]
    assert len(expected) > 0
    obj = mtcm.MetatorObject(metator_object, name, "unused.fa", CONTIG_DATA, [], min_size=0)
    obj.set_contigs()
    assert obj.contigs == expected.tolist()


def test_contactmap_default_outdir(tmp_path, monkeypatch, run_cli):
    """Without --outdir the map goes to ./contact_map_<name>."""
    assembly = tmp_path / "assembly.fa"
    shutil.copy("tests_data/subset.assembly.fa", assembly)
    contig_data = Path(CONTIG_DATA).resolve()
    pairs = Path("tests_data/subset.pairs.gz").resolve()
    monkeypatch.chdir(tmp_path)
    args = f"-a {assembly} -c {contig_data} -e DpnII -n 1 -O overlapping_bin -s 0 {pairs}"
    run_cli("contactmap", *args.split(" "))
    outdir = tmp_path / "contact_map_1"
    assert (outdir / "1.fa").stat().st_size > 0


@pytest.mark.parametrize(
    "metator_object, name, column, expected",
    [
        ("contig", "NODE_1", "Name", "NODE_1"),
        ("core_bin", 3, "Core_bin_ID", "3"),
        ("overlapping_bin", "MetaTOR_7_0", "Overlapping_bin_ID", "7"),
        ("overlapping_bin", "MetaTOR_7", "Overlapping_bin_ID", "7"),
        ("overlapping_bin", "sample_metator_00007_00002", "Overlapping_bin_ID", "00007"),
        ("recursive_bin", "MetaTOR_7_2", "Recursive_bin_ID", "2"),
        ("recursive_bin", "sample_metator_00007_00002", "Recursive_bin_ID", "00002"),
        ("recursive_bin", "MetaTOR_7_0", "Recursive_bin_ID", "0"),  # logged as an error
        ("final_bin", "MetaTOR_7_2", "Final_bin", "MetaTOR_7_2"),
        ("Other", "12", "Other", 12),
        ("Other", "group_a", "Other", "group_a"),
    ],
)
def test_metator_object_names(metator_object, name, column, expected):
    obj = mtcm.MetatorObject(metator_object, name, "unused.fa", CONTIG_DATA, [])
    assert obj.object == column
    assert obj.name == expected


@pytest.mark.parametrize(
    "metator_object, name",
    [
        ("core_bin", "MetaTOR_1"),
        ("overlapping_bin", "MetaTOR_x_0"),
        ("recursive_bin", "MetaTOR_1_x"),
        # A recursive bin ID alone: recursive bins are numbered within their overlapping bin.
        ("recursive_bin", "4"),
        ("recursive_bin", "MetaTOR_4"),
    ],
)
def test_metator_object_invalid_names(metator_object, name):
    with pytest.raises(ValueError):
        mtcm.MetatorObject(metator_object, name, "unused.fa", CONTIG_DATA, [])


def test_recursive_bin_within_its_overlapping_bin(tmp_path):
    """Recursive bins are numbered within their overlapping bin: a recursive
    bin is the contigs of its recursive bin in its overlapping bin. The
    recursive bins of the same ID of every overlapping bin were selected,
    and a prefixed name (--prefix) failed."""
    table = tmp_path / "contig_data_final.txt"
    pd.DataFrame(
        {
            "ID": [1, 2, 3, 4],
            "Name": ["c1", "c2", "c3", "c4"],
            "Size": [10, 20, 30, 40],
            "Overlapping_bin_ID": ["00001", "00001", "00002", "00002"],
            "Recursive_bin_ID": ["00001", "00002", "00001", "00001"],
        }
    ).to_csv(table, sep="\t", index=False)
    obj = mtcm.MetatorObject("recursive_bin", "my_sample_metator_00002_00001", "unused.fa", str(table), [], min_size=0)
    obj.set_contigs()
    assert (obj.contigs, obj.contigs_size) == (["c3", "c4"], [30, 40])
    obj = mtcm.MetatorObject("overlapping_bin", "my_sample_metator_00001_00000", "unused.fa", str(table), [], min_size=0)
    obj.set_contigs()
    assert obj.contigs == ["c1", "c2"]


def test_metator_object_unknown_column():
    obj = mtcm.MetatorObject("Missing_column", "1", "unused.fa", CONTIG_DATA, [])
    with pytest.raises(ValueError):
        obj.set_contigs()


# Contact maps of a bin of the mock metagenome (tests/synthetic.py).


@pytest.fixture
def final_contig_data(metagenome, metagenome_network, tmp_path):
    """contig_data_final.txt with the true genomes as final bins."""
    data = pd.read_csv(metagenome_network["contig_data"], sep="\t")
    bins = {g.name: f"metator_{k:05d}_00000" for k, g in enumerate(metagenome.genomes, start=1)}
    data["Final_bin"] = data.Name.map(metagenome.truth.Host).map(bins)
    path = tmp_path / "contig_data_final.txt"
    data.to_csv(path, sep="\t", index=False)
    return path, bins


def _pairs_between(metagenome, contigs):
    pairs = metagenome.pairs
    return int((pairs.chr1.isin(contigs) & pairs.chr2.isin(contigs)).sum())


def _matrix_total(out_dir, name, mat_fmt):
    """Sum of the contacts of the matrix written by hicstuff."""
    if mat_fmt == "graal":
        return pd.read_csv(out_dir / f"{name}.mat.tsv", sep="\t").iloc[:, 2].sum()
    if mat_fmt == "bg2":
        return pd.read_csv(out_dir / f"{name}.bg2", sep="\t", header=None).iloc[:, 6].sum()
    import cooler

    return cooler.Cooler(str(out_dir / f"{name}.cool")).pixels()[:]["count"].sum()


@pytest.mark.parametrize("mat_fmt", ["graal", "bg2", "cool"])
def test_generate_contact_map(metagenome, metagenome_inputs, metagenome_indexed_pairs, final_contig_data, tmp_path, mat_fmt):
    contig_data, bins = final_contig_data
    name = bins["genome_3"]
    contigs = [c.name for c in metagenome.molecule("genome_3").contigs]
    out, tmp = tmp_path / "out", tmp_path / "tmp"
    out.mkdir()
    tmp.mkdir()
    n_pairs = mtcm.generate_contact_map(
        str(metagenome_inputs["assembly"]),
        str(contig_data),
        "DpnII",
        name,
        [str(metagenome_indexed_pairs)],
        str(out),
        str(tmp),
        mat_fmt=mat_fmt,
        min_size=0,
    )
    assert n_pairs == _pairs_between(metagenome, contigs)
    assert sorted(synthetic.read_fasta(out / f"{name}.fa")) == sorted(contigs)
    assert _matrix_total(out, name, mat_fmt) == n_pairs


def test_contact_map_min_size(metagenome, metagenome_inputs, metagenome_indexed_pairs, final_contig_data, tmp_path):
    """Contigs below --min-size are left out of the map."""
    contig_data, bins = final_contig_data
    large = [c.name for c in metagenome.molecule("genome_3").contigs if c.length >= 15_000]
    assert 0 < len(large) < len(metagenome.molecule("genome_3").contigs)
    out, tmp = tmp_path / "out", tmp_path / "tmp"
    out.mkdir()
    tmp.mkdir()
    n_pairs = mtcm.generate_contact_map(
        str(metagenome_inputs["assembly"]),
        str(contig_data),
        "DpnII",
        bins["genome_3"],
        [str(metagenome_indexed_pairs)],
        str(out),
        str(tmp),
        min_size=15_000,
    )
    assert n_pairs == _pairs_between(metagenome, large)
    assert sorted(synthetic.read_fasta(out / f"{bins['genome_3']}.fa")) == sorted(large)
    assert _matrix_total(out, bins["genome_3"], "graal") == n_pairs


def test_contact_map_without_pairs(metagenome, metagenome_inputs, metagenome_indexed_pairs, final_contig_data, tmp_path):
    """A contig smaller than --min-size: nothing to extract, no map."""
    contig_data, _ = final_contig_data
    small = min(metagenome.contigs, key=lambda c: c.length)
    n_pairs = mtcm.generate_contact_map(
        str(metagenome_inputs["assembly"]),
        str(contig_data),
        "DpnII",
        small.name,
        [str(metagenome_indexed_pairs)],
        str(tmp_path),
        str(tmp_path),
        metator_object="contig",
        min_size=small.length + 1,
    )
    assert n_pairs == 0
    assert not list(tmp_path.glob("*.mat.tsv"))


def test_metator_object_write_fasta(tmp_path, monkeypatch):
    """The fasta of an object, the records of its contigs, is written
    without the pyfastx command, which is not on the PATH here."""
    monkeypatch.setenv("PATH", str(tmp_path))
    # A private copy of the assembly: pyfastx indexes it (.fxi) in place.
    assembly = tmp_path / "assembly.fa"
    shutil.copyfile("tests_data/assembly.fa", assembly)
    obj = mtcm.MetatorObject("contig", "NODE_522", str(assembly), CONTIG_DATA, [])
    obj.contigs = ["NODE_1404", "NODE_522"]
    obj.write_fasta(str(tmp_path))
    names = [line[1:] for line in Path(obj.fasta).read_text().splitlines() if line.startswith(">")]
    assert names == ["NODE_1404", "NODE_522"]
