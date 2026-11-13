# Tests for the contact_map module and the contactmap command.

import shutil
from pathlib import Path

import metator.commands as mtc
import metator.contact_map as mtcm
import pandas as pd
import pytest

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


def test_contactmap_default_outdir(tmp_path, monkeypatch):
    """Without --outdir the map goes to ./contact_map_<name>."""
    assembly = tmp_path / "assembly.fa"
    shutil.copy("tests_data/subset.assembly.fa", assembly)
    contig_data = Path(CONTIG_DATA).resolve()
    pairs = Path("tests_data/subset.pairs.gz").resolve()
    monkeypatch.chdir(tmp_path)
    args = f"-a {assembly} -c {contig_data} -e DpnII -n 1 -O overlapping_bin -s 0 {pairs}"
    mtc.Contactmap(args.split(" "), {}).execute()
    outdir = tmp_path / "contact_map_1"
    assert (outdir / "1.fa").stat().st_size > 0
