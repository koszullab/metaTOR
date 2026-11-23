"""The clustering rounds of a run (metator.clusters)."""

import pandas as pd
import pytest

import metator.clusters as mtc


def _rounds():
    rounds = mtc.Rounds()
    rounds.add(mtc.OVERLAPPING, "O1", ["c1", "c2", "c3"], columns={"Completeness": 0.9})
    rounds.add(mtc.OVERLAPPING, "O2", ["c4"], columns={"Completeness": 0.1})
    rounds.close_round(mtc.OVERLAPPING)
    step = mtc.recursive_round(1)
    rounds.add(step, "R1", ["c1", "c2"], parent="O1", columns={"Contigs": 2})
    rounds.add(step, "R2", ["c3"], parent="O1")
    rounds.close_round(step)
    rounds.add(mtc.FINAL, "R1", ["c1", "c2"], parent="R1")
    rounds.add(mtc.FINAL, "O2", ["c4"], parent="O2", columns={"Completeness": 0.1})
    rounds.close_round(mtc.FINAL)
    return rounds


def test_recursive_rounds_are_named_by_step():
    assert mtc.recursive_round(3) == "recursive_3"


def test_a_parent_is_the_latest_cluster_of_its_name():
    """A final cluster's parent is the cluster of the same name of the round it
    comes from: the latest round that holds one, once that round is closed."""
    clusters = pd.DataFrame(_rounds().clusters)
    final = clusters[clusters["round"] == mtc.FINAL].set_index("name")
    assert final.loc["R1", "parent_round"] == "recursive_1"
    assert final.loc["O2", "parent_round"] == "overlapping"
    first = clusters[clusters["round"] == mtc.OVERLAPPING]
    assert first["parent"].isna().all() and first["parent_round"].isna().all()


def test_a_parent_of_no_earlier_round_is_refused():
    rounds = mtc.Rounds()
    rounds.add(mtc.OVERLAPPING, "O1", ["c1"])
    # Not closed yet: O1 is not a cluster of an earlier round.
    with pytest.raises(ValueError, match="The parent of R1, O1, is not a cluster"):
        rounds.add(mtc.recursive_round(1), "R1", ["c1"], parent="O1")


def test_the_rounds_are_written_and_read_back(tmp_path):
    clusters_file, members_file = _rounds().write(tmp_path)
    assert clusters_file == str(tmp_path / mtc.CLUSTERS_FILE)
    assert members_file == str(tmp_path / mtc.CLUSTER_MEMBERS_FILE)

    lines = (tmp_path / mtc.CLUSTERS_FILE).read_text().splitlines()
    # The round, the name and the parent first, then the other columns in the
    # order they come; NA for a missing value.
    assert lines[0].split("\t") == [
        "round",
        "name",
        "parent_round",
        "parent",
        "Completeness",
        "Contigs",
    ]
    assert lines[1].split("\t") == ["overlapping", "O1", "NA", "NA", "0.9", "NA"]

    clusters, members = mtc.read_rounds(tmp_path)
    assert clusters["name"].tolist() == ["O1", "O2", "R1", "R2", "R1", "O2"]
    assert clusters["parent"].isna().tolist() == [True, True, False, False, False, False]
    assert members.columns.tolist() == ["round", "cluster", "contig"]
    assert len(members) == 10
    assert all(pd.api.types.is_string_dtype(members[c]) for c in members)


def test_names_that_read_as_numbers_stay_text(tmp_path):
    """Cluster and contig names are read as text: ``read_csv`` would make 1.0
    of a cluster called 1 in a column that also holds NA."""
    rounds = mtc.Rounds()
    rounds.add(mtc.OVERLAPPING, "1", ["10"])
    rounds.close_round(mtc.OVERLAPPING)
    rounds.add(mtc.FINAL, "1", ["10"], parent="1")
    rounds.write(tmp_path)
    clusters, members = mtc.read_rounds(tmp_path)
    assert clusters["parent"].tolist()[1] == "1"
    assert members["contig"].tolist() == ["10", "10"]


def test_a_run_without_rounds_reads_as_none(tmp_path):
    assert mtc.read_rounds(tmp_path) is None


def test_one_table_without_the_other_is_refused(tmp_path):
    _rounds().write(tmp_path)
    (tmp_path / mtc.CLUSTER_MEMBERS_FILE).unlink()
    with pytest.raises(ValueError, match=r"cluster_members\.txt is missing"):
        mtc.read_rounds(tmp_path)


def test_contig_names_that_read_as_missing_values_stay_names(tmp_path):
    """A contig called NA, None, null or nan is a contig: read_rounds keeps
    its name (pandas made them missing values, and the map's rounds lost
    them)."""
    names = ["NA", "None", "null", "nan"]
    rounds = mtc.Rounds()
    rounds.add(mtc.OVERLAPPING, "O1", names)
    rounds.close_round(mtc.OVERLAPPING)
    rounds.write(tmp_path)
    clusters, members = mtc.read_rounds(tmp_path)
    assert members["contig"].tolist() == names
    assert clusters["parent"].isna().all()


def test_the_tables_of_the_rounds_are_in_the_log_file(tmp_path):
    """The tables of the rounds are logged with the other lines of a run, in
    its log file: they were logged to a logger that no handler hears."""
    import logging

    import metator.log as mtl

    log = tmp_path / "run.log"
    mtl.set_file_handler(str(log))
    try:
        rounds = mtc.Rounds()
        rounds.add(mtc.OVERLAPPING, "metator_00001_00000", ["c1"])
        rounds.close_round(mtc.OVERLAPPING)
        rounds.write(str(tmp_path))
    finally:
        for handler in mtl.logger.handlers[:]:
            if isinstance(handler, logging.FileHandler):
                mtl.logger.removeHandler(handler)
                handler.close()
    assert "1 clusters in 1 rounds are written to" in log.read_text()
