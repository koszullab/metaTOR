"""The clustering rounds of a run (metator.clusters)."""

import numpy as np
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
    """Cluster and contig names are read as text: ``read_csv`` would make 1.0 of a
    cluster called 1 in a column that also holds NA."""
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
    """A contig called NA, None, null or nan is a contig: read_rounds keeps its
    name (pandas made them missing values, and the map's rounds lost them)."""
    names = ["NA", "None", "null", "nan"]
    rounds = mtc.Rounds()
    rounds.add(mtc.OVERLAPPING, "O1", names)
    rounds.close_round(mtc.OVERLAPPING)
    rounds.write(tmp_path)
    clusters, members = mtc.read_rounds(tmp_path)
    assert members["contig"].tolist() == names
    assert clusters["parent"].isna().all()


def test_the_tables_of_the_rounds_are_in_the_log_file(tmp_path, debug_log):
    """The tables of the rounds are logged with the other lines of a run, in its
    log file, as a detail (DEBUG), which a log at DEBUG shows: they were
    logged to a logger that no handler hears."""
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
    assert " :: DEBUG :: 1 clusters in 1 rounds are written to" in log.read_text()


def test_the_tables_of_the_rounds_are_those_written(tmp_path):
    """Rounds.tables, which a partition run from a map stores, are the tables that
    Rounds.write writes: the same rows and columns, the parents missing for
    none."""
    clusters, members = _rounds().tables()
    _rounds().write(tmp_path)
    assert clusters.to_csv(sep="\t", index=False, na_rep="NA") == (tmp_path / mtc.CLUSTERS_FILE).read_text()
    assert members.to_csv(sep="\t", index=False) == (tmp_path / mtc.CLUSTER_MEMBERS_FILE).read_text()
    assert clusters.parent.isna().sum() == 2


def test_the_names_of_the_clusters_of_a_partition():
    """The core clusters are core_<number>, the overlapping ones
    <prefix>_<number>_00000, a prefix of underscores too: their numbers, as the
    contig data tables give them (00001...), come back from their names;
    another name gives none."""
    assert mtc.core_name(1) == "core_00001" and mtc.core_name(123456) == "core_123456"
    assert mtc.overlapping_name("my_sample_metator", 2) == "my_sample_metator_00002_00000"
    assert mtc.cluster_number(mtc.CORE, "core_00001") == 1
    assert mtc.cluster_number(mtc.CORE, "core_123456") == 123456
    assert mtc.cluster_number(mtc.OVERLAPPING, "my_sample_metator_00002_00000") == 2
    for round_name, name in [
        (mtc.CORE, "core_1"),
        (mtc.CORE, "core_0000x"),
        (mtc.CORE, "x_00001"),
        (mtc.OVERLAPPING, "metator_00002_00001"),
        (mtc.OVERLAPPING, "_00002_00000"),
        (mtc.OVERLAPPING, "metator_00002"),
        # Digits other than ASCII ones (full-width digit two).
        (mtc.OVERLAPPING, "metator_" + "\uff12" * 5 + "_00000"),
    ]:
        assert mtc.cluster_number(round_name, name) is None, name


def test_the_rounds_of_a_partition():
    """The rounds core and overlapping of a partition, as partition MAP stores
    them: every cluster, in the order of its number, with its contigs (by name:
    a contig called NA is a contig), no parent and no column; the overlapping
    clusters that miComplete checked get its columns when the validation
    records them."""
    partition = mtc.Partition(
        core={2: [3], 1: [1, 2]},
        overlapping={1: [1, 2, 3]},
        names={1: "sample_metator_00001_00000"},
    )
    names = ["c1", "NA", "c3"]
    clusters, members = mtc.partition_rounds(partition, names)
    assert list(zip(clusters["round"], clusters.name, strict=True)) == [
        (mtc.CORE, "core_00001"),
        (mtc.CORE, "core_00002"),
        (mtc.OVERLAPPING, "sample_metator_00001_00000"),
    ]
    assert list(clusters.columns) == ["round", "name", "parent_round", "parent"] and clusters.parent.isna().all()
    assert list(zip(members["round"], members.cluster, members.contig, strict=True)) == [
        (mtc.CORE, "core_00001", "c1"),
        (mtc.CORE, "core_00001", "NA"),
        (mtc.CORE, "core_00002", "c3"),
        (mtc.OVERLAPPING, "sample_metator_00001_00000", "c1"),
        (mtc.OVERLAPPING, "sample_metator_00001_00000", "NA"),
        (mtc.OVERLAPPING, "sample_metator_00001_00000", "c3"),
    ]
    assert partition.sizes([100, 50, 25]) == {1: 175}
    # The validation records them with the columns of miComplete of the
    # clusters that it checked.
    rounds = mtc.record_partition(mtc.Rounds(), partition, names, quality={"sample_metator_00001_00000": {"Completeness": 0.9}})
    clusters, _ = rounds.tables()
    assert list(clusters.Completeness.isna()) == [True, True, False]
    # A partition of no cluster (a network without an edge): empty rounds.
    clusters, members = mtc.partition_rounds(mtc.Partition(core={}, overlapping={}, names={}), names)
    assert clusters.empty and members.empty


def _mge_summary(rows):
    """A summary of MGE clusters, as metator.mge.generate_bin_summary gives it."""
    columns = ["BinName", "BinLength", "GC", "Hit", "BinningScore", "ContigsNumber", "Contigs"]
    return pd.DataFrame(rows, columns=columns)


@pytest.mark.parametrize("name", [mtc.MGE, mtc.MGE_RANDOM])
def test_the_mge_round_of_a_summary(name):
    """The round mge (or mge_random) is the summary of the MGE clusters of metator
    mge, in its order: one cluster per row, named BinName, without a parent,
    with the contigs of Contigs; its other columns but BinLength and
    ContigsNumber, which the map computes from the contigs."""
    summary = _mge_summary(
        [
            ["MetaTOR_MGE_00001", 12000, 45.1234, 100, 0.9123, 2, "nan,c7"],
            ["MetaTOR_MGE_00002", 6000, 55.0, 50, np.nan, 1, "NA"],
        ]
    )
    clusters, members = mtc.mge_round(summary, name)
    assert list(clusters.columns) == ["round", "name", "parent_round", "parent", "GC", "Hit", "BinningScore"]
    assert list(clusters["round"]) == [name] * 2 and list(clusters["name"]) == ["MetaTOR_MGE_00001", "MetaTOR_MGE_00002"]
    assert clusters["parent_round"].isna().all() and clusters["parent"].isna().all()
    # The values as computed: the map stores them unrounded.
    assert list(clusters["GC"]) == [45.1234, 55.0] and list(clusters["Hit"]) == [100, 50]
    assert clusters["BinningScore"].iloc[0] == 0.9123 and pd.isna(clusters["BinningScore"].iloc[1])
    assert members.values.tolist() == [[name, "MetaTOR_MGE_00001", "nan"], [name, "MetaTOR_MGE_00001", "c7"],
                                       [name, "MetaTOR_MGE_00002", "NA"]]  # fmt: skip


def test_the_mge_round_of_an_empty_summary():
    """No MGE cluster: a round with none."""
    clusters, members = mtc.mge_round(_mge_summary([]))
    assert clusters.empty and members.empty and list(members.columns) == ["round", "cluster", "contig"]


# Decimal numbers whose last digit pandas 2's default parser reads off by one
# (24.621666666666666 was the HiC_abundance of a final cluster of the suite's
# run), and others of 17 significant digits, as the tables of a run write them.
_LONG_NUMBERS = [24.621666666666666, *np.random.default_rng(0).uniform(0, 100, 200).tolist()]


def test_the_numbers_of_the_rounds_are_read_back_exactly(tmp_path):
    """Each number of clusters.txt reads back as the number written: the map
    stores what read_rounds reads, and metator dump writes it again."""
    rounds = mtc.Rounds()
    for i, value in enumerate(_LONG_NUMBERS):
        rounds.add(mtc.FINAL, f"metator_{i + 1:05d}_00000", [f"c{i}"], columns={"HiC_abundance": value})
    rounds.write(tmp_path)
    clusters, _ = mtc.read_rounds(tmp_path)
    assert clusters["HiC_abundance"].tolist() == _LONG_NUMBERS


def test_the_numbers_of_an_mge_summary_are_kept_exactly():
    """The round mge, which metator mge stores in the map, has the numbers of its
    summary as computed: unrounded, where mges_cluster_summary.tsv rounded them
    to two decimals."""
    summary = _mge_summary([[f"MetaTOR_MGE_{i:05d}", 100, value, 5, value, 1, f"c{i}"] for i, value in enumerate(_LONG_NUMBERS)])
    clusters, _ = mtc.mge_round(summary)
    assert clusters["GC"].tolist() == _LONG_NUMBERS and clusters["BinningScore"].tolist() == _LONG_NUMBERS
