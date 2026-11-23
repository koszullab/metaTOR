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


def test_the_tables_of_the_rounds_are_those_written(tmp_path):
    """Rounds.tables, which a partition run from a map stores, are the tables
    that Rounds.write writes: the same rows and columns, the parents
    missing for none."""
    clusters, members = _rounds().tables()
    _rounds().write(tmp_path)
    assert clusters.to_csv(sep="\t", index=False, na_rep="NA") == (tmp_path / mtc.CLUSTERS_FILE).read_text()
    assert members.to_csv(sep="\t", index=False) == (tmp_path / mtc.CLUSTER_MEMBERS_FILE).read_text()
    assert clusters.parent.isna().sum() == 2


def test_the_overlapping_round_of_a_partition(tmp_path):
    """The round overlapping of a partition, from its contig data: the
    overlapping clusters that reach --size, named as the FASTA files that
    partition writes for them (<prefix>_<ID>_00000), in the order of their
    IDs (as numbers: an ID of six digits is the last), with their contigs in
    their order; a contig of no cluster ("-", no edge in the network) is in
    none, and a contig called NA is a contig."""
    table = tmp_path / "contig_data_partition.txt"
    table.write_text(
        "ID\tName\tSize\tOverlapping_bin_ID\tOverlapping_bin_contigs\tOverlapping_bin_size\n"
        "1\tc1\t600\t100000\t2\t900\n"
        "2\tNA\t300\t100000\t2\t900\n"
        "3\tc3\t500\t00002\t1\t500\n"
        "4\tc4\t499\t00003\t1\t499\n"
        "5\tc5\t50\t-\t-\t-\n"
        "6\tc6\t700\t00010\t1\t700\n"
    )
    clusters, members = mtc.overlapping_round(table, 500, "sample_metator")
    names = ["sample_metator_00002_00000", "sample_metator_00010_00000", "sample_metator_100000_00000"]
    assert list(clusters.name) == names
    assert list(clusters.columns) == ["round", "name", "parent_round", "parent"]
    assert set(clusters["round"]) == {mtc.OVERLAPPING} and clusters.parent.isna().all()
    assert list(zip(members.cluster, members.contig, strict=True)) == [
        (names[0], "c3"),
        (names[1], "c6"),
        (names[2], "c1"),
        (names[2], "NA"),
    ]
    assert set(members["round"]) == {mtc.OVERLAPPING}
    # No cluster reaches the size: the round has none.
    clusters, members = mtc.overlapping_round(table, 10_000, "metator")
    assert clusters.empty and members.empty


def test_the_mge_round_of_a_summary(tmp_path):
    """The round mge is the summary of metator mge, as written, in its
    order: one cluster per row, named BinName, without a parent, with the
    contigs of Contigs; its other columns but BinLength and ContigsNumber,
    which the map computes from the contigs. NA is a missing value, None a
    text (the host that mge leaves unknown); a contig called NA or nan is
    a contig."""
    summary = tmp_path / "mges_bin_summary.tsv"
    summary.write_text(
        "BinName\tBinLength\tGC\tHit\tBinningScore\tContigsNumber\tContigs\tMetagenomicBin\tAssociationScore\tAssociatedBins\n"
        "MetaTOR_MGE_00001\t12000.00\t45.12\t100.00\t0.91\t2.00\tnan,c7\tNone\tNA\tNA\n"
        "MetaTOR_MGE_00002\t6000.00\t55.00\t50.00\tNA\t1.00\tNA\tNone\tNA\tNA\n"
    )
    clusters, members = mtc.mge_round(str(summary))
    assert list(clusters.columns) == ["round", "name", "parent_round", "parent", "GC", "Hit", "BinningScore",
                                      "MetagenomicBin", "AssociationScore", "AssociatedBins"]  # fmt: skip
    assert list(clusters["round"]) == [mtc.MGE] * 2 and list(clusters["name"]) == ["MetaTOR_MGE_00001", "MetaTOR_MGE_00002"]
    assert clusters["parent_round"].isna().all() and clusters["parent"].isna().all()
    assert list(clusters["GC"]) == [45.12, 55.0] and list(clusters["Hit"]) == [100, 50]
    assert clusters["BinningScore"].iloc[0] == 0.91 and pd.isna(clusters["BinningScore"].iloc[1])
    assert list(clusters["MetagenomicBin"]) == ["None", "None"] and clusters["AssociatedBins"].isna().all()
    assert clusters["AssociationScore"].isna().all()
    assert members.values.tolist() == [[mtc.MGE, "MetaTOR_MGE_00001", "nan"], [mtc.MGE, "MetaTOR_MGE_00001", "c7"],
                                       [mtc.MGE, "MetaTOR_MGE_00002", "NA"]]  # fmt: skip


def test_the_mge_round_of_an_empty_summary(tmp_path):
    """No MGE cluster: a round with none."""
    summary = tmp_path / "mges_bin_summary.tsv"
    summary.write_text(
        "BinName\tBinLength\tGC\tHit\tBinningScore\tContigsNumber\tContigs\tMetagenomicBin\tAssociationScore\tAssociatedBins\n"
    )
    clusters, members = mtc.mge_round(str(summary))
    assert clusters.empty and members.empty and list(members.columns) == ["round", "cluster", "contig"]
