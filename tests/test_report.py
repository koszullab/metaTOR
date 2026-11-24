"""The report of a run (metator.report, metator report): what the map of its
clusters holds, digested, as report.txt, report_mqc.yaml and report.html, on a
pipeline --mges run of the mock community (stub miComplete)."""

import json
import os
import random
import re
import shutil
import sys

import metacooler
import numpy as np
import pytest
import yaml
from click.testing import CliRunner

import metator.dump as mtd
import metator.metacool as mtmc
import metator.network as mtn
import metator.outputs as mto
import metator.report as mtr
from metator.cli import cli
from tests import synthetic
from tests.maps import contig_map


@pytest.fixture(scope="module")
def run(metagenome, metagenome_files, tmp_path_factory):
    """A pipeline --mges run on the mock community, with a stub miComplete that
    gives a cluster a weighted completeness of 0.8951 (an MQ MAG, 0.90 once
    rounded): its output directory, its map and its log. Read only."""
    directory = tmp_path_factory.mktemp("report")

    def micomplete(fasta_dir, outfile, threads):
        synthetic.write_micomplete_summary(fasta_dir, outfile, metagenome)
        with open(outfile) as summary:
            lines = [line.split("\t") for line in summary.read().splitlines()]
        lines[4][lines[3].index("Weighted completeness")] = "0.8951"
        with open(outfile, "w") as summary:
            summary.write("".join("\t".join(line) + "\n" for line in lines))

    contig_data = mtn.create_contig_data(str(metagenome_files["assembly"]), str(metagenome_files["depth"]), "DpnII")
    metacool = contig_map(metagenome_files["pairs"], contig_data, directory / "sample.metacool", assembly_name="synthetic")
    truth = metagenome.truth
    mges = directory / "mges.txt"
    mges.write_text("".join(f"{name}\n" for name in truth.Name[truth.Kind == "mge"]))
    out = directory / "out"
    args = ["pipeline", str(metacool), "-a", str(metagenome_files["assembly"]), "-o", str(out), "-T", str(directory / "tmp"),
            "-s", "20000", "-i", "10", "-j", "5", "--mges", str(mges), "--seed", "3"]  # fmt: skip
    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr("metator.validation.micomplete_quality", micomplete)
        monkeypatch.setattr(sys, "argv", ["metator", *args])
        random.seed(0)
        result = CliRunner().invoke(cli, args, catch_exceptions=False)
    assert result.exit_code == 0, result.output
    (log,) = out.glob("metator_pipeline_*.log")
    return {"out": out, "map": out / metacool.name, "log": log.read_text(), "assembly": metagenome_files["assembly"]}


def _read(run):
    """The map of the run, as metator report reads it."""
    return mtd.read_run(str(run["map"]))


def test_the_report_is_the_run(run):
    """The numbers of the report are those of the run: its contigs, its rounds,
    its final clusters by quality class, as its log counts them (from the
    values of miComplete: the cluster of 0.8951 is an MQ MAG), its MAGs, its
    MGE clusters and their hosts, as host counted them."""
    mc, layout, stored = _read(run)
    report = mtr.collect(str(run["map"]), mc, layout, stored)
    sections, tables = report.sections, report.tables
    assert sections["report"]["sample"] == "sample" and sections["report"]["format_version"] == mtr.FORMAT_VERSION
    assert sections["run"]["steps"] == "network,partition,validation,mge,host"
    assert sections["run"]["size_bp"] == 20000 and sections["run"]["normalization"] == "empirical_hit"
    assert sections["run"]["seed"] == sections["run"]["rec_seed"] == 3
    lengths = mc.contigs["length"].to_numpy()
    assert sections["contigs"]["contigs"] == len(lengths) and sections["contigs"]["contigs_bp"] == lengths.sum()
    ordered = np.sort(lengths)[::-1]
    assert sections["contigs"]["contig_n50_bp"] == ordered[np.argmax(np.cumsum(ordered) >= ordered.sum() / 2)]
    # The quality classes, as the log counts them.
    counted = dict(
        re.findall(r"INFO :: (HQ MAGs|MQ MAGs|LQ MAGs|Contaminated potential MAGs|Other clusters): ([0-9]+)\t", run["log"])
    )
    keys = {"HQ MAGs": "hq_mags", "MQ MAGs": "mq_mags", "LQ MAGs": "lq_mags", "Contaminated potential MAGs": "contaminated_mags",
            "Other clusters": "other_clusters"}  # fmt: skip
    assert {keys[label]: int(count) for label, count in counted.items()} == {
        key: sections["clusters"][key] for key in keys.values()
    }
    assert sections["clusters"]["mq_mags"] == 1 and sections["clusters"]["final_clusters"] == len(layout.rows)
    # The rounds, the MAGs, the MGE clusters and their hosts.
    assert tables["rounds"]["round"].tolist() == [entry["name"] for entry in mc.cluster_rounds]
    assert sorted(tables["mags"]["name"]) == sorted(mto.map_mags(mc))
    assert (tables["mags"]["quality"] == "MQ").sum() == 1
    (hosts,) = re.findall(r"([0-9]+) MGE clusters and ([0-9]+) MAGs: ([0-9]+) MGE clusters with one host", run["log"])[-1:]
    assert sections["mges"]["mge_clusters"] == int(hosts[0]) == len(tables["mge_clusters"])
    assert sections["hosts"]["mge_clusters_with_one_host"] == int(hosts[2])
    assert len(tables["hosts"]) == sum(len(hosted) for hosted in mto.hosted_mges(mc).values())


def test_report_txt_reads_back(run):
    """report.txt, as the run wrote it, reads back (see metator.report.read_txt):
    its sections, key by key, and its tables, row by row, every field one
    token, as MultiQC or a script splits them."""
    sections, tables = mtr.read_txt(run["out"] / mtr.REPORT_TXT)
    assert list(sections) == ["report", "run", "contigs", "contacts", "clusters", "mges", "hosts", "commands"]
    assert sections["report"]["format_version"] == str(mtr.FORMAT_VERSION) and sections["report"]["sample"] == "sample"
    assert sections["clusters"]["mq_mags"] == "1"
    assert sections["commands"]["validation"].startswith("metator pipeline ")
    assert list(tables) == ["rounds", "mags", "mge_clusters", "hosts"]
    assert list(tables["mags"].columns) == list(mtr.MAG_COLUMNS)
    mc, layout, stored = _read(run)
    report = mtr.collect(str(run["map"]), mc, layout, stored)
    for name, table in report.tables.items():
        assert len(tables[name]) == len(table), name
    text = (run["out"] / mtr.REPORT_TXT).read_text()
    assert text.startswith("# metaTOR report of sample\n")


def test_report_mqc_yaml(run):
    """report_mqc.yaml is MultiQC's custom content: a General Statistics row for
    the run, named after its map, which multiqc aggregates over runs (the same
    id for every run)."""
    content = yaml.safe_load((run["out"] / mtr.REPORT_MQC).read_text())
    assert content["id"] == "metator" and content["plot_type"] == "generalstats"
    (row,) = content["data"].values()
    assert list(content["data"]) == ["sample"]
    assert [next(iter(column)) for column in content["pconfig"]] == list(row)
    assert row["mq_mags"] == 1 and row["mge_clusters"] == 2 and 0 < row["trans_contig_ratio"] < 1


def test_report_html(run):
    """report.html, from its template: the page of the run, its plots as Plotly
    figures in JSON (without NaN, which JSON.parse refuses), drawn by
    plotly.js from its CDN, a row per MAG in its table, and links to the
    outputs of the run that are beside it."""
    page = (run["out"] / mtr.REPORT_HTML).read_text()
    assert f'<script src="{mtr.PLOTLY_JS}"' in page and "{{" not in page and "{%" not in page
    (figures,) = re.findall(r'<script id="metator-figures" type="application/json">(.*?)</script>', page, re.S)
    figures = json.loads(figures, parse_constant=lambda name: pytest.fail(f"{name} in the figures"))
    assert set(figures) == {"plot-quality", "plot-classes", "plot-rounds", "plot-lengths", "plot-contigs", "plot-contacts",
                            "plot-heatmap", "plot-hosts"}  # fmt: skip
    assert {trace["type"] for figure in figures.values() for trace in figure["data"]} == {
        "scatter", "bar", "scattergl", "heatmap", "sankey"}  # fmt: skip
    # The lengths of the contigs, binned: a count per bin, of every contig.
    lengths = figures["plot-lengths"]["data"]
    assert sum(count or 0 for trace in lengths for count in trace["y"]) == len(_read(run)[0].contigs)
    for element in re.findall(r'<div class="plot" id="([a-z-]+)" data-figure>', page):
        assert element in figures
    mc = _read(run)[0]
    for mag in mto.map_mags(mc):
        assert f"<td>{mag}</td>" in page and f'href="{os.path.join(mto.MAGS_DIR, mag)}.fa"' in page
    for link in re.findall(r'href="([^"#][^"]*)"', page):
        assert link.startswith("https://") or (run["out"] / link).exists(), link


def test_the_report_of_a_map_alone(run, tmp_path, run_cli):
    """metator report MAP -o DIR writes the report in DIR from the map alone,
    and the log; again, refused without -F, the report left as it was, and
    replaced with it. Without the outputs of the run beside it, its page
    links to none of them."""
    out = tmp_path / "report"
    run_cli("report", run["map"], "-o", out)
    assert sorted(path.name for path in out.iterdir() if not path.name.endswith(".log")) == sorted(mtr.REPORT_FILES)
    assert "fasta" not in (out / mtr.REPORT_HTML).read_text().split('<section id="files">')[1]
    before = (out / mtr.REPORT_TXT).read_text()
    with pytest.raises(ValueError):
        run_cli("report", run["map"], "-o", out)
    assert (out / mtr.REPORT_TXT).read_text() == before
    run_cli("report", run["map"], "-o", out, "-F")
    assert mtr.read_txt(out / mtr.REPORT_TXT)[0]["clusters"] == mtr.read_txt(run["out"] / mtr.REPORT_TXT)[0]["clusters"]


@pytest.mark.parametrize("case", ["no clusters", "a directory as report.html"])
def test_report_refusals(case, run, metagenome_files, tmp_path, run_cli):
    """What metator report refuses, before it writes anything: a map without the
    clusters of a validation, as metator dump refuses it; a directory where a
    file of the report goes, with or without -F."""
    out = tmp_path / "report"
    metacool = run["map"]
    if case == "no clusters":
        contig_data = mtn.create_contig_data(str(metagenome_files["assembly"]), str(metagenome_files["depth"]), "DpnII")
        metacool = contig_map(metagenome_files["pairs"], contig_data, tmp_path / "contigs.metacool")
        message = f"{metacool} holds no clusters of a validation: it lacks"
    else:
        (out / mtr.REPORT_HTML).mkdir(parents=True)
        message = f"{out / mtr.REPORT_HTML} is a directory, where the run writes its report"
    with pytest.raises(ValueError):
        run_cli("report", metacool, "-o", out, "-F")
    (log,) = out.glob("metator_report_*.log")
    assert f"ERROR :: {message}" in log.read_text()
    assert not (out / mtr.REPORT_TXT).exists()


def test_export_writes_the_report(run, tmp_path, run_cli, logged):
    """metator export writes the report with the other outputs, last, and links
    them; the report of an earlier run there is refused without -F."""
    out = tmp_path / "exported"
    out.mkdir()
    (out / mtr.REPORT_TXT).write_text("an earlier report\n")
    with pytest.raises(ValueError):
        run_cli("export", run["map"], "-a", run["assembly"], "-o", out, "-T", tmp_path / "tmp")
    assert (out / mtr.REPORT_TXT).read_text() == "an earlier report\n"
    assert not (out / "cool").exists()
    run_cli("export", run["map"], "-a", run["assembly"], "-o", out, "-T", tmp_path / "tmp", "-F")
    assert all((out / name).exists() for name in mtr.REPORT_FILES)
    assert mtr.read_txt(out / mtr.REPORT_TXT)[1]["mags"].equals(mtr.read_txt(run["out"] / mtr.REPORT_TXT)[1]["mags"])
    assert 'href="cool/' in (out / mtr.REPORT_HTML).read_text()


def test_sample_name():
    """The name of a run in its report: its map without .metacool, or the entry
    of a collection."""
    assert mtr.sample_name("binning/sample.metacool") == "sample"
    assert mtr.sample_name("runs.mmetacool::entries/sample_2") == "sample_2"
    assert mtr.sample_name("map") == "map"


def test_a_report_without_mges(run, tmp_path):
    """A map without MGE clusters (that of the validation of the run, before its
    steps mge and host): its report has no section mges nor hosts, and its
    page none of their plots."""
    metacool = tmp_path / "clusters.metacool"
    shutil.copyfile(run["map"], metacool)
    # The map as the validation left it: the round mge, the MGE genomes and
    # the hosts dropped, as validation drops them when it runs again.
    import metator.clusters as mtcl

    mc = metacooler.MetaCooler(str(metacool))
    keep = [entry for entry in mc.cluster_rounds if entry["name"] != mtcl.MGE]
    clusters, members = mtmc.map_rounds(mc)
    metacooler.write_clusters(str(metacool), keep, clusters[clusters["round"] != mtcl.MGE], members[members["round"] != mtcl.MGE])
    mc, layout, stored = mtd.read_run(str(metacool))
    report = mtr.collect(str(metacool), mc, layout, stored)
    assert "mges" not in report.sections and "hosts" not in report.sections
    assert "mge_clusters" not in report.tables and "plot-hosts" not in mtr.figures(report)
    mtr.write_report(report, str(tmp_path / "out"))
    page = (tmp_path / "out" / mtr.REPORT_HTML).read_text()
    assert 'id="mges"' not in page and "plot-hosts" not in page
