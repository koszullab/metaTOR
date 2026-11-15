# Test for CLI toos of metator
# Commands are simply run to test crashes.


from pathlib import Path
import gzip
import metator.cli as mtc
import metator.io as mio
import metacooler
import pypairix
import pytest
import shutil
import os


# Use global variables for temporary tmp/out folders
@pytest.fixture(scope="session")
def tmp_dir(tmp_path_factory):
    p = tmp_path_factory.mktemp("tmp")
    os.makedirs(p, exist_ok=True)
    return str(p)


global_args = {
    "BAM_FOR": "tests_data/alignment_0_for.bam",
    "BAM_REV": "tests_data/alignment_0_rev.bam",
    "DEPTH": "tests_data/depth.txt",
    "FASTA": "tests_data/assembly.fa",
    "FASTQ_FOR": "tests_data/for_paired.fq.gz",
    "FASTQ_REV": "tests_data/rev_paired.fq.gz",
    "FASTA_INDEX": "tests_data/assembly",
    "PAIRS": "tests_data/alignment.pairs",
    "NETWORK": "tests_data/network.txt",
    "CONTIGS": "tests_data/contig_data_network.txt",
    "CONTIGS_PARTITION": "tests_data/contig_data_partition.txt",
    "OVERLAPPING_BINS": "tests_data/overlapping_bin/",
    "SUBSET_FASTA": "tests_data/subset.assembly.fa",
    "SUBSET_PAIRS": "tests_data/subset.pairs.gz",
    "SUBSET_MGES": "tests_data/subset.mges.txt",
}
NORMALIZE = (
    "norm",
    ["None", "length", "abundance", "RS", "empirical_hit", "theoritical_hit"],
)
ALGORITHM = ("alg", ["louvain", "leiden"])
ALIGNER = ("aligner", ["bwa", "bowtie2"])
MODE = ("mode", ["normal", "iterative", "cutsite"])


@pytest.mark.parametrize(*ALIGNER)
def test_network1(aligner, tmp_path, run_cli):
    args = "-1 {FASTQ_FOR} -2 {FASTQ_REV} -a {FASTA} -b {0} -o {OUT_TEST} -T {TMP}".format(
        aligner,
        TMP=Path(tmp_path, "tmp"),
        OUT_TEST=Path(tmp_path, "out_test"),
        **global_args,
    )
    run_cli("network", *args.split(" "))


def test_network2(tmp_path, run_cli):
    args = "-1 {BAM_FOR} -2 {BAM_REV} -a {FASTA_INDEX} -d {DEPTH} -E 500 -o {OUT_TEST} -T {TMP} -S bam -N".format(
        TMP=Path(tmp_path, "tmp"),
        OUT_TEST=Path(tmp_path, "out_test"),
        **global_args,
    )
    run_cli("network", *args.split(" "))


@pytest.mark.parametrize(*NORMALIZE)
def test_network3(norm, tmp_path, run_cli):
    args = "-1 {PAIRS} -a {FASTA} -d {DEPTH} -o {OUT_TEST} -T {TMP} -n {0} -e HindIII,DpnII -S pair".format(
        norm,
        TMP=Path(tmp_path, "tmp"),
        OUT_TEST=Path(tmp_path, "out_test"),
        **global_args,
    )
    run_cli("network", *args.split(" "))


@pytest.mark.parametrize(*MODE)
def test_network4(mode, tmp_path, run_cli):
    args = "-1 {FASTQ_FOR} -2 {FASTQ_REV} -a {FASTA} -o {OUT_TEST} -T {TMP} -B {0} -t 4 -e HindIII,DpnII".format(
        mode,
        TMP=Path(tmp_path, "tmp"),
        OUT_TEST=Path(tmp_path, "out_test"),
        **global_args,
    )
    run_cli("network", *args.split(" "))


@pytest.mark.parametrize(*ALGORITHM)
def test_partition(alg, tmp_path, run_cli):
    args = (
        "-a {SUBSET_FASTA} -c {CONTIGS} -i 5 -n {NETWORK} -o {OUT_TEST} -s 30000 -t 8 -T {TMP} -A {0} -FC --no-clean-up"
    ).format(
        alg,
        TMP=Path(tmp_path, "tmp"),
        OUT_TEST=Path(tmp_path, "out_test"),
        **global_args,
    )
    print("CLI command:\nmetator partition", *args.split(" "))
    run_cli("partition", *args.split(" "))


def test_validation(tmp_path, run_cli):
    args = (
        "-a {SUBSET_FASTA} -c {CONTIGS_PARTITION} -f {OVERLAPPING_BINS} -i 5 -n {NETWORK} -o {OUT_TEST} -t 8 -T {TMP} -F --no-clean-up"
    ).format(
        TMP=Path(tmp_path, "tmp"),
        OUT_TEST=Path(tmp_path, "out_test"),
        **global_args,
    )
    print("CLI command:\nmetator validation", *args.split(" "))
    run_cli("validation", *args.split(" "))


def test_network_partition_validation(tmp_dir, run_cli):
    tmp = Path(tmp_dir, "tmp")
    out = Path(tmp_dir, "out_test")
    args = "-1 {SUBSET_PAIRS} -a {SUBSET_FASTA} -o {OUT_TEST} -T {TMP} -e HindIII,DpnII -S pair --no-clean-up".format(
        TMP=tmp,
        OUT_TEST=out,
        **global_args,
    )
    print("CLI command:\nmetator network", *args.split(" "))
    run_cli("network", *args.split(" "))

    args = (
        "-a {SUBSET_FASTA} -c {OUT_TEST}/contig_data_network.txt -i 5 -n {OUT_TEST}/network.txt -o {OUT_TEST} -s 30000 -t 8 -T {TMP} -FC --no-clean-up"
    ).format(
        TMP=tmp,
        OUT_TEST=out,
        **global_args,
    )
    print("CLI command:\nmetator partition", *args.split(" "))
    run_cli("partition", *args.split(" "))

    args = (
        "-a {SUBSET_FASTA} -c {OUT_TEST}/contig_data_partition.txt -f {OUT_TEST}/overlapping_bin/ -i 5 -n {OUT_TEST}/network.txt -o {OUT_TEST} -t 8 -T {TMP} -F --no-clean-up"
    ).format(
        TMP=tmp,
        OUT_TEST=out,
        **global_args,
    )
    print("CLI command:\nmetator validation", *args.split(" "))
    run_cli("validation", *args.split(" "))


def test_pipeline(tmp_path, run_cli):
    """Test the metator pipeline command with subset data."""
    args = (
        "-1 {SUBSET_PAIRS} -a {SUBSET_FASTA} -F -o {OUT_TEST} -s 30000 -S pair -e HindIII,DpnII --no-clean-up -T {TMP}"
    ).format(
        TMP=Path(tmp_path, "tmp"),
        OUT_TEST=Path(tmp_path, "out_test"),
        **global_args,
    )
    print("CLI command:\n", "metator pipeline", *args.split(" "))
    run_cli("pipeline", *args.split(" "))
    # The contact map of the contigs, grouped by bin.
    metacool = Path(tmp_path, "out_test", "metator.metacool")
    metacooler.validate(str(metacool))
    genomes = metacooler.MetaCooler(str(metacool)).genomes
    assert set(genomes.name) - {"__unbinned__"} == set(mio.read_bin_summary(Path(tmp_path, "out_test", "bin_summary.txt")).index)


def test_pairs(tmp_path, run_cli, monkeypatch):
    """Plain and compressed pairs files are sorted, compressed and indexed
    next to their input, with temporary files in the working directory."""
    OUT_TEST = Path(tmp_path, "out_test")
    os.makedirs(OUT_TEST, exist_ok=True)
    pairs1 = OUT_TEST / "align1.pairs"
    pairs2 = OUT_TEST / "align2.pairs.gz"
    pairs = "{PAIRS}".format(**global_args)
    shutil.copyfile(pairs, pairs1)
    with open(pairs, "rb") as source, gzip.open(pairs2, "wb") as target:
        shutil.copyfileobj(source, target)
    sort_calls = []
    real_sort_pairs = mio.hio.sort_pairs

    def sort_pairs(*args, **kwargs):
        sort_calls.append(kwargs)
        return real_sort_pairs(*args, **kwargs)

    monkeypatch.setattr(mio.hio, "sort_pairs", sort_pairs)
    run_cli("pairs", "-t", "2", pairs1, pairs2)
    assert [(call["tmp_dir"], call["threads"]) for call in sort_calls] == [(".", 2), (".", 2)]
    outputs = [OUT_TEST / "align1_sorted.pairs.gz", OUT_TEST / "align2.pairs_sorted.pairs.gz"]
    with open(pairs) as handle:
        lines = handle.read().splitlines()
    body = sorted(line for line in lines if not line.startswith("#"))
    contig1, _, contig2 = body[0].split("\t")[1:4]
    for output in outputs:
        with gzip.open(output, "rt") as handle:
            sorted_lines = handle.read().splitlines()
        assert "#sorted: chr1-chr2-pos1-pos2" in sorted_lines
        assert sorted(line for line in sorted_lines if not line.startswith("#")) == body
        assert next(pypairix.open(str(output)).query2D(contig1, 0, 10**9, contig2, 0, 10**9, 1))
    # The outputs exist: they are kept unless --force; --remove deletes the
    # inputs.
    with pytest.raises(ValueError):
        run_cli("pairs", str(pairs1))
    run_cli("pairs", "-rF", pairs1, pairs2)
    assert not pairs1.exists() and not pairs2.exists()
    assert all(output.exists() and Path(f"{output}.px2").exists() for output in outputs)


def test_tmpdir_contents_survive(tmp_path, run_cli):
    """--tmpdir may already hold the user's files: a run removes only its own."""
    user_tmp = tmp_path / "scratch"
    user_tmp.mkdir()
    keep = user_tmp / "precious.txt"
    keep.write_text("keep me")
    args = [
        "-1", str(Path(global_args["SUBSET_PAIRS"]).resolve()),
        "-a", str(Path(global_args["SUBSET_FASTA"]).resolve()),
        "-o", str(tmp_path / "out"),
        "-T", str(user_tmp),
        "-S", "pair",
    ]  # fmt: skip
    run_cli("network", *args)
    assert keep.read_text() == "keep me"
    assert list(user_tmp.iterdir()) == [keep]


def test_default_tmpdir_is_tidied(tmp_path, monkeypatch, run_cli):
    """Without -T, ./tmp is removed if the run created it, kept if it existed."""
    args = [
        "-1", str(Path(global_args["SUBSET_PAIRS"]).resolve()),
        "-a", str(Path(global_args["SUBSET_FASTA"]).resolve()),
        "-o", "out",
        "-S", "pair",
    ]  # fmt: skip
    monkeypatch.chdir(tmp_path)
    run_cli("network", *args)
    assert not (tmp_path / "tmp").exists()

    (tmp_path / "tmp").mkdir()
    (tmp_path / "tmp" / "mine.txt").write_text("x")
    run_cli("network", *args)
    assert [p.name for p in (tmp_path / "tmp").iterdir()] == ["mine.txt"]


def test_network_default_outdir(tmp_path, monkeypatch, run_cli):
    """Without -o, the network is written to the current directory, with the
    log of the run."""
    args = [
        "-1", str(Path(global_args["SUBSET_PAIRS"]).resolve()),
        "-a", str(Path(global_args["SUBSET_FASTA"]).resolve()),
        "-S", "pair",
    ]  # fmt: skip
    monkeypatch.chdir(tmp_path)
    run_cli("network", *args)
    assert (tmp_path / "network.txt").stat().st_size > 0
    assert (tmp_path / "contig_data_network.txt").stat().st_size > 0
    (log,) = tmp_path.glob("metator_network_*.log")
    assert f"## Command used: metator network -1 {args[1]} -a {args[3]} -S pair\n" in log.read_text()


def test_invalid_prefix_is_rejected(tmp_path, run_cli):
    """miComplete drops dots from bin names, so a dotted --prefix is refused."""
    args = ("-a {SUBSET_FASTA} -c {CONTIGS} -n {NETWORK} -o {OUT} -T {TMP} -p my.sample").format(
        OUT=tmp_path / "out", TMP=tmp_path / "tmp", **global_args
    )
    with pytest.raises(ValueError, match="Invalid --prefix"):
        run_cli("partition", *args.split(" "))


def test_spinglass_is_rejected_before_any_work(tmp_path, monkeypatch, run_cli):
    """spinglass is gone: refuse it up front, not after miComplete has run."""

    def not_reached(*args, **kwargs):
        raise AssertionError("recursive decontamination should not start")

    monkeypatch.setattr(mtc.mtv, "recursive_decontamination", not_reached)
    args = ("-a {SUBSET_FASTA} -c {CONTIGS_PARTITION} -f {OVERLAPPING_BINS} -n {NETWORK} -o {OUT} -T {TMP} -A spinglass").format(
        OUT=tmp_path / "out", TMP=tmp_path / "tmp", **global_args
    )
    result = run_cli("validation", *args.split(" "), exit_code=2)
    assert "'spinglass' is not one of 'louvain', 'leiden'" in result.output
    assert not (tmp_path / "out").exists()
    assert not (tmp_path / "tmp").exists()


# Argument checks: each command refuses inconsistent options before working.

NETWORK_ERRORS = {
    "fastq without reverse reads": ["-1", "r1.fq"],
    "bowtie2 bam without reverse alignment": ["-1", "r1.bam", "-S", "bam"],
    "RS without enzyme": ["-1", "r1.fq", "-2", "r2.fq", "-n", "RS"],
    "theoritical_hit without enzyme": ["-1", "r1.fq", "-2", "r2.fq", "-n", "theoritical_hit", "-d", "depth.txt"],
    "abundance without depth": ["-1", "r1.fq", "-2", "r2.fq", "-n", "abundance"],
    "theoritical_hit without depth": ["-1", "r1.fq", "-2", "r2.fq", "-n", "theoritical_hit", "-e", "DpnII"],
    "cutsite without enzyme, then RS": ["-1", "r1.fq", "-2", "r2.fq", "-B", "cutsite", "-n", "RS"],
    "assembly neither fasta nor index": ["-1", "x.pairs", "-S", "pair", "-a", "missing.fa"],
}


@pytest.mark.parametrize("case", list(NETWORK_ERRORS))
def test_network_argument_errors(case, tmp_path, run_cli):
    args = ["-o", str(tmp_path / "out"), "-T", str(tmp_path / "tmp"), *NETWORK_ERRORS[case]]
    if "-a" not in args:
        args += ["-a", global_args["FASTA"]]
    with pytest.raises(ValueError):
        run_cli("network", *args)


# Unknown values are refused while the command line is read: nothing runs.
NETWORK_BAD_CHOICES = {
    "unknown aligner": ["-b", "minimap2"],
    "unknown aligner mode": ["-B", "fast"],
    "unknown normalization": ["-B", "cutsite", "-n", "cosine"],
    "unknown start": ["-S", "sam"],
    "network start": ["-S", "network"],
}


@pytest.mark.parametrize("case", list(NETWORK_BAD_CHOICES))
def test_network_bad_choices(case, tmp_path, run_cli):
    args = ["-1", "r1.fq", "-2", "r2.fq", "-a", global_args["FASTA"], "-o", str(tmp_path / "out"),
            "-T", str(tmp_path / "tmp"), *NETWORK_BAD_CHOICES[case]]  # fmt: skip
    result = run_cli("network", *args, exit_code=2)
    assert f"'{NETWORK_BAD_CHOICES[case][-1]}' is not one of" in result.output
    assert not (tmp_path / "out").exists()
    assert not (tmp_path / "tmp").exists()


PIPELINE_ERRORS = {
    "RS without enzyme": ["-m", "RS"],
    "abundance without depth": ["-m", "abundance"],
    "fastq without reverse reads": ["-1", "r1.fq", "-S", "fastq"],
    "assembly neither fasta nor index": ["-1", "x.pairs", "-S", "pair", "-a", "missing.fa"],
}


@pytest.mark.parametrize("case", list(PIPELINE_ERRORS))
def test_pipeline_argument_errors(case, tmp_path, run_cli):
    args = ["-o", str(tmp_path / "out"), "-T", str(tmp_path / "tmp"), *PIPELINE_ERRORS[case]]
    if "-a" not in args:
        args += ["-a", global_args["FASTA"]]
    with pytest.raises(ValueError):
        run_cli("pipeline", *args)


PIPELINE_BAD_CHOICES = {
    "unknown aligner": ["-b", "minimap2"],
    "unknown aligner mode": ["-B", "fast"],
    "unknown algorithm": ["-B", "cutsite", "-A", "spinglass"],
    "unknown normalization": ["-m", "cosine"],
    "unknown start": ["-S", "network"],
}


@pytest.mark.parametrize("case", list(PIPELINE_BAD_CHOICES))
def test_pipeline_bad_choices(case, tmp_path, run_cli):
    args = ["-a", global_args["FASTA"], "-o", str(tmp_path / "out"), "-T", str(tmp_path / "tmp"),
            *PIPELINE_BAD_CHOICES[case]]  # fmt: skip
    result = run_cli("pipeline", *args, exit_code=2)
    assert f"'{PIPELINE_BAD_CHOICES[case][-1]}' is not one of" in result.output
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("existing", ["overlapping_bin", "recursive_bin", "final_bin_unscaffold", "scaffold_bin"])
def test_pipeline_refuses_existing_outputs(existing, tmp_path, run_cli):
    (tmp_path / "out" / existing).mkdir(parents=True)
    args = ["-a", global_args["FASTA"], "-1", "x.pairs", "-S", "pair", "--scaffold", "-o", str(tmp_path / "out"),
            "-T", str(tmp_path / "tmp")]  # fmt: skip
    with pytest.raises(ValueError):
        run_cli("pipeline", *args)


def test_pipeline_force_replaces_outputs(tmp_path, run_cli):
    """With --force the output directories are emptied, then the run goes on
    (here up to the missing reverse reads)."""
    for existing in ["overlapping_bin", "recursive_bin", "final_bin_unscaffold", "scaffold_bin"]:
        (tmp_path / "out" / existing).mkdir(parents=True)
        (tmp_path / "out" / existing / "stale.fa").write_text(">old\nA\n")
    args = ["-a", global_args["FASTA"], "-1", "r1.fq", "--scaffold", "-F", "-o", str(tmp_path / "out"),
            "-T", str(tmp_path / "tmp")]  # fmt: skip
    with pytest.raises(ValueError):
        run_cli("pipeline", *args)
    assert not list((tmp_path / "out").glob("*/stale.fa"))


@pytest.mark.parametrize("existing", ["recursive_bin", "final_bin"])
def test_validation_refuses_existing_outputs(existing, tmp_path, run_cli):
    (tmp_path / "out" / existing).mkdir(parents=True)
    args = ("-a {SUBSET_FASTA} -c {CONTIGS_PARTITION} -f {OVERLAPPING_BINS} -n {NETWORK} -o {OUT} -T {TMP}").format(
        OUT=tmp_path / "out", TMP=tmp_path / "tmp", **global_args
    )
    with pytest.raises(ValueError):
        run_cli("validation", *args.split(" "))


def test_partition_existing_outputs(metagenome_inputs, metagenome_network, tmp_path, run_cli):
    """An existing overlapping_bin is refused, or replaced with --force; the
    --prefix names the bins."""
    fasta_dir = tmp_path / "out" / "overlapping_bin"
    fasta_dir.mkdir(parents=True)
    (fasta_dir / "stale.fa").write_text(">old\nA\n")
    args = ["-a", str(metagenome_inputs["assembly"]), "-c", str(metagenome_network["contig_data"]),
            "-n", str(metagenome_network["network"]), "-o", str(tmp_path / "out"), "-T", str(tmp_path / "tmp"),
            "-i", "5", "-s", "50000", "-p", "sample-1"]  # fmt: skip
    with pytest.raises(ValueError):
        run_cli("partition", *args)
    # An unknown algorithm is refused before anything is replaced.
    run_cli("partition", *args, "-F", "-A", "spinglass", exit_code=2)
    assert (fasta_dir / "stale.fa").exists()
    run_cli("partition", *args, "-F")
    bins = sorted(p.name for p in fasta_dir.iterdir())
    assert "stale.fa" not in bins
    assert bins and all(name.startswith("sample-1_metator_") for name in bins)


def test_check_output_path(tmp_path):
    target = tmp_path / "new" / "dir" / "file.txt"
    mtc.check_output_path(str(target))
    assert target.parent.is_dir()
    target.write_text("x")
    with pytest.raises(IOError, match="--force"):
        mtc.check_output_path(str(target))
    mtc.check_output_path(str(target), force=True)
