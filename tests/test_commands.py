# Test for CLI toos of metator
# Commands are simply run to test crashes.


from pathlib import Path
import metator.commands as mtc
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
def test_network1(aligner, tmp_path):
    args = "-1 {FASTQ_FOR} -2 {FASTQ_REV} -a {FASTA} -b {0} -o {OUT_TEST} -T {TMP}".format(
        aligner,
        TMP=Path(tmp_path, "tmp"),
        OUT_TEST=Path(tmp_path, "out_test"),
        **global_args,
    )
    proc = mtc.Network(args.split(" "), {})
    proc.execute()


def test_network2(tmp_path):
    args = "-1 {BAM_FOR} -2 {BAM_REV} -a {FASTA_INDEX} -d {DEPTH} -E 500 -o {OUT_TEST} -T {TMP} -S bam -N".format(
        TMP=Path(tmp_path, "tmp"),
        OUT_TEST=Path(tmp_path, "out_test"),
        **global_args,
    )
    proc = mtc.Network(args.split(" "), {})
    proc.execute()


@pytest.mark.parametrize(*NORMALIZE)
def test_network3(norm, tmp_path):
    args = "-1 {PAIRS} -a {FASTA} -d {DEPTH} -o {OUT_TEST} -T {TMP} -n {0} -e HindIII,DpnII -S pair".format(
        norm,
        TMP=Path(tmp_path, "tmp"),
        OUT_TEST=Path(tmp_path, "out_test"),
        **global_args,
    )
    proc = mtc.Network(args.split(" "), {})
    proc.execute()


@pytest.mark.parametrize(*MODE)
def test_network4(mode, tmp_path):
    args = "-1 {FASTQ_FOR} -2 {FASTQ_REV} -a {FASTA} -o {OUT_TEST} -T {TMP} -B {0} -t 4 -e HindIII,DpnII".format(
        mode,
        TMP=Path(tmp_path, "tmp"),
        OUT_TEST=Path(tmp_path, "out_test"),
        **global_args,
    )
    proc = mtc.Network(args.split(" "), {})
    proc.execute()


@pytest.mark.parametrize(*ALGORITHM)
def test_partition(alg, tmp_path):
    args = (
        "-a {SUBSET_FASTA} -c {CONTIGS} -i 5 -n {NETWORK} -o {OUT_TEST} -s 30000 -t 8 -T {TMP} -A {0} -FC --no-clean-up"
    ).format(
        alg,
        TMP=Path(tmp_path, "tmp"),
        OUT_TEST=Path(tmp_path, "out_test"),
        **global_args,
    )
    print("CLI command:\nmetator partition", *args.split(" "))

    proc = mtc.Partition(args.split(" "), {})
    proc.execute()


def test_validation(tmp_path):
    args = (
        "-a {SUBSET_FASTA} -c {CONTIGS_PARTITION} -f {OVERLAPPING_BINS} -i 5 -n {NETWORK} -o {OUT_TEST} -t 8 -T {TMP} -F --no-clean-up"
    ).format(
        TMP=Path(tmp_path, "tmp"),
        OUT_TEST=Path(tmp_path, "out_test"),
        **global_args,
    )
    print("CLI command:\nmetator validation", *args.split(" "))

    proc = mtc.Validation(args.split(" "), {})
    proc.execute()


def test_network_partition_validation(tmp_dir):
    tmp = Path(tmp_dir, "tmp")
    out = Path(tmp_dir, "out_test")
    args = "-1 {SUBSET_PAIRS} -a {SUBSET_FASTA} -o {OUT_TEST} -T {TMP} -e HindIII,DpnII -S pair --no-clean-up".format(
        TMP=tmp,
        OUT_TEST=out,
        **global_args,
    )
    print("CLI command:\nmetator network", *args.split(" "))
    proc = mtc.Network(args.split(" "), {})
    proc.execute()

    args = (
        "-a {SUBSET_FASTA} -c {OUT_TEST}/contig_data_network.txt -i 5 -n {OUT_TEST}/network.txt -o {OUT_TEST} -s 30000 -t 8 -T {TMP} -FC --no-clean-up"
    ).format(
        TMP=tmp,
        OUT_TEST=out,
        **global_args,
    )
    print("CLI command:\nmetator partition", *args.split(" "))
    proc = mtc.Partition(args.split(" "), {})
    proc.execute()

    args = (
        "-a {SUBSET_FASTA} -c {OUT_TEST}/contig_data_partition.txt -f {OUT_TEST}/overlapping_bin/ -i 5 -n {OUT_TEST}/network.txt -o {OUT_TEST} -t 8 -T {TMP} -F --no-clean-up"
    ).format(
        TMP=tmp,
        OUT_TEST=out,
        **global_args,
    )
    print("CLI command:\nmetator validation", *args.split(" "))
    proc = mtc.Validation(args.split(" "), {})
    proc.execute()


def test_pipeline(tmp_path):
    """Test the metator pipeline command with subset data."""
    args = (
        "-1 {SUBSET_PAIRS} -a {SUBSET_FASTA} -F -o {OUT_TEST} -s 30000 -S pair -e HindIII,DpnII --no-clean-up -T {TMP}"
    ).format(
        TMP=Path(tmp_path, "tmp"),
        OUT_TEST=Path(tmp_path, "out_test"),
        **global_args,
    )
    print("CLI command:\n", "metator pipeline", *args.split(" "))

    proc = mtc.Pipeline(args.split(" "), {})
    proc.execute()


def test_pairs(tmp_path):
    OUT_TEST = Path(tmp_path, "out_test")
    os.makedirs(OUT_TEST, exist_ok=True)
    pairs1 = f"{OUT_TEST}/align1.pairs"
    pairs2 = f"{OUT_TEST}/align2.pairs"
    pairs = "{PAIRS}".format(**global_args)
    shutil.copyfile(pairs, pairs1)
    shutil.copyfile(pairs, pairs2)
    args = f"-rF {pairs1} {pairs2}"
    proc = mtc.Pairs(args.split(" "), {})
    proc.execute()


def test_tmpdir_contents_survive(tmp_path):
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
    mtc.Network(args, {}).execute()
    assert keep.read_text() == "keep me"
    assert list(user_tmp.iterdir()) == [keep]


def test_default_tmpdir_is_tidied(tmp_path, monkeypatch):
    """Without -T, ./tmp is removed if the run created it, kept if it existed."""
    args = [
        "-1", str(Path(global_args["SUBSET_PAIRS"]).resolve()),
        "-a", str(Path(global_args["SUBSET_FASTA"]).resolve()),
        "-o", "out",
        "-S", "pair",
    ]  # fmt: skip
    monkeypatch.chdir(tmp_path)
    mtc.Network(args, {}).execute()
    assert not (tmp_path / "tmp").exists()

    (tmp_path / "tmp").mkdir()
    (tmp_path / "tmp" / "mine.txt").write_text("x")
    mtc.Network(args, {}).execute()
    assert [p.name for p in (tmp_path / "tmp").iterdir()] == ["mine.txt"]


def test_invalid_prefix_is_rejected(tmp_path):
    """miComplete drops dots from bin names, so a dotted --prefix is refused."""
    args = ("-a {SUBSET_FASTA} -c {CONTIGS} -n {NETWORK} -o {OUT} -T {TMP} -p my.sample").format(
        OUT=tmp_path / "out", TMP=tmp_path / "tmp", **global_args
    )
    with pytest.raises(ValueError, match="Invalid --prefix"):
        mtc.Partition(args.split(" "), {}).execute()


def test_spinglass_is_rejected_before_any_work(tmp_path, monkeypatch):
    """spinglass is gone: refuse it up front, not after miComplete has run."""

    def not_reached(*args, **kwargs):
        raise AssertionError("recursive decontamination should not start")

    monkeypatch.setattr(mtc.mtv, "recursive_decontamination", not_reached)
    args = ("-a {SUBSET_FASTA} -c {CONTIGS_PARTITION} -f {OVERLAPPING_BINS} -n {NETWORK} -o {OUT} -T {TMP} -A spinglass").format(
        OUT=tmp_path / "out", TMP=tmp_path / "tmp", **global_args
    )
    with pytest.raises(ValueError):
        mtc.Validation(args.split(" "), {}).execute()


# Argument checks: each command refuses inconsistent options before working.

NETWORK_ERRORS = {
    "fastq without reverse reads": ["-1", "r1.fq"],
    "unknown aligner": ["-1", "r1.fq", "-2", "r2.fq", "-b", "minimap2"],
    "unknown aligner mode": ["-1", "r1.fq", "-2", "r2.fq", "-B", "fast"],
    "unknown normalization": ["-1", "r1.fq", "-2", "r2.fq", "-B", "cutsite", "-n", "cosine"],
    "RS without enzyme": ["-1", "r1.fq", "-2", "r2.fq", "-n", "RS"],
    "abundance without depth": ["-1", "r1.fq", "-2", "r2.fq", "-n", "abundance"],
    "unknown start": ["-1", "r1.fq", "-2", "r2.fq", "-S", "sam"],
    "assembly neither fasta nor index": ["-1", "x.pairs", "-S", "pair", "-a", "missing.fa"],
}


@pytest.mark.parametrize("case", list(NETWORK_ERRORS))
def test_network_argument_errors(case, tmp_path):
    args = ["-o", str(tmp_path / "out"), "-T", str(tmp_path / "tmp"), *NETWORK_ERRORS[case]]
    if "-a" not in args:
        args += ["-a", global_args["FASTA"]]
    with pytest.raises(ValueError):
        mtc.Network(args, {}).execute()


PIPELINE_ERRORS = {
    "unknown aligner": ["-b", "minimap2"],
    "unknown aligner mode": ["-B", "fast"],
    "unknown algorithm": ["-B", "cutsite", "-A", "spinglass"],
    "unknown normalization": ["-m", "cosine"],
    "RS without enzyme": ["-m", "RS"],
    "abundance without depth": ["-m", "abundance"],
    "unknown start": ["-S", "network"],
    "fastq without reverse reads": ["-1", "r1.fq", "-S", "fastq"],
    "assembly neither fasta nor index": ["-1", "x.pairs", "-S", "pair", "-a", "missing.fa"],
}


@pytest.mark.parametrize("case", list(PIPELINE_ERRORS))
def test_pipeline_argument_errors(case, tmp_path):
    args = ["-o", str(tmp_path / "out"), "-T", str(tmp_path / "tmp"), *PIPELINE_ERRORS[case]]
    if "-a" not in args:
        args += ["-a", global_args["FASTA"]]
    with pytest.raises(ValueError):
        mtc.Pipeline(args, {}).execute()


@pytest.mark.parametrize("existing", ["overlapping_bin", "recursive_bin", "final_bin_unscaffold", "scaffold_bin"])
def test_pipeline_refuses_existing_outputs(existing, tmp_path):
    (tmp_path / "out" / existing).mkdir(parents=True)
    args = ["-a", global_args["FASTA"], "-1", "x.pairs", "-S", "pair", "--scaffold", "-o", str(tmp_path / "out"),
            "-T", str(tmp_path / "tmp")]  # fmt: skip
    with pytest.raises(ValueError):
        mtc.Pipeline(args, {}).execute()


def test_pipeline_force_replaces_outputs(tmp_path):
    """With --force the output directories are emptied, then the run goes on
    (here up to the missing reverse reads)."""
    for existing in ["overlapping_bin", "recursive_bin", "final_bin_unscaffold", "scaffold_bin"]:
        (tmp_path / "out" / existing).mkdir(parents=True)
        (tmp_path / "out" / existing / "stale.fa").write_text(">old\nA\n")
    args = ["-a", global_args["FASTA"], "-1", "r1.fq", "--scaffold", "-F", "-o", str(tmp_path / "out"),
            "-T", str(tmp_path / "tmp")]  # fmt: skip
    with pytest.raises(ValueError):
        mtc.Pipeline(args, {}).execute()
    assert not list((tmp_path / "out").glob("*/stale.fa"))


@pytest.mark.parametrize("existing", ["recursive_bin", "final_bin"])
def test_validation_refuses_existing_outputs(existing, tmp_path):
    (tmp_path / "out" / existing).mkdir(parents=True)
    args = ("-a {SUBSET_FASTA} -c {CONTIGS_PARTITION} -f {OVERLAPPING_BINS} -n {NETWORK} -o {OUT} -T {TMP}").format(
        OUT=tmp_path / "out", TMP=tmp_path / "tmp", **global_args
    )
    with pytest.raises(ValueError):
        mtc.Validation(args.split(" "), {}).execute()


def test_partition_existing_outputs(metagenome_inputs, metagenome_network, tmp_path):
    """An existing overlapping_bin is refused, or replaced with --force; the
    --prefix names the bins."""
    fasta_dir = tmp_path / "out" / "overlapping_bin"
    fasta_dir.mkdir(parents=True)
    (fasta_dir / "stale.fa").write_text(">old\nA\n")
    args = ["-a", str(metagenome_inputs["assembly"]), "-c", str(metagenome_network["contig_data"]),
            "-n", str(metagenome_network["network"]), "-o", str(tmp_path / "out"), "-T", str(tmp_path / "tmp"),
            "-i", "5", "-s", "50000", "-p", "sample-1"]  # fmt: skip
    with pytest.raises(ValueError):
        mtc.Partition(args, {}).execute()
    with pytest.raises(ValueError):
        mtc.Partition([*args, "-F", "-A", "spinglass"], {}).execute()
    mtc.Partition([*args, "-F"], {}).execute()
    bins = sorted(p.name for p in fasta_dir.iterdir())
    assert "stale.fa" not in bins
    assert bins and all(name.startswith("sample-1_metator_") for name in bins)


def test_check_output_path(tmp_path):
    command = mtc.Pairs(["x.pairs"], {})
    target = tmp_path / "new" / "dir" / "file.txt"
    command.check_output_path(str(target))
    assert target.parent.is_dir()
    target.write_text("x")
    with pytest.raises(IOError, match="--force"):
        command.check_output_path(str(target))
    command.check_output_path(str(target), force=True)
