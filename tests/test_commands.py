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
    "OUT_PIPELINE": "tests_data/outdir_pipeline",
    "HOSTMAG_FASTA": "tests_data/outdir_pipeline/metator_00005_00000.fa",
    "HOST_DATA": "tests_data/host_test_data",
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


def test_qc(tmp_path):
    """Test the metator qc command with subset data."""
    args = (
        "--assembly {SUBSET_FASTA} --enzyme HindIII,DpnII --bin-summary {OUT_PIPELINE}/bin_summary.txt --contig-data {OUT_PIPELINE}/contig_data_final.txt  -o {OUT_TEST} -P -R 1.5 {SUBSET_PAIRS}"
    ).format(
        OUT_TEST=Path(tmp_path, "out_test"),
        **global_args,
    )
    print("CLI command:\n", "metator qc", *args.split(" "))

    proc = mtc.Qc(args.split(" "), {})
    proc.execute()


def test_contactmap(tmp_path):
    """Test the metator contactmap command with subset data."""
    os.makedirs(Path(tmp_path, "out_test"), exist_ok=True)
    args = (
        "--assembly {SUBSET_FASTA} --enzyme HindIII,DpnII --contig-data {OUT_PIPELINE}/contig_data_final.txt --name metator_00005_00000 --filter --mat-fmt cool -o {OUT_TEST} --pcr-dup {SUBSET_PAIRS}"
    ).format(
        OUT_TEST=Path(tmp_path, "out_test"),
        **global_args,
    )
    print("CLI command:\n", "metator contactmap", *args.split(" "))

    proc = mtc.Contactmap(args.split(" "), {})
    proc.execute()


def test_scaffold(tmp_path):
    """Test the metator scaffold command with subset data."""
    os.makedirs(Path(tmp_path, "out_test"), exist_ok=True)
    args = (
        "--bin-name metator_00005_00000 --input-fasta {HOSTMAG_FASTA} --out-fasta {OUT_TEST}/metator_00005_00000_scaffolded.fa --out-frags {OUT_TEST}/metator_00005_00000_scaffolded.frags {SUBSET_PAIRS}"
    ).format(
        OUT_TEST=Path(tmp_path, "out_test"),
        **global_args,
    )
    print("CLI command:\n", "metator scaffold", *args.split(" "))

    proc = mtc.Scaffold(args.split(" "), {})
    proc.execute()


def test_host(tmp_path):
    """Test the metator host command with subset data."""
    os.makedirs(Path(tmp_path, "out_test"), exist_ok=True)
    args = (
        "--network {HOST_DATA}/network_0.txt --binning {HOST_DATA}/bin_summary.txt --contigs-data {HOST_DATA}/contig_data_final.txt --mges-bin-data {HOST_DATA}/mges_bin_summary.tsv --outdir {OUT_TEST}"
    ).format(
        OUT_TEST=Path(tmp_path, "out_test"),
        **global_args,
    )
    print("CLI command:\n", "metator host", *args.split(" "))

    proc = mtc.Host(args.split(" "), {})
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
