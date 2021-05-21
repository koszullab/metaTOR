# Test for CLI toos of metator
# Commands are simply run to test crashes.


import metator.commands as mtc
import pytest
import shutil


# Use global variables for input files
global_args = {
    "BAM_FOR": "tests_data/outdir/alignment_0_for.bam",
    "BAM_REV": "tests_data/outdir/alignment_0_rev.bam",
    "DEPTH": "tests_data/depth.txt",
    "FASTA": "tests_data/assembly.fa",
    "FASTQ_FOR": "tests_data/for_paired.fq",
    "FASTQ_REV": "tests_data/rev_paired.fq",
    "FASTA_INDEX": "tests_data/assembly",
    "NETWORK": "tests_data/outdir/network.txt",
    "CONTIGS": "tests_data/outdir/contig_data_network.txt",
    "OUT_TEST": "tests_data/out_test",
    "PAIRS": "tests_data/outdir/alignment.pairs",
    "TMP": "tests_data/tmp/",
}
NORMALIZE = (
    "norm",
    ["None", "length", "abundance", "RS", "empirical_hit", "theoritical_hit"],
)
ALGORITHM = ("alg", ["louvain", "leiden"])


def test_network():
    args = "-1 {FASTQ_FOR} -2 {FASTQ_REV} -a {FASTA} -o {OUT_TEST} -T {TMP}".format(
        **global_args
    )
    proc = mtc.Network(args.split(" "), {})
    proc.execute()


def test_network2():
    args = "-1 {BAM_FOR} -2 {BAM_REV} -a {FASTA_INDEX} -d {DEPTH} -o {OUT_TEST} -T {TMP} -S bam".format(
        **global_args
    )
    proc = mtc.Network(args.split(" "), {})
    proc.execute()


@pytest.mark.parametrize(*NORMALIZE)
def test_network3(norm):
    args = "-1 {PAIRS} -a {FASTA} -d {DEPTH} -o {OUT_TEST} -T {TMP} -n {0} -e HindIII,DpnII -S pair".format(
        norm, **global_args
    )
    proc = mtc.Network(args.split(" "), {})
    proc.execute()


@pytest.mark.parametrize(*ALGORITHM)
def test_partition(alg):
    args = (
        "-a {FASTA} -c {CONTIGS} -i 5 -n {NETWORK} -o {OUT_TEST} -s 30000 -t 8 -T {TMP} -A {0} -FC"
    ).format(alg, **global_args)
    proc = mtc.Partition(args.split(" "), {})
    proc.execute()


# Do not do test for validation as checkM is too slow and ask two much memory.
# def test_validation():
#     args = (
#         "-a {FASTA} -c {CONTIGS2} -f {OUT_FASTA} -i 5 -n {NETWORK} -o {OUT_VAL} -s 10000 -t 8 -T {TMP}"
#     ).format(**global_args)
#     proc = mtc.Validation(args.split(" "), {})
#     proc.execute()


def test_pipeline():
    args = (
        "-1 {FASTQ_FOR} -2 {FASTQ_REV} -a {FASTA_INDEX} -v -F -o {OUT_TEST} -s 30000 -C"
    ).format(**global_args)
    proc = mtc.Pipeline(args.split(" "), {})
    proc.execute()

    shutil.rmtree("tests_data/out_test")
    shutil.rmtree("tests_data/tmp/")
