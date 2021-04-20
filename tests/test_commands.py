# Test for CLI toos of metator
# Commands are simply run to test crashes.


import metator.commands as mtc
import pytest
import os


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
    "CONTIGS2": "tests_data/outdir/contig_data_partition.txt",
    "ALIGNMENT": "tests_data/outdir/alignment.bed",
    "OUT": "tests_data/outdir/",
    "OUT_TEST": "tests_data/out_test",
    # "OUT_VAL": "tests_data/outdir_validation/",
    # "OUT_FASTA": "tests_data/outdir_fasta/",
    "TMP": "tests_data/tmp/",
}

os.makedirs(global_args["OUT_TEST"], exist_ok=True)


def test_network():
    args = "-1 {FASTQ_FOR} -2 {FASTQ_REV} -a {FASTA} -o {OUT_TEST} -T {TMP}".format(
        **global_args
    )
    proc = mtc.Network(args.split(" "), {})
    proc.execute()


NORMALIZE = (
    "norm",
    ["None", "length", "abundance", "RS", "empirical_hit", "theoritical_hit"],
)


@pytest.mark.parametrize(*NORMALIZE)
def test_network2(norm):
    args = "-1 {BAM_FOR} -2 {BAM_REV} -a {FASTA_INDEX} -d {DEPTH} -o {OUT_TEST} -T {TMP} -n {0} -e DpnII,HinfI -S bam".format(
        norm, **global_args
    )
    proc = mtc.Network(args.split(" "), {})


ALGORITHM = ("alg", ["louvain", "leiden"])


@pytest.mark.parametrize(*ALGORITHM)
def test_partition(alg):
    args = (
        "-a {FASTA} -c {CONTIGS} -i 5 -n {NETWORK} -o {OUT_TEST} -s 30000 -t 8 -T {TMP} -A {0} -F"
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
        "-1 {FASTQ_FOR} -2 {FASTQ_REV} -a {FASTA_INDEX} -v -F -o {OUT_TEST} -s 30000"
    ).format(**global_args)
    proc = mtc.Pipeline(args.split(" "), {})
    proc.execute()


# shutil.rmtree("tests_data/out_test")
# shutil.rmtree("tests_data/tmp/")
