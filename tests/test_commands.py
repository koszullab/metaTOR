# Test for CLI toos of metator
# Commands are simply run to test crashes.


import metator.commands as mtc
import pytest
import shutil
import os


# Use global variables for input files
global_args = {
    "BAM_FOR": "tests_data/outdir/alignment_0_for.bam",
    "BAM_REV": "tests_data/outdir/alignment_0_rev.bam",
    "DEPTH": "tests_data/depth.txt",
    "FASTA": "tests_data/assembly.fa",
    "FASTA_VAL": "tests_data/outdir_validation/assembly_val.fa",
    "FASTQ_FOR": "tests_data/for_paired.fq.gz",
    "FASTQ_REV": "tests_data/rev_paired.fq.gz",
    "FASTA_INDEX": "tests_data/assembly",
    "NETWORK": "tests_data/outdir/network.txt",
    "NETWORK_VAL": "tests_data/outdir_validation/network.txt",
    "CONTIGS": "tests_data/outdir/contig_data_network.txt",
    "CONTIGS_VAL": "tests_data/outdir_validation/contig_data_final.txt",
    "OUT_FASTA": "tests_data/outdir_validation/overlapping_bin",
    "OUT_TEST": "tests_data/out_test",
    "OUT_DIR": "tests_data/outdir",
    "PAIRS": "tests_data/outdir/alignment.pairs",
    "TMP": "tests_data/tmp/",
}
NORMALIZE = (
    "norm",
    ["None", "length", "abundance", "RS", "empirical_hit", "theoritical_hit"],
)
ALGORITHM = ("alg", ["louvain", "leiden"])
ALIGNER = ("aligner", ["bwa", "bowtie2"])
MODE = ("mode", ["normal", "iterative", "cutsite"])


@pytest.mark.parametrize(*ALIGNER)
def test_network(aligner):
    args = "-1 {FASTQ_FOR} -2 {FASTQ_REV} -a {FASTA} -b {0} -o {OUT_TEST} -T {TMP}".format(
        aligner,
        **global_args,
    )
    proc = mtc.Network(args.split(" "), {})
    proc.execute()


def test_network2():
    args = "-1 {BAM_FOR} -2 {BAM_REV} -a {FASTA_INDEX} -d {DEPTH} -E 500 -o {OUT_TEST} -T {TMP} -S bam -N".format(
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


@pytest.mark.parametrize(*MODE)
def test_network4(mode):
    args = "-1 {FASTQ_FOR} -2 {FASTQ_REV} -a {FASTA} -o {OUT_TEST} -T {TMP} -B {0} -t 4 -e HindIII,DpnII".format(
        mode, **global_args
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


# def test_validation():
#     args = (
#         "-a {FASTA_VAL} -c {CONTIGS_VAL} -f {OUT_FASTA} -i 5 -n {NETWORK_VAL} -o {OUT_TEST} -t 8 -T {TMP} -F"
#     ).format(**global_args)
#     proc = mtc.Validation(args.split(" "), {})
#     proc.execute()


# def test_pipeline():
#     args = (
#         "-1 {FASTQ_FOR} -2 {FASTQ_REV} -a {FASTA} -F -o {OUT_TEST} -s 30000 -C"
#     ).format(**global_args)
#     proc = mtc.Pipeline(args.split(" "), {})
#     proc.execute()


# def qc():
#     args = (
#         "-a {FASTA} -F -O {OUT_DIR} -o {OUT_TEST}/QC -p test -e HindIII,DpnII -T {TMP} -P"
#     ).format(**global_args)
#     proc = mtc.Qc(args.split(" "), {})
#     proc.execute()


def test_pairs():
    pairs1 = "{OUT_TEST}/align1.pairs".format(**global_args)
    pairs2 = "{OUT_TEST}/align2.pairs".format(**global_args)
    pairs = "{PAIRS}".format(**global_args)
    os.makedirs("{OUT_TEST}".format(**global_args), exist_ok=True)
    shutil.copyfile(pairs, pairs1)
    shutil.copyfile(pairs, pairs2)
    args = f"-rF {pairs1} {pairs2}"
    proc = mtc.Pairs(args.split(" "), {})
    proc.execute()


# shutil.rmtree("tests_data/out_test")
