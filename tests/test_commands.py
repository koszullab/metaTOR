# Test for CLI toos of metator
# Commands are simply run to test crashes.


import metator.commands as mtc
import pytest


# Use global variables for input files
global_args = {
    "FASTA": "tests_data/assembly.fa",
    "FASTQ_FOR": "tests_data/for_paired.fq",
    "FASTQ_REV": "tests_data/rev_paired.fq",
    "FASTA_INDEX": "tests_data/assembly",
    "NETWORK": "tests_data/outdir/network.txt",
    "CONTIGS": "tests_data/outdir/contig_data_network.txt",
    "CONTIGS2": "tests_data/outdir/contig_data_partition.txt",
    "ALIGNMENT": "tests_data/outdir/alignment.bed",
    "OUT": "tests_data/outdir/",
    "OUT_VAL": "tests_data/outdir_validation/",
    "OUT_FASTA": "tests_data/outdir_fasta/",
    "TMP": "tests_data/tmp/",
    "MODE": ("mode", ["for_vs_rev", "all", "pile"]),
}


def test_align():
    args = (
        "-1 {FASTQ_FOR} -2 {FASTQ_REV} -g {FASTA_INDEX} -o {OUT} -q 30 -t 8 -T {TMP}"
    ).format(**global_args)
    proc = mtc.Align(args.split(" "), {})
    proc.execute()


@pytest.mark.parametrize("mode", ["for_vs_rev", "all", "pile"])
def test_cutsite(mode):
    args = (
        "-1 {FASTQ_FOR} -2 {FASTQ_REV} -e DpnII,HinfI -m {0} -o {OUT} -t 8"
    ).format(mode, **global_args)
    proc = mtc.Cutsite(args.split(" "), {})
    proc.execute()


def test_network():
    args = ("-g {FASTA} -i {ALIGNMENT} -n -o {OUT} -t 8 -T {TMP}").format(
        **global_args
    )
    proc = mtc.Network(args.split(" "), {})
    proc.execute()


def test_partition():
    args = (
        "-a {FASTA} -c {CONTIGS} -i 5 -n {NETWORK} -o {OUT} -O 90 -s 30000 -t 8 -T {TMP}"
    ).format(**global_args)
    proc = mtc.Partition(args.split(" "), {})
    proc.execute()


def test_validation():
    args = (
        "-a {FASTA} -c {CONTIGS2} -f {OUT_FASTA} -i 5 -n {NETWORK} -o {OUT_VAL} -s 10000 -t 8 -T {TMP}"
    ).format(**global_args)
    proc = mtc.Validation(args.split(" "), {})
    proc.execute()


# def test_pipeline():
#     args = (
#         "-1 {FASTQ_FOR} -2 {FASTQ_REV} -g {FASTA_INDEX} -i 5 -n -o {OUT} -O 90 -q 30 -s 30000 -t 8 -T {TMP}"
#     ).format(**global_args)
#     proc = mtc.Pipeline(args.split(" "), {})
#     proc.execute()
