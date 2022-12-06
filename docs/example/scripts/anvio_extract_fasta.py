

"""Small script to extract fasta files from refined anvio bin file. The input 
file have as first column the contig name and the in the second the bin name.

Usage:
python3 ./anvio_extract_fasta.py input_file tmp_dir outdir assembly.fa
"""

import csv
import os
import pyfastx
import sys
from os.path import join, exists
import subprocess as sp 

bin_file = sys.argv[1]
tmpdir = sys.argv[2]
outdir = sys.argv[3]
assembly = sys.argv[4]

# Defined the output and temp directory if they do not already exist.
if not exists(outdir):
    os.makedirs(outdir)
if not exists(tmpdir):
    os.makedirs(tmpdir)

bin_name =None
with open(bin_file, "r") as bin_list:
    reader = csv.reader(bin_list, delimiter="\t")
    line = next(reader)
    contig_name = line[0].split("_split")[0]
    bin_name = line[1]
    fasta_file = join(outdir, "{0}.fa".format(bin_name))
    contigs_file = join(tmpdir, "{0}.txt".format(bin_name))
    f = open(contigs_file, "w")
    f.write(contig_name + "\n")
    for line in reader:
        if bin_name != line[1]:
            # Close the previous bin and extract the fasta.
            f.close()
            cmd = "pyfastx extract {0} -l {1} > {2}".format(
                assembly, contigs_file, fasta_file
            )
            process = sp.Popen(cmd, shell=True)

            # Open a new file for the new bin and write the line
            bin_name = line[1]
            fasta_file = join(outdir, "{0}.fa".format(bin_name))
            contigs_file = join(tmpdir, "{0}.txt".format(bin_name))
            f = open(contigs_file, "w")
            contig_name = line[0].split("_split")[0]
            f.write(contig_name + "\n")
        else:
            # write the contig_id
            contig_name = line[0].split("_split")[0]
            f.write(contig_name + "\n")
    f.close()
