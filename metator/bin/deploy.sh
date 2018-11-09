#!/usr/bin/env bash
# Meta3C team
# Martial Marbouty
# Lyam Baudry
# Théo Foutel-Rodier

#Install all dependencies on Ubuntu. Should be run as root.
#This script caters to users who just want to run things quickly
#and sets up the appropriate environment for them.
#If you don't like the automated steps it takes, it is assumed
#that you know what you are doing and can setup the environment
#yourself.

current_dir="$(cd "$(dirname "$0")" && pwd)"

# shellcheck source=config.sh
. "$current_dir"/config.sh

# shellcheck source=environment.sh
. "$current_dir"/environment.sh

set -o pipefail
set -e

if [ "$EUID" -ne 0 ]; then
  echo "Please run this script as root."
  exit 1
fi

echo "Installing third-party software dependencies..."
apt install python3-pip bowtie2 samtools hmmer prodigal
echo "OK."

echo "Installing python dependencies..."
pip3 install -U numpy scipy pysam matplotlib seaborn biopython
echo "OK."

./metator.sh dependencies

echo "You should be good to go! Run './meta3c.sh pipeline -1 reads_for.fastq -2 reads_rev.fastq -a assembly.fa' to proceed."
