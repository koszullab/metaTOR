#!/usr/bin/env bash
# Meta3C team
# Martial Marbouty
# Lyam Baudry
# Théo Foutel-Rodier

#This is a global pipeline script for the impatient. 
#It basically checks that everything it needs to run
#is properly installed, then runs all steps of the 
#pipeline sequentially.

current_dir="$( cd "$( dirname "$0" )" && pwd )"

source $current_dir/config.sh
source $current_dir/environment.sh

printf "Checking python is there..."
there_is python
echo "OK."

for lib in numpy scipy pysam matplotlib; do
    printf "Checking $lib is there..."
    python -c "import $lib" >/dev/null 2>&1 || { echo "Error! $lib is missing from your python libraries. Please install it (using either your package manager or pip)"; exit 1; }
    echo "OK."
done

printf "Checking biopython is there..."
python -c "import Bio" >/dev/null 2>&1 || { echo "Error! Biopython is missing from your python libraries. Please install it (using either your package manager or pip)"; exit 1; }
echo "OK."

if [ $minimap -eq 1 ]; then
    aligner="minimap2"
else
    aligner="bowtie2"
fi
printf "Checking $aligner is there..."
check_for $aligner
echo "OK."

printf "Checking samtools is there..."
check_for samtools
echo "OK."

printf "Checking hmmer is there..."
check_for hmmsearch
echo "OK."

printf "Checking prodigal is there..."
locate_and_set_executable prodigal_executable prodigal
echo "OK."

printf "Checking Louvain is there..."
locate_and_set_executable louvain_executable louvain
echo "OK."

printf "Checking permission for everything..."
chmod +x $current_dir/* || { echo "Error! Couldn't get permission for running scripts in ${current_dir}. Please change the location of the directory accordingly."; exit 1; }
echo "OK."

echo "Everything is set! Launching the pipeline."

echo "**** Starting alignment ****"
. $current_dir/alignment.sh
echo "**** Alignment done ****"

echo "**** Starting partition ****"
. $current_dir/partition.sh
echo "**** Partition done ****"

echo "**** Starting annotation ****"
. $current_dir/annotation.sh
echo "**** Annotation done ****"

echo "**** Starting binning ****"
. $current_dir/binning.sh
echo "**** Binning done ****"
