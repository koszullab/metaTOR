#!/bin/sh

#SBATCH --qos=normal
#SBATCH --mem=96G
#SBATCH -o checkM_big.out.txt -e checkM_big.err.txt
#SBATCH -c 16

set -eu

#For the cluster uncomment these reads. 
#source /local/gensoft2/adm/etc/profile.d/modules.sh
#module purge
#module load hmmer/3.1b2 pplacer/v1.1-alpha18-2 prodigal/2.6.3
#module load CheckM/1.0.7

#Usage checkM.sh /data/TARA fasta_recursive output_recursive 
#project=/data/TARA
#input=/data/TARA/fasta_recursive
#output=/data/TARA/checkM_data/output_recursive

projects=$1
input=$2
output=$3
end=$4
threads=$5

echo "data in progress"

mkdir -p "$projects"/checkM_data/"$output"
mkdir -p "$projects"/checkM_data/"$output"_2/
mkdir -p "$projects"/checkM_data/"$output"_lineage/
mkdir -p "$projects"/checkM_data/"$output"/tmp/

checkm tree -t $threads -x $end "$projects"/"$input" "$projects"/checkM_data/"$output"

checkm tree_qa "$projects"/checkM_data/"$output" -o 1 -f "$projects"/checkM_data/"$output"/1.txt
#checkm tree_qa "$projects"/checkM_data/"$output" -o 2 -f "$projects"/checkM_data/"$output"/2.txt
#checkm tree_qa "$projects"/checkM_data/"$output" -o 3 -f "$projects"/checkM_data/"$output"/3.txt
#checkm tree_qa "$projects"/checkM_data/"$output" -o 4 -f "$projects"/checkM_data/"$output"/4.txt
#checkm tree_qa "$projects"/checkM_data/"$output" -o 5 -f "$projects"/checkM_data/"$output"/5.txt

echo "checkM lineage marker set"

checkm lineage_set "$projects"/checkM_data/"$output" "$projects"/checkM_data/"$output"/checkM_output_marker.txt

echo "checkM analyse and qa"

checkm analyze -t $threads --tmpdir "$projects"/checkM_data/"$output"/tmp -x $end "$projects"/checkM_data/"$output"/checkM_output_marker.txt "$projects"/"$input" "$projects"/checkM_data/"$output"_2/

checkm qa -t $threads --tmpdir "$projects"/checkM_data/"$output"/tmp "$projects"/checkM_data/"$output"/checkM_output_marker.txt  "$projects"/checkM_data/"$output"_2/ -o 1 > "$projects"/checkM_data/"$output"/checkM_results_complete_1.txt
checkm qa -t $threads --tmpdir "$projects"/checkM_data/"$output"/tmp "$projects"/checkM_data/"$output"/checkM_output_marker.txt  "$projects"/checkM_data/"$output"_2/ -o 2 >  "$projects"/checkM_data/"$output"/checkM_results_complete_2.txt

#checkm lineage_wf -t $threads -x $end "$projects"/"$input" "$projects"/checkM_data/"$output"_lineage

sed 's/ \+/ /g' "$projects"/checkM_data/"$output"/checkM_results_complete_1.txt > "$projects"/checkM_data/"$output"/checkM_results_complete_simplify.txt
