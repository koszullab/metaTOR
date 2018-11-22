#!/usr/bin/env bash
# Meta3C team
# Martial Marbouty
# Lyam Baudry
# Théo Foutel-Rodier

#Merge annotations with core communities, extract subnetworks
#and FASTA files.

current_dir="$(cd "$(dirname "$0")" && pwd)"
scripts_dir="$current_dir"/../scripts

# shellcheck source=config.sh
# . "$current_dir"/config.sh

# shellcheck source=environment.sh
. "$current_dir"/environment.sh

checkm_executable="checkm"

locate_and_set_executable checkm_executable checkm

mkdir -p $validation_dir
mkdir -p ${validation_dir}_V2

checkm_executable tree -t $threads -x fa $fasta_dir $validation_dir
checkm_executable tree_qa $validation_dir -o 1 -f $chekm_dir/checkM_results.txt
checkm_executable lineage_set $validation_dir $validation_dir/checkM_output_marker.txt
checkm_executable analyze -t $threads -x fa $validation_dir/checkM_output_marker.txt $fasta_dir ${validation_dir}
checkm_executable qa -t $threads $validation_dir/checkM_output_marker.txt ${validation_dir}_V2 -o 1 >${validation_dir}/checkM_results_complete.txt
