#!/usr/bin/env bash

#Run gene prediction software on the assembly (prodigal),
#then annotation software on the predicted genes (hmmer).
#The output is an annotation, i.e. a two column file where each
#chunk is assigned a number of 'hits'.
#Any annotation file can be processed, but the ones generated
#by this script come from banks of:
# *conjugative elements
# *essential genes
# *single gene copies (SGCs)
# *viral orthologous groups (VOGs)

current_dir="$(cd "$(dirname "$0")" && pwd)"

# shellcheck source=config.sh
. "$current_dir"/config.sh

# shellcheck source=environment.sh
. "$current_dir"/environment.sh

#Default values
prodigal_executable="prodigal"
hmmsearch_executable="hmmsearch"

locate_and_set_executable prodigal_executable prodigal
locate_and_set_executable hmmsearch_executable hmmsearch hmmer

mkdir -p "${annotation_dir}"/HMM
mkdir -p "${assembly_dir}"

new_assembly_name=${project}_${size_contig_threshold}.fa

if [ ! -f "$new_assembly_name" ]; then
  python "$current_dir"/fasta_utils.py --input "$assembly" --output "${assembly_dir}"/"${new_assembly_name}" --threshold "${size_contig_threshold}"
fi

if [ ! -f "${annotation_dir}"/"${project}"_prot.fa ]; then
  "$prodigal_executable" -a "${annotation_dir}"/"${project}"_prot.fa -p meta -o "${annotation_dir}"/"${project}".gene -d "${annotation_dir}"/"${project}"_nucl.fa -i "${assembly_dir}"/"${new_assembly_name}" \
    >"${annotation_dir}"/"${project}"_prodigal.log
fi

python "$current_dir"/fasta_utils.py --proteins -i "${annotation_dir}"/"${project}"_prot.fa -o "${annotation_dir}"/"${project}"_prot_renamed.fa
grep ">" "${annotation_dir}"/"${project}"_prot_renamed.fa >"${annotation_dir}"/prot_header_renamed.fa

function annotate() {

  local modele=$1

  "$hmmsearch_executable" "${model_dir}"/"${modele}".hmm "${annotation_dir}"/"${project}"_prot_renamed.fa >"${annotation_dir}"/HMM/"${modele}".out

  grep "__gene" "${annotation_dir}"/HMM/"${modele}".out |
    awk -v evalue="$evalue" 'NF > 7 && $1 < evalue { print $9 }' |
    sed 's/\(__gene.*\)$//' |
    sort -d >"${annotation_dir}"/HMM/"${modele}"_hit_sorted.txt

  uniq -c "${annotation_dir}"/HMM/"${modele}"_hit_sorted.txt | awk '{ print $2,$1 }' | sort -d -k1,1 >"${annotation_dir}"/HMM/"$modele"_hit_weighted.txt

}

for modele in $hmm_databases; do
  annotate "$modele" &
done

wait

echo "Done."
