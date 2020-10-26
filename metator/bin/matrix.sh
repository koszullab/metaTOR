#!/usr/bin/env bash
# Jacques Serizay

#Merge annotations with core communities, extract subnetworks
#and FASTA files.

current_dir="$(cd "$(dirname "$0")" && pwd)"
#current_dir=/pasteur/homes/jaseriza/.local/bin/conda/lib/python3.7/site-packages/metator/bin
#current_dir=/pasteur/homes/jaseriza/scratch/Projects/test/metaTOR/metator/bin
#current_dir=/opt/conda/lib/python3.7/site-packages/metator/bin
#project=after0h
scripts_dir="$current_dir"/../scripts

# shellcheck source=config.sh
# . "$current_dir"/config.sh
# . config.sh

# shellcheck source=environment.sh
. "$current_dir"/environment.sh

# Create directory for network matrices
if [ ! -d "${partition_dir}"/network_matrices/ ]
then
    mkdir "${partition_dir}"/network_matrices/
fi

# Compute the matrix of pairwise Hamming distances between all pairs of core communities and CC infos
echo "Computing matrix of Hamming distances for ${iter} interations..."
#
communities="${tmp_dir}"/"${project}"_iterations_"${iter}".txt
paste $(printf "${partition_dir}/iteration/%s.community " $(seq 1 "${iter}")) > "${communities}"
#
indices="${tmp_dir}"/"${project}"_indices_"${iter}".txt
input="${partition_dir}"/partition/core_size_indices_"${iter}".txt
sed 's,^[ ]*,,' "${input}" | sed 's, ,\t,' | cut -f3 > "${indices}"
#
contigs="${partition_dir}"/network_matrices/"${project}"_"${iter}"-iters_contigID-coreID-contigSize.txt
contigs_info="${partition_dir}"/partition/chunkid_core_size_"${iter}".txt
paste \
  <(cat ${contigs_info} | cut -f1,2) \
  <(cat "${assembly_dir}"/"${project}"_"${size_contig_threshold}".fa | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' | cut -f2) \
  > "${contigs}"
#
output="${partition_dir}"/network_matrices/"${project}"_"${iter}"-iters_hamming_distance.npz
#
python3 "${scripts_dir}"/hamming.py \
  --communities "${communities}" \
  --indices "${indices}" \
  --contigs "${contigs}" \
  --output "${output}" \
  --cores "${threads}"



echo "Done."
