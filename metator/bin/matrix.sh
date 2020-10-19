#!/usr/bin/env bash
# Jacques Serizay

#Merge annotations with core communities, extract subnetworks
#and FASTA files.

current_dir="$(cd "$(dirname "$0")" && pwd)"
#current_dir=/opt/conda/lib/python3.7/site-packages/metator/bin
#project=after0h
scripts_dir="$current_dir"/../scripts

# shellcheck source=config.sh
# . "$current_dir"/config.sh
# . config.sh

# shellcheck source=environment.sh
. "$current_dir"/environment.sh

# Compute the matrix of pairwise Hamming distances between all pairs of core communities
paste $(printf "${partition_dir}/iteration/%s.community " $(seq 1 "$iter")) >"${tmp_dir}"/"${project}"_iterations_"${iter}".txt




echo "Computing matrix of Hamming distances for ${oter} interations..."

input="${partition_dir}"/partition/core_size_indices_"${iter}".txt
communities="${tmp_dir}"/"${project}"_iterations_"${iter}".txt
indices="${tmp_dir}"/"${project}"_indices_"${iter}".txt
sed 's,^[ ]*,,' "${input}" | sed 's, ,\t,' | cut -f3 > "${indices}"
output="${partition_dir}"/partition/hamming_distance_"${iter}".npz
python3 "${scripts_dir}"/hamming.py \
  --communities "${communities}" \
  --indices "${indices}" \
  --output "${output}" \
  --cores "${threads}"


echo "Done."
