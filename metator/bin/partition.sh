#!/usr/bin/env bash
# Meta3C team
# Martial Marbouty
# Lyam Baudry

#This script runs the Louvain software many times to partition the network, then
#looks for 'cores' that are easily found by identifying identical lines on the
#global Louvain output.
#The script also draws some figures to get an idea of how cores and their size
#distribution evolve with the number of Louvain iterations.

#Note that the Louvain software is not, in the strictest sense, necessary:
#any program that assigns a node to a community, does so non-deterministically
#and solely outputs a list in the form: 'node_id community_id' could be plugged
#instead.

current_dir="$(cd "$(dirname "$0")" && pwd)"
scripts_dir="$current_dir"/../scripts

# shellcheck source=config.sh
# . "$current_dir"/config.sh

# shellcheck source=environment.sh
. "$current_dir"/environment.sh

# OS detection to choose 'distance' executable
if [[ "$OSTYPE" == "linux"* ]]; then
  distance="distance"
elif [[ "$OSTYPE" == "darwin"* ]]; then
  distance="distance_OSX"
elif [[ "$OSTYPE" == "cygwin" || "$OSTYPE" == "msys" || "$OSTYPE" == "win32" ]]; then
  distance="distance_win"
else
  echo "Warning, your OS is not supported. This means the prebuilt distance binary will likely fail to run."
  echo "The original source file is found at $current_dir/distance.go, please build it yourself at the same location."
  distance="distance"
fi

#Default locations
louvain_executable="$tools_dir"/louvain/louvain
convert_executable="$tools_dir"/louvain/convert
hierarchy_executable="$tools_dir"/louvain/hierarchy

locate_and_set_executable louvain_executable louvain
locate_and_set_executable convert_executable convert louvain
locate_and_set_executable hierarchy_executable hierarchy louvain

#Sanity check: input network file. Takes the form $network_dir/network.txt if generated from alignment.sh
if [ ! -f "$network_file" ]; then
  echo >&2 "Network file not found. Aborting."
  exit 1
fi

mkdir -p "$partition_dir"
mkdir -p "$partition_dir"/iteration
mkdir -p "$partition_dir"/partition
mkdir -p "$tmp_dir"

#Convert network into binary files for Louvain to work on
if [ ! -f "${tmp_dir}"/tmp_"${project}".bin ] && [ ! -f "${tmp_dir}"/tmp_"${project}".weights ]; then
  "$convert_executable" -i "$network_file" -o "${tmp_dir}"/tmp_"${project}".bin -w "${tmp_dir}"/tmp_"${project}".weights
fi

#First, perform Louvain iterations on graphs
function perform_iteration() {
  local current_iteration=$1

  "$louvain_executable" "${tmp_dir}"/tmp_"${project}".bin -l -1 -w "${tmp_dir}"/tmp_"${project}".weights >"${tmp_dir}"/tmp_"${project}"_"${current_iteration}".tree

  "$hierarchy_executable" "${tmp_dir}"/tmp_"${project}"_"${current_iteration}".tree >"${tmp_dir}"/output_louvain_"${project}"_"${current_iteration}".txt

  level="$(tail -1 "${tmp_dir}"/output_louvain_"${project}"_"${current_iteration}".txt | gawk '{print $2}' | sed 's/://g')"

  "$hierarchy_executable" "${tmp_dir}"/tmp_"${project}"_"${current_iteration}".tree -l "${level}" | cut -f 2 -d ' ' >"${partition_dir}"/iteration/"${current_iteration}".community
}

#Then, collect all the results to try and identify bins (identified by identical lines in the pasted community file)
function resolve_partition() {

  local repet=$1

  #Merge all Louvain outputs
  paste $(printf "${partition_dir}/iteration/%s.community " $(seq 1 "$repet")) >"${tmp_dir}"/"${project}"_iterations_"${repet}".txt

  #Identify bins and sort them
  "$working_dir"/$distance -m cores -i "${tmp_dir}"/"${project}"_iterations_"${repet}".txt |
    sort -nr -k1,1 --parallel="$threads" |
    cat -n \
      >"${partition_dir}"/partition/core_size_indices_"${repet}".txt

  #We use this one-liner because awk can't process more than 32767 fields
  perl -nae 'print join("\t", $_, $F[0], $F[1]), "\n" for @F[2..$#F];' "${partition_dir}"/partition/core_size_indices_"${repet}".txt |
    sort -n -k1,1 --parallel="$threads" >"${partition_dir}"/partition/chunkid_core_size_"${repet}".txt

  #Slower but more reliable than paste: sometimes it's handy to have chunks listed by ids or by names, so we generate both
  gawk '
    NR == FNR {
      names[$1] = $2
      next
    }

    $1 in names {
      name = names[$1]
      $1 = ""
      print name,$0
    }
  ' "${network_dir}"/idx_contig_hit_size_cov.txt "${partition_dir}"/partition/chunkid_core_size_"${repet}".txt \
    >"${partition_dir}"/partition/chunkname_core_size_"${repet}".txt

  gawk '{ print $2 }' "${partition_dir}"/partition/core_size_indices_"${repet}".txt >"${tmp_dir}"/"${project}"_sizes_"${repet}".txt

  #Draw some figures to have an idea of how bin sizes evolve as more Louvain iterations are computed
  gawk -v repet="$repet" '
    BEGIN {
      count_100 = 0
      count_500 = 0
      count_1000 = 0
    }

    $1 > 100 && $1 < 500 {count_100 += 1} 
    $1 >= 500 && $1 < 1000 {count_500 += 1}
    $1 >= 1000 {count_1000 += 1}

    END {
      print repet,count_100 >> "'"${partition_dir}"/partition/regression_louvain_100.txt'"
      print repet,count_500 >> "'"${partition_dir}"/partition/regression_louvain_500.txt'"
      print repet,count_1000 >> "'"${partition_dir}"/partition/regression_louvain_1000.txt'"
    }
  ' "${tmp_dir}"/"${project}"_sizes_"${repet}".txt

  python3 "$scripts_dir"/figures.py --barplots "${tmp_dir}"/"${project}"_sizes_"${repet}".txt -o "${partition_dir}"/partition/repartition_"${repet}".pdf

}

rm -f "${partition_dir}"/partition/regression_louvain_100.txt
rm -f "${partition_dir}"/partition/regression_louvain_500.txt
rm -f "${partition_dir}"/partition/regression_louvain_1000.txt

echo "Performing iterations..."
for iteration in $(seq "$iterations"); do
  perform_iteration "$iteration"
done
wait

#We resolve partitions at different points in order to get an idea of how stable cores can get.
echo "Resolving partitions..."
for column in 1 5 10 20 30 40 50 60 70 80 90 100 150 200 500 1000 $iterations; do
  if [ "$iterations" -ge "$column" ]; then
    resolve_partition "$column"
  fi
done

wait

# Compute the matrix of pairwise Hamming distances between all pairs of core communities
echo "Computing matrix of Hamming distances for ${iterations} interations..."
input="${partition_dir}"/partition/core_size_indices_"${iterations}".txt
communities="${tmp_dir}"/"${project}"_iterations_"${iterations}".txt
indices="${tmp_dir}"/"${project}"_indices_"${iterations}".txt
sed 's,^[ ]*,,' "${input}" | sed 's, ,\t,' | cut -f3 > "${indices}"
output="${partition_dir}"/partition/hamming_distance_"${iterations}".npz
python3 "${scripts_dir}"/hamming.py \
  --communities "${communities}" \
  --indices "${indices}" \
  --output "${output}" \
  --cores "${threads}"

wait

# echo "Drawing some figures..."
# 
# for u in 100 500 1000; do
#   python3 "$scripts_dir"/figures.py --plots "${partition_dir}"/partition/regression_louvain_"${u}".txt -o "${partition_dir}"/partition/regression_"${u}".pdf &
# done
# 
# wait

echo "Cleaning up..."

if [ "$clean_up" -eq 1 ]; then
  rm "${tmp_dir}"/tmp_"${project}".weights
  rm "${tmp_dir}"/tmp_"${project}"_*.tree
  rm "${tmp_dir}"/output_louvain_"${project}"_*.txt
  rm "${tmp_dir}"/"${project}"_iterations_*.txt
  rm "${tmp_dir}"/"${project}"_sizes_*.txt
fi
