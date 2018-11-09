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
mkdir -p "${partition_dir}"/fasta/iteration"$iter"
mkdir -p "${partition_dir}"/fasta_merged/iteration"$iter"
mkdir -p "${partition_dir}"/bin_matrices/iteration"$iter"
mkdir -p "${partition_dir}"/subnetworks/iteration"$iter"
mkdir -p "${partition_dir}"/figures/iteration"$iter"/data

function annotation_distribution() {

  local modele=$1
  local iter=$2

  awk '
    NR == FNR {
        hits[$1] = $2
        next
    }
    $1 in hits {
        print $2, $3
    }
    ' "${annotation_dir}"/HMM/"$modele"_hit_weighted.txt "${partition_dir}"/partition/chunkname_core_size_"${iter}".txt |
    sort | uniq -c | awk '{ print $2, $3, $1 }' >"${partition_dir}"/figures/iteration"$iter"/data/"$modele"_core_size_hit.txt

  awk '
    NR == FNR {
        core_names[$1] = $3
        next
    }
    $1 in core_names {
        print $1, $2, $3/core_names[$1]
    }
    ' "${partition_dir}"/figures/iteration"$iter"/data/total_core_size_hit.txt "${partition_dir}"/figures/iteration"$iter"/data/"$modele"_core_size_hit.txt >"${partition_dir}"/figures/iteration"$iter"/data/"$modele"_core_size_proportional_hit.txt

  python3 "$scripts_dir"/figures.py --logplots "${partition_dir}"/figures/iteration"$iter"/data/"$modele"_core_size_proportional_hit.txt --output "${partition_dir}"/figures/iteration"$iter"/"$modele"_proportional
  python3 "$scripts_dir"/figures.py --logplots "${partition_dir}"/figures/iteration"$iter"/data/"$modele"_core_size_hit.txt --output "${partition_dir}"/figures/iteration"$iter"/"$modele"

}

sed 's/\(__gene.*\)$//' "${annotation_dir}"/prot_header_renamed.fa | cut -c2- |
  awk '
    NR == FNR { 
        core_names[$1] = $2" "$3
        next
    }
    $1 in core_names {
        total_hits[core_names[$1]] += 1
    }
    END {
        for (core in total_hits) {
            print core, total_hits[core]
        }
    }
' "${partition_dir}"/partition/chunkname_core_size_"${iter}".txt - | sort -k1,1n >"${partition_dir}"/figures/iteration"$iter"/data/total_core_size_hit.txt

echo "Drawing enrichment vs. size plots..."

model_list=""
for model in $hmm_databases; do
  model_basename="$(basename "${model%.hmm}")"
  annotation_distribution "$model_basename" "$iter"
  model_list="$model_list $model_basename"
done

echo "Drawing distribution violinplot..."
python3 "$scripts_dir"/figures.py --violin "${partition_dir}"/figures/iteration"$iter"/data/total_core_size_hit.txt $(printf "${partition_dir}/figures/iteration${iter}/data/%s_core_size_hit.txt " $model_list) --output "${partition_dir}"/figures/iteration"$iter"/distrib_annot_violinplot.pdf

echo "Extracting bin subnetworks and matrices..."
python3 "$scripts_dir"/bins.py --input "${partition_dir}"/partition/chunkid_core_size_"${iter}".txt --network "${network_dir}"/network.txt --output "${partition_dir}"/subnetworks/iteration"$iter" --chunk-size "$chunk_size" --n-bins "$n_bins"

echo "Extracting bin FASTA files..."
python3 "$scripts_dir"/bins.py --input "${partition_dir}"/partition/chunkname_core_size_"${iter}".txt --fasta "$assembly" --output "${partition_dir}"/fasta/iteration"$iter" --n-bins "$n_bins" --chunk-size "$chunk_size"

echo "Merging bin chunks..."
for f in $(seq 1 "$n_bins"); do
  if [ -f "${partition_dir}"/fasta/iteration"$iter"/core_"$f".fa ]; then
    python3 "$scripts_dir"/bins.py --merge "${partition_dir}"/fasta/iteration"$iter"/core_"$f".fa --output "${partition_dir}"/fasta_merged/iteration"$iter"
  fi
done

echo "Drawing global matrix..."
sort -k2,2n -k1,1n "${partition_dir}"/partition/chunkid_core_size_"${iter}".txt | cat -n | unexpand -a >"${partition_dir}"/partition/matid_chunkid_core_size_"${iter}".txt

awk 'NR==FNR && $3 < '"$n_bins"' {    
        idx[$2] = $1
        next
    }

    NR > FNR && $1 in idx && $2 in idx {
        print idx[$1], idx[$2], $3
    }
    ' "${partition_dir}"/partition/matid_chunkid_core_size_"${iter}".txt "${network_dir}"/network.txt >"${partition_dir}"/partition/sparse_mat_"${iter}".txt

python3 "$scripts_dir"/figures.py --sparse "${partition_dir}"/partition/sparse_mat_"${iter}".txt --output "${partition_dir}"/bin_matrices/iteration"$iter"/sparse_mat_"${iter}".eps

echo "Drawing core matrix..."
awk '$3 > 4 && NR == FNR {
        core[$1] = $2
        next
    }
    $1 in core && $2 in core {
        core1 = core[$1]
        core2 = core[$2]
        if (core1 < core2) {
            contacts[core1"\t"core2] += $3
        }
        else {
            contacts[core2"\t"core2] += $3
        }
    }
    END {
        for (edge in contacts) {
            print edge, contacts[edge]
        }
    }
' "${partition_dir}"/partition/chunkid_core_size_"${iter}".txt "${network_dir}"/network.txt >"${partition_dir}"/partition/core_network_"${iter}".txt

python3 "$scripts_dir"/figures.py --sparse "${partition_dir}"/partition/core_network_"${iter}".txt --output "${partition_dir}"/bin_matrices/iteration"$iter"/core_matrix_"${iter}".eps

echo "Done."
