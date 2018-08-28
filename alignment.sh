#!/usr/bin/env bash
# Meta3C team
# Martial Marbouty
# Lyam Baudry
# Théo Foutel-Rodier

#Maps read files onto the assembly and converts the alignments into a network.
#Reads are mapped separately, sorted by names, then interleaved (rather than
#mapped in paired-end mode). This is because bowtie2 has a tendency to leave
#out far-off matches when mapping in paired-end mode.

#An interleaved, sorted alignment file can be directly supplied to
#the network.py script if needed. Both sam and bam are accepted, and the input
#can be compressed.

current_dir="$(cd "$(dirname "$0")" && pwd)"

# shellcheck source=config.sh
. "$current_dir"/config.sh

# shellcheck source=environment.sh
. "$current_dir"/environment.sh

# Default locations
bowtie2_build_executable="bowtie2-build"
bowtie2_executable="bowtie2"
minimap2_executable="minimap2"

if [ "$minimap2" -eq 1 ]; then
  locate_and_set_executable minimap2_executable minimap2
else
  locate_and_set_executable bowtie2_executable bowtie2
  locate_and_set_executable bowtie2_build_executable bowtie2-build bowtie2
  if [[ ! "$bowtie2_executable" ]]; then
    echo "I can't find any way to load bowtie2. Aborting."
    exit 1
  fi
fi

function align() {
  reads=$1
  output=$2
  if [ "$minimap2" -eq 1 ]; then
    "$minimap2_executable" -2 -t "$t" -ax sr "$assembly_dir"/"${project}"_assembly_raw.fa "$reads" >"$output"
  else
    "$bowtie2_executable" --very-sensitive-local -p "$t" -x "$index_name" -U "$reads" >"$output"
  fi
}

locate_and_set_executable pigz_executable pigz

mkdir -p "$assembly_dir"
mkdir -p "$alignment_dir"
mkdir -p "$network_dir"
mkdir -p "$tmp_dir"

t=$((threads / 6 < 1 ? 1 : threads / 6))

#We don't want to modify the original assembly
if [ ! -f "$assembly_dir"/"${project}"_assembly_raw.fa ]; then
  cp "$assembly" "$assembly_dir"/"${project}"_assembly_raw.fa || {
    echo >&2 "I can't find fasta file to map against. Aborting."
    exit 1
  }
  echo "Moved source fasta to $assembly_dir/${project}_assembly_raw.fa"
fi

if [ "$minimap2" -eq 0 ]; then
  if [ ! -f "${index_name}".1.bt2 ]; then
    echo "Bowtie index building..."
    "$bowtie2_build_executable" -f "$assembly_dir"/"${project}"_assembly_raw.fa "$index_name" ||
      {
        echo >&2 "Something wrong happened when building index. Aborting."
        exit 1
      }
    echo "Done."
  fi
fi

# We don't pipe bam files to samtools' sort because it was found to be a huge time bottleneck
echo "Mapping reads..." &

if [ ! -f "$alignment_dir"/"${project}"_forward.bam ]; then

  align "$fastq_for" "$tmp_dir"/"${project}"_forward_raw.sam &&
    samtools view -bS -F 4 -q "$mapping_quality_threshold" -@ "$t" "$tmp_dir"/"${project}"_forward_raw.sam \
      >"$tmp_dir"/"${project}"_forward_unsorted.bam &&
    samtools sort -n -@ "$t" "$tmp_dir"/"${project}"_forward_unsorted.bam \
      >"$alignment_dir"/"${project}"_forward.bam &

fi

if [ ! -f "$alignment_dir"/"${project}"_reverse.bam ]; then

  align "$fastq_rev" "$tmp_dir"/"${project}"_reverse_raw.sam &&
    samtools view -bS -F 4 -q "$mapping_quality_threshold" -@ "$t" "$tmp_dir"/"${project}"_reverse_raw.sam \
      >"$tmp_dir"/"${project}"_reverse_unsorted.bam &&
    samtools sort -n -@ "$t" "$tmp_dir"/"${project}"_reverse_unsorted.bam \
      >"$alignment_dir"/"${project}"_reverse.bam &

fi

wait

echo "Alignment done."
echo "Merging bam files..."

if [ ! -f "$alignment_dir"/"${project}"_merge.bam ]; then
  samtools merge -n "$alignment_dir"/"${project}"_merge.bam "$alignment_dir"/"${project}"_forward.bam "$alignment_dir"/"${project}"_reverse.bam
fi

echo "Done."
echo "Generating network from alignments..."

normalize=''
if [ "$norm" -eq 1 ]; then
  normalize='--normalize'
fi

python "$current_dir"/network.py --input "$alignment_dir"/"${project}"_merge.bam --output "${network_dir}" --map-quality "$mapping_quality_threshold" --chunk-size "$chunk_size" --read-size $((read_size / 2)) --size-chunk-threshold "$size_chunk_threshold" --reference "$assembly_dir"/"${project}"_assembly_raw.fa "$normalize"

echo "Done."

if [ $clean_up -eq 1 ]; then
  echo "Cleaning up..."
  rm "${tmp_dir}"/"${project}"_forward_raw.sam
  rm "${tmp_dir}"/"${project}"_reverse_raw.sam
  rm "${tmp_dir}"/"${project}"_forward_unsorted.bam
  rm "${tmp_dir}"/"${project}"_reverse_unsorted.bam
  echo "Done."
fi
