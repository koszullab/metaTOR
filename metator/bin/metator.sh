#!/usr/bin/env bash
# Meta3C team
# Martial Marbouty
# Théo Foutel-Rodier
# Lyam Baudry

current_version="0.1.1"
last_update_date="[April 2018]"

function display_help() {
  echo ""
  echo "   metaTOR - a pipeline for binning metagenomic datasets from 3C data."
  echo ""
  echo "   Written by L. Baudry, T. Foutel-Rodier and M. Marbouty (with contributions from A. Cournac and V. Scolari)"
  echo "   Spatial Genome Regulation laboratory (Institut Pasteur, Paris)"
  echo ""
  echo "   Version $current_version $last_update_date"
  echo ""
  echo "      Usage: metator {align|partition|annotation|binning|pipeline} [parameters]"
  echo ""
  echo "      A metaTOR command takes the form 'metator action --param1 arg1 --param2 arg2' etc."
  echo ""
  echo "      There are four actions or steps in the meta3c pipeline. They must be run in this order:"
  echo ""
  echo "          -align: map paired-end reads on a preliminary assembly, then generate a network from"
  echo "                  detected contacts between DNA chunks."
  echo "          -partition: perform the Louvain community detection algorithm many times to isolate"
  echo "                      chunks that consistently cluster together for binning purposes."
  echo "          -annotation: run standard annotation software on the assembly (namely gene prediction"
  echo "                       and database comparison) to match with the bins"
  echo "          -binning: match annotations to bins, extract bin genomes and subnetworks, build bin-local"
  echo "                    and global contact maps"
  echo ""
  echo "       After the last step is completed, you should have at your disposal a set of bins, their relative"
  echo "       enrichments in various gene categories as well as the contact map of each bin."
  echo ""
  echo "       In addition, there are a number of optional, miscellaneous actions:"
  echo ""
  echo "          -pipeline: check the environment is right, then run all four of the above sequentially."
  echo "                     This can take a while."
  echo "          -deploy: set up the environment on Ubuntu 14.04 and higher"
  echo "          -dependencies: download some third party dependencies"
  echo "          -version: display current version number"
  echo "          -help: display this (hopefully useful) help message"
  echo ""
  echo "      Please refer to the metaTOR manual for detailed explanations on the parameters."
  echo ""
  echo "   Contact: lyam.baudry@pasteur.fr or romain.koszul@pasteur.fr"
  echo ""
  exit 1
}

function fetch_dependencies() {

  set -e
  hmm_url="http://dl.pasteur.fr/fop/LItxiFe9/hmm_databases.tgz"
  prodigal_base="https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal"

  cd "$tools_dir" || {
    echo "Could not access ${tools_dir}. Aborting."
    exit 1
  }

  if [[ "$OSTYPE" == "linux"* ]]; then
    prodigal_url="$prodigal_base.linux"
    louvain_url="https://lip6.github.io/Louvain-BinaryBuild/louvain_linux.tar.gz"
  elif [[ "$OSTYPE" == "darwin"* ]]; then
    louvain_url="https://github.com/lip6/Louvain-BinaryBuilds/raw/osx/louvain_osx.tar.gz"
    prodigal_url="$prodigal_base.osx.10.9.5"
  elif [[ "$OSTYPE" == "cygwin" || "$OSTYPE" == "msys" || "$OSTYPE" == "win32" ]]; then
    louvain_url="https://ci.appveyor.com/api/projects/yanntm/Louvain-BinaryBuild/artifacts/website/louvain_windows.tar.gz"
    prodigal_url="$prodigal_base.windows.exe"
  else
    echo "Your OS is not supported."
    echo "Please download Louvain and prodigal at the following urls:"
    echo "   -https://sourceforge.net/projects/louvain/files/louvain-generic.tar.gz"
    echo "   -https://github.com/hyattpd/Prodigal"
    echo "and compile them before moving them into your $tools_dir directory."
    exit 1
  fi

  echo "Fetching prodigal..."
  wget -q $prodigal_url -O "$tools_dir/prodigal"
  chmod +x "$tools_dir/prodigal"

  echo "Fetching louvain..."
  mkdir -p "$tools_dir/louvain"
  wget -q $louvain_url -O "$tools_dir/louvain/louvain.tar.gz"
  tar -xzf "$tools_dir/louvain/louvain.tar.gz" -C "$tools_dir/louvain"
  chmod +x $tools_dir/louvain/*
  rm -f "$tools_dir/louvain/louvain.tar.gz"

  echo "Fetching HMMs..."
  mkdir -p "$model_dir"
  wget -q $hmm_url -O "$model_dir"/hmm_databases.tgz
  tar -xzf "$model_dir"/hmm_databases.tgz -C "$model_dir"
  rm -f "$model_dir"/hmm_databases.tar.gz

  echo "All dependencies successfully installed in $tools_dir."
  set +e
}

#Parse a bunch of arguments for the entirety of the pipeline.
#The idea is that they can be set only once and remain unchanged
#after other scripts are called, unless they get changed manually
#or the 'reset' option is enabled.
reset=0

arguments=()
while [[ $# -gt 0 ]]; do

  key="$1"

  case $key in
  -p | --project)
    project="$2"
    shift
    shift
    ;;
  -a | --assembly)
    assembly="$2"
    shift
    shift
    ;;
  -1 | --fastq-for)
    fastq_for="$2"
    shift
    shift
    ;;
  -2 | --fastq-rev)
    fastq_rev="$2"
    shift
    shift
    ;;
  -o | --output-dir)
    output_dir="$2"
    shift
    shift
    ;;
  -i | --iterations)
    iterations="$2"
    shift
    shift
    ;;
  --iter)
    iter="$2"
    shift
    shift
    ;;
  --working-dir)
    working_dir="$2"
    shift
    shift
    ;;
  --tools-dir)
    tools_dir="$2"
    shift
    shift
    ;;
  --index-name)
    index_name="$2"
    shift
    shift
    ;;
  --model-dir)
    model_dir="$2"
    shift
    shift
    ;;
  -A | --assembly-dir)
    assembly_dir="$2"
    shift
    shift
    ;;
  -k | --annotation-dir)
    annotation_dir="$2"
    shift
    shift
    ;;
  -M | --alignment-dir)
    alignment_dir="$2"
    shift
    shift
    ;;
  -T | --tmp-dir)
    tmp_dir="$2"
    shift
    shift
    ;;
  -N | --network-dir)
    network_dir="$2"
    shift
    shift
    ;;
  -n | --network-file)
    network_file="$2"
    shift
    shift
    ;;
  -P | --partition-dir)
    partition_dir="$2"
    shift
    shift
    ;;
  -R | --read-size)
    read_size="$2"
    shift
    shift
    ;;
  -s | --size-contig-threshold)
    size_contig_threshold="$2"
    shift
    shift
    ;;
  -C | --chunk-size)
    chunk_size="$2"
    shift
    shift
    ;;
  -c | --size-chunk-threshold)
    size_chunk_threshold="$2"
    shift
    shift
    ;;
  -Q | --mapping-quality-threshold)
    mapping_quality_threshold="$2"
    shift
    shift
    ;;
  -t | --threads)
    threads="$2"
    shift
    shift
    ;;
  -b | --n-bins)
    n_bins="$2"
    shift
    shift
    ;;
  -e | --e-value)
    evalue="$2"
    shift
    shift
    ;;
  -O | --overlapping-score)
    overlapping_score="$2"
    shift
    shift
    ;;
  -H | --hmm-databases)
    hmm_databases="$2"
    shift
    shift
    ;;
  --minimap2 | minimap)
    minimap2=1
    shift
    ;;
  --norm)
    norm=1
    shift
    shift
    ;;
  --clean-up)
    clean_up=1
    shift
    ;;
  --reset)
    reset=1
    shift
    ;;
  --compress)
    compress=1
    shift
    ;;
  *)
    arguments+=("$1")
    shift
    ;;
  esac
done

set -- "${arguments[@]}" #Positional parameters (i.e. actions)

current_dir="$(cd "$(dirname "$0")" && pwd)"

#This sets the appropriate script to launch with all these parameters, e.g. ./meta3c.sh action --blabla 4 --blibli 6 --etc
mode=$1

if [ ! -f config.sh ]; then
  echo "Config file not detected, generating one from template."
  cp "$current_dir"/config_template.sh config_current.sh
elif [ "$reset" -eq 1 ]; then
  echo "Resetting config file from template."
  cp "$current_dir"/config_template.sh config_current.sh
else
  cp config.sh config_current.sh
fi

#We write all these custom parameters into a config file that's called at the beginning of every script you run. This is because some parameters are dependent on others (e.g. if you specify output_dir you probably want network_dir, partition_dir etc. to all depend on output_dir unless explicited otherwise). Parameters that aren't specified are left alone and untouched from config_template.sh (which contains the defaults). It's a bit dirty but it works.

awk -v date="$(date +"%Y-%m-%d %r")" '

BEGIN {

    FS = "="

    parameters["project"] = "'"$project"'"
    parameters["assembly"] = "'"$assembly"'"
    parameters["fastq_for"] = "'"$fastq_for"'"
    parameters["fastq_rev"] = "'"$fastq_rev"'"
    parameters["working_dir"] = "'"$working_dir"'"
    parameters["tools_dir"] = "'"$tools_dir"'"
    parameters["index_name"] = "'"$index_name"'"
    parameters["model_dir"] = "'"$model_dir"'"
    parameters["output_dir"] = "'"$output_dir"'"
    parameters["assembly_dir"] = "'"$assembly_dir"'"
    parameters["alignment_dir"] = "'"$alignment_dir"'"
    parameters["tmp_dir"] = "'"$tmp_dir"'"
    parameters["network_dir"] = "'"$network_dir"'"
    parameters["network_file"] = "'"$network_file"'"
    parameters["partition_dir"] = "'"$partition_dir"'"
    parameters["annotation_dir"] = "'"$annotation_dir"'"
    parameters["read_size"] = "'"$read_size"'"
    parameters["size_contig_threshold"] = "'"$size_contig_threshold"'"
    parameters["chunk_size"] = "'"$chunk_size"'"
    parameters["size_chunk_threshold"] = "'"$size_chunk_threshold"'"
    parameters["mapping_quality_threshold"] = "'"$mapping_quality_threshold"'"
    parameters["evalue"] = "'"$evalue"'"
    parameters["overlapping_score"] = "'"$overlapping_score"'"
    parameters["hmm_databases"] = "'"$hmm_databases"'"
    parameters["norm"] = "'"$norm"'"
    parameters["n_bins"] = "'"$n_bins"'"
    parameters["iterations"] = "'"$iterations"'"
    parameters["iter"] = "'"$iter"'"
    parameters["clean_up"] = "'"$clean_up"'"
    parameters["compress"] = "'"$compress"'"
    parameters["threads"] = "'"$threads"'"
    parameters["minimap2"] = "'"$minimap2"'"

}

#Empty lines or comments

/^#Auto-generated/ {
    next
}

/^#/ || NF == 0 || /^export/ {
    print $0
    next
}

#Specified parameters
$1 in parameters && parameters[$1] != "" {
    print $1"="parameters[$1]
}

#Unspecified parameters
$1 in parameters && parameters[$1] == "" {
    print $0
}

END {

    print "#Auto-generated config file used to launch '"$mode"'.sh at " date

}

' config_current.sh | sed 's/\/$//g' >config.sh

. config.sh

#A mode (or action) is usually just a specific script to run with the above config file.
case $mode in
m | align | align.sh)
  # shellcheck source=alignment.sh
  . "$current_dir"/alignment.sh
  ;;
p | partition | partition.sh)
  # shellcheck source=partition.sh
  . "$current_dir"/partition.sh
  ;;
a | annotation | annotation.sh)
  # shellcheck source=annotation.sh
  . "$current_dir"/annotation.sh
  ;;
b | binning | binning.sh)
  # shellcheck source=binning.sh
  . "$current_dir"/binning.sh
  ;;
A | pipeline | pipeline.sh)
  # shellcheck source=pipeline.sh
  . "$current_dir"/pipeline.sh
  ;;
D | deploy | deploy.sh)
  # shellcheck source=deploy.sh
  . "$current_dir"/deploy.sh
  ;;
d | dependencies)
  fetch_dependencies
  ;;
h | help)
  display_help
  ;;
v | version)
  echo $current_version
  ;;
*)
  #Fallback if the user typed in something wrong. May add a custom error message later if need be.
  display_help
  ;;
esac
