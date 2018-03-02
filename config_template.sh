#!/usr/bin/env bash
# Meta3C team
# Martial Marbouty
# Théo Foutel-Rodier
# Lyam Baudry

#This file contains information related to file paths and parameters for various parts of the pipeline.
#Don't modfiy this template unless you really know what you're doing. 
#However you may do so if you wish to alter the defaults for a Meta3C experiment on a standard machine.

#project name: all outputs will be contained in a global folder with that name
project=my_project



#### Input files

#Initial files
fastq_for=""
fastq_rev=""

#Initial preliminary assembly
assembly=${assembly_dir}/${project}.fa

#### Directories

#Base folder for everything. 
working_dir="$( cd "$( dirname "$0" )" && pwd )"

## Third-party folders

#Third-party tools folder. Normally contains louvain, pplacer and skewer. The pipeline will also look there for any tools that it can't find by other methods.
tools_dir=${working_dir}/tools

#Folder containing model files used for gene prediction
model_dir=${working_dir}/modele_HMM

#Folder containing general outputs. By default all outputs will be contained in that folder.
output_dir=${working_dir}/output

#Folder containing outputs related to the alignment 
alignment_dir=${output_dir}/${project}/alignement

#Temporary folder containing all intermediary files (labeled after $project so they don't erase each other)
tmp_dir=${output_dir}/${project}/temp

#Folder containing outputs related to the 3C contact network
network_dir=${output_dir}/${project}/network

#Folder containing outputs related to the partition algorithm's output
partition_dir=${output_dir}/${project}/partition

#Folder containing outputs related to the assembly's annotation
annotation_dir=${output_dir}/${project}/annotation

#Folder containing assemblies
assembly_dir=${output_dir}/${project}/assembly




#### Meta3C parameters

## Read processing

#Read size after trimming and filtering (both paired ends combined)
read_size=130


## Assembly processing parameters

#Contigs below this threshold are discarded. This should not be normally smaller than $size_chunk_threshold.
size_contig_threshold=500


## Network generation parameters

#Default network name when partitioning
network_file=${network_dir}/network.txt

#Contigs are virtually 'split' into fixed-length chunks: each chunk gets a separate node and label in the contact network. 
chunk_size=1000

#Contig tails below this threshold are discarded, while those above this threshold are treated as separate chunks.
size_chunk_threshold=500

#Alignments below this mapping quality score are discarded.
mapping_quality_threshold=10

#Normalize contacts by their respective coverages
norm=0


## Louvain parameters

#Number of iterations of the Louvain algorithm to perform. This is one of the time and memory critical steps of the pipeline, so don't set this number higher than necessary.
iterations=300

#Number of iterations to consider if for some reason you wish it to be lower than the total number of iterations
iter=${iterations}

#Number of bins to extract in size order (starting with the biggest)
n_bins=100

#Score threshold on the Louvain matrix to consider communities to be overlapping
overlapping_score=1


## Annotation parameters

#Minimum E-value for discarding matches of predicted proteins onto HMMs
evalue=1e-4



#### Miscellaneous parameters

#Remove temporary files after each pipeline step
clean_up=0

#Number of threads to use
#threads=$(grep -c ^processor /proc/cpuinfo) #Only works on Linux
threads=1

#Use minimap2 instead of bowtie2
minimap2=0

#Basename for bowtie2 index structures
index_name=${working_dir}/$(basename ${assembly%.*})
