#!/bin/bash

## Purge existing artifacts (required for local rebuild)
rm -rf artifacts/ gen-louvain/ bowtie2/ networkanalysis/

## Install louvain
tar -k -xzf louvain-generic.tar.gz
cd gen-louvain
make
cd ..

## Install leiden (https://github.com/CWTSLeiden/networkanalysis)
# wget "https://github-registry-files.githubusercontent.com/153760626/0f40f180-3ed3-11ee-916e-23eb9928c186?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAVCODYLSA53PQK4ZA%2F20250218%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20250218T134441Z&X-Amz-Expires=300&X-Amz-Signature=065d9dec2e31375461b33db29c5e6261b54eca2ec2e44c6460a3138586f86b89&X-Amz-SignedHeaders=host&response-content-disposition=filename%3Dnetworkanalysis-1.3.0.jar&response-content-type=application%2Foctet-stream" -O networkanalysis-1.3.0.jar
mkdir -p networkanalysis/build/libs/
cp networkanalysis-1.3.0.jar networkanalysis/build/libs/

# ## Install bowtie2
# wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.5.1/bowtie2-2.5.1-source.zip/download -O bowtie2-2.5.1-source.zip
# zip -d bowtie2-2.5.1-source.zip "bowtie2-2.5.1/example/*"
# unzip bowtie2-2.5.1-source.zip
# mv bowtie2-2.5.1/ bowtie2
# cd bowtie2
# make
# cd ..

## Move artifacts to the correct location
mkdir -p artifacts/networkanalysis/build artifacts/bowtie2/bin
mv gen-louvain/ artifacts/
mv networkanalysis/build artifacts/networkanalysis/
# mv bowtie2/bowtie2* artifacts/bowtie2/bin/

rm -rf gen-louvain/ bowtie2/ networkanalysis/
