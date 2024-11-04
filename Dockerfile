# syntax=docker/dockerfile:1

FROM mambaorg/micromamba:latest

LABEL Name=metator Version=1.3.3

COPY --chown=$MAMBA_USER:$MAMBA_USER . ./

# Install 3rd party packages
USER root
RUN apt update && \
    apt install -y --no-install-recommends git make g++ curl default-jre default-jdk zlib1g-dev

# Install Louvain 
RUN cd ./external && \
    tar -xzf louvain-generic.tar.gz && \
    cd gen-louvain && \
    make && \
    cd ../
ENV LOUVAIN_PATH=./external/gen-louvain

# Install Leiden through Network analysis repo
RUN git clone https://github.com/vtraag/networkanalysis.git && \
    cd ./networkanalysis && \
    ./gradlew build && \
    cd ../
ENV LEIDEN_PATH=$(pwd)/networkanalysis/build/libs/networkanalysis-1.2.0.jar

## Install dependencies
USER mambauser
RUN micromamba install -y -n base --file metator.yaml && \
    micromamba install -y -n base pip && \
    micromamba clean --all --yes

# Install metator
RUN micromamba run -n base python3 -m pip install -e .

WORKDIR /home/mambauser/
ENTRYPOINT [ "/bin/bash" ]
