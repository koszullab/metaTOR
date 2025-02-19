FROM mambaorg/micromamba:latest

LABEL Name=metator Version=1.3.3

COPY --chown=$MAMBA_USER:$MAMBA_USER . ./

# Install 3rd party packages
USER root
RUN apt update && \
    apt install -y --no-install-recommends git make g++ curl default-jre default-jdk zlib1g-dev unzip

# ## Install dependencies
USER mambauser
RUN micromamba install -y -n base --file metator.yaml && \
    micromamba clean --all --yes

# Install metator
RUN micromamba run -n base pip install -e .[dev]

WORKDIR /home/mambauser/
ENTRYPOINT [ "/bin/bash" ]
