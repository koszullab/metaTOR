Bootstrap: docker
From: ubuntu:18.04

%labels
  Maintainer tpall

%post
  # Get dependencies
  apt-get update
  apt-get install -y --no-install-recommends \
  git \
  python \
  python-pip \
  python3-setuptools \
  python3 \
  python3-dev \
  python3-pip \
  python3-virtualenv \
  bowtie2 \
  samtools \
  hmmer \
  prodigal \
  libfreetype6-dev \
  libpng-dev \
  pkg-config \
  wget
  
  # Install
  pip3 install metator

  # Install louvain
  metator dependencies

  # Clean up
  rm -rf /var/lib/apt/lists/*

%runscript
  exec metator "$@"
