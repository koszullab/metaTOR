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
  
  mkdir -p /tools
  cd /tools
  
  # Fetching prodigal
  wget -q https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux -O /tools/prodigal
  chmod +x /tools/prodigal
  
  # Fetching louvain
  mkdir -p /tools/louvain
  wget -q https://lip6.github.io/Louvain-BinaryBuild/louvain_linux.tar.gz -O /tools/louvain/louvain.tar.gz
  cd /tools/louvain
  tar -xzf louvain.tar.gz
  chmod +x /tools/louvain/*
  rm -f /tools/louvain/louvain.tar.gz
  
  # Fetching HMMs
  mkdir -p /HMM_databases
  cd /HMM_databases
  # Original database link was dead let's use something else for now
  wget -q https://raw.githubusercontent.com/MadsAlbertsen/multi-metagenome/master/R.data.generation/essential.hmm
  cd /
  
  # Install metator and requests
  pip3 install requests metator

  # Clean up
  rm -rf /var/lib/apt/lists/*

%runscript
  exec metator "$@"
