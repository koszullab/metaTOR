FROM ubuntu:18.04

RUN \
    apt-get update && \
    apt-get install -y --no-install-recommends python python-pip python-setuptools python3 python3-dev python3-pip python3-virtualenv bowtie2 samtools hmmer prodigal libfreetype6-dev libpng-dev pkg-config wget && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /data

ADD *.* /data/

RUN pip3 install metator

RUN metator dependencies

ENTRYPOINT ["metator"]