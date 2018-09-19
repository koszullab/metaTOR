FROM ubuntu:18.04

RUN \
    apt-get update && \
    apt-get install -y python3 python3-dev python3-pip python3-virtualenv bowtie2 samtools hmmer prodigal libfreetype6-dev libpng-dev pkg-config wget && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /data

ADD *.* /data/

RUN pip3 install -Ur requirements.txt

RUN ./meta3c.sh dependencies

ENTRYPOINT ["./meta3c.sh"]