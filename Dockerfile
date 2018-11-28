FROM ubuntu:18.04

RUN \
    apt-get update && \
    apt-get install -y --no-install-recommends git python python-pip python3-setuptools python3 python3-dev python3-pip python3-virtualenv bowtie2 samtools hmmer prodigal libfreetype6-dev libpng-dev pkg-config wget && \
    rm -rf /var/lib/apt/lists

COPY metator /app/metator

COPY *.* /app/

WORKDIR /app

RUN cd /app && \
    pip3 install .

USER root

RUN chmod 777 -R /usr/local/lib/python3.6/dist-packages/metator/bin/

RUN metator dependencies

ENTRYPOINT ["metator"]