# Distributed under the terms of the Modified BSD License.
ARG OWNER=jupyter
ARG BASE_CONTAINER=$OWNER/minimal-notebook
FROM $BASE_CONTAINER as base

# Fix: https://github.com/hadolint/hadolint/wiki/DL4006
# Fix: https://github.com/koalaman/shellcheck/wiki/SC3014
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

USER root

RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends \
    # for cython: https://cython.readthedocs.io/en/latest/src/quickstart/install.html
    build-essential \
    # for latex labels
    cm-super \
    dvipng \
    # for matplotlib anim
    ffmpeg \
    time && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

RUN  apt update && \
     apt-get install -y --no-install-recommends \
     man-db \
     g++ \
     less \
     zlib1g-dev \
     && \
     apt-get clean && rm -rf /var/lib/apt/lists/*


# Install dependencies for bam2hints and filterBam
RUN apt update && \
    apt-get install -y libbamtools-dev augustus

# compile VARUS
RUN cd /opt && \
    git clone https://github.com/Gaius-Augustus/VARUS.git && \
    cd VARUS/Implementation && \
    make

USER root

ENV PATH=${PATH}:/opt/VARUS

RUN apt update && \
    apt install --yes --no-install-recommends samtools bamtools hisat2 parallel && \
    apt clean all

#RUN apt update && \
#    apt install -yq libyaml-perl \
#                    libhash-merge-perl \
#                    libparallel-forkmanager-perl \
#                    libscalar-util-numeric-perl \
#                    libclass-data-inheritable-perl \
#                    libexception-class-perl \
#                    libtest-pod-perl \
#                    libfile-which-perl \
#                    libmce-perl \
#                    libthread-queue-perl \
#                    libmath-utils-perl \
#                    libscalar-list-utils-perl && \
#    apt clean all

RUN cd /opt && \
    wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.1.0/sratoolkit.3.1.0-ubuntu64.tar.gz && \
    tar -xf sratoolkit.3.1.0-ubuntu64.tar.gz && \
    mv sratoolkit.3.1.0-ubuntu64 sratoolkit


RUN cd /opt && \
    mkdir datasets && \
    cd datasets && \
    wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets && \
    chmod a+x datasets



ENV PATH=${PATH}:/opt/datasets:/opt/VARUS/Implementation/bin:/opt/VARUS/Implementation/scripts:/opt/sratoolkit/bin:/opt/Augustus/bin/

USER ${NB_UID}
