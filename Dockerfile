# Dockerfile for FENRIR RNA-seq pipeline
FROM ubuntu:16.04
MAINTAINER Farrah Roy

RUN apt-get update && apt-get install -y software-properties-common && add-apt-repository -y ppa:openjdk-r/ppa && \
    apt-get update && apt-get install -y \
        build-essential \
        cmake \
        curl \
        libboost-all-dev \
        libbz2-dev \
        libcurl3-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        openjdk-8-jdk \
        python3 \
        python3-pip \
        unzip \
        tar \
        vim-common \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*


#-----------------------------
# Pipeline components
#-----------------------------

#SRA Tools
RUN wget "http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz" && \
    tar -xvf sratoolkit.current-centos_linux64.tar.gz && rm sratoolkit.current-centos_linux64.tar.gz 

# samtools
RUN  wget  https://github.com/samtools/samtools/releases/download/1.5/samtools-1.5.tar.bz2 && \
    tar -xf samtools-1.5.tar.bz2 && rm samtools-1.5.tar.bz2 && cd samtools-1.5 && \
    make && make install && make clean

#bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz && \
    tar -xf bedtools-2.26.0.tar.gz && rm bedtools-2.26.0.tar.gz &&\
    cd bedtools2 && \
    make

# HISAT2
RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip && \
    unzip hisat2-2.1.0-Linux_x86_64.zip
RUN cp -p hisat2-2.1.0/hisat2 hisat2-2.1.0/hisat2-* /usr/bin
ENV PATH /usr/bin/hisat2-2.1.0/hisat2 hisat2-2.1.0/hisat2-*:$PATH

# python modules
RUN pip3 install --upgrade pip
RUN pip install HTSeq

# Install Python 3.6.5
ENV PYTHON_VERSION="3.6.5"
RUN wget https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tar.xz \
    && tar xvf Python-${PYTHON_VERSION}.tar.xz \
    && rm Python-${PYTHON_VERSION}.tar.xz \
    && cd Python-${PYTHON_VERSION} \
    && ./configure \
    && make altinstall \
    && cd / \
    && rm -rf Python-${PYTHON_VERSION}

#Trimmomatic
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip && \
 unzip Trimmomatic-0.36.zip && \
 cp -r Trimmomatic-0.36 /usr/bin 

#install FastQC 0.11.5
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && unzip fastqc_v0.11.5.zip && chmod +x FastQC/fastqc 

# Install OpenJDK-8
RUN apt-get update && \
    apt-get install -y openjdk-8-jdk && \
    apt-get install -y ant && \
    apt-get clean;

# Fix certificate issues for OpenJDK
RUN apt-get update && \
    apt-get install ca-certificates-java && \
    apt-get clean && \
    update-ca-certificates -f;

# Setup JAVA_HOME -- useful for docker commandline
ENV JAVA_HOME /usr/bin/java-8-openjdk-amd64/:$PATH
RUN export JAVA_HOME

#Picard 
RUN wget -r https://github.com/broadinstitute/picard/releases/download/2.18.16/picard.jar

#featureCounts
RUN wget https://excellmedia.dl.sourceforge.net/project/subread/subread-1.6.3/subread-1.6.3-Linux-x86_64.tar.gz && \
   tar -xvf subread-1.6.3-Linux-x86_64.tar.gz && rm subread-1.6.3-Linux-x86_64.tar.gz
   
# clean up
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

#set programs as executable    
ADD ./scripts/Trimming_protocol.sh /home/Trimming_protocol.sh
ADD ./scripts/pseudovirome_alignment_second_pass.sh /home/pseudovirome_alignment_second_pass.sh
ADD ./scripts/PostAlignment_hg38ToPseudo.sh /home/PostAlignment_hg38ToPseudo.sh
ADD ./scripts/Index_Generation.sh /home/Index_Generation.sh
ADD ./scripts/PostAlignment_ERV.sh /home/PostAlignment_ERV.sh
ADD ./scripts/BLAST_prep.sh /home/BLAST_prep.sh

RUN chmod +x /home/Trimming_protocol.sh /home/PostAlignment_hg38ToPseudo.sh \
  /home/pseudovirome_alignment_second_pass.sh  /home/Index_Generation.sh /home/PostAlignment_ERV.sh /home/BLAST_prep.sh 
RUN rm *.zip
