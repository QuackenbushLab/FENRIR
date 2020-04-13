#Description: Downloads necessary programs for the FENRIR pipeline
#Input: None

#Output: Downloaded software for FENRIR

#Author:
  #Farrah Roy

#!/bin/bash
#programs needed for FENRIR

sudo mv /var/lib/apt/lists
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install -y software-properties-common
sudo add-apt-repository -y ppa:openjdk-r/ppa
sudo apt-get install -y build-essential  cmake curl libboost-all-dev libbz2-dev libcurl3-dev liblzma-dev libncurses5-dev \
   libssl-dev openjdk-8-jdk python3 python3-pip unzip tar vim-common wget zlib1g-dev samtools bamtools trimmomatic \
   fastqc r-base git bedtools sra-toolkit python2.7-dev python-numpy python-matplotlib python-pysam python-pip python-htseq

#install HISAT2
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip
unzip hisat2-2.1.0-Linux_x86_64.zip
rm hisat2-2.1.0-Linux_x86_64.zip
sudo cp -p hisat2-2.1.0/hisat2 hisat2-2.1.0/hisat2-* /usr/bin

#install Trimmomatic
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
unzip Trimmomatic-0.36.zip 
rm Trimmomatic-0.36.zip
sudo cp -r Trimmomatic-0.36 /usr/bin 

#install gffread
git clone https://github.com/gpertea/gclib
  git clone https://github.com/gpertea/gffread
  cd gffread
  make release

#install htseq
pip install 'matplotlib>=1.4'
pip install Cython
pip install 'pysam>=0.9'
pip install 'HTSeq==0.11.0'

#install fastq pair
git clone https://github.com/linsalrob/fastq-pair.git
cd fastq-pair/
mkdir build && cd build
cmake ..
make && sudo make install
sudo cp -r fastq_pair /usr/bin 

#install Picard
wget -r https://github.com/broadinstitute/picard/releases/download/2.18.16/picard.jar
sudo cp -r ./github.com/broadinstitute/picard/releases/download/2.18.16/picard.jar /usr/bin

#install RStudio
sudo apt-get install gdebi-core
wget https://download2.rstudio.org/server/trusty/amd64/rstudio-server-1.2.5033-amd64.deb
sudo gdebi rstudio-server-1.2.5033-amd64.deb
