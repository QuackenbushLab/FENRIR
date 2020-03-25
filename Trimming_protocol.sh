#!/bin/bash

# grab out filename from list and trim

Trimming_Path=/Trimmomatic-0.36/adapters
FASTQC_Path=/FastQC
Input=/root/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/input
for f in $(find $Input -name "*_1.fastq"); do
  FILENAME=$(echo $f | rev | cut -c 9- | rev | uniq)
  $FASTQC_Path/fastqc ${FILENAME}_1.fastq ${FILENAME}_2.fastq 
  java -jar /Trimmomatic-0.36/trimmomatic-0.36.jar  PE \
  -phred33 ${FILENAME}_1.fastq ${FILENAME}_2.fastq ${FILENAME}_forward_paired.fq ${FILENAME}_forward_unpaired.fq ${FILENAME}_reverse_paired.fq ${FILENAME}_reverse_unpaired.fq \
  ILLUMINACLIP:$Trimming_Path/merged.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:75
  $FASTQC_Path/fastqc ${FILENAME}_forward_paired.fq ${FILENAME}_forward_unpaired.fq ${FILENAME}_reverse_paired.fq ${FILENAME}_reverse_unpaired.fq
done
