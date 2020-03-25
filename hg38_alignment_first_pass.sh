#!/bin/bash

#variables
File_Dir=/Users/far122/Desktop/Docker/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/input
Ref_Dir=/Users/far122/Desktop/Docker/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/ref

for f in $(find $File_Dir -name "*forward_paired.fq"); do
  FILENAME=$(echo $f |  rev | cut -c 19- | rev | uniq)
  echo ${FILENAME}
  hisat2 --phred33 --known-splicesite-infile ${Ref_Dir}/hg38_ss.txt --dta-cufflinks \
    -x ${Ref_Dir}/hg38_genome_index -1 ${FILENAME}_forward_paired.fq -2 ${FILENAME}_reverse_paired.fq \
    -U ${FILENAME}_forward_unpaired.fq,${FILENAME}_reverse_unpaired.fq -S ${FILENAME}_HISAT2_hg38_HML2_alignment.sam --un ${FILENAME}_HISAT2_hg38_HML2_alignment_un-seqs \
    --un-conc ${FILENAME}_HISAT2_hg38_HML2_alignment_un-conc-mate
done


