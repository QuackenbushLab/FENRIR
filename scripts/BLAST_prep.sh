#!/bin/bash

#variables
File_Dir=/root/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/input

for f in $(find $File_Dir/ -name *HISAT2_virome_alignment_Cufflinks.sam); do
  FILENAME=$(echo $f | rev | cut -c 39- | rev | uniq)

  # sort and extract mapped reads from bam
  samtools sort -n ${FILENAME}_HISAT2_virome_alignment_Cufflinks_sorted_MrkDup_unique_noTanRep.bam | samtools view -b -F 4 -  > \
    ${FILENAME}_HISAT2_virome_alignment_Cufflinks_sorted_MrkDup_unique_noTanRep_mapped_sorted.bam

  # bam to fasta
  samtools fasta ${FILENAME}_HISAT2_virome_alignment_Cufflinks_sorted_MrkDup_unique_noTanRep_mapped_sorted.bam > \
    ${FILENAME}_HISAT2_virome_alignment_Cufflinks_sorted_MrkDup_unique_noTanRep_mapped.fa

  # rename fasta to txt for the sake of BLAST
  mv ${FILENAME}_HISAT2_virome_alignment_Cufflinks_sorted_MrkDup_unique_noTanRep_mapped.fa ${FILENAME}_HISAT2_virome_alignment_Cufflinks_sorted_MrkDup_unique_noTanRep_mapped.txt
done
