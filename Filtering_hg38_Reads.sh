#!/bin/bash
#variables
File_Dir=root/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/input
Picard_Path=/github.com/broadinstitute/picard/releases/download/2.18.16

for f in $(find $File_Dir -name *HISAT2_hg38_HML2_alignment_second_pass.sam); do
  FILENAME=$(echo $f | cut -d'_' -f1,2)
  #get read names from sam file
  cut -f1 ${FILENAME}_HISAT2_hg38_HML2_alignment_second_pass.sam > ${FILENAME}_HISAT2_hg38_HML2_alignment_second_pass_readnames.txt
  #remove reads in text file from filtered virome bam file
  java -jar $Picard_Path/picard.jar FilterSamReads I=${FILENAME}_HISAT2_virome_alignment_first_pass_Cufflinks_sorted_MrkDup_unique_noTanRep.bam \
    O=${FILENAME}_HISAT2_virome_alignment_first_pass_Cufflinks_sorted_MrkDup_unique_noTanRep_hg38readsremoved.bam \
    READ_LIST_FILE=${FILENAME}_HISAT2_hg38_HML2_alignment_second_pass_readnames.txt FILTER=excludeReadList 
  #index
  samtools index ${FILENAME}_HISAT2_virome_alignment_first_pass_Cufflinks_sorted_MrkDup_unique_noTanRep_hg38readsremoved.bam
  #count reads
  samtools idxstats ${FILENAME}_HISAT2_virome_alignment_first_pass_Cufflinks_sorted_MrkDup_unique_noTanRep_hg38readsremoved.bam > \
  ${FILENAME}virome_counts_MrkDup_unique_noTanRep_hg38readsremoved.txt
done
