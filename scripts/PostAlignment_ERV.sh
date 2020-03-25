#!/bin/bash

#variables
File_Dir=/root/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/input
Picard_Path=/github.com/broadinstitute/picard/releases/download/2.18.16
Ref_Dir=/root/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/ref
Bedtools_Path=/bedtools2/bin
FeatureCounts_Path=/subread-1.6.3-Linux-x86_64/bin

for f in $(find $File_Dir/ -name *_HISAT2_hg38_HML2_alignment.sam); do
  FILENAME=$(echo $f | rev | cut -c 32- | rev | uniq)

  #convert .sam to .bam and sort
  samtools view -Su ${FILENAME}_HISAT2_hg38_HML2_alignment.sam | samtools sort - > ${FILENAME}_HISAT2_hg38_HML2_alignment_sorted.bam

  #MarkDuplicates with Picard
  java -jar $Picard_Path/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=${FILENAME}_HISAT2_hg38_HML2_alignment_sorted.bam \
    O=${FILENAME}_HISAT2_hg38_HML2_alignment_sorted_MrkDup.bam M=${FILENAME}_marked_dup_hg38_metrics.txt 

  #capture unique reads 
  samtools view -bq 50 ${FILENAME}_HISAT2_hg38_HML2_alignment_sorted_MrkDup.bam > ${FILENAME}_HISAT2_hg38_HML2_alignment_sorted_MrkDup_unique.bam

  #index
  samtools index ${FILENAME}_HISAT2_hg38_HML2_alignment_sorted_MrkDup_unique.bam
  samtools index ${FILENAME}_HISAT2_hg38_HML2_alignment_sorted_MrkDup.bam
  samtools index ${FILENAME}_HISAT2_hg38_HML2_alignment_sorted.bam
  
  #start calculating read depth
  samtools depth -a ${FILENAME}_HISAT2_hg38_HML2_alignment_sorted_MrkDup_unique.bam > \
    ${FILENAME}_HISAT2_hg38_HML2_alignment_sorted_MrkDup_unique.coverage

  #featureCounts
  ${FeatureCounts_Path}/featureCounts -p -O -M -t exon -g gene_id -a ${Ref_Dir}/hg38.gtf -o \
     ${FILENAME}_featureCounts_hg38_HML2_sorted_MrkDup_unique.txt ${FILENAME}_HISAT2_hg38_HML2_alignment_sorted_MrkDup_unique.bam
  ${FeatureCounts_Path}/featureCounts -p -O -M -t exon -g gene_id -a ${Ref_Dir}/hg38.gtf -o \
     ${FILENAME}_featureCounts_hg38_HML2_sorted_MrkDup.txt ${FILENAME}_HISAT2_hg38_HML2_alignment_sorted_MrkDup.bam
  ${FeatureCounts_Path}/featureCounts -p -O -M -t exon -g gene_id -a ${Ref_Dir}/hg38.gtf -o \
     ${FILENAME}_featureCounts_hg38_HML2_sorted.txt ${FILENAME}_HISAT2_hg38_HML2_alignment_sorted.bam
done
