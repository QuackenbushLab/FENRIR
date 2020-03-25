#!/bin/bash

#variables
File_Dir=/root/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/input
Ref_Dir=/root/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/ref
Picard_Path=/github.com/broadinstitute/picard/releases/download/2.18.16
Bedtools_Path=/bedtools2/bin
FeatureCounts_Path=/subread-1.6.3-Linux-x86_64/bin

for f in $(find $File_Dir/ -name *HISAT2_virome_alignment_Cufflinks.sam); do
  FILENAME=$(echo $f | rev | cut -c 39- | rev | uniq)
  #convert .sam to .bam and sort
  samtools view -Su ${FILENAME}_HISAT2_virome_alignment_Cufflinks.sam | samtools sort - > ${FILENAME}_HISAT2_virome_alignment_Cufflinks_sorted.bam

  #MarkDuplicates with Picard

  java -jar $Picard_Path/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=${FILENAME}_HISAT2_virome_alignment_Cufflinks_sorted.bam \
  O=${FILENAME}_HISAT2_virome_alignment_Cufflinks_sorted_MrkDup.bam M=${FILENAME}_marked_dup_metrics.txt

  #capture unique reads
  samtools view -bq 50 ${FILENAME}_HISAT2_virome_alignment_Cufflinks_sorted_MrkDup.bam > ${FILENAME}_HISAT2_virome_alignment_Cufflinks_sorted_MrkDup_unique.bam

  #remove tandem repeats
  ${Bedtools_Path}/bedtools intersect -abam ${FILENAME}_HISAT2_virome_alignment_Cufflinks_sorted_MrkDup_unique.bam -b $Ref_Dir/Pseudovirome_2.7.7.80.10.50.2000.bed \
    -v > ${FILENAME}_HISAT2_virome_alignment_Cufflinks_sorted_MrkDup_unique_noTanRep.bam

  #start calculating read depth
  samtools depth -a ${FILENAME}_HISAT2_virome_alignment_Cufflinks_sorted_MrkDup_unique_noTanRep.bam > \
    ${FILENAME}_HISAT2_virome_alignment_Cufflinks_sorted_MrkDup_unique_noTanRep.coverage

  #index
  samtools index ${FILENAME}_HISAT2_virome_alignment_Cufflinks_sorted_MrkDup_unique_noTanRep.bam
  samtools index ${FILENAME}_HISAT2_virome_alignment_Cufflinks_sorted_MrkDup_unique.bam
  samtools index ${FILENAME}_HISAT2_virome_alignment_Cufflinks_sorted_MrkDup.bam
  samtools index ${FILENAME}_HISAT2_virome_alignment_Cufflinks_sorted.bam

  #count reads
  ${FeatureCounts_Path}/featureCounts -p -O -M -t exon -g transcript_id -a ${Ref_Dir}/Pseudovirome.gtf -o \
     ${FILENAME}_featureCounts_virome_sorted_MrkDup_unique_noTanRep.txt ${FILENAME}_HISAT2_virome_alignment_Cufflinks_sorted_MrkDup_unique_noTanRep.bam
  ${FeatureCounts_Path}/featureCounts -p -O -M -t exon -g transcript_id -a ${Ref_Dir}/Pseudovirome.gtf -o \
     ${FILENAME}_featureCounts_virome_sorted_MrkDup_unique.txt ${FILENAME}_HISAT2_virome_alignment_Cufflinks_sorted_MrkDup_unique.bam
  ${FeatureCounts_Path}/featureCounts -p -O -M -t exon -g transcript_id -a ${Ref_Dir}/Pseudovirome.gtf -o \
     ${FILENAME}_featureCounts_virome_sorted_MrkDup.txt ${FILENAME}_HISAT2_virome_alignment_Cufflinks_sorted_MrkDup.bam
  ${FeatureCounts_Path}/featureCounts -p -O -M -t exon -g transcript_id -a ${Ref_Dir}/Pseudovirome.gtf -o \
     ${FILENAME}_featureCounts_virome_sorted.txt ${FILENAME}_HISAT2_virome_alignment_Cufflinks_sorted.bam
done
