#!/bin/bash
source Variables.txt
echo $READING
echo


echo 
echo IMPORTING LIST OF EXOGNEOUS VIRUSES AND ENDOGENOUS VIRUSES
echo
Exo_input=${File_Dir}/"list_of_viruses.txt"
ERV_input=${File_Dir}/"HERV_list.txt"
echo
echo DONE IMPORTING
echo

echo
echo EXTRACTING COORDINATES FOR ENDOGENOUS VIRUSES
echo
ERV_coordinates=$(grep -o ',"[^"]\+"' ${ERV_input} | sed 's/"//g' | sed 's/,//g' | awk '{print $2}' -)
echo
echo DONE EXTRACTING COORDINATES
echo

echo
echo BLAST PSEUDOVIROME RESULTS
echo

#extract reads for exogenous viruses
for f in $(find ${File_Dir} -name "*_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep.bam"); do
  FILENAME=$(echo $f | rev | cut -c 59- | rev | uniq)
  
  echo
  echo SORTING AND EXTRACTING MAPPED READS FROM ${FILENAME}
  echo

  # sort and extract mapped reads from bam
  samtools sort -n ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep.bam | samtools view -b -f 3 -  > \
    ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep_mapped_sorted_flag3.bam

  samtools sort -n ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep.bam | samtools view -b -f 73 -  > \
    ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep_mapped_sorted_flag73.bam

  samtools sort -n ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep.bam | samtools view -b -f 89 -  > \
    ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep_mapped_sorted_flag89.bam

  samtools sort -n ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep.bam | samtools view -b -f 153 -  > \
    ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep_mapped_sorted_flag153.bam

  samtools sort -n ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep.bam | samtools view -b -f 137 -  > \
    ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep_mapped_sorted_flag137.bam
  
  echo
  echo MERGING BAM FILES FROM ${FILENAME}
  echo
  
  # merge bam files
  samtools merge ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep_mapped_sorted_flag_merged.bam  \
  ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep_mapped_sorted_flag3.bam \
  ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep_mapped_sorted_flag73.bam \
  ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep_mapped_sorted_flag89.bam \
  ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep_mapped_sorted_flag153.bam \
  ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep_mapped_sorted_flag137.bam
  
  samtools sort ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep_mapped_sorted_flag_merged.bam -o \
     ${FILENAME}_HISAT2_virome_alignment_BLAST_input.bam
  samtools index ${FILENAME}_HISAT2_virome_alignment_BLAST_input.bam

  while read -r line
  do
    echo
    echo ISOLATING READS FOR EXOGENOUS VIRUS "$line"
    echo
    
    samtools view -h ${FILENAME}_HISAT2_virome_alignment_BLAST_input.bam "${line}" > ${FILENAME}_HISAT2_virome_alignment_${line}_BLAST_input.bam
  
    # bam to fasta
    samtools fasta ${FILENAME}_HISAT2_virome_alignment_${line}_BLAST_input.bam > ${FILENAME}_HISAT2_virome_alignment_${line}_BLAST_input.fa
    echo
    echo DONE EXTRACTING READS FROM ${FILENAME} FOR EXOGENOUS VIRUS ${line}
    echo
  done < "$Exo_input"
  
done

#extract reads for ERVs
for f in $(find ${File_Dir} -name "*_HISAT2_reference_alignment_sorted_MrkDup_unique.bam"); do
  FILENAME=$(echo $f | rev | cut -c 53- | rev | uniq)
  
  echo
  echo SORTING AND EXTRACTING MAPPED READS FROM ${FILENAME}
  echo

  # sort and extract mapped reads from bam
  samtools sort -n ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique.bam  | samtools view -b -f 3 -  > \
    ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique_mapped_sorted_flag3.bam

  samtools sort -n ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique.bam | samtools view -b -f 73 -  > \
    ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique_mapped_sorted_flag73.bam

  samtools sort -n ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique.bam | samtools view -b -f 89 -  > \
    ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique_mapped_sorted_flag89.bam

  samtools sort -n ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique.bam | samtools view -b -f 153 -  > \
    ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique_mapped_sorted_flag153.bam

  samtools sort -n ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique.bam | samtools view -b -f 137 -  > \
    ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique_mapped_sorted_flag137.bam
  
  echo
  echo MERGING BAM FILES FROM ${FILENAME}
  echo
  
  # merge bam files
  samtools merge ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique_mapped_sorted_flag_merged.bam  \
  ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique_mapped_sorted_flag3.bam \
  ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique_mapped_sorted_flag73.bam \
  ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique_mapped_sorted_flag89.bam \
  ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique_mapped_sorted_flag153.bam \
  ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique_mapped_sorted_flag137.bam
  
  samtools sort ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique_mapped_sorted_flag_merged.bam -o \
     ${FILENAME}_HISAT2_reference_alignment_BLAST_input.bam
  samtools index ${FILENAME}_HISAT2_reference_alignment_BLAST_input.bam

  while IFS= read -r line
  do
    echo
    echo ISOLATING READS FOR ENDOGENOUS VIRUS "$line"
    echo
    
    samtools view -h ${FILENAME}_HISAT2_reference_alignment_BLAST_input.bam "${line}" > ${FILENAME}_HISAT2_reference_alignment_${line}_BLAST_input.bam
  
    # bam to fasta
    samtools fasta ${FILENAME}_HISAT2_reference_alignment_${line}_BLAST_input.bam > ${FILENAME}_HISAT2_reference_alignment_${line}_BLAST_input.fa
    echo
    echo DONE EXTRACTING READS FROM ${FILENAME} FOR ENDOGENOUS VIRUS ${line}
    echo
  done < <(printf '%s\n' "$ERV_coordinates")
  
done

#move all files to their final directories
mkdir ${Working_Dir}/BLAST_results
mkdir ${Working_Dir}/JPEG_Results
mkdir ${Working_Dir}/htseq_counts_data
mkdir ${Working_Dir}/individual_counts_files
mkdir ${Working_Dir}/coverage_data
mkdir ${Working_Dir}/Picard_metrics

mv ${File_Dir}/*BLAST_input.fa ${Working_Dir}/BLAST_results
mv ${File_Dir}/*HISAT2_reference_alignment* ${Working_Dir}/Reference_Alignment
mv ${File_Dir}/*HISAT2_virome_alignment* ${Working_Dir}/Pseudovirome_Alignment
mv ${File_Dir}/*jpg ${Working_Dir}/JPEG_Results
mv ${File_Dir}/*counts.txt ${Working_Dir}/individual_counts_files
mv ${File_Dir}/*coverage_merged.csv ${Working_Dir}/coverage_data
mv ${File_Dir}/*csv ${Working_Dir}/htseq_counts_data
mv${File_Dir}/*metrics.txt ${Working_Dir}/Picard_metrics

echo
echo READY FOR BLAST
echo
