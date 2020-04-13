#Description: This script is responsible for trimming, aligning, and counting aligned reads to the reference genome (hg38)
  #and to our pseudovirome.
    #1) Reads Variables.txt
    #2) Check for reference and pseudovirome .fasta and .gtf
    #3) Checks for the presence of reference and pseudovirome index files for HISAT2. If missing, it builds index files
    #4) FastQC on input fastq files
    #5) Trims fastq files
    #6) FastQC on trimmed fastq files
    #7) Aligns trimmed fastq files to reference genome (ex. hg38)
    #8) Aligns unmapped reads to the Pseudovirome
    #9) Convert .sam to .bam and filter reference and pseudovirome alignments
    #10) Calculate counts for reference and pseudovirome alignments
    #11) Move original .fastq, trimmed .fastq, FastQC output, and alignment files into their own respective directories

#Input:
  #1) Fastq files
  #2) Reference genome and pseudovirome
  #3) Depending on size of the reference genome, index files for HISAT2
  #4) .bed file containing tandem repeats in pseudovirome
  #5) Variables.txt - filled out to match your experimental design

#Output:
  #1) Counts matricies of endogenous and exogenous viruses for individual samples and individual SAM flags
  #2) Samtools coverage files for each sample to calculate read depth for endogenous and exogenous viruses

#Author:
  #Farrah Roy

#!/bin/bash
source Variables.txt
echo $READING
echo

#generate index for HISAT2 aligner
if [ ! -s $Ref_Dir/${Pseudovirome} ];
then
  echo
  echo MISSING OR EMPTY PSEUDOVIROME FASTA. CHECK THAT FILE EXISTS, OR THAT THE NAME IS CORRECT
  echo
  exit
else
  echo FOUND PSEUDOVIROME FASTA, CONTINUE
  echo
fi

if [ ! -s $Ref_Dir/${Pseudovirome_GTF} ];
then
  echo
  echo MISSING OR EMPTY PSEUDOVIROME GTF. CHECK THAT FILE EXISTS, THAT IT IS IN THE RIGHT PLACE, OR THAT THE NAME IS CORRECT
  echo
  exit
else
  echo FOUND PSEUDOVIROME GTF, CONTINUE
  echo
fi

if [ ! -s $Ref_Dir/${Ref_Anno} ];
then
  echo
  echo MISSING OR EMPTY REFERENCE FASTA. CHECK THAT FILE EXISTS, THAT IT IS IN THE RIGHT PLACE, OR THAT THE NAME IS CORRECT
  echo
  exit
else
  echo FOUND REFERENCE FASTA, CONTINUE
  echo
fi

if [ ! -s $Ref_Dir/${Ref_GTF} ];
then
  echo
  echo MISSING OR EMPTY REFERENCE GTF. CHECK THAT FILE EXISTS, THAT IT IS IN THE RIGHT PLACE, OR THAT THE NAME IS CORRECT
  echo
  exit
else
  echo FOUND REFERENCE GTF, CONTINUE
  echo
fi

if [[ -s $Ref_Dir/${Ref_Name}_ss.txt && $Ref_Dir/${Ref_Name}_exons.txt ]];
then
  echo
  echo REFERENCE EXON AND SS FILES EXIST, MOVING ALONG
  echo
else
  echo
  echo GENERATING REFERENCE EXON AND SS FILES
  echo
  hisat2-2.1.0/hisat2_extract_splice_sites.py $Ref_Dir/${Ref_GTF} > $Ref_Dir/${Ref_Name}_ss.txt
  hisat2-2.1.0/hisat2_extract_exons.py $Ref_Dir/${Ref_GTF} > $Ref_Dir/${Ref_Name}_exons.txt
fi

if [[ -s $Ref_Dir/genome.1.ht2 && $Ref_Dir/genome.2.ht2 && $Ref_Dir/genome.3.ht2 && $Ref_Dir/genome.4.ht2 && $Ref_Dir/genome.5.ht2 && $Ref_Dir/genome.6.ht2 && $Ref_Dir/genome.7.ht2 && $Ref_Dir/genome.8.ht2 ]]; 
then
  echo
  echo REFERENCE INDEX FILES EXIST, MOVING ALONG
  echo
else
  echo GENERATE INDEX REQUIREMENTS FOR HISAT2 HG38
  hisat2-2.1.0/hisat2-build --ss $Ref_Dir/${Ref_Name}_ss.txt --exon $Ref_Dir/${Ref_Name}_exons.txt -f $Ref_Dir/${Ref_Anno} $Ref_Dir/${Ref_Name}
  echo
  echo DONE WITH INDEX GENERATION FOR REFERENCE
  echo
fi

if [[ -s $Ref_Dir/${Pseudovirome_Name}_ss.txt && $Ref_Dir/${Pseudovirome_Name}_exons.txt ]];
then
  echo
  echo PSEUDOVIROME EXON AND SS FILES EXIST, MOVING ALONG
  echo
else
  echo
  echo GENERATING PSEUDOVIROME EXON AND SS FILES
  echo
  hisat2-2.1.0/hisat2_extract_splice_sites.py $Ref_Dir/${Pseudovirome_GTF} > $Ref_Dir/${Pseudovirome_Name}_ss.txt
  hisat2-2.1.0/hisat2_extract_exons.py $Ref_Dir/${Pseudovirome_GTF} > $Ref_Dir/${Pseudovirome_Name}_exons.txt
fi

if [[ -s $Ref_Dir/Pseudovirome_index.1.ht2 && $Ref_Dir/Pseudovirome_index.2.ht2 && $Ref_Dir/Pseudovirome_index.3.ht2 && $Ref_Dir/Pseudovirome_index.4.ht2 && $Ref_Dir/Pseudovirome_index.5.ht2 && $Ref_Dir/Pseudovirome_index.6.ht2 && $Ref_Dir/Pseudovirome_index.7.ht2 && $Ref_Dir/Pseudovirome_index.8.ht2 ]];
then
  echo
  echo PSEUDOVIROME INDEX FILES EXIST, MOVING ALONG
  echo
else
  echo
  echo GENERATE INDEX REQUIREMENTS FOR HISAT2 PSEUDOVIROME
  echo
  hisat2-2.1.0/hisat2-build --ss $Ref_Dir/${Pseudovirome_Name}_ss.txt --exon $Ref_Dir/${Pseudovirome_Name}_exons.txt -f $Ref_Dir/${Pseudovirome} $Ref_Dir/${Pseudovirome_Name}_index
  echo
  echo DONE WITH INDEX GENERATION PSEUDOVIROME
  echo
fi

# trim input files with Trimmomatic and print quality report with FastQC
echo
echo START TRIMMING
echo

for f in $(find $File_Dir -name "*.fastq"); do
  FILENAME=$(echo $f | rev | cut -c 9- | rev | uniq)
  if [[ ! -s ${FILENAME}${PE1_FileTag} && ${FILENAME}${PE2_FileTag} ]];
  then
    echo
    echo MISSING OR EMTPY INPUT FASTQ FILES. CHECK THAT FILES EXIST OR THAT NAMES ARE CORRECT
    echo
    exit
  else  
    if [[ -s ${FILENAME}_1_fastqc.html && ${FILENAME}_2_fastqc.html ]];
    then
      echo
      echo ORIGINAL FASTQ FILES ALREADY RUN THROUGH FASTQC. SKIPPING.
      echo
    else
      echo
      echo RUNNING FASTQC ON ORIGINAL ${FILENAME} FASTQ FILES
      echo
      $FASTQC_Path/fastqc ${FILENAME}${PE1_FileTag} ${FILENAME}${PE2_FileTag}
      echo
      echo DONE WITH FASTQC ON ORIGINAL ${FILENAME} FASTQ FILES
      echo
    fi
    if [[ -s ${FILENAME}_forward_paired.fq && ${FILENAME}_forward_unpaired.fq && ${FILENAME}_reverse_paired.fq && ${FILENAME}_reverse_unpaired.fq ]];
    then
      echo
      echo FILES HAVE ALREADY BEEN TRIMMED. SKIP.
      echo
    else
      if [ ! -s ${Adapter_Dir}/${Adapter_File} ];
      then
       echo
       echo HAVING A HARD TIME FINDING THE ADAPTER FILE, OR IT IS EMPTY. QUIT.
       exit
      else
       echo
       echo TRIMMING ${FILENAME}
       echo
       java -jar ${Trimmomatic_Path}/Trimmomatic-0.36/trimmomatic-0.36.jar  PE \
        -phred33 ${FILENAME}${PE1_FileTag} ${FILENAME}${PE2_FileTag} ${FILENAME}_forward_paired.fq ${FILENAME}_forward_unpaired.fq ${FILENAME}_reverse_paired.fq ${FILENAME}_reverse_unpaired.fq \
       ILLUMINACLIP:${Adapter_Dir}/${Adapter_File}:${CLIP} LEADING:${LEADING} TRAILING:${TRAILING} SLIDINGWINDOW:${SLIDINGWINDOW} MINLEN:${MINLEN}
       echo
       echo DONE TRIMMING ${FILENAME}
       echo
      fi
    fi 
    if [[ -s ${FILENAME}_forward_paired_fastqc.html && ${FILENAME}_forward_unpaired_fastqc.html && ${FILENAME}_reverse_paired_fastqc.html && ${FILENAME}_reverse_unpaired_fastqc.html ]];
    then
      echo
      echo TRIMMED FILES HAVE ALREADY RUN THROUGH FASTQC. SKIP.
      echo
    else
      echo
      echo RUNNING FASTQC ON TRIMMED ${FILENAME} FASTQ FILES
      echo
      $FASTQC_Path/fastqc ${FILENAME}_forward_paired.fq ${FILENAME}_forward_unpaired.fq ${FILENAME}_reverse_paired.fq ${FILENAME}_reverse_unpaired.fq 
      echo
      echo DONE WITH FASTQC ON TRIMMED ${FILENAME} FASTQ FILES
      echo
    fi
  fi
done
echo
echo DONE TRIMMING
echo

#align reads to reference genome
for f in $(find $File_Dir -name "*forward_paired.fq"); do
  FILENAME=$(echo $f |  rev | cut -c 19- | rev | uniq)
  echo ${FILENAME}
  if [[ ! -s ${FILENAME}_forward_paired.fq && ${FILENAME}_forward_unpaired.fq && ${FILENAME}_reverse_paired.fq && ${FILENAME}_reverse_unpaired.fq ]];
    then
      echo
      echo TRIMMED READS ARE MISSING OR EMPTY. CHECK NAME AND INPUT FILE. 
      echo
      exit
    else
    if [ -s ${FILENAME}_HISAT2_reference_alignment.sam ];
      then
        echo
        echo REFERENCE ALIGNMENT FOR ${FILENAME} ALREADY EXISTS. SKIPPING.
        echo
      else
        echo
        echo START ALIGNING TO REFERENCE
        echo
        hisat2 --phred33 --known-splicesite-infile ${Ref_Dir}/${Ref_Name}_ss.txt --dta-cufflinks -x ${Ref_Dir}/genome -1 ${FILENAME}_forward_paired.fq -2 ${FILENAME}_reverse_paired.fq \
        -U ${FILENAME}_forward_unpaired.fq,${FILENAME}_reverse_unpaired.fq -S ${FILENAME}_HISAT2_reference_alignment.sam --un ${FILENAME}_HISAT2_reference_alignment_un-seqs \
        --un-conc ${FILENAME}_HISAT2_reference_alignment_un-conc-mate
        echo
        echo DONE ALIGNING ${FILENAME}
        echo
      fi
    fi
done
echo DONE ALIGNING ALL SAMPLES TO REFERENCE
echo

#align reads to pseudovirome
echo
echo START ALIGNING TO PSEUDOVIROME
echo
for f in $(find $File_Dir -name "*_HISAT2_reference_alignment.sam"); do
  FILENAME=$(echo $f |  rev | cut -c 32- | rev | uniq)
  echo ${FILENAME}
  if [[ ! -s ${FILENAME}_HISAT2_reference_alignment_un-conc-mate.1 && ${FILENAME}_HISAT2_reference_alignment_un-conc-mate.2 && ${FILENAME}_HISAT2_reference_alignment.sam ]];
    then
      echo
      echo UNALIGNED REFERENCE FILES NOT PRESENT FOR ${FILENAME} OR IS EMPTY. QUITTING.
      echo
      exit
    else
    if [ -s ${FILENAME}_HISAT2_virome_alignment.sam ];
      then
        echo
        echo PSEUDOVIROME ALIGNMENT FOR ${FILENAME} ALREADY EXISTS. SKIPPING.
        echo
      else
        echo
        echo START ALIGNING TO PSEUDOVIROME
        echo
        #make sure unmapped reads are paired, and if there are any unpaired reads, combine them
        ${Fastq_Pair_Path}/fastq_pair ${FILENAME}_HISAT2_reference_alignment_un-conc-mate.1 ${FILENAME}_HISAT2_reference_alignment_un-conc-mate.2
        cat ${FILENAME}_HISAT2_reference_alignment_un-conc-mate.1.single.fq ${FILENAME}_HISAT2_reference_alignment_un-conc-mate.2.single.fq \
        ${FILENAME}_HISAT2_reference_alignment_un-seqs > ${FILENAME}_HISAT2_reference_alignment_unpaired.fq
        # run!
        hisat2 --phred33 --known-splicesite-infile $Ref_Dir/${Pseudovirome_Name}_ss.txt --dta-cufflinks \
       -x $Ref_Dir/${Pseudovirome_Name}_index  -1 ${FILENAME}_HISAT2_reference_alignment_un-conc-mate.1.paired.fq \
       -2 ${FILENAME}_HISAT2_reference_alignment_un-conc-mate.2.paired.fq -U ${FILENAME}_HISAT2_reference_alignment_unpaired.fq \
       -S ${FILENAME}_HISAT2_virome_alignment.sam --un ${FILENAME}_virome_un-seqs --un-conc ${FILENAME}_virome_un-conc-mate
       echo
       echo DONE ALIGNING ${FILENAME}
       echo
      fi
    fi
done
echo
echo DONE ALIGNING ALL SAMPLES TO PSEUDOVIROME
echo

#post-alignment processing for reference
echo STARTED FILTERING ALIGNMENT FILES
echo
for f in $(find $File_Dir -name "*_HISAT2_reference_alignment.sam"); do
  FILENAME=$(echo $f | rev | cut -c 32- | rev | uniq)
   if [[ -s ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique_HML2.coverage && ${FILENAME}_HISAT2_reference_alignment_sorted.bam && ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup.bam && ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique.bam ]];
    then
      echo
      echo FILTERED REFERENCE ALIGNMENT FILES FOR ${FILENAME} ALREADY EXIST. SKIPPING.
      echo
    else
      echo
      echo START FILTERING REFERENCE ALIGNMENT FILES FOR ${FILENAME}
      echo
      #convert .sam to .bam and sort
      samtools view -Su ${FILENAME}_HISAT2_reference_alignment.sam | samtools sort - > ${FILENAME}_HISAT2_reference_alignment_sorted.bam

      #MarkDuplicates with Picard
      java -jar $Picard_Path/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=${FILENAME}_HISAT2_reference_alignment_sorted.bam \
      O=${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup.bam M=${FILENAME}_marked_dup_reference_metrics.txt 
    
      #capture unique reads 
       samtools view -bq 50 ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup.bam > ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique.bam

     #index
      samtools index ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique.bam
      samtools index ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup.bam
      samtools index ${FILENAME}_HISAT2_reference_alignment_sorted.bam

      #start calculating ERV read depth
      samtools depth -b ${Ref_Dir}/${ERV_bed} ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique.bam > \
      ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique_HML2.coverage
    fi
done
  echo DONE FILTERING ${FILENAME} REFERENCE ALIGNMENT
  echo

#post-alignment processing for Pseudovirome
for f in $(find $File_Dir -name *_HISAT2_virome_alignment.sam); do
  FILENAME=$(echo $f | rev | cut -c 29- | rev | uniq)
    if [ ! -s ${FILENAME}_HISAT2_virome_alignment.sam  ];
    then
      echo
      echo PSEUDOVIROME ALIGNMENT FILE NOT PRESENT FOR ${FILENAME} OR IS EMPTY. QUITTING.
      echo
      exit
    else
     if [[ -s ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep.coverage && ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep.bam && ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique.bam && ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup.bam && ${FILENAME}_HISAT2_virome_alignment_sorted.bam ]];
       then
         echo
         echo FILTERED PSEUDOVIROME ALIGNMENT FILES FOR ${FILENAME} ALREADY EXIST. SKIPPING.
         echo
       else
        #convert .sam to .bam and sort
        samtools view -Su ${FILENAME}_HISAT2_virome_alignment.sam | samtools sort - > ${FILENAME}_HISAT2_virome_alignment_sorted.bam

        #MarkDuplicates with Picard
        java -jar $Picard_Path/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=${FILENAME}_HISAT2_virome_alignment_sorted.bam \
        O=${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup.bam M=${FILENAME}_marked_dup_pseudovirome_metrics.txt 

        #capture unique reads 
        samtools view -bq 50 ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup.bam > ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique.bam

        #remove tandem repeats
        bedtools intersect -abam ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique.bam -b $Ref_Dir/${Pseudovirome_TanRep} \
        -v > ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep.bam

        #start calculating read depth
        samtools depth -a ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep.bam > \
        ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep.coverage

        #index
        samtools index ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep.bam
        samtools index ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique.bam
        samtools index ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup.bam
        samtools index ${FILENAME}_HISAT2_virome_alignment_sorted.bam
      fi
    fi
  echo DONE FILTERING ${FILENAME}
  echo
done
echo
echo DONE FILTERING ALL FILES
echo


#estimate counts with htseq-counts
echo
echo START HTSEQ-COUNTS FOR REFERENCE
echo

if [[ ! -s ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique.bam && ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup.bam && ${FILENAME}_HISAT2_reference_alignment_sorted.bam ]];
    then
      echo
      echo FILTERED REFERENCE FILES DO NOT EXIST OR ARE EMPTY FOR ${FILENAME}. QUITTING.
      echo
      exit
    else
      for f in $(find $File_Dir -name *_HISAT2_reference_alignment_sorted_MrkDup_unique.bam); do
      FILENAME=$(echo $f | rev | cut -c 53- | rev | uniq)
      if [ -s ${FILENAME}_HISAT2_reference_htseq_sorted_single_end_137_counts.txt ];
      then
        echo
        echo COUNTING HAS ALREADY BEEN COMPLETED FOR ${FILENAME} REFERENCE ALIGNMENT
        echo
      else
        echo
        echo START HTSEQ-COUNTS CALCULATIONS FOR ${FILENAME}
        echo
        samtools view -bf 3 ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique.bam | \
        htseq-count -f ${F} --nonunique ${nonunique} -r ${R} -s ${S} -t ${T} -i ${I} - ${Ref_Dir}/${Ref_GTF} > ${FILENAME}_HISAT2_reference_htseq_sorted_MrkDup_unique_paired_end_3_counts.txt

        samtools view -bf 73 ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique.bam | \
        htseq-count -f ${F} --nonunique ${nonunique} -r ${R} -s ${S} -t ${T} -i ${I} - ${Ref_Dir}/${Ref_GTF} > ${FILENAME}_HISAT2_reference_htseq_sorted_MrkDup_unique_single_end_73_counts.txt
   
        samtools view -bf 89 ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique.bam | \
        htseq-count -f ${F} --nonunique ${nonunique} -r ${R} -s ${S} -t ${T} -i ${I} - ${Ref_Dir}/${Ref_GTF} > ${FILENAME}_HISAT2_reference_htseq_sorted_MrkDup_unique_single_end_89_counts.txt

        samtools view -bf 153 ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique.bam | \
        htseq-count -f ${F} --nonunique ${nonunique} -r ${R} -s ${S} -t ${T} -i ${I} - ${Ref_Dir}/${Ref_GTF} > ${FILENAME}_HISAT2_reference_htseq_sorted_MrkDup_unique_single_end_153_counts.txt

        samtools view -bf 137 ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup_unique.bam | \
        htseq-count -f ${F} --nonunique ${nonunique} -r ${R} -s ${S} -t ${T} -i ${I} - ${Ref_Dir}/${Ref_GTF} > ${FILENAME}_HISAT2_reference_htseq_sorted_MrkDup_unique_single_end_137_counts.txt

        samtools view -bf 3 ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup.bam | \
        htseq-count -f ${F} --nonunique ${nonunique} -r ${R} -s ${S} -t ${T} -i ${I} - ${Ref_Dir}/${Ref_GTF} > ${FILENAME}_HISAT2_reference_htseq_sorted_MrkDup_paired_end_3_counts.txt

        samtools view -bf 73 ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup.bam | \
        htseq-count -f ${F} --nonunique ${nonunique} -r ${R} -s ${S} -t ${T} -i ${I} - ${Ref_Dir}/${Ref_GTF} > ${FILENAME}_HISAT2_reference_htseq_sorted_MrkDup_single_end_73_counts.txt
   
        samtools view -bf 89 ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup.bam | \
        htseq-count -f ${F} --nonunique ${nonunique} -r ${R} -s ${S} -t ${T} -i ${I} - ${Ref_Dir}/${Ref_GTF} > ${FILENAME}_HISAT2_reference_htseq_sorted_MrkDup_single_end_89_counts.txt

        samtools view -bf 153 ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup.bam | \
        htseq-count -f ${F} --nonunique ${nonunique} -r ${R} -s ${S} -t ${T} -i ${I} - ${Ref_Dir}/${Ref_GTF} > ${FILENAME}_HISAT2_reference_htseq_sorted_MrkDup_single_end_153_counts.txt

        samtools view -bf 137 ${FILENAME}_HISAT2_reference_alignment_sorted_MrkDup.bam | \
        htseq-count -f ${F} --nonunique ${nonunique} -r ${R} -s ${S} -t ${T} -i ${I} - ${Ref_Dir}/${Ref_GTF} > ${FILENAME}_HISAT2_reference_htseq_sorted_MrkDup_single_end_137_counts.txt

        samtools view -bf 3 ${FILENAME}_HISAT2_reference_alignment_sorted.bam | \
        htseq-count -f ${F} --nonunique ${nonunique} -r ${R} -s ${S} -t ${T} -i ${I} - ${Ref_Dir}/${Ref_GTF} > ${FILENAME}_HISAT2_reference_htseq_sorted_paired_end_3_counts.txt

        samtools view -bf 73 ${FILENAME}_HISAT2_reference_alignment_sorted.bam | \
        htseq-count -f ${F} --nonunique ${nonunique} -r ${R} -s ${S} -t ${T} -i ${I} - ${Ref_Dir}/${Ref_GTF} > ${FILENAME}_HISAT2_reference_htseq_sorted_single_end_73_counts.txt
   
        samtools view -bf 89 ${FILENAME}_HISAT2_reference_alignment_sorted.bam | \
        htseq-count -f ${F} --nonunique ${nonunique} -r ${R} -s ${S} -t ${T} -i ${I} - ${Ref_Dir}/${Ref_GTF} > ${FILENAME}_HISAT2_reference_htseq_sorted_single_end_89_counts.txt

        samtools view -bf 153 ${FILENAME}_HISAT2_reference_alignment_sorted.bam | \
        htseq-count -f ${F} --nonunique ${nonunique} -r ${R} -s ${S} -t ${T} -i ${I} - ${Ref_Dir}/${Ref_GTF} > ${FILENAME}_HISAT2_reference_htseq_sorted_single_end_153_counts.txt

        samtools view -bf 137 ${FILENAME}_HISAT2_reference_alignment_sorted.bam | \
        htseq-count -f ${F} --nonunique ${nonunique} -r ${R} -s ${S} -t ${T} -i ${I} - ${Ref_Dir}/${Ref_GTF} > ${FILENAME}_HISAT2_reference_htseq_sorted_single_end_137_counts.txt
        echo
        echo DONE HTSEQ-COUNTS CALCULATIONS FOR ${FILENAME} REFERENCE
        echo
      fi
     done
     fi

if [[ ! -s ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep.bam && ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique.bam && ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup.bam && ${FILENAME}_HISAT2_virome_alignment_sorted.bam ]];
    then
      echo
      echo FILTERED PSEUDOVIROME FILES DO NOT EXIST OR ARE EMPTY FOR ${FILENAME}. QUITTING.
      echo
      exit
    else
     for f in $(find $File_Dir  -name *_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep.bam); do
     FILENAME=$(echo $f | rev | cut -c 59- | rev | uniq)
     if [ -s ${FILENAME}_HISAT2_virome_htseq_sorted_single_end_137_counts.txt ];
       then
        echo
        echo COUNTING HAS ALREADY BEEN COMPLETED FOR ${FILENAME} PSEUDOVIROME ALIGNMENT
        echo
       else
        echo
        echo START HTSEQ-COUNTS CALCULATIONS FOR ${FILENAME}  PSEUDOVIROME
        echo
        samtools view -bf 3 ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep.bam | \
        htseq-count -f ${F} -r ${R} -s ${S} -t ${T} -i ${IP} --nonunique ${nonunique} - $Ref_Dir/${Pseudovirome_GTF} > ${FILENAME}_HISAT2_virome_htseq_sorted_MrkDup_unique_noTanRep_paired_end_3_counts.txt

        samtools view -bf 73 ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep.bam | \
        htseq-count -f ${F}  -r ${R} -s ${S} -t ${T} -i ${IP} --nonunique ${nonunique} - $Ref_Dir/${Pseudovirome_GTF} > ${FILENAME}_HISAT2_virome_htseq_sorted_MrkDup_unique_noTanRep_single_end_73_counts.txt
   
        samtools view -bf 89 ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep.bam | \
        htseq-count -f ${F} -r ${R} -s ${S} -t ${T} -i ${IP} --nonunique ${nonunique} - $Ref_Dir/${Pseudovirome_GTF} > ${FILENAME}_HISAT2_virome_htseq_sorted_MrkDup_unique_noTanRep_single_end_89_counts.txt

        samtools view -bf 153 ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep.bam | \
        htseq-count -f ${F}  -r ${R} -s ${S} -t ${T} -i ${IP} --nonunique ${nonunique} - $Ref_Dir/${Pseudovirome_GTF} > ${FILENAME}_HISAT2_virome_htseq_sorted_MrkDup_unique_noTanRep_single_end_153_counts.txt

        samtools view -bf 137 ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep.bam | \
        htseq-count -f ${F}  -r ${R} -s ${S} -t ${T} -i ${IP} --nonunique ${nonunique} - $Ref_Dir/${Pseudovirome_GTF} > ${FILENAME}_HISAT2_virome_htseq_sorted_MrkDup_unique_noTanRep_single_end_137_counts.txt

        samtools view -bf 3 ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique.bam | \
        htseq-count -f ${F}  -r ${R} -s ${S} -t ${T} -i ${IP} --nonunique ${nonunique} - $Ref_Dir/${Pseudovirome_GTF} > ${FILENAME}_HISAT2_virome_htseq_sorted_MrkDup_unique_paired_end_3_counts.txt

        samtools view -bf 73 ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique.bam | \
        htseq-count -f ${F}  -r ${R} -s ${S} -t ${T} -i ${IP} --nonunique ${nonunique} - $Ref_Dir/${Pseudovirome_GTF} > ${FILENAME}_HISAT2_virome_htseq_sorted_MrkDup_unique_single_end_73_counts.txt
   
        samtools view -bf 89 ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique.bam | \
        htseq-count -f ${F}  -r ${R} -s ${S} -t ${T} -i ${IP} --nonunique ${nonunique} - $Ref_Dir/${Pseudovirome_GTF} > ${FILENAME}_HISAT2_virome_htseq_sorted_MrkDup_unique_single_end_89_counts.txt

        samtools view -bf 153 ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique.bam | \
        htseq-count -f ${F}  -r ${R} -s ${S} -t ${T} -i ${IP} --nonunique ${nonunique} - $Ref_Dir/${Pseudovirome_GTF} > ${FILENAME}_HISAT2_virome_htseq_sorted_MrkDup_unique_single_end_153_counts.txt

        samtools view -bf 137 ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup_unique.bam | \
        htseq-count -f ${F} -r ${R} -s ${S} -t ${T} -i ${IP} --nonunique ${nonunique} - $Ref_Dir/${Pseudovirome_GTF} > ${FILENAME}_HISAT2_virome_htseq_sorted_MrkDup_unique_single_end_137_counts.txt

        samtools view -bf 3 ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup.bam | \
        htseq-count -f ${F}  -r ${R} -s ${S} -t ${T} -i ${IP} --nonunique ${nonunique} - $Ref_Dir/${Pseudovirome_GTF} > ${FILENAME}_HISAT2_virome_htseq_sorted_MrkDup_paired_end_3_counts.txt

        samtools view -bf 73 ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup.bam | \
        htseq-count -f ${F}  -r ${R} -s ${S} -t ${T} -i ${IP} --nonunique ${nonunique} - $Ref_Dir/${Pseudovirome_GTF} > ${FILENAME}_HISAT2_virome_htseq_sorted_MrkDup_single_end_73_counts.txt
   
        samtools view -bf 89 ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup.bam | \
        htseq-count -f ${F}  -r ${R} -s ${S} -t ${T} -i ${IP} --nonunique ${nonunique} - $Ref_Dir/${Pseudovirome_GTF} > ${FILENAME}_HISAT2_virome_htseq_sorted_MrkDup_single_end_89_counts.txt

        samtools view -bf 153 ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup.bam | \
        htseq-count -f ${F}  -r ${R} -s ${S} -t ${T} -i ${IP} --nonunique ${nonunique} - $Ref_Dir/${Pseudovirome_GTF} > ${FILENAME}_HISAT2_virome_htseq_sorted_MrkDup_single_end_153_counts.txt
  
        samtools view -bf 137 ${FILENAME}_HISAT2_virome_alignment_sorted_MrkDup.bam | \
        htseq-count -f ${F}  -r ${R} -s ${S} -t ${T} -i ${IP} --nonunique ${nonunique} - $Ref_Dir/${Pseudovirome_GTF} > ${FILENAME}_HISAT2_virome_htseq_sorted_MrkDup_single_end_137_counts.txt

        samtools view -bf 3 ${FILENAME}_HISAT2_virome_alignment_sorted.bam | \
        htseq-count -f ${F}  -r ${R} -s ${S} -t ${T} -i ${IP} --nonunique ${nonunique} - $Ref_Dir/${Pseudovirome_GTF} > ${FILENAME}_HISAT2_virome_htseq_sorted_paired_end_3_counts.txt

        samtools view -bf 73 ${FILENAME}_HISAT2_virome_alignment_sorted.bam | \
        htseq-count -f ${F}  -r ${R} -s ${S} -t ${T} -i ${IP} --nonunique ${nonunique} - $Ref_Dir/${Pseudovirome_GTF} > ${FILENAME}_HISAT2_virome_htseq_sorted_single_end_73_counts.txt
     
        samtools view -bf 89 ${FILENAME}_HISAT2_virome_alignment_sorted.bam | \
        htseq-count -f ${F}  -r ${R} -s ${S} -t ${T} -i ${IP} --nonunique ${nonunique} - $Ref_Dir/${Pseudovirome_GTF} > ${FILENAME}_HISAT2_virome_htseq_sorted_single_end_89_counts.txt

        samtools view -bf 153 ${FILENAME}_HISAT2_virome_alignment_sorted.bam | \
        htseq-count -f ${F}  -r ${R} -s ${S} -t ${T} -i ${IP} --nonunique ${nonunique} - $Ref_Dir/${Pseudovirome_GTF} > ${FILENAME}_HISAT2_virome_htseq_sorted_single_end_153_counts.txt

        samtools view -bf 137 ${FILENAME}_HISAT2_virome_alignment_sorted.bam | \
        htseq-count -f ${F}  -r ${R} -s ${S} -t ${T} -i ${IP} --nonunique ${nonunique} - $Ref_Dir/${Pseudovirome_GTF} > ${FILENAME}_HISAT2_virome_htseq_sorted_single_end_137_counts.txt
        echo
        echo DONE HTSEQ-COUNTS CALCULATIONS FOR ${FILENAME} PSEUDOVIROME
        echo
        fi
     done
     fi

mkdir ${Working_Dir}/Trimmed_Reads ${Working_Dir}/FASTQC_output ${Working_Dir}/fastq_original ${Working_Dir}/Pseudovirome_Alignment ${Working_Dir}/Reference_Alignment
mv ${File_Dir}/*forward_paired.fq ${File_Dir}/*forward_unpaired.fq ${File_Dir}/*reverse_paired.fq ${File_Dir}/*reverse_unpaired.fq ${Working_Dir}/Trimmed_Reads
mv ${File_Dir}/*html ${File_Dir}/*zip ${Working_Dir}/FASTQC_output
mv ${File_Dir}/*fastq ${Working_Dir}/fastq_original
mv ${File_Dir}/*_HISAT2_virome_alignment.sam ${File_Dir}/*virome_un-conc-mate* ${File_Dir}/*virome_un-seqs ${Working_Dir}/Pseudovirome_Alignment
mv ${File_Dir}/*HISAT2_reference_alignment.sam ${File_Dir}/*HISAT2_reference_alignment_un-seqs ${File_Dir}/*HISAT2_reference_alignment_un-conc* ${Working_Dir}/Reference_Alignment

echo
echo DONE WITH HTSEQ-COUNTS
echo

echo
echo READY FOR R
echo

