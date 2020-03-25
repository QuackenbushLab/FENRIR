#!/bin/bash
#variables
Ref_Dir=/root/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/ref

# run!
for f in $(find $File_Dir/ -name "*HISAT2_hg38_HML2_alignment.sam"); do
  FILENAME=$(echo $f | rev | cut -c 32- | rev | uniq)
  echo ${FILENAME}

   hisat2 --phred33 --known-splicesite-infile $Ref_Dir/Pseudovirome_ss.txt --dta-cufflinks -x $Ref_Dir/Pseudovirome_index \
      -1 ${FILENAME}_HISAT2_hg38_HML2_alignment_un-conc-mate.1.paired.fq -2 ${FILENAME}_HISAT2_hg38_HML2_alignment_un-conc-mate.2.paired.fq \
      -U ${FILENAME}_HISAT2_hg38_HML2_alignment_unpaired.fq -S ${FILENAME}_HISAT2_virome_alignment_Cufflinks.sam \
      --un ${FILENAME}_HISAT2_virome_alignment_Cufflinks_un-seqs --un-conc ${FILENAME}_HISAT2_virome_alignment_Cufflinks_un-conc-mate
done
