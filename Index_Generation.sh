#!/bin/bash
Ref_Dir=/root/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/ref
#generate .ss and .exon files
/hisat2-2.1.0/extract_splice_sites.py $Ref_Dir/Pseudovirome.gtf > $Ref_Dir/Pseudovirome_ss.txt
/hisat2-2.1.0/extract_exons.py $Ref_Dir/Pseudovirome.gtf > $Ref_Dir/Pseudovirome_exons.txt
/hisat2-2.1.0/extract_splice_sites.py $Ref_Dir/hg38.gtf > $Ref_Dir/hg38_ss.txt
/hisat2-2.1.0/extract_exons.py $Ref_Dir/hg38.gtf > $Ref_Dir/hg38_exons.txt


#generate  index for pseudovirome
/hisat2-2.1.0/hisat2-build --ss $Ref_Dir/Pseudovirome_ss.txt --exon $Ref_Dir/Pseudovirome_exons.txt -f $Ref_Dir/Pseudovirome.fa $Ref_Dir/Pseudovirome_index
