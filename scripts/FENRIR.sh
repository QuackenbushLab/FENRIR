#Description: Overworld script in charge of executing individual steps within the FENRIR pipeline

#Input:
  #1) fastq files
  #2) Reference genome and pseudovirome
  #3) list of proviruses of interest, otherwise use the one provided
  #4) Depending on size of the reference genome, index files for HISAT2
  #5) .bed file containing tandem repeats in pseudovirome
  #6) .bed of ERVs, otherwise use the one provided
  #7) Variables.txt - filled out to match your experimental design

#Output:
  #1) combined counts matrix for endogenous and exogenous viruses detected in RNA-seq, respectively
  #2) graphed viral counts for endogenous and exogenous viruses
  #3) percent coverage matricies for endogenous and exogenous viruses
  #4) graphed coverage of endogenous and exogenous viruses
  #5) reads aligning to a given virus stored in a .fasta for BLAST verification

#Author:
  #Farrah Roy

#!/bin/bash
source Variables.txt
echo $READING
echo


${Working_Dir}/FENRIR_Alignment.sh

Rscript ${Working_Dir}/FENRIR_R_packages.R

Rscript --vanilla ${Working_Dir}/FENRIR_Counts.R ${Working_Dir}/Variables.txt

${Working_Dir}/FENRIR_BLAST_Prep.sh
