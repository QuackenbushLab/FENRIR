#!/bin/bash
source Variables.txt
echo $READING
echo


${Working_Dir}/FENRIR_Alignment.sh

Rscript ${Working_Dir}/FENRIR_R_packages.R

Rscript --vanilla ${Working_Dir}/FENRIR_Counts.R ${Working_Dir}/Variables.txt

${Working_Dir}/FENRIR_BLAST_Prep.sh
