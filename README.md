# Focused Exploration of NGS Reads for Infectious RNAs (FENRIR): an RNA-Seq pipeline

### Last updated: April 9th, 2020
- *more compact script with two execution commands instead of many*


We present a pipeline for the detection of exogenous human-tropic and endogenous virus-specific reads within a given dataset, with the goal of providing an easily executable pipeline for those who are unfamiliar with bioinformatics and wish to collect virus expression data from their paired-end, unstranded RNA-seq data. FENRIR includes our trimming, alignment, quality filtering, and viral counts protocols for both endogenous and exogenous viruses. This pipeline also graphs virus expression data and RNA-seq coverage over detected viruses, exporting matricies of virus fragment counts and percent coverage over a viral genome, and jpegs of expression data and virus genome coverage. 

For details on what is in this repo as well as a guide for how to download and run this package, please keep reading.

Here is an overview of the FENRIR pipeline before we begin.

<p align="center">
  <img src="https://github.com/QuackenbushLab/FENRIR/blob/master/FENRIR_pipeline.jpg" height="500">
</p>


## Pipeline components

The following tools are required for FENRIR. A script within the FENRIR pipeline is responsible for downloading and installing all of the necessary tools, so installing these tools on ones own isn't needed:
**For shell scripts**:
- [**FastQC**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): *to check quality of reads pre- and post-trimming*
- [**Trimmomatic**](http://www.usadellab.org/cms/?page=trimmomatic): *trimming software to remove short and/or low-quality reads*
- [**HISAT2**](https://ccb.jhu.edu/software/hisat2/index.shtml): *aligner used to remove human reads and align unmapped reads to the pseudovirome reference*
- [**Samtools**](http://samtools.sourceforge.net/): *filtering pseudovirome alignment output to get viral read counts*
- [**Bedtools**](https://bedtools.readthedocs.io/en/latest/index.html): *help remove tandem repeats from the pseudovirome alignment output*
- [**Python**](https://www.python.org/downloads/): *HISAT2 needs this*
- [**Java**](https://packages.ubuntu.com/hu/xenial/openjdk-8-jdk): *Trimmomatic and Picard need this*
- [**Picard**](https://broadinstitute.github.io/picard/): *used to remove PCR duplicates from alignment output*
- [**HTSeq-Count**](https://htseq.readthedocs.io/en/release_0.11.1/count.html): *used to count exogenous and endogenous reads*
- [**R**](https://www.r-project.org/): *used to merge read counts into one .csv file and used to calculate coverage of both a virus and ORF of interest*
- [**TandemRepeatsFinder**](https://tandem.bu.edu/trf/trf.html): *removes tandem repeats from alignment output. There is  an online tool that one can use to generate a list of tandem repeats found within a desired genome. The list of tandem repeats for our pseudovirome is included in this repo (Pseudovirome_2.7.7.80.10.50.2000.bed).*
- [**fastq-pair**](https://github.com/linsalrob/fastq-pair): *pairs unmapped reads prior to pseudovirome alignment*

**For R**:
- [**Reshape**](https://cran.r-project.org/web/packages/reshape/index.html): *reorganizes count matricies so they may be graphed*
- [**Pheatmap**](https://cran.r-project.org/web/packages/pheatmap/index.html): *makes heatmaps out of endogenous virus counts data*
- [**data.table**](https://cran.r-project.org/web/packages/data.table/index.html): *to join columns together to generate an aggregated counts matrix across all samples*
- [**tidyverse**](https://cran.r-project.org/web/packages/tidyverse/index.html): *to create tidy data*
- [**BiocManager**](https://www.bioconductor.org/install/): *to download packages from Bioconductor*
- [**GenomicAlignments**](https://bioconductor.org/packages/release/bioc/html/GenomicAlignments.html): *to work with RNA-seq data in R*
- [**refGenome**](https://cran.r-project.org/web/packages/refGenome/index.html):*to work with RNA-seq data in R*
- [**RSQLite**](https://cran.r-project.org/web/packages/RSQLite/index.html)
- [**doBy**](https://cran.r-project.org/web/packages/doBy/index.html)


**The following programs will also be installed, but are optional depending on what one needs to do**:
- [**SRA Toolkit**](https://ncbi.github.io/sra-tools/install_config.html): *to download test data and split .sra files into .fastq files*
- [**gffread**](http://ccb.jhu.edu/software/stringtie/gff.shtml): *used to convert gff3 to gtf, in the event one wishes to make their own reference and/or pseudovirome*
- [**RStudio**](https://rstudio.com/products/rstudio/download/):*in case one wishes a more visual approach to running the R scripts*
- [**BEDOPS**](https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/convert2bed.html):*to create a .bed file of ones own endogenous viruses of interest*


**The following will not be automatically installed, but since these are optional one can download them at the given locations**:
- [**Tandem Repeat Finder**](https://tandem.bu.edu/trf/trf.download.html): *used to identify tandem repeats from the given pseudovirome, generating a .dat file that can be converted into a .bed. Can be downloaded from the web browser.*
- [**Tandem Repeat Finder dat to bed/txt**](https://github.com/hdashnow/TandemRepeatFinder_scripts): *used to convert the Tandem Repeat Finder output .dat file into a .bed file. This .bed file is used to remove tandem repeats from RNA-seq data. Can be cloned from the github repository with the command git clone https://github.com/hdashnow/TandemRepeatFinder_scripts.git .*



# Assembling reference files before running FENRIR

There are several items one needs to have assembled before beginning the pipeline: 
- the host reference genome 
- a reference gtf
- the pseudovirome
- a pseudovirome gtf
- a .bed file of all your proviruses of interest
- a .bed file containing all the detectable tandem repeats in your pseudovirome reference
- a list of all your proviruses of interest in a .txt
- in the case of large genomes, such as hg38, index files for your reference genome 

Within the ref folder, we provide all of these files. Below is an explaination on either where we acquired the files or how we generated them in the event one needs to generate their own.

## Host reference genome 

The human genome reference build (hg38) can be downloaded [here](https://support.illumina.com/sequencing/sequencing_software/igenome.html) from the UCSC source. This contains both the fasta sequence and the gtf annotation. Since the .gtf typically does not annotate all of the known endogenous proviruses, one should append the necessary information for the provirus clade one is interested in to the end of the host.gtf. The format for these added proviruses should follow [these guidelines](https://genome.ucsc.edu/FAQ/FAQformat.html#format4) and can be generated using [Genbank](https://www.ncbi.nlm.nih.gov/genbank/) to obtain the coordinates as a .gff3 and [gffread](http://ccb.jhu.edu/software/stringtie/gff.shtml) to convert .gff3 to .gtf. These new entries can then be appended either via command line (example: cat provirus.gtf >> reference.gtf) or with a common text editor (TextWrangler, for example). In the case of HERV-K(HML-2), we have already appended their sequences to the end of the provided hg38.gtf. Considering the size of hg38 and the resource intensive process of building an index with HISAT2, HISAT2 provides an [index for hg38](https://ccb.jhu.edu/software/hisat2/index.shtml). 

## Pseudovirome reference genome

The pseudovirome reference was built through combining the sequences of 3,152 virus sequences from an array of virus families and isolates - as was the case with Ebolavirus totaling 1,981 isolate entries of our pseudovirome reference. Virus sequences and their accompanying .gff3 were collected from the [Viral Genomes Resource](https://www.ncbi.nlm.nih.gov/genome/viruses/) and [Genbank](https://www.ncbi.nlm.nih.gov/genbank/), the latter being used to extract the sequences and annotations for all of the Ebola virus isolates. The appended .gff3 was converted into a gtf with gffread, and both the genome and .gtf are provided. 

## .bed file of proviruses

While we also provide this file for HML-2s, one can generate this from their appended reference.gtf using BEDOPS convert2bed followed by grep to save only proviruses to a .bed file.

## .bed file of tandem repeats in Pseudovirome

To get this, one can run Tandem Repeat Finder on the pseudovirome.fasta to generate a .dat file containing all of the detected tandem repeats. This .dat file can be converted to a .bed using python and the above mentioned Tandem Repeat Finder scripts for .dat to .bed conversion. 

## List of proviruses of interest in a .txt

This should be easy enough to generate. Use whichever method you prefer as long as the text file one provirus name per row. As with all of the other variables, we provide a list for HERV-K(HML-2).


# Before beginning the pipeline

## File structure

Before running the pipeline, assemble all files into the right file structure. 

In the very least, the main directory should contain the following and look something like this:

- Programs.sh: the script responsible for downloading all of the programs needed for this pipeline (with the exception of R libraries)
- FENRIR.sh: the overworld script in charge of executing all of the other scripts with exception to Programs.sh
- FENRIR_Alignment.sh: the part of the pipeline responsible for the QC, trimming, aligning, and counting fragments
- FENRIR_R_packages.R: the R script reponsible for downloading all of the necessary R libraries
- FENRIR_Counts.R: merges the counts data into one matrix per alignment type (i.e. reference or pseudovirome), graphs the expression data, calculates coverage over the exogenous or endogenous viral genome, and graphs coverage over said viral genome.
- FENRIR_BLAST_Prep.sh: extracts reads that map to a specified exogenous or endogenous virus and saves them into their own file for BLAST
- Variables.txt: a text file containing all of the variables one wishes to change from the standard pipeline. **this must be checked and modified before running the pipeline to fit with ones experiment**
- input: where all of the input .fastq files should be stored
- ref: where all reference files (reference genome, pseudovirome, .gtfs, .beds, etc) should be stored

<p align="center">
  <img src="https://github.com/QuackenbushLab/FENRIR/blob/master/working_directory_structure.png" height="25">
</p>           

# Running the pipeline

## Installing programs

Once all files are properly organized, run the Programs.sh script:

```
./Programs.sh
```

Don't wander too far because the script will require interactive input as it's running. Once it's completed, the overworld directory should look like this:

<p align="center">
  <img src="https://github.com/QuackenbushLab/FENRIR/blob/master/working_directory_structure_programs.png" height="25">
</p>           

As for the reference directory, it should look something like this before you begin the alignment software:

<p align="center">
  <img src="https://github.com/QuackenbushLab/FENRIR/blob/master/reference.png" height="50 width="75">
</p>

Of note, there is one additional file that is not necessary to use but we use as a catch all: merged.fa. This file contains all of the paired end adapter sequences Trimmomatic provides within its software. This is also provided with the pipeline.

As an aside, there are a few extra files/programs in this directory that do not need to be present. These files/programs are: Pseudovirome.fasta.2.7.7.80.10.50.2000.dat (can be deleted once finished with Tandem Repeat Finder), TandemRepeatFinder_scripts (the directory containing the dat to bed script, which can be stored in any directory), hg38 (the directory the genome index files were downloaded in and then moved out of. This can be deleted), hg38.bed (used to get HML2.bed, so once HML2.bed is generated this file can be deleted), hg38.tar.gz (the hg38 index file original download, which can be deleted), igenomes.illumina.com.s3-website-us-east-1.amazonaws.com (the hg38 genome download, which can be deleted), trf409.linux64 (Tandem Repeat Finder software that can be stored elsewhere), and sequence.gff3 (the original Pseudovirome .gff3 download, which can be deleted once converted to .gtf).  

## _Optional: downloading test data_

The data we used to build this pipeline was obtained from [Peng et al's 2014 paper on HIV-infected cells](https://www.ncbi.nlm.nih.gov/pubmed/24850744), with their [RNA-seq data stored in the SRA database](https://www.ncbi.nlm.nih.gov/sra?term=SRP035316). To download these files and extract fastq data, one can use the SRA toolkit with the following commands:

```
prefetch SRR1106189 SRR1106190 SRR1106191 SRR1106192 SRR1106193 SRR1106194 SRR1106195 SRR1106196 SRR1106197 SRR1106198 SRR1106199 SRR1106200 SRR1106201 SRR1106202 SRR1106203 SRR1106204

cd ncbi/public/sra/

fastq-dump --split-files SRR1106189.sra SRR1106190.sra SRR1106191.sra SRR1106192.sra SRR1106193.sra SRR1106194.sra SRR1106195.sra SRR1106196.sra SRR1106197.sra SRR1106198.sra SRR1106199.sra SRR1106200.sra SRR1106201.sra SRR1106202.sra SRR1106203.sra SRR1106204.sra

mv *fastq path/to/input
```

The above command also moves the .fastq files into the input directory.

## Running FENERIR

Once all files and programs are present and accounted for, one executes this command to run the rest of the pipeline:

```
./FENRIR.sh
```

If you wish to save standard out as a log for later, you can instead execute the following command:

```
./FENRIR.sh &> FENRIR_stdout.txt
```

**Do not forget to edit the Variables.txt file so it is in agreement with the input data (fastq files, different reference files if one choses to use a different set of reference files, etc) prior to running FENRIR.sh.**

# Output

Once the pipeline is done running, the overworld directory should look something like this:

<p align="center">
  <img src="https://github.com/QuackenbushLab/FENRIR/blob/master/overworld_output.png" height="75" width="2000">
</p>

You should find the following within each directory:
- input: empty
- ref: reference files, including fasta sequences for the host and pseudovirome, gtf annotations for both genomes, HISAT2 index files for both genomes, and exon and splice site reports for HISAT2 (see below)
<p align="center">
  <img src="https://github.com/QuackenbushLab/FENRIR/blob/master/Ref_output.png" height="75" width="2000">
</p>

- Reference_Alignment: reference alignment files (.sam, .bam, .bai, and unmapped reads in fastq format) (see below an example for one sample)
<p align="center">
  <img src="https://github.com/QuackenbushLab/FENRIR/blob/master/Ref_Alignment_output.png" height="250" width="600">
</p>

- Pseudovirome_Alignment: pseudovirome alignment files(.sam, .bam, .bai, and unmapped reads in fastq format) (see below an example for one sample)
<p align="center">
  <img src="https://github.com/QuackenbushLab/FENRIR/blob/master/Pseudovirome_alignment_output.png" height="250" width="600">
</p>

- FASTQC_output: the FASTQC results for the original fastq files and the post-trimming fastq files

- Trimmed_Reads: all of the trimmed fastq files

- coverage_data: percent coverage matricies for both exogenous and endogenous viruses

- htseq_counts_data: fragment counts matricies for both exogenous and endogenous viruses

- JPEG_Results: where the graphed exogenous and endogenous counts data along with graphed coverage data is stored

- fastq_original: where the originial fastq files are moved to after FENRIR_Alignment.sh is done

- individual_counts_files: where all of the individual sample counts matricies are stored on case one wishes to look at them

- BLAST_results: files (.fa) containing extracted reads that mapped to a detected virus 

- Picard_metrics: contains the Picard report from the mark duplicates step

Also of note is the FENRIR_R.log report. This contains the standard out from running FENRIR_Counts.R in case one wishes to look at the output.               
                      
  
Happy hunting.
