# Focused Exploration of NGS Reads for Infectious RNAs (FENRIR): an RNA-Seq pipeline

### Last updated: April 9th, 2020
- *more compact script with two execution commands instead of many*


We present a pipeline for the detection of exogenous human-tropic and endogenous virus-specific reads within a given dataset, with the goal of providing an easily executable pipeline for those who are unfamiliar with bioinformatics and wish to collect virus expression data from their paired-end, unstranded RNA-seq data. FENRIR includes our trimming, alignment, quality filtering, and viral counts protocols for both endogenous and exogenous viruses. This pipeline also graphs virus expression data and RNA-seq coverage over detected viruses, exporting matricies of virus fragment counts and percent coverage over a viral genome, and jpegs of expression data and virus genome coverage. 

For details on what is in this repo as well as a guide for how to download and run this package, please keep reading.

<img src="https://github.com/QuackenbushLab/FENRIR/blob/master/FENRIR_pipeline.jpg" height="2000">

## Image contents and pipeline components

The following tools are required for FENRIR. A script within the FENRIR pipeline is responsible for downloading and installing all of the necessary tools, so installing these tools on ones own isn't needed:
- [**FastQC**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): *to check quality of reads pre- and post-trimming*
- [**Trimmomatic**](http://www.usadellab.org/cms/?page=trimmomatic): *trimming software to remove short and/or low-quality reads*
- [**HISAT2**](https://ccb.jhu.edu/software/hisat2/index.shtml): *aligner used to remove human reads and align unmapped reads to the pseudovirome reference*
- [**Samtools**](http://samtools.sourceforge.net/): *filtering pseudovirome alignment output to get viral read counts*
- [**Bedtools**](https://bedtools.readthedocs.io/en/latest/index.html): *help remove tandem repeats from the pseudovirome alignment output*
- [**Python**](https://www.python.org/downloads/): *HISAT2 needs this*
- [**Java**](https://packages.ubuntu.com/hu/xenial/openjdk-8-jdk): *Trimmomatic and Picard need this*
- [**Picard**](https://broadinstitute.github.io/picard/): *used to remove PCR duplicates from alignment output*
- [**featureCounts**](http://bioinf.wehi.edu.au/featureCounts/): *used to count exogenous and endogenous reads*
- [**R**](https://www.r-project.org/): *used to merge read counts into one .csv file and used to calculate coverage of both a virus and ORF of interest*
- [**TandemRepeatsFinder**](https://tandem.bu.edu/trf/trf.html): *removes tandem repeats from alignment output. There is  an online tool that one can use to generate a list of tandem repeats found within a desired genome. The list of tandem repeats for our pseudovirome is included in this repo (Pseudovirome_2.7.7.80.10.50.2000.bed).*
- [**fastq-pair**](https://github.com/linsalrob/fastq-pair): *pairs unmapped reads prior to pseudovirome alignment*


- [**BLAST**](https://blast.ncbi.nlm.nih.gov/Blast.cgi): *used to screen the output for reads that do not align to exogenous human-tropic viruses so they may be removed, and to look up the function of ORFs with 50% or more coverage.*
- [**IGV**](http://software.broadinstitute.org/software/igv/): *used to view alignment output. Provides an image for virus and ORF coverage.*

# Installation

## Reference genome 

The human genome reference build (hg38) can be downloaded [here](https://support.illumina.com/sequencing/sequencing_software/igenome.html) from the UCSC source. Considering the size of hg38 and the resource intensive process of building an index with HISAT2, HISAT2 provides an [index for hg38](https://ccb.jhu.edu/software/hisat2/index.shtml). We renamed the hg38 index and fasta files to distinguish them from the pseudovirome, so once the hg38 reference files are downloaded rename them to “hg38_genome.fa” and “hg38_genome_index” otherwise alignments to hg38 will not work.  The pseudovirome reference was built through combining the sequences of 3,152 virus sequences from an array of virus families and isolates - as was the case with Ebolavirus totaling 1,981 isolate entries of our pseudovirome reference. The pseudovirome reference was built through combining the sequences of 3,152 virus sequences from an array of virus families and isolates - as was the case with Ebolavirus totaling 1,981 isolate entries of our pseudovirome reference. Virus sequences were collected from the [Viral Genomes Resource](https://www.ncbi.nlm.nih.gov/genome/viruses/), [PaVE](https://pave.niaid.nih.gov/), and the [NCI](https://home.ccr.cancer.gov/lco/PyVE.asp). The pseudovirome reference will be updated routinely. 

## Running the pipeline

### Building the index

Generating the index uses three steps: creating two text files containing a list of splice sites and exons, respectively, and the genomic indexed files. The way the script is currently formatted, it will build the list and index files subsequently. 

To initiate index build:

`docker run -v /path/to/your/files/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/ref:/root/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/ref yourdockerimage /home/Index_Generation.sh /gtf`

While this script will build the splice site and exon files for both hg38 and the pseudovirome, this script will only build the pseudovirome index. As mentioned above, HISAT2 has pre-built indexes one can download. For full disclosure, our hg38 index build was executed in a computing cluster with the following command: 

 `hisat2-build --ss hg38_ss.txt --exon hg38_exons.txt -f hg38_genome.fa hg38_genome_index`
 

### Pre-alignment processing

Prior to the hg38 alignment, one should check on the quality of their reads and remove low quality reads. We have done this using FastQC and Trimmomatic. This step of the pipeline will check on the quality of the reads before and after trimming as well as remove/trim reads that fall below a phred score of 30 and a length of 75bp. User must specify which adapter file they wish to use prior to running this script. These variables can be changed within the body of the script, if desired. 

To initiate this pre-alignment step, create an empty text file to save your trimming statistics to and then run:

`docker run -v /path/to/your/files/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/input:/root/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/input yourdockerimage /home/Trimming_protocol.sh /fastq 2>&1 | tee -a /path/to/your/files/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/input/Trimming_output.txt`

### hg38 alignment to remove human-specific reads

This is the recommended next step of our pipeline. The purpose of this step is to remove any human-specific sequences. This step will also remove endogenous viral sequences found within the human reference, so this pipeline is focused on detecting the presence of exogenous viral sequences within a given database. The output of this step provides the alignment output to the human reference (.sam and unaligned reads). Unfortunately, this step cannot be run from within the Docker application due to what appears to be how memory is allocated between Docker, HISAT2, and hg38. When attempted, the process is killed with an unspecific error message (value 137). Therefore, to align trimmed reads to the human reference, set a locally downloaded HISAT2 to your path, make an empty file called "hg38_alignment_output.txt" and run the following command from within the “input” folder:

`/path/to/your/files/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/scripts/hg38_alignment_first_pass.sh 2>&1 | tee -a /path/to/your/files/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/input/hg38_alignment_output.txt`

This will align your reads to hg38 while saving the alignment efficiency results to a local file. 

### Aligning unmapped reads to the pseudovirome

This step picks up immediately after the hg38 alignment. The input for this step is the unmapped fastq files from the hg38 alignment step. The output of this step will yield a .sam file to the pseudovirome, and three files containing unmapped reads (two of which contained paired mates and one contains unpaired reads). The .sam output from this step will be used to determine the presence and expression of virus reads. 

Prior to alignment, make sure that your unmapped reads are properly paired by running fastq-pair. To do this, set a locally downloaded fastq-pair to your path, generate an empty text file called "fastqpair_step.txt" to save the pairing results to, and run the following commands from within the "input" folder:

`for FILENAME in $(ls *HISAT2_hg38_HML2_alignment.sam |rev | cut -c 32- | rev | uniq); do fastq_pair ${FILENAME}_HISAT2_hg38_HML2_alignment_un-conc-mate.1 ${FILENAME}_HISAT2_hg38_HML2_alignment_un-conc-mate.2; done 2>&1 | tee -a /path/to/your/files/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/input/fastqpair_step.txt`

`for FILENAME in $(ls *HISAT2_hg38_HML2_alignment.sam |rev | cut -c 32- | rev | uniq); do cat ${FILENAME}_HISAT2_hg38_HML2_alignment_un-conc-mate.1.single.fq ${FILENAME}_HISAT2_hg38_HML2_alignment_un-conc-mate.2.single.fq ${FILENAME}_HISAT2_hg38_HML2_alignment_un-seqs > ${FILENAME}_HISAT2_hg38_HML2_alignment_unpaired.fq; done`

These commands will run fastq-pair on your unmapped output and combined any unpaired reads into one file per sample.

To align unmapped reads to the pseudovirome, create an empty text file called "pseudovirome_alignment_output.txt" to save your alignment statistics to, then run:

`docker run -v /path/to/your/files/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/input:/root/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/input -v /path/to/your/files/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/ref:/root/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/ref yourdockerimage /home/pseudovirome_alignment_second_pass.sh /fq 2>&1 | tee -a /path/to/your/files/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/input/pseudovirome_alignment_output.txt`

### Post-alignment steps: Filtering the alignment output and calculating the number of virus-specific mapped reads

After alignment to the pseudovirome, one should filter the output to identify viruses of interest to pursue further in a wet lab environment. This step removes PCR duplicates, multimapping reads, and tandem repeats (which could be the result of sequencing error rather than proper virus expression). The final two steps of this file indexes the alignment and filtering output as well as calculates the number of virus-specific reads using featureCounts. 

To filter the alignment output and calculate read counts:

`docker run -v /path/to/your/files/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/input:/root/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/input -v /path/to/your/files/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/ref:/root/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/ref yourdockerimage /home/PostAlignment_hg38ToPseudo.sh /sam`

###  Calculating virus counts: Merge_Virus_Counts.R

One can combine the counts output from the Post-Alignment filtering into one file. Therefore, if one wanted to see which virus families are present – and to what degree of expression - in x number of samples after completing the filtration process this script will combine the counts output from those x samples into one .csv. This script can be modified so the user can look at how the different steps of the filtration process can affect read counts (ex. removing PCR duplicates vs not removing PCR duplicates). 

As of now, our application cannot run R. We plan on fixing this in future updates, but for now, one could follow [this](https://ropenscilabs.github.io/r-docker-tutorial/02-Launching-Docker.html) tutorial for our R scripts. 

###  Identification of Interesting ORFs: ORF_Coverage.R

Viruses can affect host health in a multitude of ways, one of which being through viral protein expression. One way to screen for interesting viruses that pop up in an alignment for further investigation is to look for viruses that are well covered by reads. This R script will report coverage of a virus of interest (user must specify within the script) as well as report coverage of an ORF of interest (user must specify within the script). The input requires a .bam file, and the name of your favorite virus (Genbank accession number), which will generate a list of ORFs and their coverage (%). This script will also graph the depth and coverage over your virus of interest. 

### Counting endogenous virus reads

To calculate expression of endogenous retroviruses, we added one new script: PostAlignment_ERV.sh. This will filter the hg38 alignment output by removing duplicates and isolating unique reads, index, and count the number of hg38-mapped reads using featureCounts. Since this script relies heavily upon the hg38 annotation, your ERVs of interest should either be annotated in the standard or provided hg38.gtf OR should be added to the annotation file prior to index generation. 

To filter the alignment output and calculate read counts:

`docker run -v /path/to/your/files/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/input:/root/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/input -v /path/to/your/files/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/ref:/root/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/ref yourdockerimage /home/PostAlignment_ERV.sh /sam`

### Port to BLAST

To submit mapped reads to BLAST, we added an additional script that extracts mapped reads from the pseudovirome alignment and saves them to a .txt. This .txt can be uploaded to BLAST to search for reads that map to hg38 so they can be removed from the pseudovirome alignment. 

To extract mapped reads and save them:

`docker run -v /path/to/your/files/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/input:/root/Formidable-Exploration-of-humaN-Reads-for-vIRal-sequences/input  yourdockerimage /home/BLAST_prep.sh /sam`
