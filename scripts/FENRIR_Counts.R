#Description: Responsible for data processing and graphing
  #1) Loads in needed libraries
  #2) Read in Variables.txt
  #3) Set variables based on Variables.txt
  #4) Merges counts matricies to create a master counts matrix for endogenous and exogenous viruses, respectively
  #5) Graphs endogenous and exogenous counts 
  #6) Calculates percent coverage for endogenous and exogenous viruses
  #7) Graphs coverage for endogenous and exogenous viruses
  #8) Merges individual coverage files to create a master coverage matrix for endogenous and exogenoous viruses, respectively

#Input:
  #1) Individual counts matricies for each sample for each SAM flag for both endogenous and exogenous viruses
  #2) Samtools-generaged coverage file
  #3) .gtf for reference and pseudovirome
  #4) List of all proviruses you're interested in. Otherwise, use the one provided.
  #5) .bam alignment files
  #6) Variables.txt - filled out to match your experimental design

#Output:
  #1) combined counts matrix for endogenous and exogenous viruses detected in RNA-seq, respectively
  #2) graphed viral counts for endogenous and exogenous viruses
  #3) percent coverage matricies for endogenous and exogenous viruses
  #4) graphed coverage of endogenous and exogenous viruses

#Author:
  #Farrah Roy



#FENRIR RNA-Seq pipeline for detection of virus-specific fragments
#R script merges counts files, graphs their expression, calculates coverage, and plots coverage
args = commandArgs(trailingOnly=TRUE)
print(args)

sink(file="FENRIR_R.log", append=TRUE, type="output")
print("load needed libraries")
#needed librariese
library(ggplot2)
library(readr)
library(stringr)
library(plyr)
library(reshape)
library(pheatmap)
library(tidyr)
library(GenomicAlignments)
library(GenomicRanges)
library(Rsamtools)
library(refGenome)
library(data.table)

list.files()

print("read in Variable.txt file to get WD and counts data file extensions")
#read in variable text file

mydata <- read.table((args[1]), colClasses = c("character"))
#extract WD path
as.character(mydata[2,1])
WD=gsub('"', '', regmatches(as.character(mydata[2,1]), gregexpr('"([^"]*)"', as.character(mydata[2,1])))[[1]])
setwd(WD)
getwd()
#extract file extensions for counts data
R_FileList=as.character(mydata[38,1])
R_FileList1=as.character(mydata[39,1])
R_FileList2=as.character(mydata[40,1])
R_FileList3=as.character(mydata[41,1])
R_FileList4=as.character(mydata[42,1])
R_FileList5=as.character(mydata[43,1])
R_FileList6=as.character(mydata[44,1])


File=gsub('"', '', regmatches(R_FileList, gregexpr('"([^"]*)"', R_FileList))[[1]])
File1=gsub('"', '', regmatches(R_FileList1, gregexpr('"([^"]*)"', R_FileList1))[[1]])
File2=gsub('"', '', regmatches(R_FileList2, gregexpr('"([^"]*)"', R_FileList2))[[1]])
File3=gsub('"', '', regmatches(R_FileList3, gregexpr('"([^"]*)"', R_FileList3))[[1]])
File4=gsub('"', '', regmatches(R_FileList4, gregexpr('"([^"]*)"', R_FileList4))[[1]])
File5=gsub('"', '', regmatches(R_FileList5, gregexpr('"([^"]*)"', R_FileList5))[[1]])
File6=gsub('"', '', regmatches(R_FileList6, gregexpr('"([^"]*)"', R_FileList6))[[1]])

file_list=c(File,File1,File2,File3,File4,File5,File6)
file_list

print("finished generating list of counts data file extensions")

print("merging counts matricies for each alignment type")
for(value in file_list){
  #combine counts data for each alignment type
  paste(value)
  #assign SAM flag counts and alignment type to variable
  Input_Text_PE=paste(value, "paired_end_3_counts.txt", sep="_") 
  Input_Text_SE_73=paste(value, "single_end_73_counts.txt", sep="_")
  Input_Text_SE_89=paste(value, "single_end_89_counts.txt", sep="_")
  Input_Text_SE_153=paste(value, "single_end_153_counts.txt", sep="_")
  Input_Text_SE_137=paste(value, "single_end_137_counts.txt", sep="_")
  
  #list all files specified above that ARE NOT empty files. This is done for each SAM flag.
  filelist_PE_3 = list.files(pattern = Input_Text_PE)
  filelist_PE_3
  list.of.files_3 = file.info(filelist_PE_3)
  sizes = list.of.files_3$size
  list.of.non.empty.files_3 = rownames(list.of.files_3)[which(sizes != 0)]
  list.of.non.empty.files_3
  
  filelist_SE_73 = list.files(pattern = Input_Text_SE_73)
  filelist_SE_73
  list.of.files_73 = file.info(filelist_SE_73)
  sizes = list.of.files_73$size
  list.of.non.empty.files_73 = rownames(list.of.files_73)[which(sizes != 0)]
  list.of.non.empty.files_73
  
  filelist_SE_89 = list.files(pattern = Input_Text_SE_89)
  filelist_SE_89
  list.of.files_89 = file.info(filelist_SE_89)
  sizes = list.of.files_89$size
  list.of.non.empty.files_89 = rownames(list.of.files_89)[which(sizes != 0)]
  list.of.non.empty.files_89
  
  filelist_SE_153 = list.files(pattern = Input_Text_SE_153)
  filelist_SE_153
  list.of.files_153 = file.info(filelist_SE_153)
  sizes = list.of.files_153$size
  list.of.non.empty.files_153 = rownames(list.of.files_153)[which(sizes != 0)]
  list.of.non.empty.files_153
  
  filelist_SE_137 = list.files(pattern = Input_Text_SE_137)
  filelist_SE_137
  list.of.files_137 = file.info(filelist_SE_137)
  sizes = list.of.files_137$size
  list.of.non.empty.files_137 = rownames(list.of.files_137)[which(sizes != 0)]
  list.of.non.empty.files_137
  
  #Read in the counts files from the above lists
  datalist_PE_3 = lapply(list.of.non.empty.files_3, read.csv, header=FALSE, sep ="\t", stringsAsFactors=FALSE)
  datalist_SE_73 = lapply(list.of.non.empty.files_73, read.csv, header=FALSE, sep ="\t", stringsAsFactors=FALSE)
  datalist_SE_89 = lapply(list.of.non.empty.files_89, read.csv, header=FALSE, sep ="\t", stringsAsFactors=FALSE)
  datalist_SE_153 = lapply(list.of.non.empty.files_153, read.csv, header=FALSE, sep ="\t", stringsAsFactors=FALSE)
  datalist_SE_137 = lapply(list.of.non.empty.files_137, read.csv, header=FALSE, sep ="\t", stringsAsFactors=FALSE)
  
  #get list of gene names in each count matrix
  gene_names_PE_3 = read.table(list.of.non.empty.files_3[1], header=FALSE, sep="\t")[,1]
  lst2_pseudo_PE_3 = lapply(datalist_PE_3, "[", c("V2")) 
  head(lst2_pseudo_PE_3)
  
  gene_names_SE_73 = read.table(list.of.non.empty.files_73[1], header=FALSE, sep="\t")[,1]
  lst2_pseudo_SE_73 = lapply(datalist_SE_73, "[", c("V2")) 
  head(lst2_pseudo_SE_73)
  
  gene_names_SE_89 = read.table(list.of.non.empty.files_89[1], header=FALSE, sep="\t")[,1]
  lst2_pseudo_SE_89 = lapply(datalist_SE_89, "[", c("V2")) 
  head(lst2_pseudo_SE_89)
  
  gene_names_SE_153 = read.table(list.of.non.empty.files_153[1], header=FALSE, sep="\t")[,1]
  lst2_pseudo_SE_153 = lapply(datalist_SE_153, "[", c("V2")) 
  head(lst2_pseudo_SE_153)
  
  gene_names_SE_137 = read.table(list.of.non.empty.files_137[1], header=FALSE, sep="\t")[,1]
  lst2_pseudo_SE_137 = lapply(datalist_SE_137, "[", c("V2")) 
  head(lst2_pseudo_SE_137)
  
  #merge data
  lst2_pseudo_PE_3 <- do.call(data.frame,lst2_pseudo_PE_3)
  head(lst2_pseudo_PE_3)
  
  lst2_pseudo_SE_73 <- do.call(data.frame,lst2_pseudo_SE_73)
  head(lst2_pseudo_SE_73)
  
  lst2_pseudo_SE_89 <- do.call(data.frame,lst2_pseudo_SE_89)
  head(lst2_pseudo_SE_89)
  
  lst2_pseudo_SE_153 <- do.call(data.frame,lst2_pseudo_SE_153)
  head(lst2_pseudo_SE_153)
  
  lst2_pseudo_SE_137 <- do.call(data.frame,lst2_pseudo_SE_137)
  head(lst2_pseudo_SE_137)
  
  #add sample names to each column
  colnames(lst2_pseudo_PE_3) <- unlist(lapply(strsplit(list.of.non.empty.files_3, "_HISAT2"), function(x) x[[1]])) 
  head(lst2_pseudo_PE_3)
  
  colnames(lst2_pseudo_SE_73) <- unlist(lapply(strsplit(list.of.non.empty.files_73, "_HISAT2"), function(x) x[[1]])) 
  head(lst2_pseudo_SE_73)
  
  colnames(lst2_pseudo_SE_89) <- unlist(lapply(strsplit(list.of.non.empty.files_89, "_HISAT2"), function(x) x[[1]])) 
  head(lst2_pseudo_SE_89)
  
  colnames(lst2_pseudo_SE_153) <- unlist(lapply(strsplit(list.of.non.empty.files_153, "_HISAT2"), function(x) x[[1]])) 
  head(lst2_pseudo_SE_153)
  
  colnames(lst2_pseudo_SE_137) <- unlist(lapply(strsplit(list.of.non.empty.files_137, "_HISAT2"), function(x) x[[1]])) 
  head(lst2_pseudo_SE_137)
  
  #bind gene names to above combined data frame of counts
  res_PE_3 <- cbind(gene_names_PE_3,lst2_pseudo_PE_3)
  head(res_PE_3)
  
  res_SE_73 <- cbind(gene_names_SE_73,lst2_pseudo_SE_73)
  head(res_SE_73)
  
  res_SE_89 <- cbind(gene_names_SE_89,lst2_pseudo_SE_89)
  head(res_SE_89)
  
  res_SE_153 <- cbind(gene_names_SE_153,lst2_pseudo_SE_153)
  head(res_SE_153)
  
  res_SE_137 <- cbind(gene_names_SE_137,lst2_pseudo_SE_137)
  head(res_SE_137)
  
  #write the data frames to .csv
  write.csv(res_SE_73, file = paste(value, '73_single_end_HISAT2.csv', sep="_")) 
  
  write.csv(res_SE_89, file = paste(value, '89_single_end_HISAT2.csv', sep="_")) 
  
  write.csv(res_SE_153, file = paste(value, '153_single_end_HISAT2.csv', sep="_")) 
  
  write.csv(res_SE_137, file = paste(value, '137_single_end_HISAT2.csv', sep="_")) 
  
  write.csv(res_PE_3, file = paste(value, '3_paired_end_HISAT2.csv', sep="_")) 
  
  #combine all the differet counts matrices for each SAM flag into one
  csv1 <- res_PE_3[,-1]
  rownames(csv1) <- res_PE_3[,1]
  csv2 <- res_SE_73[,-1]
  rownames(csv2) <- res_SE_73[,1]
  csv3 <- res_SE_89[,-1]
  rownames(csv3) <- res_SE_89[,1]
  csv4 <- res_SE_137[,-1]
  rownames(csv4) <- res_SE_137[,1]
  csv5 <- res_SE_153[,-1]
  rownames(csv5) <- res_SE_153[,1]
  
  z <- cbind(csv1, csv2, csv3, csv4, csv5)
  res <- t(rowsum(t(z), colnames(z)))
  
  #remove any rows that have a sum of 0 and write it
  greater_than_0 = res[rowSums(res[,-1])>0, ]
  write.csv(greater_than_0, file = paste(value, '.csv', sep="_"))
}

print("finished merging")
print("removing individual SAM flag counts files")

#remove all SAM flag counts files to keep only the merged counts files
junk=dir(path=WD,pattern="*_end_HISAT2.csv")
file.remove(junk)
#junk=dir(path=WD,pattern="*paired_end_3_counts.txt")
#file.remove(junk)
#junk=dir(path=WD,pattern="*single_end_73_counts.txt")
#file.remove(junk)
#junk=dir(path=WD,pattern="*single_end_89_counts.txt")
#file.remove(junk)
#junk=dir(path=WD,pattern="*single_end_153_counts.txt")
#file.remove(junk)
#junk=dir(path=WD,pattern="*single_end_137_counts.txt")
#file.remove(junk)


print("finished removing individual SAM flag counts files")
print("start graphing endogenous virus counts")

#graph counts: endogenous viruses
#load in individual expression files for each tissue and assign to a list
reference_data = read.csv("HISAT2_reference_htseq_sorted_MrkDup_unique_.csv", header=TRUE, sep =",", stringsAsFactors=FALSE)
reference_data = rename(reference_data, c("X"="gene_name"))
head(reference_data)
class(reference_data)

as.character(mydata[45,1])
All_Proviruses=gsub('"', '', regmatches(as.character(mydata[45,1]), gregexpr('"([^"]*)"', as.character(mydata[45,1])))[[1]])
All_Proviruses
#load in a list of proviruses that you're interested in to later assign as row names
gene_names = read.table(All_Proviruses, header=FALSE)
colnames(gene_names) = c("gene_name")
head(gene_names)
dim(gene_names)

#make empty data frame and extract counts data matching the first column
HERV_individual_Counts = as.data.frame(gene_names)
head(HERV_individual_Counts)
HERV_individual_Counts <- reference_data[reference_data$gene_name %in% HERV_individual_Counts$gene_name,]
head(HERV_individual_Counts)
dim(HERV_individual_Counts)

#assign provirus names as the row names and remove the old column containing just provirus names
rownames(HERV_individual_Counts)=HERV_individual_Counts$gene_name
HERV_individual_Counts$gene_name <- NULL
head(HERV_individual_Counts)

#save rownames of HERV individual counts for BLAST input later on
HERV_list = data.frame(stringsAsFactors = FALSE)
for(value in rownames(HERV_individual_Counts)){
  print(value)
  Virus_Identifier=value
  print(Virus_Identifier)
  Chromosome_Identifier=paste("chr", vapply(strsplit(sapply(strsplit(Virus_Identifier,"_"), `[`, 2), "p|q"), `[`, 1, FUN.VALUE=character(1)), sep = "")
  print(Chromosome_Identifier)
  GTF="../ref/hg38.gtf"
  
  #read in hg38.gtf for ORF coordinates
  ens = ensemblGenome()
  read.gtf(ens, GTF)
  class(ens)
  Virus_Annotations = extractByGeneId(ens, Virus_Identifier)#isolates ERVs of interest from gtf
  print(Virus_Annotations)
  Virus_Annotation_df = data.frame(start=getGtf(Virus_Annotations)$start,end=getGtf(Virus_Annotations)$end,gene_name=getGtf(Virus_Annotations)$gene_id )
  print(Virus_Annotation_df)#converts to dataframe and extracts start, end, and name for each extracted ERV
  
  #get coordinates for ERV of interest
  my_ERV = subset(Virus_Annotation_df, gene_name==Virus_Identifier)
  print(my_ERV)
  coordinates=paste(Chromosome_Identifier,":",my_ERV$start,"-",my_ERV$end,sep="")
  print(coordinates)
  entries=paste(value,coordinates,sep=" ")
  print(entries)
  HERV_list = rbind(HERV_list,entries,stringsAsFactors = FALSE)
}
write.csv(HERV_list,file="HERV_list.txt")

#convert all NA designations to 0
HERV_individual_Counts[is.na(HERV_individual_Counts)] <- 0

#save data 
write.csv(HERV_individual_Counts, file = 'HERV_Expression_sorted_MrkDup_unique.csv')

# plot a basic heatmap
paletteLength <- 50
myColor <- colorRampPalette(c("white", "blue"))(paletteLength)
map = pheatmap(HERV_individual_Counts, color = myColor, cluster_rows = FALSE, cluster_cols = FALSE, 
               main = "ERV expression", fontsize = 14, fontsize_row = 4, 
               fontsize_col = 8, breaks = seq(0,8,by=0.2))
mypath <- file.path("HERV_Expression_sorted_MrkDup_unique_heatmap.jpg")
jpeg(file=mypath)
map
dev.off()

print("finished graphing endogenous virus counts")
print("start plotting exogenous virus counts")

#graph counts: exogenous viruses  
Exogenous_Virus_Expression = read.csv("HISAT2_virome_htseq_sorted_MrkDup_unique_noTanRep_.csv", header=TRUE, sep=",")
Exogenous_Virus_Expression = rename(Exogenous_Virus_Expression, c("X"="virus_ID"))
head(Exogenous_Virus_Expression)
Exogenous_Virus_Expression<-Exogenous_Virus_Expression[!(Exogenous_Virus_Expression$virus_ID=="__ambiguous"),]
Exogenous_Virus_Expression<-Exogenous_Virus_Expression[!(Exogenous_Virus_Expression$virus_ID=="__no_feature"),]
Exogenous_Virus_Expression
cat(as.character(Exogenous_Virus_Expression$virus_ID), file = "list_of_viruses.txt", sep ="\n")

# Basic bar plot
Exogenous_Virus_Expression_melted = reshape2::melt(Exogenous_Virus_Expression, id.vars = "virus_ID")
head(Exogenous_Virus_Expression_melted)  
plot = ggplot(Exogenous_Virus_Expression_melted, aes(fill=virus_ID, x=variable, y=value)) + 
  geom_bar(position="stack", stat="identity") + theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Exogenous virus expression") +
  xlab("Sample ID") + ylab("Fragment counts")
mypath <- file.path("Exogenous_Virus_Expression_sorted_MrkDup_unique_noTanRep_stacked_barplot.jpg")
jpeg(file=mypath)
plot
dev.off()

print("finished plotting exogenous virus counts")
print("calculate coverage for ERVs")

#coverage for ERVs
#get a list of bam files and proviruses detected in RNA-seq data
Reference_bam_list = list.files(pattern = "HISAT2_reference_alignment_sorted_MrkDup_unique.bam$")
Reference_bam_list
Detected_Proviruses=rownames(HERV_individual_Counts)
Detected_Proviruses
#for each bam file in character vector
for(value in Reference_bam_list){
  print(value)
  Bam_File=value
  Bam_File
  Sample_ID=sapply(strsplit(Bam_File,"_"), `[`, 1)
  Sample_ID
  #for each provirus id in counts matrix
  for(value in Detected_Proviruses){
    print(value)
    Virus_Identifier=value
    Virus_Identifier
    Chromosome_Identifier=paste("chr", vapply(strsplit(sapply(strsplit(Virus_Identifier,"_"), `[`, 2), "p|q"), `[`, 1, FUN.VALUE=character(1)), sep = "")
    Chromosome_Identifier
    peak_threshold=1 #find areas of your favorite virus that are covered equal to or greater than this threshold
    Samtools_Coverage_File=paste(Sample_ID,"HISAT2_reference_alignment_sorted_MrkDup_unique_HML2.coverage",sep="_")
    GTF="../ref/hg38.gtf"
    
    #read in bam 
    gal <- readGAlignments(Bam_File)
    cvg <- coverage(gal)
    cvg
    
    #read in hg38.gtf for ORF coordinates
    ens = ensemblGenome()
    read.gtf(ens, GTF)
    class(ens)
    Virus_Annotations = extractByGeneId(ens, Virus_Identifier)#isolates ERVs of interest from gtf
    Virus_Annotations
    Virus_Annotation_df = data.frame(start=getGtf(Virus_Annotations)$start,end=getGtf(Virus_Annotations)$end,gene_name=getGtf(Virus_Annotations)$gene_id )
    Virus_Annotation_df#converts to dataframe and extracts start, end, and name for each extracted ERV
    
    #get coordinates for ERV of interest
    my_ERV = subset(Virus_Annotation_df, gene_name==Virus_Identifier)
    my_ERV
    start = my_ERV[, "start"]#assign for GRanges
    end = my_ERV[,"end"]
    start
    end
    ranges = GRanges(Chromosome_Identifier, IRanges(start,end))
    ranges
    grl=GRangesList(Chromosome_Identifier=ranges)
    names(grl)=Chromosome_Identifier
    grl
    
    #calculate coverage for whole ERV
    rangeCoverage <- function(ranges, cvg)
    {
      ranges <- reduce(ranges)
      unlist(revElements(cvg[ranges], strand(ranges) == "-"), use.names=FALSE)
    }
    Virus_Annotation = RleList(lapply(grl, rangeCoverage, cvg))
    Virus_Annotation
    Total_Bases <- sum(runLength(Virus_Annotation))#saves the length of your favorite ORF as a variable
    Total_Bases
    
    #converts the data from the RleList containing total virus coverage into a data frame
    x <- runLength(Virus_Annotation)
    x
    y <- runValue(Virus_Annotation)
    y
    Coverage_Total_Virus_df <- data.frame(x,y)
    Coverage_Total_Virus_df
    colnames(Coverage_Total_Virus_df)
    Coverage_Total_Virus_df[,c("group", "group_name", "group.1", "group_name.1")] <- list(NULL)
    Coverage_Total_Virus_df
    colnames(Coverage_Total_Virus_df) <- c("Lengths","Values")
    Coverage_Total_Virus_df
    
    #selects portions of your favorite virus that are covered by at least one read and saves this data to a variable
    Total_Virus_Reads_subdf <- subset(Coverage_Total_Virus_df, Values > 0)
    Total_Virus_Reads_subdf
    Full_Length_Reads_Bases <- sum(Total_Virus_Reads_subdf$Lengths)#calcuates # of covered bases
    Full_Length_Reads_Bases
    
    #calculates coverage of your favorite virus as a percentage
    percentage=(Full_Length_Reads_Bases/Total_Bases)*100
    Virus_Coverage_Percentage <- as.data.frame(percentage)
    Virus_Coverage_Percentage$provirus_name = Virus_Identifier
    head(Virus_Coverage_Percentage)
    write.csv(Virus_Coverage_Percentage, file=paste(Sample_ID,Virus_Identifier,"ERV_coverage.csv", sep="_"))
    
    
    ############################ graphing read depth from samtools output
    
    Samtools_Depth = fread(Samtools_Coverage_File, header=FALSE, sep="\t")
    head(Samtools_Depth)
    tail(Samtools_Depth)
    Samtools_Depth_df = data.frame(Samtools_Depth)
    colnames(Samtools_Depth_df) = c("Gene_ID", "locus", "depth")
    head(Samtools_Depth_df)
    ERV_Samtools_Depth = subset(Samtools_Depth_df, locus >= start & locus <= end)
    #if ERV_Samtools_Depth is empty skip the rest, if not then proceed
    if(dim(ERV_Samtools_Depth)[1] == 0){
      print(paste("The provirus",Virus_Identifier,"was not detected in",Sample_ID))
    }else{
      tail(ERV_Samtools_Depth)
      head(ERV_Samtools_Depth) 
      mypath = file.path(paste(Sample_ID,Virus_Identifier,"genome_coverage.jpg", sep="_"))
      jpeg(file=mypath)
      plot(ERV_Samtools_Depth$locus, ERV_Samtools_Depth$depth, type ="l", 
           main = paste("depth over locus", Virus_Identifier, "for sample", Sample_ID, sep = " "), 
           xlab = "locus", ylab = "depth", pch=1)
      dev.off()
    }
  }
}

print("finished calculating coverage for ERVs")
print("start merging ERV coverage data into one file")
#merge all ERV coverage plots into one csv by sample ID, then delete the multiple *ERV_coverage.csv files
#get list of sample IDs by making a list of total coverage files
ERV_Coverage_Files = list.files(pattern = "_ERV_coverage.csv")
ERV_Coverage_Files
ERV_Sample_ID=unique(sapply(strsplit(ERV_Coverage_Files,"_"), `[`, 1))
ERV_Sample_ID

#combine covereage files based on sample ID
for(value in ERV_Sample_ID){
  print(value)
  #list all coverage files for sample ID
  ERV_Coverage_Files_List = intersect(list.files(pattern = value), list.files(pattern = "_ERV_coverage.csv$"))
  ERV_Coverage_Files_List
  
  #Read in the counts files from the above lists
  ERV_Coverage_datalist = lapply(ERV_Coverage_Files_List, read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
  head(ERV_Coverage_datalist)
  
  #merge data
  ERV_Coverage_df <- do.call(rbind,ERV_Coverage_datalist)
  colnames(ERV_Coverage_df) = c("chr","percentage","provirus_name")
  rownames(ERV_Coverage_df)=ERV_Coverage_df$provirus_name
  ERV_Coverage_df$provirus_name = NULL
  ERV_Coverage_df
  write.csv(ERV_Coverage_df, file=paste(value,"Total_ERV_coverage_merged.csv", sep="_"))
}

print("finished merging ERV coverage files")
print("remove individual coverage files: ERVs")

#remove individual coverage files
junk=dir(path=WD,pattern="*ERV_coverage.csv")
file.remove(junk)

print("finished removing ERV coverage files")
print("calculate coverage for exogenous viruses")

#coverage for exogenous viruses
#get a list of bam files and proviruses detected in RNA-seq data
Pseudovirome_bam_list = list.files(pattern = "HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep.bam$")
Pseudovirome_bam_list
Detected_Viruses=Exogenous_Virus_Expression$virus_ID
Detected_Viruses
#for each bam file in character vector
for(value in Pseudovirome_bam_list){
  print(value)
  Bam_File=value
  Bam_File
  Sample_ID=sapply(strsplit(Bam_File,"_"), `[`, 1)
  Sample_ID
  #for each virus id in counts matrix
  for(value in Detected_Viruses){
    print(value)
    Virus_Identifier=value
    Virus_Identifier
    
    peak_threshold=1 #find areas of your favorite virus that are covered equal to or greater than this threshold
    Samtools_Coverage_File = paste(Sample_ID,"HISAT2_virome_alignment_sorted_MrkDup_unique_noTanRep.coverage",sep = "_")
    GTF="../ref/Pseudovirome.gtf"
    
    #read in bam and calculate coverage of all viruses in the set sample
    gal <- readGAlignments(Bam_File)
    cvg <- coverage(gal)
    cvg
    
    #calculate coverage of your favorite virus
    cvg[[Virus_Identifier]]
    mean(cvg[[Virus_Identifier]])
    max(cvg[[Virus_Identifier]])
    Total_Virus_Bases <- sum(runLength(cvg[[Virus_Identifier]]))#saves the full length of your favorite virus as a variable
    Total_Virus_Bases
    
    #converts the data from the RleList containing total virus coverage into a data frame
    x <- runLength(cvg[[Virus_Identifier]])
    x
    y <- runValue(cvg[[Virus_Identifier]])
    y
    Coverage_Total_Virus_df <- data.frame(x,y)
    colnames(Coverage_Total_Virus_df) <- c("Lengths","Values")
    Coverage_Total_Virus_df
    
    #selects portions of your favorite virus that are covered by at least one read and saves this data to a variable
    Total_Virus_Reads_subdf <- subset(Coverage_Total_Virus_df, Values > 0)
    Total_Virus_Reads_subdf
    Full_Length_Reads_Bases <- sum(Total_Virus_Reads_subdf$Lengths)#calcuates # of covered bases
    Full_Length_Reads_Bases
    
    #calculates coverage of your favorite virus as a percentage
    Total_Virus_Coverage_Percentage <- (Full_Length_Reads_Bases/Total_Virus_Bases)*100
    Total_Virus_Coverage_Percentage
    
    #read in pseudovirome.gtf for ORF coordinates
    ens = ensemblGenome()
    read.gtf(ens, GTF)
    class(ens)
    ens
    Virus_Annotations = extractSeqids(ens, Virus_Identifier)
    Virus_Annotations
    Virus_Annotation_df = data.frame(start=getGtf(Virus_Annotations)$start,end=getGtf(Virus_Annotations)$end, gene_id=getGtf(Virus_Annotations)$gene_id )
    Virus_Annotation_df
    
    ############## create an empty data frame with 0 column and 1 row
    temp_df <- data.frame(matrix(ncol=0,nrow=1))
    
    #calculates coverage over your favorite ORF
    for (row in 1:nrow(Virus_Annotation_df)) {
      start = Virus_Annotation_df[row, "start"]
      end = Virus_Annotation_df[row,"end"]
      start
      end
      
      ranges = IRanges(start,end) 
      gr=GRanges(Virus_Identifier, ranges)
      grl=GRangesList(Virus_Identifier = gr)
      names(grl)=Virus_Identifier
      gr
      grl
      rangeCoverage <- function(gr, cvg)
      {
        gr <- reduce(gr)
        unlist(revElements(cvg[gr], strand(gr) == "-"), use.names=FALSE)
      }
      Virus_ORF_Annotation = RleList(lapply(grl, rangeCoverage, cvg))
      Virus_ORF_Annotation
      Total_ORF_Bases <- sum(runLength(Virus_ORF_Annotation))#saves the length of your favorite ORF as a variable
      Total_ORF_Bases
      
      #converts the data from the RleList containing ORF coverage into a data frame
      x <- runLength(Virus_ORF_Annotation)
      x
      y <- runValue(Virus_ORF_Annotation)
      y
      Coverage_ORF_df <- data.frame(x,y)
      Coverage_ORF_df
      
      ############################### try this line to replace the above codes to remove empty columns in Coverage_ORF_df
      Coverage_ORF_df <- Coverage_ORF_df[,-c(1,2,4,5)]
      colnames(Coverage_ORF_df) <- c("Lengths","Values")
      Coverage_ORF_df
      
      #selects portions of your favorite ORF that are covered by at least one read and saves this data to a variable
      Reads_ORF_subdf <- subset(Coverage_ORF_df, Values > 0)
      Reads_ORF_subdf
      Reads_ORF_Bases <- sum(Reads_ORF_subdf$Lengths)#calcuates # of covered bases
      Reads_ORF_Bases
      
      #calculates coverage of your favorite ORF as a percentage
      ORF_Coverage_Percentage <- (Reads_ORF_Bases/Total_ORF_Bases)*100
      ORF_Coverage_Percentage
      ############################# ORF_Coverage_Percentage conbinds temp_df (line 66) by column.
      temp_df <- cbind(temp_df,ORF_Coverage_Percentage)
    }
    
    ############################# rename the column of temp_df before writting out.
    colnames(temp_df) <- Virus_Annotation_df[,3]
    temp_df$Total_Virus_Coverage_Percentage=Total_Virus_Coverage_Percentage
    temp_df$SampleID=Sample_ID
    head(temp_df)
    
    # write out coverage for all orfs
    str <- paste(Virus_Identifier,Sample_ID,"ORF_Coverage_output.csv",sep = "_")
    str
    write.table(temp_df, file=str,  sep=",", append=F)
    
    ############################ graphing read depth from samtools output
    
    Samtools_Depth = read.csv(Samtools_Coverage_File, header=FALSE, sep="\t")
    head(Samtools_Depth)
    tail(Samtools_Depth)
    Virus_Samtools_Depth = subset(Samtools_Depth, V1 == Virus_Identifier)
    tail(Virus_Samtools_Depth)
    colnames(Virus_Samtools_Depth)= c("Virus", "locus", "depth")
    head(Virus_Samtools_Depth)
    #if Samtools_Depth is empty skip the rest, if not then proceed
    if(dim(Virus_Samtools_Depth)[1] == 0){
      print(paste("The virus",Virus_Identifier,"was not detected in",Sample_ID))
    }else{
      mypath = file.path(paste(Sample_ID,Virus_Identifier,"genome_coverage.jpg", sep="_"))
      jpeg(file=mypath)
      plot(Virus_Samtools_Depth$locus, Virus_Samtools_Depth$depth, type ="l", 
           main = paste(Virus_Identifier, "depth by locus for", Sample_ID,sep=" "), xlab = "locus", ylab = "depth", pch=16)
      dev.off()
    }
  }
}

print("finished calculating coverage for exogenous viruses")
print("merge all coverage files into one: exogenous viruses")

#merge all Exogenous virus coverage plots into one csv by sample ID, then delete the multiple *ORF_Coverage_output.csv files
#get list of sample IDs by making a list of total coverage files
Exo_Coverage_Files = list.files(pattern = "ORF_Coverage_output.csv")
Exo_Coverage_Files
Exo_Virus_ID=paste(unique(sapply(strsplit(Exo_Coverage_Files,"_"), `[`, 1)),unique(sapply(strsplit(Exo_Coverage_Files,"_"), `[`, 2)), sep="_")
Exo_Virus_ID

#combine covereage files based on sample ID
for(value in Exo_Virus_ID){
  print(value)
  #list all coverage files for sample ID
  Exo_Coverage_Files_List = intersect(list.files(pattern = value), list.files(pattern = "ORF_Coverage_output.csv$"))
  Exo_Coverage_Files_List
  
  #Read in the counts files from the above lists
  Exo_Coverage_datalist = lapply(Exo_Coverage_Files_List, read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
  head(Exo_Coverage_datalist)
  
  #merge data
  Exo_Coverage_df <- do.call(rbind,Exo_Coverage_datalist)
  rownames(Exo_Coverage_df)=Exo_Coverage_df$SampleID
  Exo_Coverage_df$SampleID = NULL
  Exo_Coverage_df
  write.csv(Exo_Coverage_df, file=paste(value,"Total_Exogenous_Virus_coverage_merged.csv", sep="_"))
}

print("finished merging coverage data for exogenous viruses")
print("removing individual coverage files")

#remove individual coverage files
junk=dir(path=WD,pattern="*ORF_Coverage_output.csv")
file.remove(junk)

print("finished removing individual coverage files")
print("ready for BLAST")


sink()     
