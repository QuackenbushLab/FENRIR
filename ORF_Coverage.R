#Variables and libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomicAlignments", version = "3.8")
BiocManager::install("GenomicRanges", version = "3.8")
BiocManager::install("Rsamtools", version = "3.8")
install.packages('refGenome', dependencies=TRUE, repos='http://cran.r-project.org/')
library(GenomicAlignments)
library(GenomicRanges)
library(Rsamtools)
library(refGenome)
WD="~/input"
Bam_File="IPM0904_S20_HISAT2_virome_alignment_first_pass_Cufflinks_sorted_MrkDup_unique_noTanRep.bam"
Sample_ID="IPM0904_S20"
Virus_Identifier="NC_006273.2"
Virus_Identifier
Virus_Name="CMV"
Virus_Name
peak_threshold=1 #find areas of your favorite virus that are covered equal to or greater than this threshold
setwd(WD)
getwd()

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
Virus_Coverage_Percentage <- (Full_Length_Reads_Bases/Total_Virus_Bases)*100
Virus_Coverage_Percentage

#looks for regions of your favorite virus that are covered to the set threshold
#Virus_coverage <- slice(cvg[[Virus_Identifier]], peak_threshold)
#Virus_coverage

#read in pseudovirome.gtf for ORF coordinates
ens = ensemblGenome()
read.gtf(ens, "Pseudovirome.gtf")
class(ens)
Virus_Annotations = extractSeqids(ens, Virus_Identifier)
Virus_Annotations
Virus_Annotation_df = data.frame(start=getGtf(Virus_Annotations)$start,end=getGtf(Virus_Annotations)$end,gene_name=getGtf(Virus_Annotations)$gene_name )
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
  grl=GRangesList(Virus_Name = gr)
  names(grl)=Virus_Name
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

# write out coverage for all orfs
str <- paste(Virus_Name,Sample_ID,"ORF_Coverage_output.csv",sep = "_")
str
write.table(temp_df, file=str,  sep=",", append=F)

#output ORFs to blast so the end user can see what their function is?
#output ORFs to blastn/blastx so the end user can compare mapped reads to nr and human? (it's own script?)