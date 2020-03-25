#Variables
Input_Text="*_featureCounts_virome_sorted.txt" #change this to reflect what input file you're using (i.e. hg38 or pseudovirome)
WD="<set this>"

setwd(WD)
getwd()

filelist = list.files(pattern = Input_Text)
filelist

datalist = lapply(filelist, read.csv, header=FALSE, sep ="\t", stringsAsFactors=FALSE)

lst1_pseudo = read.table(filelist[1], header=TRUE, sep="\t")[,1:6]
head(lst1_pseudo)
lst2_pseudo = lapply(datalist, "[", c("V7"))
head(lst2_pseudo)

x <- do.call(data.frame,lst2_pseudo)
x <- x[-1,]
colnames(x) = x[1,]
x = x[-1,]
colnames(x)

tmp <- unlist(lapply(strsplit(colnames(x),"/"), function(x) x[[5]]))
colnames(x) <- unlist(lapply(strsplit(tmp, "_"), function(x) x[[1]]))

res <- cbind(lst1_pseudo,x)
head(res)

write.csv(res, file = 'Pseudo_counts.csv') # change output file name for exogenous or endogenous viruses

