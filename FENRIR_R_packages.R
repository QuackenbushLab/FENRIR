print("Installing necessary libraries")

install.packages("reshape")
print("y")
print("y")
install.packages("pheatmap")
install.packages("data.table")
install.packages("tidyverse")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomicAlignments")

# Download refgenome package tarball from CRAN archive and install
url <- "https://cran.r-project.org/src/contrib/Archive/refGenome/refGenome_1.7.3.tar.gz"
pkgFile <- "refGenome_1.7.3.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages("RSQLite")
install.packages(pkgs=pkgFile, type="source", repos=NULL)
unlink(pkgFile)

url1="https://cran.r-project.org/src/contrib/Archive/doBy/doBy_4.6-1.tar.gz"
pkgFile1 <- "doBy_4.6-1.tar.gz"
download.file(url = url1, destfile = pkgFile1)
install.packages(pkgs=pkgFile1, type="source", repos=NULL)
unlink(pkgFile1)

