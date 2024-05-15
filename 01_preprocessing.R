rm(list = ls())
ay<- gc()
options(scipen = 5)
# load packages ####
sessioninfo::session_info()
#set the working directory
setwd("/Aakash/temp_dir/Mat_Pap")
for (i in c(
  "data.table",
  "kknn",
  "tidyr",
  "ggplot2",
  "dplyr",
  "GEOquery",
  "limma",
  "umap",
  "edgeR",
  "tidyr",
  "here",
  "tibble",
  "car",
  "DESeq2"
)
) {
  suppressPackageStartupMessages(
    library(i, character.only = TRUE
    ))
}

## set root and parallel settings ####
# Knitr should use the project root and not the script location as root
# base.dir refers to the plot location, that should remain with the script
knitr::opts_knit$set(
  root.dir = here()
)

# Give data.table enough threads
writeLines(paste0("Threads available: ", parallel::detectCores()))
writeLines(paste0("Threads given to data.table: ", parallel::detectCores() / 2))
setDTthreads(parallel::detectCores() / 2)

# load the required data ####
# Load the counts file into R
#counts <- read.table("countsall.txt", header = TRUE, row.names = 1, check.names = FALSE)
# Load the counts file into R
counts1 <- read.table("counttweenR2.txt", header = TRUE, row.names = 1, check.names = FALSE)
# Load the counts file into R
counts2 <- read.table("countssaltR2.txt", header = TRUE, row.names = 1, check.names = FALSE)
# Load the counts file into R
counts3 <- read.table("countsoilR2.txt", header = TRUE, row.names = 1, check.names = FALSE)
counts4 <- read.table("countstween.txt", header = TRUE, row.names = 1, check.names = FALSE)
counts5 <- read.table("countssalt.txt", header = TRUE, row.names = 1, check.names = FALSE)
#counts6 <- read.table("countstween.txt", header = TRUE, row.names = 1, check.names = FALSE)
counts6 <- read.table("countsoil.txt", header = TRUE, row.names = 1, check.names = FALSE)


counts1$rn <- rownames(counts1)
counts2$rn <- rownames(counts2)
counts3$rn <- rownames(counts3)
counts4$rn <- rownames(counts4)
counts5$rn <- rownames(counts5)
counts6$rn <- rownames(counts6)

#countsall<-merge(counts1,counts2, by = 'row.names')


counts <- join_all(list(counts1,counts2,counts3,counts4,counts5,counts6), by = 'rn', type = 'full')

#counts <- as.matrix(counts[, -c(2:5)])
counts <- as.matrix(counts[, -c(1:5)])
counts<- counts[,c(2,7,4,6,3,5,1)]
colnames(counts)<- c("geneid","oil_R1","oil_R2","salt_R1","salt_R2","tween_R1","tween_R2")
counts<- as.data.frame(counts)
counts<-counts %>%
  remove_rownames() %>%
  column_to_rownames(var = 'geneid')

fwrite(counts,"all_counts.txt",row.names = TRUE)










