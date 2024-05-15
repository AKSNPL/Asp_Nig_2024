rm(list = ls())
ay<- gc()
options(scipen = 5)
# load packages ####
sessioninfo::session_info()
#set the working directory
setwd("/Aakash/temp_dir/Mat_Pap")
requiredPackages <- c("pacman")
for (package in requiredPackages) { #Installs packages if not yet installed
  if (!requireNamespace(package, quietly = TRUE))
    install.packages(package)
}

for (i in c(
  "data.table",
  "kknn",
  "tidyr",
  "ggplot2",
  "dplyr",
  "GEOquery",
  "limma",
  "umap",
  "tidyr",
  "here",
  "tibble",
  "car",
  "DESeq2",
  "readxl",
  "clusterProfiler",
  "biomaRt",
  "pheatmap",
  "AnnotationHub"
)
) {
  suppressPackageStartupMessages(
    pacman::p_load(i, character.only = TRUE
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

# Oil as Treatment vs Salt as Control ####
#load the annotation file:
annotation <- read.table("annotation.txt", header = TRUE, row.names = 1, check.names = FALSE)

count1s <- read.table("all_counts.txt", header = TRUE, row.names = 1, check.names = FALSE,sep = ",")
#make the counts file ready:
countsOvS<-count1s[,c(1,2,3,4)]

genes_all <- as.data.frame(read_excel("DAVID_LOCUS_TAG.Aspergillus_niger.xlsx"))
colnames(genes_all)<- c("entrez_ID","From","To")
rownames(genes_all)<- genes_all$From
genes_all$From<- NULL
merged_test<- merge(countsOvS,genes_all,by=0)
rownames(merged_test) <- merged_test$To

merged_test<- merged_test %>%
  dplyr::select(2,3,4,5)

countsOvS<- merged_test
#check the dimensions
dim(countsOvS)
head(countsOvS)


## samples information
head(annotation)
mat_fun <- function(m){
  m2 <- apply(m,  2,  function(x) as.numeric(paste(x)))
  colnames(m2) <- colnames(m)
  rownames(m2) <- rownames(m)
  return(m2)
}
countsOvS<-mat_fun(countsOvS)

class(countsOvS)
## transform counts into matrix
countsOvS.matrix <-as.matrix(countsOvS)


countsOvS0 <- countsOvS[, c("oil_R1", "oil_R2", "salt_R1", "salt_R2")]

all(colnames(countsOvS0) == rownames(annotation))


## group membership for all samples
gsms <- "1100"
sml <- strsplit(gsms, split="")[[1]]



## log2 transformation
ex <- countsOvS0
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
countsOvS0 <- log2(ex) }

## assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("Oil_Treatment","Salt_Control"))
levels(gs) <- groups
design <- model.matrix(~group + 0, data.frame(group = gs))
colnames(design) <- levels(gs)

## remove missings
countsOvS0 <- countsOvS0[complete.cases(countsOvS0), ] 

## fit a linear model:
fit <- lmFit(countsOvS0, design)

# set up contrasts of interest and recalculate model coefficients
cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.05)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
tT <- subset(tT, select=c("adj.P.Val","P.Value","t","B","logFC","AveExpr"))
tT_Ordered <- tT[order(tT$P.Value),]
tT_Sig <- subset(tT_Ordered, adj.P.Val < 0.05)

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=0.5)

dtf<- as.data.frame((dT))
## ############Deseq2#####################------------------------------------########
countsOvS.matrix <-as.matrix(countsOvS)
str(countsOvS.matrix)

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countsOvS.matrix, colData = annotation, design = ~ condition)

# keep genes with at least 10 reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# differentially expressed genes (and other calculations)
dds <- DESeq(dds)
res <- results(dds,alpha=0.05)

#show the results
res

#summarize some basic tallies using the summary function
summary(res)

# order the results table by the smallest p-value
resOrdered <- res[order(res$pvalue),]
resOrdered <- as.data.frame(resOrdered)

# export only the results which pass an adjusted p-value threshold
resSig <- subset(resOrdered, padj < 0.05)

## merged significant genes of limma and deseq2##############################################################
merged_df<- merge(resSig,tT_Sig,by=0)
rownames(merged_df)<- merged_df$Row.names
merged_df$Row.names<- NULL
merged_dtf<- merge(merged_df,dtf, by =0)


# select 100 genes with largest log fold change scores of limma:
select <- order(abs(merged_dtf$logFC), decreasing=TRUE)[1:100]

library(pheatmap)
reso <- 1200
length <- 100*reso/72
png("02_heatmaps/OvS_heatmap.png",width     = length, height    = length)
# plot simple heatmap with clustering for merged data
pheatmap(countsOvS0[select,], cluster_rows=TRUE, cluster_cols=TRUE,
         show_rownames=TRUE, show_colnames=TRUE, fontsize_col=8,annotation_col=annotation)
dev.off()

# Build histogram of P-values for all genes. Normal test
png("03_histograms/OvS_histogram.png")
# assumption is that most genes are not differentially expressed.
hist(merged_df$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")
dev.off()

# summarize test results as "up", "down" or "not expressed"
rownames(merged_dtf)<- merged_dtf$Row.names
merged_dtf$Row.names<- NULL

dt_up<- merged_dtf %>%
  filter(`Oil_Treatment-Salt_Control`=="1")

dt_up_row<- rownames(dt_up)

write.table(dt_up_row, "01_up_down_genes/oilT_saltC_upreg.txt",row.names = F, quote = F)

dt_down<- merged_dtf %>%
  filter(`Oil_Treatment-Salt_Control`=="-1")

dt_down_row<- rownames(dt_down)

write.table(dt_down_row, "01_up_down_genes/oilT_saltC_downreg.txt",row.names = F, quote = F)

dt_no_sig<- merged_dtf %>%
  filter(`Oil_Treatment-Salt_Control`=="0")

png("04_venndiagram/OvS_venndiagram.png")
# Venndiagram showing expressed and non expressed genes:
vennDiagram(merged_dtf$`Oil_Treatment-Salt_Control`, circle.col=palette(),cex=c(4,0.5,0.5),names = c("Expressed", "Not expressed"),main="Comparison of Oil Treatment and Salt Control")
dev.off()

# add a column of NAs
merged_dtf$diffexpressed <- "NO"

# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
merged_dtf$diffexpressed[merged_dtf$logFC > 0.5 & merged_dtf$P.Value < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
merged_dtf$diffexpressed[merged_dtf$logFC < -0.5 & merged_dtf$P.Value < 0.05] <- "DOWN"

# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
merged_dtf$delabel <- NA
merged_dtf$delabel[merged_dtf$diffexpressed != "NO"] <- row.names(merged_dtf)[merged_dtf$diffexpressed != "NO"]


# Assuming 'merged_dtf' contains your data
top_20 <- merged_dtf[order(merged_dtf$logFC),][1:20,]

png("05_volcano_plots/OvS_volcanoplots.png",width     = 800, height    = 800)
#ggplot showing volcano plot
ggplot(data=merged_dtf, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label=ifelse(rank(-abs(logFC)) <= 20, delabel, ""))) + 
  geom_point() + 
  theme_minimal() +
  geom_text(size = 2.5) +
  scale_color_manual(values = c("DOWN" = "blue", "UP" = "red")) 
dev.off()

## KEGG ENRICHMENT ###############################################################
gene_up <- rownames(merged_dtf %>%
                   filter(diffexpressed=="UP"))


kk <- enrichKEGG(gene         = unique(gene_up),
                 organism     = 'ang')

kk_df<- as.data.frame(kk)

# Replace the text in the specified column
kk@result[["Description"]] <- gsub(" - Aspergillus niger.*", "", kk@result[["Description"]])

png("06_kegg_enrich_plots/OvS_dotplots.png")
dotplot(kk, showCategory=1000) + ggtitle("dotplot ")
dev.off()

gene_down <- rownames(merged_dtf %>%
                   filter(diffexpressed=="DOWN"))


kk <- enrichKEGG(gene         = gene_down,
                 organism     = 'ang')

kk_df<- as.data.frame(kk)

# Replace the text in the specified column
kk@result[["Description"]] <- gsub(" - Aspergillus niger.*", "", kk@result[["Description"]])

png("06_kegg_enrich_plots/OvS_dotplots2.png")
dotplot(kk, showCategory=1000) + ggtitle("dotplot ")
dev.off()

## functional annotations ###########################################################
library(biomaRt)
ensembl_fungi <- useMart(host="https://fungi.ensembl.org", 
                         biomart="fungi_mart", 
                         port = 443)

head(listDatasets(ensembl_fungi))

searchDatasets(mart = ensembl_fungi, pattern = "niger")


ensembl <- useDataset(dataset = "aniger_eg_gene", mart = ensembl_fungi)

aay<-listAttributes(ensembl)

ensembl_ids <- getBM(attributes = c("ensembl_gene_id","start_position","end_position", "external_gene_name"),
                     filters = "ensembl_gene_id",
                     values = gene_up,
                     mart = ensembl)


###################################################
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.16")
#install.packages("devtools")
#devtools::install_version("dbplyr", version = "2.3.4")
#BiocManager::install("BiocFileCache", version = "3.15")

library(AnnotationHub) 
library(AnnotationDbi)
NigerDb <- loadDb("filenew1.sqlite")
keytypes(NigerDb)
str(NigerDb)
egid <- head(keys(NigerDb, "PMID"))

genes_all <- as.data.frame(read_excel("DAVID_LOCUS_TAG.Aspergillus_niger.xlsx"))
colnames(genes_all)<- c("entrez_ID","From","To")
rownames(genes_all)<- genes_all$To
genes_all$To<- NULL

gene_up_df<- data.frame(gene_up)
rownames(gene_up_df)<- gene_up_df$gene_up
#gene_up_df$gene_up<- NULL
merged_test<- merge(gene_up_df,genes_all, by=0)
rownames(merged_test) <- merged_test$To
merged_test$gene_up<- NULL


genes_to_test<- merged_test$entrez_ID

# Perform Gene Ontology (GO) enrichment analysis for Biological Process (BP) terms using "enrichGO" function.
GO_results <- enrichGO(gene = genes_to_test, OrgDb = NigerDb, keyType = "ENTREZID", ont = "BP")
as.data.frame(GO_results)

png("07_GO_plots/OvS_barplot_BP.png")
# Create a barplot visualizing the top 15 enriched GO terms for Biological Process.
fit <- plot(barplot(GO_results)) #showCategory = 15))
dev.off()


# Perform GO enrichment analysis for Molecular Function (MF) terms using "enrichGO" function.
GO_results2 <- enrichGO(gene = genes_to_test, OrgDb = NigerDb, keyType = "ENTREZID", ont = "MF")
as.data.frame(GO_results2)

png("07_GO_plots/OvS_barplot_MF.png")
# Create a barplot visualizing the top 15 enriched GO terms for Molecular Function.
fit2 <- plot(barplot(GO_results2))
dev.off()

# Perform GO enrichment analysis for Cellular Component (CC) terms using "enrichGO" function.
GO_results3 <- enrichGO(gene = genes_to_test, OrgDb = NigerDb, keyType = "ENTREZID", ont = "CC")
as.data.frame(GO_results3)

png("07_GO_plots/OvS_barplot_CC.png")
# Create a barplot visualizing the top 15 enriched GO terms for Cellular Component.
fit3 <- plot(barplot(GO_results3))
dev.off()


## Tween_control_vs_Oil_treatment:                        #### ############## ####

rm(list = ls())
ay<- gc()
options(scipen = 5)
# load packages ####
sessioninfo::session_info()
#set the working directory
setwd("/Aakash/temp_dir/Mat_Pap")
requiredPackages <- c("pacman")
for (package in requiredPackages) { #Installs packages if not yet installed
  if (!requireNamespace(package, quietly = TRUE))
    install.packages(package)
}

for (i in c(
  "data.table",
  "kknn",
  "tidyr",
  "ggplot2",
  "dplyr",
  "GEOquery",
  "limma",
  "umap",
  "tidyr",
  "here",
  "tibble",
  "car",
  "DESeq2",
  "readxl",
  "clusterProfiler",
  "biomaRt",
  "pheatmap",
  "AnnotationHub"
)
) {
  suppressPackageStartupMessages(
    pacman::p_load(i, character.only = TRUE
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

# Oil as Treatment vs Salt as Control ####
#load the annotation file:
annotation <- read.table("annotationTO.txt", header = TRUE, row.names = 1, check.names = FALSE)

count1s <- read.table("all_counts.txt", header = TRUE, row.names = 1, check.names = FALSE,sep = ",")
#make the counts file ready:
countsOvS<-count1s[,c(1,2,5,6)]

genes_all <- as.data.frame(read_excel("DAVID_LOCUS_TAG.Aspergillus_niger.xlsx"))
colnames(genes_all)<- c("entrez_ID","From","To")
rownames(genes_all)<- genes_all$From
genes_all$From<- NULL
merged_test<- merge(countsOvS,genes_all,by=0)
rownames(merged_test) <- merged_test$To

merged_test<- merged_test %>%
  dplyr::select(2,3,4,5)

countsOvS<- merged_test
#check the dimensions
dim(countsOvS)
head(countsOvS)


## samples information
head(annotation)
mat_fun <- function(m){
  m2 <- apply(m,  2,  function(x) as.numeric(paste(x)))
  colnames(m2) <- colnames(m)
  rownames(m2) <- rownames(m)
  return(m2)
}
countsOvS<-mat_fun(countsOvS)

class(countsOvS)
## transform counts into matrix
countsOvS.matrix <-as.matrix(countsOvS)


countsOvS0 <- countsOvS[, c("oil_R1", "oil_R2", "tween_R1", "tween_R2")]

all(colnames(countsOvS0) == rownames(annotation))


## group membership for all samples
gsms <- "1100"
sml <- strsplit(gsms, split="")[[1]]



## log2 transformation
ex <- countsOvS0
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
countsOvS0 <- log2(ex) }

## assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("Oil_Treatment","Tween_Control"))
levels(gs) <- groups
design <- model.matrix(~group + 0, data.frame(group = gs))
colnames(design) <- levels(gs)

## remove missings
countsOvS0 <- countsOvS0[complete.cases(countsOvS0), ] 

## fit a linear model:
fit <- lmFit(countsOvS0, design)

# set up contrasts of interest and recalculate model coefficients
cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.05)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
tT <- subset(tT, select=c("adj.P.Val","P.Value","t","B","logFC","AveExpr"))
tT_Ordered <- tT[order(tT$P.Value),]
tT_Sig <- subset(tT_Ordered, adj.P.Val < 0.05)

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=0.5)

dtf<- as.data.frame((dT))
## ############Deseq2#####################------------------------------------########
countsOvS.matrix <-as.matrix(countsOvS)
str(countsOvS.matrix)

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countsOvS.matrix, colData = annotation, design = ~ condition)

# keep genes with at least 10 reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# differentially expressed genes (and other calculations)
dds <- DESeq(dds)
res <- results(dds,alpha=0.05)

#show the results
res

#summarize some basic tallies using the summary function
summary(res)

# order the results table by the smallest p-value
resOrdered <- res[order(res$pvalue),]
resOrdered <- as.data.frame(resOrdered)

# export only the results which pass an adjusted p-value threshold
resSig <- subset(resOrdered, padj < 0.05)

## merged significant genes of limma and deseq2##############################################################
merged_df<- merge(resSig,tT_Sig,by=0)
rownames(merged_df)<- merged_df$Row.names
merged_df$Row.names<- NULL
merged_dtf<- merge(merged_df,dtf, by =0)


# select 100 genes with largest log fold change scores of limma:
select <- order(abs(merged_dtf$logFC), decreasing=TRUE)[1:100]

library(pheatmap)
reso <- 1200
length <- 100*reso/72
png("02_heatmaps/OvT_heatmap.png",width     = length, height    = length)
# plot simple heatmap with clustering for merged data
pheatmap(countsOvS0[select,], cluster_rows=TRUE, cluster_cols=TRUE,
         show_rownames=TRUE, show_colnames=TRUE, fontsize_col=8,annotation_col=annotation)
dev.off()

# Build histogram of P-values for all genes. Normal test
png("03_histograms/OvT_histogram.png")
# assumption is that most genes are not differentially expressed.
hist(merged_df$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")
dev.off()

# summarize test results as "up", "down" or "not expressed"
rownames(merged_dtf)<- merged_dtf$Row.names
merged_dtf$Row.names<- NULL

dt_up<- merged_dtf %>%
  filter(merged_dtf[,13]=="1")

dt_up_row<- rownames(dt_up)

write.table(dt_up_row, "01_up_down_genes/oilT_tweenC_upreg.txt",row.names = F, quote = F)

dt_down<- merged_dtf %>%
  filter(merged_dtf[,13]=="-1")

dt_down_row<- rownames(dt_down)

write.table(dt_down_row, "01_up_down_genes/oilT_tweenC_downreg.txt",row.names = F, quote = F)

dt_no_sig<- merged_dtf %>%
  filter(merged_dtf[,13]=="0")

png("04_venndiagram/OvT_venndiagram.png")
# Venndiagram showing expressed and non expressed genes:
vennDiagram(merged_dtf[,13], circle.col=palette(),cex=c(4,0.5,0.5),names = c("Expressed", "Not expressed"),main="Comparison of Oil Treatment and Tween Control")
dev.off()

# add a column of NAs
merged_dtf$diffexpressed <- "NO"

# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
merged_dtf$diffexpressed[merged_dtf$logFC > 0.5 & merged_dtf$P.Value < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
merged_dtf$diffexpressed[merged_dtf$logFC < -0.5 & merged_dtf$P.Value < 0.05] <- "DOWN"

# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
merged_dtf$delabel <- NA
merged_dtf$delabel[merged_dtf$diffexpressed != "NO"] <- row.names(merged_dtf)[merged_dtf$diffexpressed != "NO"]


# Assuming 'merged_dtf' contains your data
top_20 <- merged_dtf[order(merged_dtf$logFC),][1:20,]

png("05_volcano_plots/OvT_volcanoplots.png",width     = 800, height    = 800)
#ggplot showing volcano plot
ggplot(data=merged_dtf, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label=ifelse(rank(-abs(logFC)) <= 20, delabel, ""))) + 
  geom_point() + 
  theme_minimal() +
  geom_text(size = 2.5) +
  scale_color_manual(values = c("DOWN" = "blue", "UP" = "red"))
dev.off()

## KEGG ENRICHMENT ###############################################################
gene_up <- rownames(merged_dtf %>%
                      filter(diffexpressed=="UP"))


kk <- enrichKEGG(gene         = unique(gene_up),
                 organism     = 'ang')

kk_df<- as.data.frame(kk)

# Replace the text in the specified column
kk@result[["Description"]] <- gsub(" - Aspergillus niger.*", "", kk@result[["Description"]])

png("06_kegg_enrich_plots/OvT_dotplots.png")
dotplot(kk, showCategory=1000) + ggtitle("dotplot ")
dev.off()

gene_down <- rownames(merged_dtf %>%
                        filter(diffexpressed=="DOWN"))


kk <- enrichKEGG(gene         = gene_down,
                 organism     = 'ang')

kk_df<- as.data.frame(kk)

# Replace the text in the specified column
kk@result[["Description"]] <- gsub(" - Aspergillus niger.*", "", kk@result[["Description"]])

png("06_kegg_enrich_plots/OvT_dotplots2.png")
dotplot(kk, showCategory=1000) + ggtitle("dotplot ")
dev.off()

## functional annotations ###########################################################
library(biomaRt)
ensembl_fungi <- useMart(host="https://fungi.ensembl.org", 
                         biomart="fungi_mart", 
                         port = 443)

head(listDatasets(ensembl_fungi))

searchDatasets(mart = ensembl_fungi, pattern = "niger")


ensembl <- useDataset(dataset = "aniger_eg_gene", mart = ensembl_fungi)

aay<-listAttributes(ensembl)

ensembl_ids <- getBM(attributes = c("ensembl_gene_id","start_position","end_position", "external_gene_name"),
                     filters = "ensembl_gene_id",
                     values = gene_up,
                     mart = ensembl)


###################################################
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.16")
#install.packages("devtools")
#devtools::install_version("dbplyr", version = "2.3.4")
#BiocManager::install("BiocFileCache", version = "3.15")

library(AnnotationHub) 
library(AnnotationDbi)
NigerDb <- loadDb("filenew1.sqlite")
keytypes(NigerDb)
str(NigerDb)
egid <- head(keys(NigerDb, "PMID"))

genes_all <- as.data.frame(read_excel("DAVID_LOCUS_TAG.Aspergillus_niger.xlsx"))
colnames(genes_all)<- c("entrez_ID","From","To")
rownames(genes_all)<- genes_all$To
genes_all$To<- NULL

gene_up_df<- data.frame(gene_up)
rownames(gene_up_df)<- gene_up_df$gene_up
#gene_up_df$gene_up<- NULL
merged_test<- merge(gene_up_df,genes_all, by=0)
rownames(merged_test) <- merged_test$To
merged_test$gene_up<- NULL


genes_to_test<- merged_test$entrez_ID

# Perform Gene Ontology (GO) enrichment analysis for Biological Process (BP) terms using "enrichGO" function.
GO_results <- enrichGO(gene = genes_to_test, OrgDb = NigerDb, keyType = "ENTREZID", ont = "BP")
as.data.frame(GO_results)

png("07_GO_plots/OvT_barplot_BP.png")
# Create a barplot visualizing the top 15 enriched GO terms for Biological Process.
fit <- plot(barplot(GO_results)) #showCategory = 15))
dev.off()


# Perform GO enrichment analysis for Molecular Function (MF) terms using "enrichGO" function.
GO_results2 <- enrichGO(gene = genes_to_test, OrgDb = NigerDb, keyType = "ENTREZID", ont = "MF")
as.data.frame(GO_results2)

png("07_GO_plots/OvT_barplot_MF.png")
# Create a barplot visualizing the top 15 enriched GO terms for Molecular Function.
fit2 <- plot(barplot(GO_results2))
dev.off()

# Perform GO enrichment analysis for Cellular Component (CC) terms using "enrichGO" function.
GO_results3 <- enrichGO(gene = genes_to_test, OrgDb = NigerDb, keyType = "ENTREZID", ont = "CC")
as.data.frame(GO_results3)

png("07_GO_plots/OvT_barplot_CC.png")
# Create a barplot visualizing the top 15 enriched GO terms for Cellular Component.
fit3 <- plot(barplot(GO_results3))
dev.off()

## Tween_treatment_vs_Salt_control:                        #### ############## ####

rm(list = ls())
ay<- gc()
options(scipen = 5)
# load packages ####
sessioninfo::session_info()
#set the working directory
setwd("/Aakash/temp_dir/Mat_Pap")
requiredPackages <- c("pacman")
for (package in requiredPackages) { #Installs packages if not yet installed
  if (!requireNamespace(package, quietly = TRUE))
    install.packages(package)
}

for (i in c(
  "data.table",
  "kknn",
  "tidyr",
  "ggplot2",
  "dplyr",
  "GEOquery",
  "limma",
  "umap",
  "tidyr",
  "here",
  "tibble",
  "car",
  "DESeq2",
  "readxl",
  "clusterProfiler",
  "biomaRt",
  "pheatmap",
  "AnnotationHub"
)
) {
  suppressPackageStartupMessages(
    pacman::p_load(i, character.only = TRUE
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

# Oil as Treatment vs Salt as Control ####
#load the annotation file:
annotation <- read.table("annotationST.txt", header = TRUE, row.names = 1, check.names = FALSE)

count1s <- read.table("all_counts.txt", header = TRUE, row.names = 1, check.names = FALSE,sep = ",")
#make the counts file ready:
countsOvS<-count1s[,c(5,6,3,4)]

genes_all <- as.data.frame(read_excel("DAVID_LOCUS_TAG.Aspergillus_niger.xlsx"))
colnames(genes_all)<- c("entrez_ID","From","To")
rownames(genes_all)<- genes_all$From
genes_all$From<- NULL
merged_test<- merge(countsOvS,genes_all,by=0)
rownames(merged_test) <- merged_test$To

merged_test<- merged_test %>%
  dplyr::select(2,3,4,5)

countsOvS<- merged_test
#check the dimensions
dim(countsOvS)
head(countsOvS)


## samples information
head(annotation)
mat_fun <- function(m){
  m2 <- apply(m,  2,  function(x) as.numeric(paste(x)))
  colnames(m2) <- colnames(m)
  rownames(m2) <- rownames(m)
  return(m2)
}
countsOvS<-mat_fun(countsOvS)

class(countsOvS)
## transform counts into matrix
countsOvS.matrix <-as.matrix(countsOvS)


countsOvS0 <- countsOvS[, c("tween_R1", "tween_R2","salt_R1","salt_R2")]

all(colnames(countsOvS0) == rownames(annotation))


## group membership for all samples
gsms <- "1100"
sml <- strsplit(gsms, split="")[[1]]



## log2 transformation
ex <- countsOvS0
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
countsOvS0 <- log2(ex) }

## assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("Tween_Treatment","Salt_Control"))
levels(gs) <- groups
design <- model.matrix(~group + 0, data.frame(group = gs))
colnames(design) <- levels(gs)

## remove missings
countsOvS0 <- countsOvS0[complete.cases(countsOvS0), ] 

## fit a linear model:
fit <- lmFit(countsOvS0, design)

# set up contrasts of interest and recalculate model coefficients
cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.05)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
tT <- subset(tT, select=c("adj.P.Val","P.Value","t","B","logFC","AveExpr"))
tT_Ordered <- tT[order(tT$P.Value),]
tT_Sig <- subset(tT_Ordered, adj.P.Val < 0.05)

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

dtf<- as.data.frame((dT))
## ############Deseq2#####################------------------------------------########
countsOvS.matrix <-as.matrix(countsOvS)
str(countsOvS.matrix)

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countsOvS.matrix, colData = annotation, design = ~ condition)

# keep genes with at least 10 reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# differentially expressed genes (and other calculations)
dds <- DESeq(dds)
res <- results(dds,alpha=0.05)

#show the results
res

#summarize some basic tallies using the summary function
summary(res)

# order the results table by the smallest p-value
resOrdered <- res[order(res$pvalue),]
resOrdered <- as.data.frame(resOrdered)

# export only the results which pass an adjusted p-value threshold
resSig <- subset(resOrdered, padj < 0.05)

## merged significant genes of limma and deseq2##############################################################
merged_df<- merge(resSig,tT_Sig,by=0)
rownames(merged_df)<- merged_df$Row.names
merged_df$Row.names<- NULL
merged_dtf<- merge(merged_df,dtf, by =0)


# select 100 genes with largest log fold change scores of limma:
select <- order(abs(merged_dtf$logFC), decreasing=TRUE)[1:100]

library(pheatmap)
reso <- 1200
length <- 100*reso/72
png("02_heatmaps/SvT_heatmap.png",width     = length, height    = length)
# plot simple heatmap with clustering for merged data
pheatmap(countsOvS0[select,], cluster_rows=TRUE, cluster_cols=TRUE,
         show_rownames=TRUE, show_colnames=TRUE, fontsize_col=8,annotation_col=annotation)
dev.off()

# Build histogram of P-values for all genes. Normal test
png("03_histograms/SvT_histogram.png")
# assumption is that most genes are not differentially expressed.
hist(merged_df$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")
dev.off()

# summarize test results as "up", "down" or "not expressed"
rownames(merged_dtf)<- merged_dtf$Row.names
merged_dtf$Row.names<- NULL

dt_up<- merged_dtf %>%
  filter(merged_dtf[,13]=="1")

dt_up_row<- rownames(dt_up)

write.table(dt_up_row, "01_up_down_genes/TweenT_SaltC_upreg.txt",row.names = F, quote = F)

dt_down<- merged_dtf %>%
  filter(merged_dtf[,13]=="-1")

dt_down_row<- rownames(dt_down)

write.table(dt_down_row, "01_up_down_genes/TweenT_SaltC_downreg.txt",row.names = F, quote = F)

dt_no_sig<- merged_dtf %>%
  filter(merged_dtf[,13]=="0")

png("04_venndiagram/SvT_venndiagram.png")
# Venndiagram showing expressed and non expressed genes:
vennDiagram(merged_dtf[,13], circle.col=palette(),cex=c(4,0.5,0.5),names = c("Expressed", "Not expressed"),main="Comparison of Oil Treatment and Tween Control")
dev.off()

# add a column of NAs
merged_dtf$diffexpressed <- "NO"

# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
merged_dtf$diffexpressed[merged_dtf$logFC > 0.2 & merged_dtf$P.Value < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
merged_dtf$diffexpressed[merged_dtf$logFC < -0.2 & merged_dtf$P.Value < 0.05] <- "DOWN"

# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
merged_dtf$delabel <- NA
merged_dtf$delabel[merged_dtf$diffexpressed != "NO"] <- row.names(merged_dtf)[merged_dtf$diffexpressed != "NO"]


# Assuming 'merged_dtf' contains your data
top_20 <- merged_dtf[order(merged_dtf$logFC),][1:20,]

png("05_volcano_plots/SvT_volcanoplots.png",width     = 800, height    = 800)
#ggplot showing volcano plot
ggplot(data=merged_dtf, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label=ifelse(rank(-abs(logFC)) <= 20, delabel, ""))) + 
  geom_point() + 
  theme_minimal() +
  geom_text(size = 2.5) +
  scale_color_manual(values = c("DOWN" = "blue", "UP" = "red"))
dev.off()

## KEGG ENRICHMENT ###############################################################
gene_up <- rownames(merged_dtf %>%
                      filter(diffexpressed=="UP"))


kk <- enrichKEGG(gene         = unique(gene_up),
                 organism     = 'ang')

kk_df<- as.data.frame(kk)

# Replace the text in the specified column
kk@result[["Description"]] <- gsub(" - Aspergillus niger.*", "", kk@result[["Description"]])

#png("06_kegg_enrich_plots/SvT_dotplots.png")
dotplot(kk) + ggtitle("dotplot ")
dev.off()

gene_down <- rownames(merged_dtf %>%
                        filter(diffexpressed=="DOWN"))


kk <- enrichKEGG(gene         = gene_down,
                 organism     = 'ang')

kk_df<- as.data.frame(kk)

# Replace the text in the specified column
kk@result[["Description"]] <- gsub(" - Aspergillus niger.*", "", kk@result[["Description"]])

png("06_kegg_enrich_plots/SvT_dotplots2.png")
dotplot(kk, showCategory=1000) + ggtitle("dotplot ")
dev.off()

## functional annotations ###########################################################
library(biomaRt)
ensembl_fungi <- useMart(host="https://fungi.ensembl.org", 
                         biomart="fungi_mart", 
                         port = 443)

head(listDatasets(ensembl_fungi))

searchDatasets(mart = ensembl_fungi, pattern = "niger")


ensembl <- useDataset(dataset = "aniger_eg_gene", mart = ensembl_fungi)

aay<-listAttributes(ensembl)

ensembl_ids <- getBM(attributes = c("ensembl_gene_id","start_position","end_position", "external_gene_name"),
                     filters = "ensembl_gene_id",
                     values = gene_up,
                     mart = ensembl)


###################################################
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.16")
#install.packages("devtools")
#devtools::install_version("dbplyr", version = "2.3.4")
#BiocManager::install("BiocFileCache", version = "3.15")

library(AnnotationHub) 
library(AnnotationDbi)
NigerDb <- loadDb("filenew1.sqlite")
keytypes(NigerDb)
str(NigerDb)
egid <- head(keys(NigerDb, "PMID"))

genes_all <- as.data.frame(read_excel("DAVID_LOCUS_TAG.Aspergillus_niger.xlsx"))
colnames(genes_all)<- c("entrez_ID","From","To")
rownames(genes_all)<- genes_all$To
genes_all$To<- NULL

gene_up_df<- data.frame(gene_up)
rownames(gene_up_df)<- gene_up_df$gene_up
#gene_up_df$gene_up<- NULL
merged_test<- merge(gene_up_df,genes_all, by=0)
rownames(merged_test) <- merged_test$To
merged_test$gene_up<- NULL


genes_to_test<- merged_test$entrez_ID

# Perform Gene Ontology (GO) enrichment analysis for Biological Process (BP) terms using "enrichGO" function.
GO_results <- enrichGO(gene = genes_to_test, OrgDb = NigerDb, keyType = "ENTREZID", ont = "BP")
as.data.frame(GO_results)

png("07_GO_plots/SvT_barplot_BP.png")
# Create a barplot visualizing the top 15 enriched GO terms for Biological Process.
fit <- plot(barplot(GO_results)) #showCategory = 15))
dev.off()


# Perform GO enrichment analysis for Molecular Function (MF) terms using "enrichGO" function.
GO_results2 <- enrichGO(gene = genes_to_test, OrgDb = NigerDb, keyType = "ENTREZID", ont = "MF")
as.data.frame(GO_results2)

png("07_GO_plots/SvT_barplot_MF.png")
# Create a barplot visualizing the top 15 enriched GO terms for Molecular Function.
fit2 <- plot(barplot(GO_results2))
dev.off()

# Perform GO enrichment analysis for Cellular Component (CC) terms using "enrichGO" function.
GO_results3 <- enrichGO(gene = genes_to_test, OrgDb = NigerDb, keyType = "ENTREZID", ont = "CC")
as.data.frame(GO_results3)

png("07_GO_plots/SvT_barplot_CC.png")
# Create a barplot visualizing the top 15 enriched GO terms for Cellular Component.
fit3 <- plot(barplot(GO_results3))
dev.off()

