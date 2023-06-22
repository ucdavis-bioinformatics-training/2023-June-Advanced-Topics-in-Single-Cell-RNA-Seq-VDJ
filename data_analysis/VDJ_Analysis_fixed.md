---
title: "Single Cell V(D)J Analysis"
author: "UCD Bioinformatics Core"
output:
    html_document:
      keep_md: TRUE
---

Last Updated March 23, 2022

# Single Cell V(D)J Analysis

## Packages
```{r, libraries, warning=FALSE,error=FALSE,message=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}

if (!any(rownames(installed.packages()) == "ggplot2")){
  BiocManager::install("ggplot2")
}

if (!any(rownames(installed.packages()) == "knitr")){
  BiocManager::install("knitr")
}

if (!any(rownames(installed.packages()) == "kableExtra")){
  BiocManager::install("kableExtra")
}

if (!any(rownames(installed.packages()) == "dplyr")){
  BiocManager::install("dplyr")
}

if (!any(rownames(installed.packages()) == "tidyr")){
  BiocManager::install("tidyr")
}
if (!any(rownames(installed.packages()) == "magrittr")){
  BiocManager::install("magrittr")
}

if (!any(rownames(installed.packages()) == "scRepertoire")){
  BiocManager::install("scRepertoire")
}

library(ggplot2)
library(tidyr)
library(magrittr)
library(knitr)
library(kableExtra)
library(dplyr)
library(scRepertoire)
<div class='r_output'>
 Download Cell Ranger results
```{r, eval=FALSE}
options(timeout=1200)
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2022-March-Advanced-Topics-in-Single-Cell-RNA-Seq-VDJ/main/data_analysis/cellranger_vdj_results.zip", "cellranger_vdj_results.zip")
system("unzip cellranger_vdj_results.zip")
</div>
## Set-up
```{r}
experiment_name = "VDJ Example"
dataset_loc <- "./cellranger_vdj_results"
ids <- c("Pool1")
<div class='r_output'>
 Sequencing metrics
```{r}
metrics <- paste(dataset_loc, ids, "metrics_summary.csv", sep = "/")
metrics_table <- do.call("cbind", lapply(metrics, function(x) {
  as.data.frame(t(read.csv(x)))
  }))
colnames(metrics_table) <- ids
rownames(metrics_table) <- gsub(".", " ", rownames(metrics_table), fixed = TRUE)
metrics_table  %>%
  kable(caption = 'Cell Ranger Results') %>%
  pack_rows("Overview", 1, 3, label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Sequencing Characteristics", 4, 8, label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Mapping Characteristics", 9, 26, label_row_css = "background-color: #666; color: #fff;") %>%
  kable_styling("striped", fixed_thead = TRUE)
</div>
The majority of the following functions and figures come from [scRepertoire](https://ncborcherding.github.io/vignettes/vignette.html). We will be exploring and making changes to the code as we go, so please take notes and don't be afraid to experiment and ask questions!

## Read in Cell Ranger VDJ Data
```{r}
clonotypes <- paste(dataset_loc, ids, "filtered_contig_annotations.csv", sep = "/")
vdj <- combineTCR(lapply(clonotypes, read.csv),
                  samples = ids,
                  ID = ids,
                  cells = "T-AB",
                  removeMulti = TRUE)
class(vdj)
str(vdj)
class(vdj[[1]])
head(vdj[[1]])
<div class='r_output'>
 Number of unique clonotypes
```{r}
quantContig(vdj, cloneCall="aa", group = "sample", scale = FALSE)
</div>## Distribution of clonotypes by abundance
```{r}
abundanceContig(vdj, cloneCall = "gene", group = "sample", scale = FALSE)
<div class='r_output'>
 Contig length distribution
```{r}
lengthContig(vdj, cloneCall="nt", scale=TRUE, chains = "combined", group="sample")
lengthContig(vdj, cloneCall="aa", chains = "single", group = "sample")
</div>
### Shared clonotypes

This function does not make much sense with only one sample, and is included just for the purposes of demonstration.

```{r}
compareClonotypes(vdj, numbers = 10, cloneCall = "aa", graph = "alluvial") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(caption = "Results of compareClonotypes() with numbers = 10.")
<div class='r_output'>
# Relative abundance of clones by frequency
```{r}
clonalHomeostasis(vdj, cloneCall = "aa") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
</div>
### Relative abundance of clones by index

Clonal index 1 represents the most frequent clone in a given sample, while index 1000 represents the 1000th most frequent clone.
```{r}
clonalProportion(vdj, cloneCall = "aa", split = c(10, 50, 100, 500, 1000)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
<div class='r_output'>
# Overlap analysis

Clonal overlap is scaled to the number of unique clonotypes in the smaller sample. This code errors on fewer than two samples.
```{r, eval=FALSE}
clonalOverlap(vdj, cloneCall = "gene+nt")
</div>
## Clonal diversity
```{r}
clonalDiversity(vdj, cloneCall = "aa", group = "samples")
<div class='r_output'>
 TCR clustering

This is slow. I suggest skipping it for now so that you don't get stuck at this point.
```{r, eval=FALSE}
tcr.clusters <- clusterTCR(vdj[[1]], chain = "TCRA", sequence = "aa", threshold = 0.9)
</div>
## Combine V(D)J and expression data

In the terminal, from your project directory, run `scp username@tadpole.genomecenter.ucdavis.edu:/share/workshop/vdj_workshop/R_objects/seurat_object.rds .`.

```{r}
library(Seurat)
expression <- readRDS("seurat_object.rds")
expression$barcode <- sapply(strsplit(colnames(expression), split = "-"), "[[", 1)
expression <- RenameCells(expression, new.names = expression$barcode)
vdj <- lapply(vdj, function (x) {
  b = sapply(strsplit(sapply(strsplit(x$barcode, split = "_"), "[[", 3), split = "-"), "[[", 1)
  x$barcode = b
  x
  })
expression <- combineExpression(vdj, expression, cloneCall="gene")
head(expression@meta.data)
<div class='r_output'># Find markers for "large" clones
```{r}
DimPlot(expression, group.by = "cloneType")
Idents(expression) <- "cloneType"
large.markers <-FindMarkers(expression, ident.1 = "Large (0.01 < X <= 0.1)")
head(large.markers)
Idents(expression) <- "orig.ident"
</div>
### Split V(D)J data into separate samples and re-run scRepertoire functions
```{r}
vdj <- expression2List(expression, group = "orig.ident")
quantContig(vdj, cloneCall="aa", scale = FALSE)
abundanceContig(vdj, cloneCall = "gene", scale = FALSE)
compareClonotypes(vdj, numbers = 20, cloneCall = "aa", graph = "alluvial") +
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(caption = "Results of compareClonotypes() with numbers = 20.")
clonalHomeostasis(vdj, cloneCall = "aa") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
clonalProportion(vdj, cloneCall = "aa", split = c(10, 50, 100, 500, 1000)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
clonalOverlap(vdj, cloneCall = "gene+nt")
clonalDiversity(vdj, cloneCall = "aa")
<div class='r_output'>
 Circos plots

This code is under construction! If we have time, I will attempt to update this and run it.

```{r, eval=FALSE}
?getCirclize
circles <- getCirclize(vdj, groupBy = "orig.ident" )

#Just assigning the normal colors to each cluster
grid.cols <- hue_pal()(length(unique(seurat$orig.ident)))
names(grid.cols) <- levels(seurat$orig.ident)

#Graphing the chord diagram
chordDiagram(circles, self.link = 1, grid.col = grid.cols)

data_to_circlize <- experiment.aggregate[[]][experiment.aggregate$RNA_snn_res.0.75 %in% b_cells & !is.na(experiment.aggregate$CTgene),]
dim(data_to_circlize)
head(data_to_circlize)

aa_seqs <- strsplit(as.character(unlist(data_to_circlize$CTaa)),split="_")
table(sapply(aa_seqs, length))
data_to_circlize$A_chain = sapply(aa_seqs, "[[", 1L)
data_to_circlize$B_chain = sapply(aa_seqs, "[[", 2L)

data_to_circlize$IGH = sapply(strsplit(data_to_circlize$CTstrict, split="_"), function(x) paste(unique(x[c(1)]),collapse="_"))
data_to_circlize$IGL = sapply(strsplit(data_to_circlize$CTstrict, split="_"), function(x) paste(unique(x[c(3)]),collapse="_"))
                              
# get optimal sequence order from trivial plot
chordDiagram(data.frame(data_to_circlize$IGH[1:15], data_to_circlize$IGL[1:15], times = 1), annotationTrack = "grid" )
seq.order <- get.all.sector.index()
circos.clear()
</div>
## Session Information
```{r}
sessionInfo()
<div class='r_output'>
