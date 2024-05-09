library(Seurat)
library(dplyr)
library(patchwork)
library(scales)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(SoupX)
library(conflicted)

# On error, save debugging info to file last.dump.rda
dump_and_quit <- function() {
  dump.frames(to.file = TRUE)
  quit(status = 1)
}
options(error = dump_and_quit)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Sample acession number must be provided", call.=FALSE)
}

PROJECT_DIR <- Sys.getenv("PROJECT_DIR")

SAMPLE <- args[1]
IN_DIR <- file.path(PROJECT_DIR, "out", "phase2", SAMPLE, "outs")
OUT_DIR <- file.path(PROJECT_DIR, "out", "phase3")

if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR)
}

# SoupX - Automatic mode using cellranger outputs 
# Obtain better results - Previous basic clustering information
# Using 10X data mapped with cellranger 
# Default clustering produced by cellranger is automatically loaded and used
# Load cellranger outputs, to create SoupChannel object, and to estime soup profile

sample_SpX <- load10X(IN_DIR)

# Estime the "rho" parameter that represents the fraction of contamination
# rho = 0 means no contamination,  rho = 1 means 100% of UMIs into a droplet are "soup"

sample_SpX <- autoEstCont(sample_SpX, doPlot=FALSE)

#Generate the corrected matrix without contamination.

out_sample_SpX <- adjustCounts(sample_SpX)

# Using corrected matrix from SoupX directly in the Seurat through the CreateSeuratObject function

sample <- CreateSeuratObject(counts = out_sample_SpX, project = SAMPLE, min.cells=1)

#Calculate the percentage of mitochondrial RNA reads
sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^MT-")

# Metrics
feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

png(file.path(OUT_DIR, paste("before_filters_", SAMPLE, ".png", sep="")))

VlnPlot(sample, layer = "counts", features = feats, pt.size = 0.1, ncol = 3) + NoLegend()
dev.off()

# Step of quality control filters defined in previous analyzes based on Wauters et al, 2021
# Cell filtering parameters for selecting cells based on number of genes/cell, UMI counts/cell, and percent mitochondrial genes according to Wauters et al, 2021.
sample_filtered <- subset(sample, subset=nFeature_RNA>150 & nFeature_RNA<3000 & percent.mt<20 & nCount_RNA>301)

png(file.path(OUT_DIR, paste("after_filters_", SAMPLE, ".png", sep="")))

VlnPlot(sample_filtered, layer = "counts", features = feats, pt.size = 0.1, ncol = 3) + NoLegend()
dev.off()

# The next steps involve looking for viral RNA in the samples - it just kepts the commands from the previous phase 3 script
# Export .tsv and .csv files

features_file <- file.path(OUT_DIR, paste("features_", SAMPLE, ".tsv", sep=""))
metadata_file <- file.path(OUT_DIR, paste("metadata_", SAMPLE, ".csv", sep=""))

out_test <- as.matrix(sample_filtered@assays$RNA$counts)

write.table(out_test, file=features_file, quote=FALSE, sep='\t', col.names = TRUE)
write.csv(sample_filtered@meta.data, file=metadata_file)

# Import of seurat filter file
sample <- read.table(features_file, sep='\t', header=TRUE)

# Select rows containing sarscov2

#covid_rows <- grep("virus-v6", row.names(sample))
covid_rows <- grep("SARS", row.names(sample)) #TODO: verificar nome da feature

sample_sarscov2 <- sample[covid_rows,,drop=FALSE]

# Transpose data frame
sample_sarscov2 <- t(sample_sarscov2)

#print(sample_sarscov2)

#print rows where all columns are zero
not_infected <- row.names(sample_sarscov2)[which(rowSums(sample_sarscov2)==0)]

#print rows where some columns are different from zero
infected <- rownames(sample_sarscov2)[which(rowSums(sample_sarscov2)>0)]

features_not_infected_file <- file.path(OUT_DIR, paste("features_not_infected_", SAMPLE, ".tsv", sep=""))
features_infected_file <- file.path(OUT_DIR, paste("features_infected_", SAMPLE, ".tsv", sep=""))

# Dataframe not infected 
write.table(sample[,not_infected,drop=FALSE], file=features_not_infected_file, quote=FALSE, sep='\t', col.names=TRUE)

# Dataframe infected 
write.table(sample[,infected,drop=FALSE], file=features_infected_file, quote=FALSE, sep='\t', col.names=TRUE)
