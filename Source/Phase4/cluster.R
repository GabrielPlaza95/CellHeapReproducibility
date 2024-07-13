library(cellassign)
library(DropletUtils)
library(Seurat)
library(dplyr)
library(patchwork)
library(scales)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(MAST)
library(SingleCellExperiment)
library(Matrix)
library(reshape2)
library(org.Hs.eg.db)
library(edgeR)
library(scater)
library(loomR)
library(reticulate)

# On error, save debugging info to file last.dump.rda
dump_and_quit <- function() {
  print(reticulate::py_last_error())
  dump.frames(to.file = TRUE)
  quit(status = 1)
}
options(error = dump_and_quit)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Sample acession number must be provided", call.=FALSE)
}

PROJECT_DIR <- Sys.getenv("PROJECT_DIR")

IN_DIR <- file.path(PROJECT_DIR, "out", "phase3")
OUT_DIR <- file.path(PROJECT_DIR, "out", "phase4")

if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR)
}

########## Seurat objects ##########

alldata.list = list()

for (sample in args) {
  counts_file <- file.path(IN_DIR, paste("features_not_infected_", sample, ".tsv", sep=""))
  counts_table <- read.table(counts_file, header=TRUE, sep = "\t")
  so <- CreateSeuratObject(counts = counts_table, project = sample)
  
  meta_file <- file.path(IN_DIR, paste("metadata_not_infected_", sample, ".csv", sep=""))
  meta_table <- read.csv(meta_file, header=TRUE, row.names=1)
  so[[]] <- meta_table
  
  alldata.list <- append(alldata.list, list(so))
}

########## MERGE #########

normalize_and_find_features <- function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x
}

alldata.list <- lapply(X = alldata.list, FUN = normalize_and_find_features)

features <- SelectIntegrationFeatures(object.list = alldata.list)

immune.anchors <- FindIntegrationAnchors(object.list = alldata.list, anchor.features = features)

##########

#Create the integrated assays
immune.combined <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(immune.combined) <- "integrated"

# ScaleData() pattern - to dimension high variable genes (previously identified - 2000). It does not affect the PCA and clustering results.
immune.combined <- ScaleData(immune.combined)

# PCA clustering and visualization - pattern of the npcs = 50 from Seurat vignette
immune.combined <- RunPCA(immune.combined, npcs = 50)

# Dimension test - above 16 it is stable, to use dim = 20 
# Number of replicas - based on Seurat vignette
immune.combined <- JackStraw(immune.combined, num.replicate = 100)
immune.combined <- ScoreJackStraw(immune.combined, dims = 1:20)

# Graphs related to dimension test
dpi=300
png(file=file.path(OUT_DIR, "dimension_test_jackstraw.png"), width = dpi*20, height = dpi*14, units = "px",res = dpi,type='cairo')
JackStrawPlot(immune.combined, dims = 1:20)
dev.off()

dpi=300
png(file=file.path(OUT_DIR, "dimension_test_elbow.png"), width = dpi*20, height = dpi*14, units = "px",res = dpi,type='cairo')
ElbowPlot(immune.combined)
dev.off()

# UMAP - visualization and exploration of datasets. UMAP input - Seurat suggestion, same PCs as input to the clustering analysis
# UMAP clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:16)

#Calculate k.param nearest neighbors for a dataset.
# FindNeighbors() function - the input is the previously defined dataset dimension to construct the KNN graph
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:16)

# Identify clusters of cells by SNN modularity.
# Resolution parameter is based on Liao et al 2020's paper
# To group cells - Seurat uses modularity optimization techniques (Louvain - pattern, opctions: Louvain algorithm with multilevel refinement; SLM algorithm e Leiden)
# Resolution parameter defines the clustering granularity - greater values greater the number of clusters. Seurat vignette suggest values 0.4-1.2 for 3k cells.
# In general, ideal resolution increases for bigger datasets. Silvin et al 2020's paper uses 0.3 (to test 0.5).
immune.combined <- FindClusters(immune.combined, resolution = 1.2)

DefaultAssay(immune.combined) <- "RNA"

# To convert Seurat object into SingleCellExperiment object
immune.combined <- JoinLayers(immune.combined)

#Find markers for every cluster compared to all remaining cells
immune.combined.markers <- FindAllMarkers(immune.combined, assay = 'RNA',logfc.threshold = 0.25, only.pos = TRUE, test.use = 'MAST')

test.sce <- as.SingleCellExperiment(immune.combined)

# Matrix construction
marker_gene_list <- list(
  InflammatoryMacrophages = c("STAT1", "TNF", "IL6", "CD68", "CD14", "IL1B", "IRF5"),
  NonInflammatoryMacrophages = c("CD14", "CD68","IRF4", "IL10", "TGFB1","ARG1", "TGFBR2", "CD163"),
  Neutrophils = c("FCGR3B", "PI3", "G0S2", "ELANE", "LCN2", "ORM1", "MMP8"),
  MastCells = c("CPA3","FCER1A", "TPSAB1", "RGS13", "KIT", "TPSG1", "SLC18A2", "TPSB2"),
  Basophils = c("ANPEP", "CD22", "FCGR2B", "FCER1A","CD33", "IL3RA", "ENPP3"),
  #innate lymphoid
  NKCells = c("NCR1", "NCAM1", "KIR3DL1", "ITGAE", "KLRC1", "NKG7"),
  ILC1sCitotoxic= c("IFNG","TBX21", "EOMES", "NCR2", "ITGAE"),
  ILC1sNonCitotoxic= c("IFNG","IL7R", "TBX21"),
  ILC2s = c("GATA3", "PTGDR2", "KLRB1", "IL33", "IL1RL1", "KLRG1"),
  ILC3ssubpopNCRneg = c("RORGT", "IL23R", "IL17A"),
  ILC3ssubpopNCRpos = c("NCR2", "NCR1", "RORGT", "TBX21"),
  LTi  = c("CD4", "IL7R", "CCR6", "RORGT"),
  # adaptiive lymphoid
  TCD4Th1 = c("CD3E", "CD4","TBX21", "IFNG", "TNF"),
  TCD4Th2 = c("CD3E", "CD4", "GATA3", "IL4", "IL5", "IL13"),
  TCD4Th17 = c("CD3E", "CD4","IL17A", "IL17F", "IL21"),
  TCD4TReg = c("CD3E", "CD4", "FOXP3", "TGFB1", "IL10"),
  TCD8Citotoxic = c("GNLY", "PRF1", "CD3E", "CD8A"),
  TGammaDeltaCells = c("TRGC1", "TRDC", "CD3E"),
  PlasmaCells = c("CD79A", "TNFRSF13C", "KRT20", "IGHM", "IGHD", "IGKC","IGLC2", "JCHAIN", "XBP1", "MZB1"),
  BregCells = c("CD19", "CD24", "CD27", "GZMB", "IL2RA", "TFRC","CD274"),
  # epithelial
  Secretory = c("SCGB1A1", "SCGB3A1", "MSMB"),
  Basal = c("KRT5", "AQP3", "TP63"),
  Ciliated = c("CAPS", "TPPP3", "RSPH1"),
  Squamous = c("KRT13", "KRT4", "SPRR3"),
  Inflammatory = c("KRT8", "KRT18", "MMP7"),
  AT2 = c("SFTPC", "SFTPA1", "SFTPB")
)

marcadores <- marker_list_to_mat(marker_gene_list, include_other = FALSE)

# Considera apenas marcadores que existem na amostra
shared <- intersect(rownames(marcadores), rownames(test.sce))
test.sce <- test.sce[shared,]

# Garante que não nenhuma coluna ou linha esteja vazia
test.sce <- test.sce[which(rowSums(counts(test.sce)) > 0),]
test.sce <- test.sce[,which(colSums(counts(test.sce)) > 0)]

# Ajusta lista de marcadores após filtro anterior
shared <- intersect(rownames(marcadores), rownames(test.sce))
marcadores <- marcadores[shared,]

# Size Factors
test.sce <- scran::computeSumFactors(test.sce)
s1 <- sizeFactors(test.sce, onAbsence = "warn")

# Cellassign fit
fit <- cellassign(exprs_obj = test.sce, marker_gene_info = marcadores, s = s1, learning_rate = 1e-2, shrinkage = TRUE,  verbose = TRUE)

test.sce$cellassign_type <- celltypes(fit, assign_prob = 0.95)

# Plot heatmap
dpi = 300
png(file=file.path(OUT_DIR, "heatmap.png"), width = dpi*30, height = dpi*24, units = "px", res = dpi, type='cairo')
pheatmap::pheatmap(cellprobs(fit))

#Plot the expression of a given gene with cellassign_type coloring.
dpi = 300
png(file=file.path(OUT_DIR, "violin_tnf.png"), width = dpi*30, height = dpi*24, units = "px",res = dpi,type='cairo')
plotExpression(test.sce, features = "TNF", x = "cellassign_type", colour_by = 'cellassign_type')
