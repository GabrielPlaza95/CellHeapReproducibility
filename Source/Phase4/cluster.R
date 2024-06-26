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
  #TODO mudar para arquivos not_infected
  counts_file <- file.path(IN_DIR, paste("features_", sample, ".tsv", sep=""))
  counts_table <- read.table(counts_file, header=TRUE, sep = "\t")
  so <- CreateSeuratObject(counts = counts_table, project = sample)
  
  meta_file <- file.path(IN_DIR, paste("metadata_", sample, ".csv", sep=""))
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
immune.combined = JoinLayers(immune.combined)

test.sce <- as.SingleCellExperiment(immune.combined)

#TODO reintroduzir esta etapa
#Find markers for every cluster compared to all remaining cells
#immune.combined <- FindAllMarkers(immune.combined, assay = 'RNA',logfc.threshold = 0.25, only.pos = TRUE, test.use = 'MAST')

# Matrix construction
# Marker list - based on cellassign vignette
# https://irrationone.github.io/cellassign/articles/introduction-to-cellassign.html#constructing-a-marker-gene-matrix
marker_gene_list <- list(
  InflammatoryMacrophages = c("human----STAT1", "human----TNF", "human----IL6", "human----CD68", "human----CD14", "human----IL1B", "human----IRF5"),
  NonInflammatoryMacrophages = c("human----CD14", "human----CD68","human----IRF4", "human----IL10", "human----TGFB1","human----ARG1", "human----TGFBR2", "human----CD163"),
  Neutrophils = c("human----FCGR3B", "human----PI3", "human----G0S2", "human----ELANE", "human----LCN2", "human----ORM1", "human----MMP8"),
  MastCells = c("human----CPA3","human----FCER1A", "human----TPSAB1", "human----RGS13", "human----KIT", "human----TPSG1", "human----SLC18A2", "human----TPSB2"),
  Basophils = c("human----ANPEP", "human----CD22", "human----FCGR2B", "human----FCER1A","human----CD33", "human----IL3RA", "human----ENPP3"),
  #innate lymphoid
  NKCells = c("human----NCR1", "human----NCAM1", "human----KIR3DL1", "human----ITGAE", "human----KLRC1", "human----NKG7"),
  ILC1sCitotoxic= c("human----IFNG","human----TBX21", "human----EOMES", "human----NCR2", "human----ITGAE"),
  ILC1sNonCitotoxic= c("human----IFNG","human----IL7R", "human----TBX21"),
  ILC2s = c("human----GATA3", "human----PTGDR2", "human----KLRB1", "human----IL33", "human----IL1RL1", "human----KLRG1"),
  ILC3ssubpopNCRneg = c("human----RORGT", "human----IL23R", "human----IL17A"),
  ILC3ssubpopNCRpos = c("human----NCR2", "human----NCR1", "human----RORGT", "human----TBX21"),
  LTi  = c("human----CD4", "human----IL7R", "human----CCR6", "human----RORGT"),
  # adaptiive lymphoid
  TCD4Th1 = c("human----CD3E", "human----CD4","human----TBX21", "human----IFNG", "human----TNF"),
  TCD4Th2 = c("human----CD3E", "human----CD4", "human----GATA3", "human----IL4", "human----IL5", "human----IL13"),
  TCD4Th17 = c("human----CD3E", "human----CD4","human----IL17A", "human----IL17F", "human----IL21"),
  TCD4TReg = c("human----CD3E", "human----CD4", "human----FOXP3", "human----TGFB1", "human----IL10"),
  TCD8Citotoxic = c("human----GNLY", "human----PRF1", "human----CD3E", "human----CD8A"),
  TGammaDeltaCells = c("human----TRGC1", "human----TRDC", "human----CD3E"),
  PlasmaCells = c("human----CD79A", "human----TNFRSF13C", "human----KRT20", "human----IGHM", "human----IGHD", "human----IGKC","human----IGLC2", "human----JCHAIN", "human----XBP1", "human----MZB1"),
  BregCells = c("human----CD19", "human----CD24", "human----CD27", "human----GZMB", "human----IL2RA", "human----TFRC","human----CD274"),
  # epithelial
  Secretory = c("human----SCGB1A1", "human----SCGB3A1", "human----MSMB"),
  Basal = c("human----KRT5", "human----AQP3", "human----TP63"),
  Ciliated = c("human----CAPS", "human----TPPP3", "human----RSPH1"),
  Squamous = c("human----KRT13", "human----KRT4", "human----SPRR3"),
  Inflammatory = c("human----KRT8", "human----KRT18", "human----MMP7"),
  AT2 = c("human----SFTPC", "human----SFTPA1", "human----SFTPB")
)

#TODO voltar a usar a lista real
#marcadores <- marker_list_to_mat(marker_gene_list, include_other = FALSE)

marker_gene_list_teste <- list(
  TCD4Th1 = sample(rownames(test.sce), 3),
  TCD4Th2 = sample(rownames(test.sce), 6)
)

marcadores <- marker_list_to_mat(marker_gene_list_teste, include_other = FALSE)

marcadores_teste <- match(rownames(marcadores), rownames(test.sce))
stopifnot(all(!is.na(marcadores_teste)))

test.sce <- test.sce[marcadores_teste,]
stopifnot(all.equal(rownames(marcadores), rownames(test.sce)))

# Garante que não nenhuma coluna ou linha esteja vazia
test.sce <- test.sce[which(rowSums(counts(test.sce)) > 0),]
test.sce <- test.sce[,which(colSums(counts(test.sce)) > 0)]

# Ajusta lista de marcadores após filtro anterior
shared <- intersect(rownames(marcadores), rownames(test.sce))
marcadores <- marcadores[shared,]

# Size Factors
test.sce <- scran::computeSumFactors(test.sce)
s1 <- sizeFactors(test.sce, onAbsence = "warn")

# Cellassign fit - based on cellassign vignette
# https://irrationone.github.io/cellassign/articles/introduction-to-cellassign.html#constructing-a-marker-gene-matrix
fit <- cellassign(exprs_obj = test.sce, marker_gene_info = marcadores, s = s1, learning_rate = 1e-2, shrinkage = TRUE,  verbose = TRUE)

# Cell types
celltypes(fit,assign_prob = 0.95)

# https://bioconductor.org/packages/devel/bioc/vignettes/scater/inst/doc/overview.html
# https://bioconductor.riken.jp/packages/3.4/bioc/vignettes/scater/inst/doc/vignette.html#plots-of-expression-values
# https://www.singlecellcourse.org/single-cell-rna-seq-analysis-using-seurat.html#cell-type-annotation-using-singler
# Add a column - cell types identified by the cellassign in the test.sce object
test.sce$cellassign_type <- celltypes(fit,assign_prob = 0.95)

# Number of identified cell types
table(test.sce$cellassign_type)

# Output file names - total_violin_XX.png - it is possible to create other output files formats
# dpi=300
# png(file=file.path(OUT_DIR, "total_violin_fc.png"), width = dpi*20, height = dpi*14, units = "px",res = dpi,type='cairo')
# plotExpression(test.sce, features = c("human----FCGRT","human----FCGR1A", "human----FGL2", "human----FCGR2B", "human----FCAR","human----FCER1A","human----TRIM21", "human----FCGR3A", "human----FCGR3B","human----FCGR2A"), x = "cellassign_type",colour_by = 'cellassign_type')
# dev.off()
# 
# png(file=file.path(OUT_DIR, "total_violin_fc2.png"), width = dpi*20, height = dpi*14, units = "px",res = dpi,type='cairo')
# plotExpression(test.sce, features = c("human----FCGRT","human----FCGR1A", "human----FGL2", "human----FCGR2B"), x = "cellassign_type",colour_by = 'cellassign_type')
# dev.off()
# 
# png(file=file.path(OUT_DIR, "total_violin_fc3.png"), width = dpi*20, height = dpi*14, units = "px",res = dpi,type='cairo')
# plotExpression(test.sce, features = c("human----FCAR","human----FCER1A","human----TRIM21", "human----FCGR3A", "human----FCGR3B"), x = "cellassign_type",colour_by = 'cellassign_type')
# dev.off()
# 
# png(file=file.path(OUT_DIR, "total_violin_vs.png"), width = dpi*20, height = dpi*14, units = "px",res = dpi,type='cairo')
# plotExpression(test.sce, features = c("human----B2M","human----FCER1G", "human----HCK", "human----SYK","human----RAF1", "human----MAPK8","human----MAPK11","human----MAPK1", "human----LYN", "human----FYN", "human----LAT","human----LCP2", "human----PTPN6", "human----PTPN11","human----INPP5D", "human----BTK","human----UBE2W","human----UBE2N", "human----UBE2V2", "human----PSMD14", "human----IKBKG","human----NFKB1", "human----NFKB2", "human----VCP","human----CGAS", "human----DDX58"), x = "cellassign_type",colour_by = 'cellassign_type')
# dev.off()
# 
# png(file=file.path(OUT_DIR, "total_violin_vs2.png"), width = dpi*20, height = dpi*14, units = "px",res = dpi,type='cairo')
# plotExpression(test.sce, features = c("human----B2M","human----FCER1G", "human----HCK", "human----SYK","human----RAF1", "human----MAPK8","human----MAPK11","human----MAPK1"), x = "cellassign_type",colour_by = 'cellassign_type')
# dev.off()
# 
# png(file=file.path(OUT_DIR, "total_violin_vs3.png"), width = dpi*20, height = dpi*14, units = "px",res = dpi,type='cairo')
# plotExpression(test.sce, features = c("human----LYN", "human----FYN", "human----LAT","human----LCP2", "human----PTPN6", "human----PTPN11","human----INPP5D", "human----BTK"), x = "cellassign_type",colour_by = 'cellassign_type')
# dev.off()
# 
# png(file=file.path(OUT_DIR, "total_violin_vs4.png"), width = dpi*20, height = dpi*14, units = "px",res = dpi,type='cairo')
# plotExpression(test.sce, features = c("human----UBE2W","human----UBE2N", "human----UBE2V2", "human----PSMD14", "human----IKBKG","human----NFKB1", "human----NFKB2", "human----VCP","human----CGAS", "human----DDX58"), x = "cellassign_type",colour_by = 'cellassign_type')
# dev.off()
# 
# png(file=file.path(OUT_DIR, "total_violin_cito.png"), width = dpi*20, height = dpi*14, units = "px",res = dpi,type='cairo')
# plotExpression(test.sce, features = c("human----IFNG", "human----IFNA1", "human----TNF", "human----IL18","human----IL1B", "human----IL10","human----CCL5","human----IL6", "human----IL4", "human----CCL2", "human----CXCL10","human----CCL3", "human----CXCL8","human----CX3CL1", "human----TGFB1"), x = "cellassign_type",colour_by = 'cellassign_type')
# dev.off()
# 
# png(file=file.path(OUT_DIR, "total_violin_cito2.png"), width = dpi*20, height = dpi*14, units = "px",res = dpi,type='cairo')
# plotExpression(test.sce, features = c("human----IFNG","human----IFNA1", "human----TNF", "human----IL18","human----IL1B", "human----IL10","human----CCL5","human----IL6"), x = "cellassign_type",colour_by = 'cellassign_type')
# dev.off()
# 
# png(file=file.path(OUT_DIR, "total_violin_cito3.png"), width = dpi*20, height = dpi*14, units = "px",res = dpi,type='cairo')
# plotExpression(test.sce, features = c("human----IL4", "human----CCL2", "human----CXCL10","human----CCL3", "human----CXCL8","human----CX3CL1", "human----TGFB1"), x = "cellassign_type",colour_by = 'cellassign_type')
# dev.off()
# 
# png(file=file.path(OUT_DIR, "total_dots_fc.png"), width = dpi*20, height = dpi*14, units = "px",res = dpi,type='cairo')
# plotDots(test.sce, features = c("human----FCGRT","human----FCGR1A", "human----FGL2", "human----FCGR2B","human----FCAR","human----FCER1A","human----TRIM21", "human----FCGR3A", "human----FCGR3B"), group='cellassign_type')
# dev.off()
# 
# png(file=file.path(OUT_DIR, "total_dots_vs.png"), width = dpi*20, height = dpi*14, units = "px",res = dpi,type='cairo')
# plotDots(test.sce, features = c("human----B2M","human----FCER1G", "human----HCK", "human----SYK","human----RAF1", "human----MAPK8","human----MAPK11","human----MAPK1", "human----LYN", "human----FYN", "human----LAT","human----LCP2", "human----PTPN6", "human----PTPN11","human----INPP5D", "human----BTK","human----UBE2W","human----UBE2N", "human----UBE2V2", "human----PSMD14", "human----IKBKG","human----NFKB1", "human----NFKB2", "human----VCP","human----CGAS", "human----DDX58"), group='cellassign_type')
# dev.off()
# 
# png(file=file.path(OUT_DIR, "total_dots_cito2.png"), width = dpi*20, height = dpi*14, units = "px",res = dpi,type='cairo')
# plotDots(test.sce, features = c("human----IFNG", "human----TNF", "human----IL18","human----IL1B", "human----IL10","human----CCL5","human----IL6", "human----IL4", "human----CCL2", "human----CXCL10","human----CCL3", "human----CXCL8","human----CX3CL1", "human----TGFB1"), group='cellassign_type')
# dev.off()
# 
# png(file=file.path(OUT_DIR, "total_heatmap_fc.png"), width = dpi*20, height = dpi*14, units = "px",res = dpi,type='cairo')
# plotHeatmap(test.sce, features = c("human----FCGRT","human----FCGR1A", "human----FGL2", "human----FCGR2B","human----FCAR","human----FCER1A","human----TRIM21", "human----FCGR3A", "human----FCGR3B"))
# dev.off()
# 
# png(file=file.path(OUT_DIR, "total_heatmap_vs.png"), width = dpi*20, height = dpi*14, units = "px",res = dpi,type='cairo')
# plotHeatmap(test.sce, features = c("human----B2M","human----FCER1G", "human----HCK", "human----SYK","human----RAF1", "human----MAPK8","human----MAPK11","human----MAPK1", "human----LYN", "human----FYN", "human----LAT","human----LCP2", "human----PTPN6", "human----PTPN11","human----INPP5D", "human----BTK","human----UBE2W","human----UBE2N", "human----UBE2V2", "human----PSMD14", "human----IKBKG","human----NFKB1", "human----NFKB2", "human----VCP","human----CGAS", "human----DDX58"))
# dev.off()
# 
# png(file=file.path(OUT_DIR, "total_heatmap_cito.png"), width = dpi*20, height = dpi*14, units = "px",res = dpi,type='cairo')
# plotHeatmap(test.sce, features = c("human----IFNG","human----IFNA1", "human----TNF", "human----IL18","human----IL1B", "human----IL10","human----CCL5","human----IL6", "human----IL4", "human----CCL2", "human----CXCL10","human----CCL3", "human----CXCL8","human----CX3CL1", "human----TGFB1"))
# dev.off()

#Counting number of cellsTCD4Th1
TCD4Th1.only <- test.sce[, test.sce$cellassign_type == "TCD4Th1"]
ncol(counts(TCD4Th1.only))
print("TCD4Th1 - greater than zero expression of a gene - FCGR1A")
TCD4Th1.FCGR1A.only <- TCD4Th1.only[,which(assay(TCD4Th1.only)['human----FCGR1A',] > 0)]
ncol(counts(TCD4Th1.FCGR1A.only))
print("TCD4Th1 - greater than zero expression of a gene  - FCGR2A")
TCD4Th1.FCGR2A.only <- TCD4Th1.only[,which(assay(TCD4Th1.only)['human----FCGR2A',] > 0)]
ncol(counts(TCD4Th1.FCGR2A.only))
print("TCD4Th1 - greater than zero expression of a gene  - FCGR2B")
TCD4Th1.FCGR2B.only <- TCD4Th1.only[,which(assay(TCD4Th1.only)['human----FCGR2B',] > 0)]
ncol(counts(TCD4Th1.FCGR2B.only))
print("TCD4Th1 - greater than zero expression of a gene  - FCGR3A")
TCD4Th1.FCGR3A.only <- TCD4Th1.only[,which(assay(TCD4Th1.only)['human----FCGR3A',] > 0)]
ncol(counts(TCD4Th1.FCGR3A.only))
print("TCD4Th1 - greater than zero expression of a gene  - FCGR3B")
TCD4Th1.FCGR3B.only <- TCD4Th1.only[,which(assay(TCD4Th1.only)['human----FCGR3B',] > 0)]
ncol(counts(TCD4Th1.FCGR3B.only))
