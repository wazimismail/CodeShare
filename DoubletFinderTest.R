
library("DoubletFinder")

## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
counts <- Read10X_h5("filtered_feature_bc_matrix.h5")
obj <- CreateSeuratObject(counts = counts, assay = 'RNA', project = '10x_RNA')
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = 3000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 100)
obj <- RunTSNE(obj, dims = 1:30)
obj <- RunUMAP(obj, dims = 1:30)
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list <- paramSweep_v3(obj, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(obj@meta.data$seurat_clusters) ## ex: annotations <- obj@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(obj@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
obj <- doubletFinder_v3(obj, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)


