library("DoubletFinder")
setwd("/mforge/research/labs/microbiome/Beyond_DNA/shared/cellranger_output/scRNA")

basedir <- "/mforge/research/labs/microbiome/Beyond_DNA/shared/cellranger_output/scRNA/"
dirnames <- c("072", "108", "120", "454", "570", "789", "890", "984", "L010-1", "L010-2", "L027", "M3399", "M4666", "M5167", "PM001", "cv-A5", "cv-A6", "cv-A12", "cv-A16", "cv-A25", "cv-A26", "cv-A30", "cv-B23", "cv-C1", "cv-C2", "cv-C3", "cv-C4")

doDF <- function(samplepath){
    ## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
    counts <- Read10X_h5(paste0(samplepath, "/outs/filtered_feature_bc_matrix.h5"))
    obj <- CreateSeuratObject(counts = counts, assay = 'RNA', project = '10x_RNA')
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    obj <- subset(obj, subset = percent.mt < 50 & nFeature_RNA > 200)
    obj <- obj %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% FindNeighbors() %>% FindClusters()

    ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
    sweep.res.list <- paramSweep_v3(obj, PCs = 1:10, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    temp <- bcmvn %>% slice_max(BCmetric) %>% select(pK) %>% unlist()
    pK <- as.numeric(levels(temp))[temp]

    ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
    homotypic.prop <- modelHomotypic(obj@meta.data$seurat_clusters) ## ex: annotations <- obj@meta.data$ClusteringResults
    nExp_poi <- round(0.075*nrow(obj@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

    pN <- 0.25
    ## Run DoubletFinder
    obj <- doubletFinder_v3(obj, PCs = 1:10, pN = pN, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)

    DF_tag <- paste0("DF.classifications_", as.character(pN), "_", as.character(pK), "_", as.character(nExp_poi.adj))

    scrublet <- read.table(paste0(samplepath, "/outs/filtered_feature_bc_matrix/scrublet.tsv"), sep = "\t", col.names = c('scrublet.observed', 'scrublet.simulated', 'scrublet.prediction'))
    
    return(data.frame(TotalCells = nrow(scrublet), DoubletFinder = (obj[[DF_tag]] %>% table())["Doublet"], Scrublet = (scrublet$scrublet.prediction %>% table())["True"]))
}

res <- lapply(paste0(basedir, dirnames), doDF)
names(res) <- dirnames
res <- bind_rows(res, .id = "Sample")

write.table(res, file = "DoubletFinder_Scrublet_compare", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)










