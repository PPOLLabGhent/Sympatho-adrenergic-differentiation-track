# Demultiplexing data (noALK) and integration with public data (GSE147821) from Kameneva et al. 2021

# Load libraries
# install.packages("dplyr") - V1.0.9
library(dplyr)
# install.packages("Seurat") - V4.1
library(Seurat)
# install.packages("patchwork") - V1.1.1
library(patchwork)
# install.packages("ggplot2") - V3.3.6
library(ggplot2)
# install.packages("hdf5r") - V1.3.5
library(hdf5r)


# Create lists to loop over
seurList <- list()
scDIRS <- c("CRSC0","CRSC1","CRSC2","CRSC3")
scnames <- c("sc0","sc1","sc2","sc3")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

for (i in 1:4){

# Load in the UMI matrix
rna <- Read10X(paste(scDIRS[i],"/RNA/",sep=""))
# Load the HTO matrix
hto <- Read10X(paste(scDIRS[i],"/HTO/",sep=""))
# Select cell barcodes detected by both RNA and HTO
joint.bcs <- intersect(colnames(rna), colnames(hto))
# Subset RNA and HTO counts by joint cell barcodes
rna <- rna[, joint.bcs]
hto <- as.matrix(hto[, joint.bcs])

# Confirm that the HTO have the correct names
#rownames(hto)

# Setup Seurat object
pbmc.hashtag <- CreateSeuratObject(counts = rna,project = scnames[i])

# Add HTO data as a new assay independent from RNA
pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = hto)
# Cleanup
rm(hto)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")


# QC and filtering
pbmc.hashtag[["percent.mt"]] <- PercentageFeatureSet(pbmc.hashtag, pattern = "^MT-")

#subset seuratObject
# Set hard thresholds
# Hard thresholds after examining violinplots mt>10% , 300<feat<4000
pbmc.hashtag <- subset(pbmc.hashtag, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 10)
# Cleanup
rm(rna)
# Normalize RNA data with log normalization
pbmc.hashtag <- NormalizeData(pbmc.hashtag)
# Find and scale variable features
#pbmc.hashtag <- FindVariableFeatures(pbmc.hashtag, selection.method = "mean.var.plot")
pbmc.hashtag <- FindVariableFeatures(pbmc.hashtag,verbose=FALSE)
#pbmc.hashtag <- ScaleData(pbmc.hashtag, features = VariableFeatures(pbmc.hashtag))


# You can also play with additional parameters (see documentation for HTODemux())
# to adjust the threshold for classification Here we are using the default settings
pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.99)

print("table should be coming")
# Global classification results
print(table(pbmc.hashtag$HTO_classification.global))
print(table(pbmc.hashtag$HTO_classification))
# Group cells based on the max HTO signal
Idents(pbmc.hashtag) <- "HTO_maxID"
Idents(pbmc.hashtag) <- "HTO_classification.global"


# First, we will remove negative cells from the object
#  -> Already filtered so no negative cells in the dataset
pbmc.hashtag.subset <- pbmc.hashtag
# Extract the singlets
seurList[[i]]  <- subset(pbmc.hashtag, idents = "Singlet")
#cleanup
rm(pbmc.hashtag)
# Add CellCycleMarkers
seurList[[i]] <- CellCycleScoring(seurList[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

}

# Public data GSE147821
adaFiles <- list.files("GSE147821_RAW",full.names = T)
#adaNames <- c("w8_1","w9_63","w6_88","w14_123","w12_124","w8_125","w9_5","w11_6","w9_7","w8_16","w9_31_pGan","w12_35","w12_36_exAdr")
adaNames <- c("ada1","ada2","ada3","ada4","ada5","ada6","ada7","ada8","ada9","ada10","ada11","ada12","ada13")


for (i in 1:13){
  j=i+4
  rna <- Read10X_h5(adaFiles[i], use.names = TRUE, unique.features = TRUE)
  seurList[[j]] <- CreateSeuratObject(counts = rna,project = adaNames[i],min.cells = 3, min.features = 200)

  # QC filtering
  # Here I should make cutoffs equal for all datasets instead of choosing cutoff for each dataset separately
  seurList[[j]][["percent.mt"]] <- PercentageFeatureSet(seurList[[j]], pattern = "^MT-")
  seurList[[j]] <- subset(seurList[[j]], subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 10)

  # Normalize and identify variable features for each dataset independently
  seurList[[j]] <- NormalizeData(seurList[[j]])
  seurList[[j]] <- FindVariableFeatures(seurList[[j]])

}


# Select features that are repeatedly variable across datasets for integration run PCA on each dataset using these features
features <- SelectIntegrationFeatures(object.list = seurList)
seurList <- lapply(X = seurList, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

seur.anchors <- FindIntegrationAnchors(object.list = seurList, anchor.features = features, reduction = "rpca")
seur.combined <- IntegrateData(anchorset = seur.anchors)
DefaultAssay(seur.combined) <- "integrated"

# Add CellCycleMarkers
seur.combined <- CellCycleScoring(seur.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)


# Run the standard workflow for visualization and clustering
seur.combined <- ScaleData(seur.combined, verbose = FALSE)
seur.combined <- RunPCA(seur.combined, verbose = FALSE)

pdf("ElbowPlot_full_ada_plus_Stef.pdf")
ElbowPlot(seur.combined)
dev.off()

seur.combined <- FindNeighbors(seur.combined, reduction = "pca", dims = 1:20)
seur.combined <- FindClusters(seur.combined, resolution = 0.6)
seur.combined <- RunUMAP(seur.combined, reduction = "pca", dims = 1:20)

saveRDS(seur.combined,file="ada_and_stefNoALK_dim20_res06.rds")

seur.markers <- FindAllMarkers(seur.combined, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
write.table(seur.markers,"ada_and_stefNoALK_dim20_res06.txt",sep="\t",quote=F)

pdf("DimPlot_ada_plus_Stef_dim20.pdf")
DimPlot(seur.combined,label=TRUE) + theme(legend.position="none")
dev.off()
