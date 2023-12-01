#### load packages ####
library(tidyverse)
library(Seurat)
options(Seurat.object.assay.version = "v5")

#### load data ####
cts <- read.csv(file = "data/GSE213755_expression_matrix.csv.gz", header = TRUE, row.names = 1)
meta <- read.csv(file = "data/GSE213755_metadata.csv.gz", header = TRUE, row.names = 1)

#### check meta data ####
str(meta)
table(meta$Cluster)
table(meta$Sample_ID)

# add annotation from the paper to the meta data
anno <- data.frame(
  Cluster = c(0:9),
  celltype = c("Neck","Chief1","Pit","Isthmus1","Isthmus2","Chief2","Endocrine1","Endocrine2","Parietal","Tuft"))
meta <- inner_join(meta, anno, by = "Cluster")
meta$GFP <- sapply(strsplit(meta$Sample_ID, "_"), "[", 1)

#### create seurat object ####
seu <- CreateSeuratObject(counts = cts, meta.data = meta, min.cells = 3, min.features = 200, project = "2023_Shiokawa")

#### QC ####
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
# seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 25)
# already filtered with 1500 < nFeature < 8000 and percent.mt < 15

#### cluster without integration ####
# seu <- JoinLayers(seu)
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
DimPlot(seu, group.by = "orig.ident") + NoAxes()
DimPlot(seu, group.by = "celltype") + NoAxes()
DimPlot(seu, group.by = "GFP") + NoAxes()

Idents(seu) <- "celltype"

saveRDS(seu, file = "RDSfiles/seu.RDS")

# FeaturePlot(seu,features = "Ptprc", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
# FeaturePlot(seu,features = "Pecam1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Mki67", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Stmn1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Muc5ac", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Muc6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Cblif", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Atp4b", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Chga", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Dclk1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Muc4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Cd44", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

#### monocle3 ####
# https://rpubs.com/mahima_bose/Seurat_and_Monocle3_p
library(SeuratWrappers)
library(monocle3)
seu <- subset(seu, subset = celltype %in% c("Endocrine1", "Endocrine2", "Parietal", "Tuft"), invert = TRUE)
cds <- as.cell_data_set(seu)
fData(cds)$gene_short_name <- rownames(fData(cds))

# assign partitions
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

# assign clusters
list.cluster <- seu@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

# assign UMAP coordinates
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- seu@reductions$umap@cell.embeddings

# learn trajectory
cds <- learn_graph(cds, use_partition = F)
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = F, label_roots = T, label_leaves = F,
           group_label_size = 5, graph_label_size = 0, cell_size = 0.75) + NoAxes()

# order cells in pseudotime
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) == "Isthmus2"]))
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.75) + NoAxes()

#### CytoTRACE ####
library(CytoTRACE)
counts_matrix <- LayerData(seu, assay='RNA', layer='counts') %>% as.data.frame()
obj_cell_type_anno <- as.data.frame(seu@meta.data$celltype)
results <- CytoTRACE(counts_matrix, ncores = 4)
pheno <- as.character(seu@meta.data$celltype)
names(pheno) <- colnames(seu)
plotCytoTRACE(results, phenotype = pheno)


