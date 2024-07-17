install.packages("remotes")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)
library(DoubletFinder)
install.packages("harmony")
library(glmGamPoi)
library(harmony)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(dplyr)

qc <- function(folder_path) {
  counts <- Read10X(data.dir = folder_path)
  seurat_obj <- CreateSeuratObject(counts = counts, min.cells = 3)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 100 & nFeature_RNA < quantile(seurat_obj@meta.data$nFeature_RNA, probs = 0.98) & percent.mt < 20)
  seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE) 
  seurat_obj <- RunPCA(seurat_obj, assay = "SCT")
  return(seurat_obj)
}

folder_path <- "C:/Users/mehak/Downloads/IDC"
subdirs <- list.dirs(folder_path, recursive = FALSE)

all_obj = list()
for (subdir in subdirs) {
  seurat_obj <- qc(subdir)
  all_obj[[subdir]] = seurat_obj 
}
names(all_obj) <- basename(subdirs)
all_obj

condition_labels <- c(rep("dcis", 7), rep("idc", 6))
condition_labels
sample_labels <- 1:13
sample_labels

for (i in 1:length(all_obj)) {
  cell_names <- colnames(all_obj[[i]])
  metadata <- data.frame(sample = rep(sample_labels[i], length(cell_names)), 
                         condition = rep(condition_labels[i], length(cell_names)), 
                         row.names = cell_names)
  all_obj[[i]] <- AddMetaData(all_obj[[i]], metadata = metadata)
}

#saveRDS(all_obj, file = "C:/Users/mehak/Desktop/seurat_objects_list.rds")
all_obj <- readRDS("C:/Users/mehak/Desktop/seurat_objects_list.rds")
all_obj

find_pk_dd <- function(seurat_obj) {
  sweep.res.list <- paramSweep(seurat_obj, PCs = 1:12, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  return(bcmvn)
}

pk_dcis1 <- find_pk_dd(all_obj$DCIS_1)
pk_dcis2 <- find_pk_dd(all_obj$DCIS_2)
pk_dcis3 <- find_pk_dd(all_obj$DCIS_3)
pk_dcis4 <- find_pk_dd(all_obj$DCIS_4)
pk_dcis5 <- find_pk_dd(all_obj$DCIS_5)
pk_dcis6 <- find_pk_dd(all_obj$DCIS_6)
pk_dcis7 <- find_pk_dd(all_obj$DCIS_7)
pk_idc1 <- find_pk_dd(all_obj$IDC_1)
pk_idc2 <- find_pk_dd(all_obj$IDC_2)
pk_idc3 <- find_pk_dd(all_obj$IDC_3)
pk_idc4 <- find_pk_dd(all_obj$IDC_4)
pk_idc5 <- find_pk_dd(all_obj$IDC_5)
pk_idc6 <- find_pk_dd(all_obj$IDC_6)

pk_df <- list(pk_dcis1, pk_dcis2, pk_dcis3, pk_dcis4, pk_dcis5, pk_dcis6, pk_dcis7, 
                    pk_idc1, pk_idc2, pk_idc3, pk_idc4, pk_idc5, pk_idc6)

optimal_pk <- list()
for (i in 1:length(pk_df)) {
  df <- pk_df[[i]]
  highest_value <- max(df$BCmetric, na.rm = TRUE)
  optimal_pk[[i]] <- highest_value
}
optimal_pk
optimal_pk <- readRDS("C:/Users/mehak/Desktop/optimal_pk_list.rds")

#saveRDS(optimal_pk, file = "C:/Users/mehak/Desktop/optimal_pk_list.rds")

ind_clustering <- function(seurat_obj) {
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:12)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:12)
  return(seurat_obj)
}

for (i in 1:length(all_obj)) {
  all_obj[[i]] <- ind_clustering(all_obj[[i]])
}

annotations_list <- list()
for (i in 1:length(all_obj)) {
  annotations <- all_obj[[i]]@meta.data$seurat_clusters
  annotations_list[[i]] <- annotations
}
annotations_list

#saveRDS(annotations_list, file = "C:/Users/mehak/Desktop/annotations_list.rds")
annotations_list <- readRDS("C:/Users/mehak/Desktop/annotations_list.rds")

run_doublet_finder <- function(seurat_obj, annotations, pk) {
  homotypic.prop <- modelHomotypic(annotations)     
  nExp_poi <- round(0.075*nrow(seurat_obj@meta.data))  
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  seurat_obj <- doubletFinder(seurat_obj, PCs = 1:12, pN = 0.25, pK = pk, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)
  return(seurat_obj)
}

#for (i in 1:length(all_obj)) {
#  all_obj[[i]] <- run_doublet_finder(all_obj[[i]], annotations_list[[i]], optimal_pk[[i]])
#}    computationally very costly

all_obj$DCIS_1 <- run_doublet_finder(all_obj$DCIS_1, annotations_list[[1]], optimal_pk[[1]])
all_obj$DCIS_2 <- run_doublet_finder(all_obj$DCIS_2, annotations_list[[2]], optimal_pk[[2]])
all_obj$DCIS_3 <- run_doublet_finder(all_obj$DCIS_3, annotations_list[[3]], optimal_pk[[3]])
all_obj$DCIS_4 <- run_doublet_finder(all_obj$DCIS_4, annotations_list[[4]], optimal_pk[[4]])
all_obj$DCIS_5 <- run_doublet_finder(all_obj$DCIS_5, annotations_list[[5]], optimal_pk[[5]])
all_obj$DCIS_6 <- run_doublet_finder(all_obj$DCIS_6, annotations_list[[6]], optimal_pk[[6]])
all_obj$DCIS_7 <- run_doublet_finder(all_obj$DCIS_7, annotations_list[[7]], optimal_pk[[7]])
all_obj$IDC_1 <- run_doublet_finder(all_obj$IDC_1, annotations_list[[8]], optimal_pk[[8]])
all_obj$IDC_2 <- run_doublet_finder(all_obj$IDC_2, annotations_list[[9]], optimal_pk[[9]])
all_obj$IDC_3 <- run_doublet_finder(all_obj$IDC_3, annotations_list[[10]], optimal_pk[[10]])
all_obj$IDC_4 <- run_doublet_finder(all_obj$IDC_4, annotations_list[[11]], optimal_pk[[11]])
all_obj$IDC_5 <- run_doublet_finder(all_obj$IDC_5, annotations_list[[12]], optimal_pk[[12]])
all_obj$IDC_6 <- run_doublet_finder(all_obj$IDC_6, annotations_list[[13]], optimal_pk[[13]])

#saveRDS(all_obj, file = "C:/Users/mehak/Desktop/doublet_seurat_objects_list.rds")
#DimPlot(all_obj$DCIS_2, reduction = "umap", group.by = "DF.classifications_0.25_6206.12755384966_549")  

all_obj_filtered <- all_obj
for (i in 1:length(all_obj)) {
  old_colname <- grep("^DF.classifications", colnames(all_obj[[i]]@meta.data), value = TRUE)
  colnames(all_obj[[i]]@meta.data)[colnames(all_obj[[i]]@meta.data) == old_colname] <- "doublet_status"
  all_obj_filtered[[i]] <- subset(all_obj[[i]], subset = doublet_status == "Singlet")
}
all_obj
all_obj_filtered

saveRDS(all_obj_filtered, file = "C:/Users/mehak/Desktop/filtered_seurat_obj_list.rds")

combined_seurat_obj <- merge(x = all_obj_filtered[[1]], y = all_obj_filtered[-1])
combined_seurat_obj_wo_doublets <- merge(x = all_obj[[1]], y = all_obj[-1])
integ_features <- SelectIntegrationFeatures(object.list = all_obj_filtered, nfeatures = 3000) 
integ_features_wo_doublets <- SelectIntegrationFeatures(object.list = all_obj, nfeatures = 3000) 
VariableFeatures(combined_seurat_obj) <- integ_features
VariableFeatures(combined_seurat_obj_wo_doublets) <- integ_features_wo_doublets
combined_seurat_obj_wo_doublets
combined_seurat_obj

combined_seurat_obj <- RunPCA(combined_seurat_obj, assay = "SCT")
combined_seurat_obj_wo_doublets <- RunPCA(combined_seurat_obj_wo_doublets, assay = "SCT")
ElbowPlot(combined_seurat_obj)
ElbowPlot(combined_seurat_obj_wo_doublets)

combined_seurat_obj <- RunHarmony(combined_seurat_obj, group.by.vars = "sample", dims = 1:12)
combined_seurat_obj_wo_doublets <- RunHarmony(combined_seurat_obj_wo_doublets, group.by.vars = "sample", dims = 1:12)
combined_seurat_obj <- RunUMAP(combined_seurat_obj, reduction = "harmony", dims = 1:12)
combined_seurat_obj_wo_doublets <- RunUMAP(combined_seurat_obj_wo_doublets, reduction = "harmony", dims = 1:12)
DimPlot(combined_seurat_obj, reduction = "umap", group.by = c("condition", "sample"))
DimPlot(combined_seurat_obj_wo_doublets, reduction = "umap", group.by = c("condition", "sample"))
combined_seurat_obj_wo_doublets@meta.data$doublet_status
combined_seurat_obj@meta.data$doublet_status
DimPlot(combined_seurat_obj_wo_doublets, reduction = "umap", group.by = "doublet_status")

combined_seurat_obj <- ind_clustering(combined_seurat_obj)
combined_seurat_obj_wo_doublets <- ind_clustering(combined_seurat_obj_wo_doublets)

DimPlot(combined_seurat_obj, reduction = "umap", group.by = "seurat_clusters",label = TRUE) + NoLegend()
DimPlot(combined_seurat_obj_wo_doublets, reduction = "umap", group.by = "seurat_clusters",label = TRUE) + NoLegend()

saveRDS(combined_seurat_obj, file = "C:/Users/mehak/Desktop/filtered_combined.rds")
saveRDS(combined_seurat_obj_wo_doublets, file = "C:/Users/mehak/Desktop/filtered_combined_wo_doublets.rds")

combined_seurat_obj <- readRDS("C:/Users/mehak/Desktop/filtered_combined.rds")

DefaultAssay(combined_seurat_obj) <- "RNA"
combined_seurat_obj <- NormalizeData(combined_seurat_obj, normalization.method = "LogNormalize")
combined_seurat_obj <- FindVariableFeatures(combined_seurat_obj, selection.method = "vst", nfeatures = 3000, layer = "data")
saveRDS(combined_seurat_obj, file = "C:/Users/mehak/Desktop/rna_markers_data.rds")
combined_seurat_obj <- readRDS("C:/Users/mehak/Desktop/rna_markers_data.rds")
combined_seurat_obj

get_conserved <- function(cluster) {
  conserved_markers <- FindConservedMarkers(
    combined_seurat_obj,
    ident.1 = cluster,
    grouping.var = "condition",
    only.pos = TRUE
  )
  conserved_markers <- rownames_to_column(conserved_markers, var = "gene")
  conserved_markers <- cbind(cluster_id = cluster, conserved_markers)
  return(conserved_markers)
}

library(purrr)
BiocManager::install('multtest')
install.packages('metap')
install.packages('devtools')
devtools::install_github('immunogenomics/presto')
library(metap)
library(tibble)

conserved_markers <- map_dfr(c(0:21), get_conserved)
head(conserved_markers)

saveRDS(conserved_markers, file = "C:/Users/mehak/Desktop/conserved_markers.rds")
conserved_markers <- readRDS("C:/Users/mehak/Desktop/conserved_markers.rds")
conserved_markers

top10 <- conserved_markers %>% 
  mutate(avg_fc = (idc_avg_log2FC + dcis_avg_log2FC) /2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, wt = avg_fc) %>%
  arrange(cluster_id, desc(avg_fc))
View(top10)

cluster_14_markers <- conserved_markers %>% filter(cluster_id == 14)
head(cluster_14_markers)
cluster_21_markers <- conserved_markers %>% filter(cluster_id == 21)
head(cluster_21_markers)

#Automatic annotation using reference

library(Seurat)
library(SingleR)
library(SingleCellExperiment)
library(scater)

data_path <- "C:/Users/mehak/Downloads/GSE176078_Wu_etal_2021_BRCA_scRNASeq/Wu_etal_2021_BRCA_scRNASeq"
counts <- Read10X(data_path, gene.column=1)
metadata <- read.csv("C:/Users/mehak/Desktop/metadata.csv")

counts_sparse <- as(counts, "sparseMatrix")
reference_sparse <- SingleCellExperiment(
  assays = list(counts = counts_sparse),
  colData = DataFrame(metadata)
)
reference_sparse <- logNormCounts(reference_sparse)
saveRDS(reference_sparse, file = "C:/Users/mehak/Desktop/reference.rds")
reference_sparse <- readRDS("C:/Users/mehak/Desktop/reference.rds")
reference_sparse

combined_seurat_obj <- readRDS("C:/Users/mehak/Desktop/rna_markers_data.rds")
query_counts <- GetAssayData(combined_seurat_obj, layer = 'counts')
query_counts
saveRDS(query_counts, file = "C:/Users/mehak/Desktop/query_counts.rds")
combined_seurat_obj
query_counts_sparse <- as(combined_seurat_obj@assays$RNA$counts, "sparseMatrix")
query_sparse <- SingleCellExperiment(
  assays = list(counts = query_counts_sparse),
  colData = DataFrame(combined_seurat_obj@meta.data)
)
query_sparse
query_sparse <- logNormCounts(query_sparse)
saveRDS(query_sparse, file = "C:/Users/mehak/Desktop/query.rds")
query_sparse

#Submitted a job to run this

data <- readRDS("C:/Users/mehak/Desktop/annotations.rds")
head(data)
combined_seurat_obj <- readRDS("C:/Users/mehak/Desktop/rna_markers_data.rds")
combined_seurat_obj$cell_type <- data$labels[match(rownames(combined_seurat_obj@meta.data), rownames(data))]
View(combined_seurat_obj$cell_type)
saveRDS(combined_seurat_obj, file = "C:/Users/mehak/Desktop/annotated_obj.rds")
combined_seurat_obj <- readRDS("C:/Users/mehak/Desktop/annotated_obj.rds")
DimPlot(combined_seurat_obj, reduction = 'umap', group.by = 'cell_type', label = TRUE)
View(combined_seurat_obj@meta.data)
combined_seurat_obj$celltype.condition <- paste0(combined_seurat_obj$cell_type,'_', combined_seurat_obj$condition)
DimPlot(combined_seurat_obj, reduction = 'umap', group.by = 'celltype.condition', split.by = 'condition', label = TRUE)

cell_counts <- combined_seurat_obj@meta.data %>%
  group_by(cell_type, condition) %>%
  summarise(count = n()) %>%
  ungroup()
View(cell_counts)

total_counts <- combined_seurat_obj@meta.data %>%
  group_by(cell_type) %>%
  summarise(total = n())
View(total_counts)

cell_proportions <- cell_counts %>%
  left_join(total_counts, by = "cell_type") %>%
  mutate(proportion = count / total)
View(cell_proportions)

library(ggplot2)
ggplot(cell_proportions, aes(x = proportion, y = cell_type, fill = condition)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("dcis" = "blue", "idc" = "red")) +
  labs(
    title = "Proportion of cells in each cluster",
    x = "Proportion of cells",
    y = NULL,
  )

