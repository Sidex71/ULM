#######spatial transcriptomics#############
devtools::install_github('satijalab/seurat-data')
library(devtools)
#install_github("satijalab/seurat", ref = "release/5.1.0")
setwd("/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/ovarian cancer/spatial data")

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)

###############Ovarian cancer ST data######################################
my_files <- list.files()
my_files <- substr(my_files, 1, 14)
my_files <- unique(my_files)


spatialObj_list <- lapply(seq_along(my_files), function (x){
  temp_mat <- ReadMtx(mtx =paste0(my_files[x], '_' ,"matrix.mtx"), 
                      cells =paste0(my_files[x], '_', "barcodes.tsv") , 
                      features = paste0(my_files[x], '_', "features.tsv"))
  slide_seq <- CreateSeuratObject(temp_mat, assay = 'Spatial')
  img <- Read10X_Image(paste0(my_files[x], '_', 'spatial', '/', 'spatial')) 
  image <- img[Cells(slide_seq)]
  DefaultAssay(object = image) <- 'Spatial'
  slide_seq[[paste0('slice', x)]] <- image
  slide_seq$origin <- my_files[x]
  slide_seq <- SCTransform(slide_seq, assay = "Spatial", verbose = FALSE)
  return(slide_seq)
})

spatialObj <- spatialObj_list[[1]]
for (i in 2:length(spatialObj_list)) {
  spatialObj <- merge(spatialObj, spatialObj_list[[i]])
}

setwd("/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/ovarian cancer")
#saveRDS(spatialObj_list, 'spatialObj_list.rds')
#saveRDS(spatialObj, 'spatialObj.rds')
plot1 <- VlnPlot(spatialObj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(spatialObj, features = "nCount_Spatial") + theme(legend.position = "right")
library(patchwork)
wrap_plots(plot1, plot2)

#spatialObj <- SCTransform(spatialObj, assay = "spatial", verbose = FALSE)
SpatialFeaturePlot(spatialObj, features = c("EPCAM", "CD8A"))
p1 <- SpatialFeaturePlot(spatialObj, features = "EPCAM", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(spatialObj, features = "CD8A", alpha = c(0.1, 1))
p1 / p2

DefaultAssay(spatialObj) <- "SCT"
my_var <- c()
for (i in spatialObj_list){
  my_var <- c(my_var,VariableFeatures(i))
}
VariableFeatures(spatialObj) <- my_var
spatialObj <- RunPCA(spatialObj, verbose = FALSE)
spatialObj <- FindNeighbors(spatialObj, dims = 1:30)
spatialObj <- FindClusters(spatialObj, verbose = FALSE)
spatialObj <- RunUMAP(spatialObj, dims = 1:30)
spatialObj$origin <- str_split(spatialObj$origin, '_', simplify =T)[,2]
DimPlot(spatialObj, reduction = "umap", label = TRUE, group.by = c("ident", "origin"))
SpatialDimPlot(spatialObj)

SpatialDimPlot(spatialObj, cells.highlight = CellsByIdentities(object = spatialObj, idents = c(2, 1, 4, 3,
                                                                                               5, 8)), facet.highlight = TRUE, ncol = 3)
#SpatialDimPlot(spatialObj, interactive = TRUE)
#SpatialFeaturePlot(spatialObj, features = "CD4", interactive = TRUE)
#LinkedDimPlot(spatialObj)

#de_markers <- FindMarkers(spatialObj, ident.1 = 5, ident.2 = 6)
#SpatialFeaturePlot(object = spatialObj, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)
#brain <- FindSpatiallyVariableFeatures(spatialObj, assay = "SCT", features = VariableFeatures(spatialObj)[1:1000],
#                                      selection.method = "moransi")
#spatialObj <- brain
#top.features <- head(SpatiallyVariableFeatures(spatialObj, selection.method = "moransi"), 6)
#SpatialFeaturePlot(spatialObj, features = top.features, ncol = 3, alpha = c(0.1, 1))


#################RCTD integration
#options(timeout = 600000000)
#devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
#install.packages("remotes")
#remotes::install_github("dmcable/RCTD")
library(spacexr)
#ElenaData@meta.data <- OC_metadata_2
ref <- ElenaData
Idents(ref) <- ElenaData$cell_type
counts <- ref[["RNA3"]]$counts
cluster <- as.factor(ref$cell_type)
cluster <- str_replace(cluster, 'B/Plasma', 'B_Plasma')
names(cluster) <- colnames(ref)
cluster<- as.factor(cluster)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)

# set up query with the RCTD function SpatialRNA
#slide.seq <- SeuratData::LoadData("ssHippo")
weight_list <- lapply(spatialObj_list, function(x){
  coords <- GetTissueCoordinates(x)
  colnames(coords) <- c("x", "y")
  coords[is.na(colnames(coords))] <- NULL
  counts <- x[['Spatial']]$counts
  query <- SpatialRNA(coords, counts, colSums(counts))
  RCTD <- create.RCTD(query, reference, max_cores = 40)
  RCTD <- run.RCTD(RCTD, doublet_mode = "full")
  weights <- RCTD@results$weights
  #norm_weights <- normalize_weights(weights)
  spot_labels <- as.data.frame(weights)
  return(spot_labels)
})

saveRDS(weight_list, 'weight_list_2.rds')

spot_labels <- do.call(rbind, weight_list)

saveRDS(spot_labels, 'spot_labels_2.rds')
#saveRDS(RCTD, 'RCTD_2.rds')
all(colnames(spatialObj) %in% rownames(spot_labels))
colnames(spatialObj) <- 
  spatialObj <- spatialObj[, colnames(spatialObj) %in% rownames(spot_labels)]
spatialObj@assays$RCTD <- CreateAssay5Object(as.matrix(t(spot_labels_2)))
DefaultAssay(spatialObj) <- 'RCTD'
#install.packages("Hmisc")

##########################colocalisation analysis
library(Hmisc)
result <- rcorr(as.matrix(spot_labels))
saveRDS(result, 'result_2.rds')
# Extract correlation coefficients and p-values
cor_matrix <- result$r
pval_matrix <- result$P
corrplot::corrplot(cor_matrix, type = 'upper',method = 'color',
                   addCoef.col = T)


png('/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/ovarian cancer/plots/corrplot_RCTD_2.png', width = 9, height = 9, units = 'in', res = 600)
corrplot::corrplot(cor_matrix, type = 'upper',method = 'color',
                   addCoef.col = T)
dev.off()

##################breast cancer ST data#######################
setwd("/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/breast cancer/spatial_data")
my_files <- list.files()
my_files <- str_split(my_files, '_', simplify = T)[,1]
my_files <- unique(my_files)

brObj_list <- lapply(seq_along(my_files), function(x){
  temp_mat <- ReadMtx(mtx =paste0(my_files[x], '_filtered_count_matrix/', "matrix.mtx.gz") , 
                      cells =paste0(my_files[x], '_filtered_count_matrix/', "barcodes.tsv.gz"), 
                      features =paste0(my_files[x], '_filtered_count_matrix/', "features.tsv.gz"), feature.column = 1)
  slide_seq <- CreateSeuratObject(temp_mat, assay = 'Spatial')
  meta <- read_csv(paste0(my_files[x], "_metadata.csv"))
  slide_seq$classification <- meta$Classification
  slide_seq$subtype <- meta$subtype
  img <- Read10X_Image(paste0(my_files[x], '_', 'spatial')) 
  image <- img[Cells(slide_seq)]
  DefaultAssay(object = image) <- 'Spatial'
  slide_seq[[paste0('slice', x)]] <- image
  slide_seq$origin <- my_files[x]
  slide_seq <- SCTransform(slide_seq, assay = "Spatial", verbose = FALSE)
  return(slide_seq)
})

BrspatialObj <- brObj_list[[1]]
for (i in 2:length(brObj_list)) {
  BrspatialObj <- merge(BrspatialObj, brObj_list[[i]])
}

setwd("/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/breast cancer")
saveRDS(BrspatialObj, 'BrspatialObj.rds')
saveRDS(brObj_list, 'brObj_list.rds')
VlnPlot(BrspatialObj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
SpatialFeaturePlot(BrspatialObj, features = "nCount_Spatial") + theme(legend.position = "right")


#BrspatialObj <- SCTransform(BrspatialObj, assay = "spatial", verbose = FALSE)
SpatialFeaturePlot(BrspatialObj, features = c("EPCAM", "CD8A"))
SpatialFeaturePlot(BrspatialObj, features = "EPCAM", pt.size.factor = 1)
SpatialFeaturePlot(BrspatialObj, features = "CD8A", alpha = c(0.1, 1))

DefaultAssay(BrspatialObj) <- "SCT"
my_var <- c()
for (i in brObj_list){
  my_var <- c(my_var,VariableFeatures(i))
}
VariableFeatures(BrspatialObj) <- my_var
BrspatialObj <- RunPCA(BrspatialObj, verbose = FALSE)
BrspatialObj <- FindNeighbors(BrspatialObj, dims = 1:30)
BrspatialObj <- FindClusters(BrspatialObj, verbose = FALSE)
BrspatialObj <- RunUMAP(BrspatialObj, dims = 1:30)
BrspatialObj$origin <- str_split(BrspatialObj$origin, '_', simplify =T)[,2]
DimPlot(BrspatialObj, reduction = "umap", label = TRUE, group.by = c("ident", "origin"))

#####################RCTD integration
####set up reference
library(spacexr)
library(Seurat)
library(tidyverse)
ref <- sunnyData
Idents(ref) <- sunnyData$celltype_major
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$celltype_major)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)

# set up query with the RCTD function SpatialRNA

#counts <- BrspatialObj[["Spatial"]]$counts
#BrspatialObj[['Spatial2']] <- as(object = BrspatialObj[["Spatial"]], Class = "Assay")
#DefaultAssay(BrspatialObj) <- 'Spatial2'
#dim(counts)
#coords <- GetTissueCoordinates(BrspatialObj, image = 'slice1')
#colnames(coords) <- c("x", "y")
#query <- SpatialRNA(coords, counts, colSums(counts))
#RCTD <- create.RCTD(query, reference, max_cores = 40)
#RCTD <- run.RCTD(RCTD, doublet_mode = "full")
#weights <- RCTD@results$weights
#norm_weights <- normalize_weights(weights)
#spot_labels <- as.data.frame(norm_weights)

weight_list <- lapply(brObj_list, function(x){
  coords <- GetTissueCoordinates(x)
  colnames(coords) <- c("x", "y")
  coords[is.na(colnames(coords))] <- NULL
  counts <- x[['Spatial']]$counts
  query <- SpatialRNA(coords, counts, colSums(counts))
  RCTD <- create.RCTD(query, reference, max_cores = 40)
  RCTD <- run.RCTD(RCTD, doublet_mode = "full")
  weights <- RCTD@results$weights
  #norm_weights <- normalize_weights(weights)
  spot_labels <- as.data.frame(weights)
  return(spot_labels)
})

#saveRDS(weight_list, 'weight_list.rds')

spot_labels <- do.call(rbind, weight_list)

#saveRDS(spot_labels, 'spot_labels.rds')

#########################################coloclisation analysi
#install.packages("Hmisc")
library(Hmisc)
result <- rcorr(as.matrix(spot_labels))
saveRDS(result, 'result.rds')
# Extract correlation coefficients and p-values
cor_matrix <- result$r
pval_matrix <- result$P
corrplot::corrplot(cor_matrix, type = 'upper',method = 'color',
                   addCoef.col = T)

#dir.create('plots')

png('/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/breast cancer/plots/corrplot_RCTD.png', width = 9, height = 9, units = 'in', res = 600)
corrplot::corrplot(cor_matrix, type = 'upper',method = 'color',
                   addCoef.col = T)
dev.off()

