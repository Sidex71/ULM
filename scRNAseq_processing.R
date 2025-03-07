########################DC-T data###########################
library(tidyverse)
library(Seurat)
setwd("~/GSE135382_RAW")
#my_obj <- read.table('GSM4007219_AB6032.txt')
#my_obj <- CreateSeuratObject(my_obj)
my_files <- list.files(getwd())
my_len <- c()
obj_all <- as.list(rep(NA, length(my_files)))
for (i in 1: length(my_files)){
  my_obj <- read.table(my_files[i])
  obj_all[[i]] <- my_obj
  my_rows <- length(rownames(my_obj))
  my_len <- c(my_len, my_rows)
}

obj_final <- Reduce(cbind, obj_all)
saveRDS(obj_final, 'obj_final.rds')
obj_meta <- read_tsv('/home/sodiq/GSE135382_metadata2.txt')
saveRDS(obj_meta, 'obj_meta.rds')
Giladi_Data <- CreateSeuratObject(obj_final, meta.data = obj_meta, min.features =500)
Giladi_Data[["percent.mt"]] <-PercentageFeatureSet(Giladi_Data, pattern = "^mt-")##mouse 'mt', human 'MT'
VlnPlot(Giladi_Data,features = c("nFeature_RNA",
                                 "nCount_RNA",
                                 "percent.mt"),ncol = 3)
Giladi_Data<-subset(Giladi_Data,subset = percent.mt <25 & nFeature_RNA < 5000)
Giladi_Data <- NormalizeData(Giladi_Data)
Giladi_Data <- FindVariableFeatures(Giladi_Data)
Giladi_Data <- ScaleData(Giladi_Data)
Giladi_Data <- RunPCA(Giladi_Data)
ElbowPlot(Giladi_Data, ndims = 50)
Giladi_Data<-FindNeighbors(Giladi_Data,dims = 1:20)
Giladi_Data<-FindClusters(Giladi_Data,resolution= 0.3)
Giladi_Data<-RunUMAP(Giladi_Data, dims = 1:20)
obj_meta <- obj_meta %>% column_to_rownames('well')%>% rownames_to_column('barcodes')
my_meta <- Giladi_Data@meta.data %>% rownames_to_column('barcodes') %>%
  left_join(obj_meta, 'barcodes') %>% column_to_rownames('barcodes') 
Giladi_Data@meta.data <- my_meta
DimPlot(Giladi_Data, reduction = 'umap', group.by = 'sorting_scheme', label=T)
table(Giladi_Data$tissue)
Idents(Giladi_Data) <- Giladi_Data$tissue
Giladi_Data[['RNA3']]<- as(object = Giladi_Data[["RNA"]], Class = "Assay")

#saveRDS(Giladi_Data, 'Giladi_Data.rds')
invitro_obj <- subset(Giladi_Data, idents= 'in vitro')
invivo_obj <- subset(Giladi_Data, idents = 'LN')
p1 <-DimPlot(invitro_obj, reduction = 'umap', group.by = 'sorting_scheme', label=T) + ggtitle('')

png("/mnt/8TB/users/shameed/shameed/Doublet predictions/figures/UMAP_invitro.png", width = 6, height = 4.5, units = 'in', res = 600)
p1
dev.off()

Idents(invitro_obj) <- invitro_obj$sorting_scheme
marker_sub <- c('Cd3d', 'Cd3e', 'Cd3g','Cd4', 'Trac') 
p1 <-DoHeatmap(invitro_obj, features = marker_sub, slot = 'data')
p1
png("/mnt/8TB/users/shameed/shameed/Doublet predictions/figures/T_markers_heatmap.png", width = 10, height = 6.5, units = 'in', res = 600)
p1
dev.off()

marker_sub <- c("Cd80", 'Cd83', 'Cd80', "Clec7a", "Cd86", "H2-DMb2", "H2-Aa", "H2-Ab1", "H2-Eb1")
p1 <-DoHeatmap(invitro_obj, features = marker_sub, slot = 'data')
png("/mnt/8TB/users/shameed/shameed/Doublet predictions/figures/DC_markers_heatmap.png", width = 12, height = 7.5, units = 'in', res = 600)
p1
dev.off()

####################Liver data#################################################
Paired <- read_table("/mnt/Data/shameed/Doublet predictions/liver/GSE108561_PairedCells_umitab.txt")
####Paired <- Paired %>% column_to_rownames('gene_symbol')
library(tools)
rownames(Paired) <- toTitleCase(rownames(Paired))
#saveRDS(Paired, 'Paired.rds')
PairedData <- CreateSeuratObject(Paired)
PairedData[["percent.mt"]] <-PercentageFeatureSet(PairedData, pattern = "^Mt-")
summary(PairedData$percent.mt)
PairedData <- NormalizeData(PairedData)
PairedData <- FindVariableFeatures(PairedData)
PairedData <- ScaleData(PairedData)
PairedData <- RunPCA(PairedData)
ElbowPlot(PairedData, ndims = 50)
PairedData<-FindNeighbors(PairedData,dims = 1:30)
PairedData<-FindClusters(PairedData,resolution= 0.3)
PairedData<-RunUMAP(PairedData, dims = 1:30)
DimPlot(PairedData, reduction = 'umap', label=T)
PairedData$CellType <- 'Hep_Endo'
PairedData[['RNA3']]<- as(object = PairedData[["RNA"]], Class = "Assay")
VlnPlot(PairedData, features = c('Cyp2e1','Apoa1',  'Krt18', 'Oit3', 'Dpp4', 'Fcgr2b'), raster = F)
#saveRDS(PairedData, 'PairedData.rds')

Liver <- read_table("/mnt/Data/shameed/Doublet predictions/liver/GSE166504_cell_raw_counts.20220204.txt")
liver_meta <- read_tsv("/mnt/Data/shameed/Doublet predictions/liver/GSE166504_cell_metadata.20220204.tsv")
rownames(liver_meta) <- paste0(liver_meta$FileName, '_', liver_meta$CellID)
rownames(Liver) <- Liver$NPC_15weeks_Animal1_Capture2_AAACCTGAGACAAAGG
LiverData <- CreateSeuratObject(Liver, meta.data = liver_meta)
LiverData@assays$RNA@layers$counts <-replace(LiverData@assays$RNA@layers$counts, 
                                             is.na(LiverData@assays$RNA@layers$counts), 0)
rownames(LiverData) <- Liver$NPC_15weeks_Animal1_Capture2_AAACCTGAGACAAAGG
LiverData[["percent.mt"]] <-PercentageFeatureSet(LiverData, pattern = "^mt-")
summary(LiverData$percent.mt)
LiverData <- NormalizeData(LiverData)
LiverData <- FindVariableFeatures(LiverData)
LiverData <- ScaleData(LiverData)
LiverData <- RunPCA(LiverData)
ElbowPlot(LiverData, ndims = 50)
LiverData<-FindNeighbors(LiverData,dims = 1:30)
LiverData<-FindClusters(LiverData,resolution= 0.3)
LiverData<-RunUMAP(LiverData, dims = 1:30)
DimPlot(LiverData, reduction = 'umap', label=T, group.by = 'CellType')
VlnPlot(LiverData, features = c('Dpp4', 'Bmp3', 'Oit3'), group.by = 'CellType', raster = F)
#saveRDS(raw_object, 'LiverData.rds')
LiverData$FileName <- liver_meta$FileName
LiverData$CellType <- liver_meta$CellType
LiverData$CellID <- liver_meta$CellID
LiverData[['RNA3']]<- as(object = LiverData[["RNA"]], Class = "Assay")
saveRDS(LiverData, 'LiverData.rds')

####################Intestine datasets#####################################
#############Andrews et al singlet
int_sing <- read_tsv("/mnt/Data/shameed/Doublet predictions/new data/Newmans/raw/GSM5343586_SI_singlet_counts.tsv")
int_sing_met <-read_csv("/mnt/Data/shameed/Doublet predictions/new data/Newmans/raw/GSM5343586_GEO_SI_Singlet.csv")
head(int_sing[1:3, 1:3])
my_barcodes <- colnames(int_sing)
int_sing <- column_to_rownames(int_sing, 'AAACGAAAGAGGTCGT')
data_class <-sapply(int_sing, class)
data_class <- data_frame(data_class)
table(data_class)
data_class$cells <- colnames(int_sing)
which(data_class$data_class=='character') ##column 5278
bad <- int_sing[5278 ]
colnames(bad) <- 'badcol'
bad <- separate(bad, col =badcol, sep = '', into = c('empty', 'first', 'empty_2', 'second'))
bad <- bad %>% select(c(first, second))
bad <- as.data.frame(bad)
head(bad)
bad$first <-as.numeric(bad$first) %>% replace_na(0)
bad$second <- as.numeric(bad$second) %>% replace_na(0)
sum(is.na(bad))
int_sing <- int_sing[-5278]
int_sing <- cbind(int_sing, bad)
colnames(int_sing) <- my_barcodes
sum(is.na(int_sing))
int_sing_met <- int_sing_met[int_sing_met$Sample %in% colnames(int_sing),]
int_sing_met <- int_sing_met %>% column_to_rownames('Sample')
int_sing_met <- int_sing_met %>% select(-...1)
int_singData <- CreateSeuratObject(int_sing, meta.data = int_sing_met)
DefaultAssay(int_singData) <- 'RNA'
int_singData <- NormalizeData(int_singData)
int_singData <- FindVariableFeatures(int_singData, selection.method = 'vst', nfeatures = 2000)
int_singData <- ScaleData(int_singData)
int_singData <- RunPCA(int_singData)
ElbowPlot(int_singData, ndims = 50)
int_singData<-FindNeighbors(int_singData,dims = 1:20)
int_singData<-FindClusters(int_singData,resolution= 0.3)
int_singData<-RunUMAP(int_singData, dims = 1:20)
DimPlot(int_singData, reduction = 'umap', label=T, group.by = 'Cell_Type')
VlnPlot(int_singData, features = c('Cd3e', 'Cd79a', 'Cd3e', 'Cd19'), group.by = 'Cell_Type')


int_singData[['RNA3']]<- as(object = int_singData[["RNA"]], Class = "Assay")

saveRDS(int_singData, "/mnt/Data/shameed/Doublet predictions/new data/Newmans/int_singData.rds")

###########Andrews et al multiplet
int_mult <- read_tsv("/mnt/8TB/users/shameed/shameed/Doublet predictions/new data/Newmans/raw/GSM5343587_SI_multiplet_counts.tsv")
int_mult_met <-read_csv("/mnt/8TB/users/shameed/shameed/Doublet predictions/new data/Newmans/raw/GSM5343587_GEO_SI_Multiplet.csv")
head(int_mult[1:3, 1:3])
my_barcodes <- colnames(int_mult)
int_mult <- column_to_rownames(int_mult, 'AAACCCACAAATCGTC')
data_class <-sapply(int_mult, class)
data_class <- data_frame(data_class)
table(data_class)
data_class$cells <- colnames(int_mult)
which(data_class$data_class=='character') ###column 3670
bad <- int_mult[3670]
colnames(bad) <- 'badcol'
bad <- separate(bad, col =badcol, sep = '', into = c('empty', 'first', 'empty_2', 'second'))
bad <- bad %>% select(c(first, second))
bad <- as.data.frame(bad)
head(bad)
bad$first <-as.numeric(bad$first) %>% replace_na(0)
bad$second <- as.numeric(bad$second) %>% replace_na(0)
sum(is.na(bad))
int_mult <- int_mult[-3670]
int_mult <- cbind(int_mult, bad)
colnames(int_mult) <- my_barcodes
sum(is.na(int_mult))
int_mult_met <- int_mult_met[int_mult_met$Sample %in% colnames(int_mult),]
int_mult_met <- int_mult_met %>% column_to_rownames('Sample')
int_mult_met <- int_mult_met %>% select(-...1)
int_multData <- CreateSeuratObject(int_mult, meta.data = int_mult_met)
DefaultAssay(int_multData) <- 'RNA'
int_multData <- NormalizeData(int_multData)
int_multData <- FindVariableFeatures(int_multData, selection.method = 'vst', nfeatures = 2000)
int_multData <- ScaleData(int_multData)
int_multData <- RunPCA(int_multData)
ElbowPlot(int_multData, ndims = 50)
int_multData<-FindNeighbors(int_multData,dims = 1:20)
int_multData<-FindClusters(int_multData,resolution= 0.3)
int_multData<-RunUMAP(int_multData, dims = 1:20)
DimPlot(int_multData, reduction = 'umap', label=T)
VlnPlot(int_multData, features = c('Cd3e', 'Cd79a', 'Cd3e', 'Cd19'))

int_multData[['RNA3']]<- as(object = int_multData[["RNA"]], Class = "Assay")

saveRDS(int_multData, "/mnt/Data/shameed/Doublet predictions/new data/Newmans/int_multData.rds")

###############Manco et al multiplet data
clumps <- read.table("/mnt/Data/shameed/Doublet predictions/new data/clumpseq/GSE154714_Clumps_Seq_clumps_umitab.txt")
clumpsData <- CreateSeuratObject(clumps, min.features = 200)
clumpsData[["percent.mt"]] <-PercentageFeatureSet(clumpsData, pattern = "^mt-")
VlnPlot(clumpsData,features = c("nFeature_RNA",
                                "nCount_RNA",
                                "percent.mt"),ncol = 3)
clumpsData<-subset(clumpsData,subset = nFeature_RNA >200
                   & percent.mt < 30)
clumpsData <- NormalizeData(clumpsData)
clumpsData <- FindVariableFeatures(clumpsData)
clumpsData <- ScaleData(clumpsData)
clumpsData <- RunPCA(clumpsData)
ElbowPlot(clumpsData, ndims = 50)
clumpsData<-FindNeighbors(clumpsData,dims = 1:18)
clumpsData<-FindClusters(clumpsData,resolution= 0.3)
clumpsData<-RunUMAP(clumpsData, dims = 1:18)
DimPlot(clumpsData, reduction = 'umap', label=T)
VlnPlot(clumpsData, features = c('Dpp4', 'Bmp3', 'Oit3'),raster = F)

clumpsData[['RNA3']]<- as(object = clumpsData[["RNA"]], Class = "Assay")

clumpsData$num_cells <- substr(colnames(clumpsData), 19, 20)
clumpsData$num_cells <- ifelse(clumpsData$num_cells=='4n', '2', '> 2')
saveRDS(clumpsData,"/mnt/Data/shameed/Doublet predictions/new data/clumpseq/clumpsData.rds")

############################lung datasets##################################
##################singlet
lung_sing <- read_tsv("/mnt/8TB/users/shameed/shameed/Doublet predictions/new data/Newmans/raw/GSM5343588_Lung_singlet_counts.tsv")
lung_sing_met <-read_csv("/mnt/8TB/users/shameed/shameed/Doublet predictions/new data/Newmans/raw/GSM5343588_GEO_Lung_Singlet.csv")
head(lung_sing[1:3, 1:3])
my_barcodes <- colnames(lung_sing)
lung_sing <- column_to_rownames(lung_sing, 'AAACCCACAAAGGTTA')
data_class <-sapply(lung_sing, class)
data_class <- data_frame(data_class)
table(data_class)
data_class$cells <- colnames(lung_sing)
which(data_class$data_class=='character') ##column 6083
bad <- lung_sing[6083]
colnames(bad) <- 'badcol'
bad <- separate(bad, col =badcol, sep = '', into = c('empty', 'first', 'empty_2', 'second'))
bad <- bad %>%dplyr:: select(c(first, second))
bad <- as.data.frame(bad)
head(bad)
bad$first <-as.numeric(bad$first) %>% replace_na(0)
bad$second <- as.numeric(bad$second) %>% replace_na(0)
sum(is.na(bad))
lung_sing <- lung_sing[-6083]
lung_sing <- cbind(lung_sing, bad)
colnames(lung_sing) <- my_barcodes
sum(is.na(lung_sing))
lung_sing_met <- lung_sing_met[lung_sing_met$Sample %in% colnames(lung_sing),]
lung_sing_met <- lung_sing_met %>% column_to_rownames('Sample')
lung_sing_met <- lung_sing_met %>% select(-...1)
lung_singData <- CreateSeuratObject(lung_sing, meta.data = lung_sing_met)
DefaultAssay(lung_singData) <- 'RNA'
lung_singData <- NormalizeData(lung_singData)
lung_singData <- FindVariableFeatures(lung_singData, selection.method = 'vst', nfeatures = 2000)
lung_singData <- ScaleData(lung_singData)
lung_singData <- RunPCA(lung_singData)
ElbowPlot(lung_singData, ndims = 50)
lung_singData<-FindNeighbors(lung_singData,dims = 1:20)
lung_singData<-FindClusters(lung_singData,resolution= 0.3)
lung_singData<-RunUMAP(lung_singData, dims = 1:20)
DimPlot(lung_singData, reduction = 'umap', label=T, group.by = 'Cell_Type')
VlnPlot(lung_singData, features = c('Cd3e', 'Cd79a', 'Cd3e', 'Cd19'), group.by = 'Cell_Type')

lung_singData[['RNA3']]<- as(object = lung_singData[["RNA"]], Class = "Assay")

saveRDS(lung_singData, "/mnt/8TB/users/shameed/shameed/Doublet predictions/new data/Newmans/lung_singData.rds")

##############multiplet

lung_mult <- read_tsv("/mnt/8TB/users/shameed/shameed/Doublet predictions/new data/Newmans/raw/GSM5343589_Lung_multiplet_counts.tsv")
lung_mult_met <-read_csv("/mnt/8TB/users/shameed/shameed/Doublet predictions/new data/Newmans/raw/GSM5343589_GEO_Lung_Multiplet.csv")
head(lung_mult[1:3, 1:3])
my_barcodes <- colnames(lung_mult)
lung_mult <- column_to_rownames(lung_mult, 'AAACCCACAACGTAAA')
data_class <-sapply(lung_mult, class)
data_class <- data_frame(data_class)
table(data_class)
data_class$cells <- colnames(lung_mult)
which(data_class$data_class=='character') ##column 4728
bad <- lung_mult[4728]
colnames(bad) <- 'badcol'
bad <- separate(bad, col =badcol, sep = '', into = c('empty', 'first', 'empty_2', 'second'))
bad <- bad %>%dplyr:: select(c(first, second))
bad <- as.data.frame(bad)
head(bad)
bad$first <-as.numeric(bad$first) %>% replace_na(0)
bad$second <- as.numeric(bad$second) %>% replace_na(0)
sum(is.na(bad))
lung_mult <- lung_mult[-4728]
lung_mult <- cbind(lung_mult, bad)
colnames(lung_mult) <- my_barcodes
sum(is.na(lung_mult))
lung_mult_met <- lung_mult_met[lung_mult_met$Sample %in% colnames(lung_mult),]
lung_mult_met <- lung_mult_met %>% column_to_rownames('Sample')
lung_mult_met <- lung_mult_met %>%dplyr:: select(-...1)
lung_multData <- CreateSeuratObject(lung_mult, meta.data = lung_mult_met)
DefaultAssay(lung_multData) <- 'RNA'
lung_multData <- NormalizeData(lung_multData)
lung_multData <- FindVariableFeatures(lung_multData, selection.method = 'vst', nfeatures = 2000)
lung_multData <- ScaleData(lung_multData)
lung_multData <- RunPCA(lung_multData)
ElbowPlot(lung_multData, ndims = 50)
lung_multData<-FindNeighbors(lung_multData,dims = 1:25)
lung_multData<-FindClusters(lung_multData,resolution= 0.3)
lung_multData<-RunUMAP(lung_multData, dims = 1:25)
DimPlot(lung_multData, reduction = 'umap', label=T)
VlnPlot(lung_multData, features = c('Cd3e', 'Cd79a', 'Cd3e', 'Cd19'))

lung_multData[['RNA3']]<- as(object = lung_multData[["RNA"]], Class = "Assay")

saveRDS(lung_multData, "/mnt/8TB/users/shameed/shameed/Doublet predictions/new data/Newmans/lung_multData.rds")

##########################breast cancer data##########################################
setwd("/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/breast cancer/scRNAseq")
library(Seurat)
temp_mat <- ReadMtx(mtx = "count_matrix_sparse.mtx", cells = "count_matrix_barcodes.tsv", features = "count_matrix_genes.tsv", feature.column = 1)
metadata <- column_to_rownames(metadata, '...1')
sunnyData <- CreateSeuratObject(temp_mat, meta.data = metadata)
sunnyData[["percent.mt"]] <-PercentageFeatureSet(sunnyData, pattern = "^MT-")
sunnyData$project <- 'breast'
Idents(sunnyData) <- sunnyData$project
VlnPlot(sunnyData,features = c("nFeature_RNA",
                               "nCount_RNA",
                               "percent.mt"),ncol = 3)
#sunnyData<-subset(sunnyData,subset = nFeature_RNA >500
#                 & percent.mt < 15)
sunnyData <- NormalizeData(sunnyData)
sunnyData <- FindVariableFeatures(sunnyData)
sunnyData <- ScaleData(sunnyData)
sunnyData <- RunPCA(sunnyData)
ElbowPlot(sunnyData, ndims = 50)
sunnyData<-FindNeighbors(sunnyData,dims = 1:30)
sunnyData<-FindClusters(sunnyData,resolution= 0.35)
sunnyData<-RunUMAP(sunnyData, dims = 1:30)
DimPlot(sunnyData, reduction = 'umap', label=T, group.by = 'celltype_major')
VlnPlot(sunnyData, features = c('CD79B', 'CD14', 'CD19'),raster = F, group.by = 'celltype_major')

sunnyData[['RNA3']]<- as(object = sunnyData[["RNA"]], Class = "Assay")

#saveRDS(sunnyData, 'sunnyData.rds')
p1 <-DimPlot(sunnyData, reduction = 'umap', group.by = 'celltype_major', label=T) + ggtitle('')
png("/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/breast cancer/plots/UMAP.png", width = 8, height = 5.5, units = 'in', res = 600)
p1
dev.off()

###########################Ovarian Cancer data#######################################################
setwd("/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/ovarian cancer/ScRNA")
my_files <- list.files()

Elena_mat <- lapply(my_files, function(x){
  temp_mat <- read.table(x)
  return(temp_mat)
}) 

Elena_mat <- do.call(cbind, Elena_mat)
setwd("/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/ovarian cancer")
Elena_meta <- read_excel("Elena_meta.xlsx")
Elena_meta <- column_to_rownames(Elena_meta, 'Cell ID')
colnames(Elena_meta)[2] <- 'annotations'
ElenaData <- CreateSeuratObject(Elena_mat, min.cells = 3, min.features = 500, 
                                meta.data = Elena_meta)
ElenaData[["percent.mt"]] <-PercentageFeatureSet(ElenaData, pattern = "^MT-")
VlnPlot(ElenaData,features = c("nFeature_RNA",
                               "nCount_RNA",
                               "percent.mt"),ncol = 3)
ElenaData<-subset(ElenaData,subset = nFeature_RNA >500
                  & percent.mt < 15)
ElenaData <- NormalizeData(ElenaData)
ElenaData <- FindVariableFeatures(ElenaData)
ElenaData <- ScaleData(ElenaData)
ElenaData <- RunPCA(ElenaData)
ElbowPlot(ElenaData, ndims = 50)
ElenaData<-FindNeighbors(ElenaData,dims = 1:30)
ElenaData<-FindClusters(ElenaData,resolution= 0.35)
ElenaData<-RunUMAP(ElenaData, dims = 1:30)
DimPlot(ElenaData, reduction = 'umap', label=T, group.by = 'annotations')
VlnPlot(ElenaData, features = c('CD79B', 'CD14', 'CD19'),raster = F, group.by = 'annotations')

ElenaData[['RNA3']]<- as(object = ElenaData[["RNA"]], Class = "Assay")

temp_df <- data_frame(annotations = unique(ElenaData$annotations))
temp_df <- temp_df %>% mutate(cell_type = 
                                ifelse(str_detect(annotations, 
                                                  regex("Fibro", ignore_case = F)), 
                                       "fibroblasts", annotations))
temp_meta <- ElenaData@meta.data
temp_meta <- temp_meta %>%left_join(temp_df, by=  'annotations')
rownames(temp_meta) <- rownames(ElenaData@meta.data)
ElenaData@meta.data <- temp_meta

p1 <-DimPlot(ElenaData, reduction = 'umap', group.by = 'cell_type', label=T) + ggtitle('')
png("/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/ovarian cancer/plots/UMAP.png", width = 6, height = 4.5, units = 'in', res = 600)
p1
dev.off()
