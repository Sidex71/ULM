Idents(sunnyData) <- sunnyData$celltype_ulm
plasma_T <- subset(sunnyData, idents= c('B-cells_T-cells','Plasmablasts_T-cells'))
plasma_T$cellsub <- 'Plasma-T cells'
plasma_T <- NormalizeData(plasma_T)
plasma_T <- FindVariableFeatures(plasma_T, selection.method = "vst", nfeatures = 2000)
plasma_T <- ScaleData(plasma_T)

plasma_T <- RunPCA(plasma_T)
#ElbowPlot(plasma_T)
plasma_T <- FindNeighbors(plasma_T, dims = 1:30)
plasma_T <- FindClusters(plasma_T, resolution = .4)
ElbowPlot(plasma_T)
plasma_T <- RunUMAP(plasma_T, dims = 1:10)
DimPlot(plasma_T, reduction = 'umap', label=T, group.by = 'seurat_clusters')
DimPlot(plasma_T, reduction = 'umap', group.by = 'celltype_ulm')
plasma_T$annot <- 'T-B doublets'
Idents(plasma_T) <- plasma_T$annot
VlnPlot(plasma_T, features = c('CD79A', 'CD3G'))

p1 <-DoHeatmap(plasma_T, features = c('CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'TRAC', 
                                      'IGKC', 'IGLV3-1', 'JCHAIN', 'CD24', 'CD27',
                                      'CD79A', 'CD79B'), slot = 'data',label = F)

p1
png("/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/breast cancer/plots/plasma_T.png", width = 8, height = 5.5, units = 'in', res = 600)
p1
dev.off()

#############################LR analysis##############################################################
LR_Rowman <- read_excel("/mnt/8TB/users/shameed/shameed/LR_Rowman.xlsx", 
                        sheet = "All.Pairs")

predf <- LR_Rowman %>% filter(Ligand.ApprovedSymbol %in% rownames(plasma_T)&
                                Receptor.ApprovedSymbol %in% rownames(plasma_T))
LRP <- unique(predf$Pair.Name)
length(LRP)
LRdf <- sapply(LRP, function(x){
  LRSplit <- str_split(x,'_')
  mylig <- LRSplit[[1]][1]
  myrec <- LRSplit[[1]][2]
  ligcell <- plasma_T@assays$RNA3@data[mylig,]
  ligmean <- mean(ligcell)
  reccell <- plasma_T@assays$RNA3@data[myrec,]
  recmean <- mean(reccell)
  ligexp <-as.numeric(ligcell > ligmean)
  recexp <- as.numeric(reccell > recmean)
  LRexp <- ligexp + recexp
})
LRdf <- LRdf ==2
LRP_enriched <- as.data.frame(colSums(LRdf))
colnames(LRP_enriched) <- 'Freq'
LRP_enriched$percent <- LRP_enriched$Freq * 100/length(plasma_T$orig.ident)
LRP_filt <- LRP_enriched[LRP_enriched$Freq >= 5,]  #130 LRPs

saveRDS(LRP_enriched, 'LRP_enriched.rds')
LRP_filt <-LRP_filt[order(LRP_filt$percent, decreasing = T),]
LRP_filt <- LRP_filt %>% rownames_to_column('LRP')
saveRDS(LRP_filt, 'LRP_filt.rds')
p1<-ggplot(LRP_filt[1:50,], aes(x=LRP, y=percent, fill = LRP)) + 
  geom_bar(stat = 'identity') + 
  labs(title = ' ',
       y= 'cell proportion (%)', x=NULL)+
  theme_bw() + coord_flip() + NoLegend()
png("/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/breast cancer/plots/LRP.png", width = 15, height = 10.5, units = 'in', res = 600)
p1
dev.off()

##########################spatial LRP##########################################
spot_labels <- readRDS("/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/breast cancer/spot_labels.rds")
head(spot_labels)
spot_labels <- normalize_weights(spot_labels)
head(spot_labels)
spot_labels <- as.data.frame(spot_labels)
Can_B_spots <- spot_labels[, c('Plasmablasts', 'B-cells', 'T-cells')]
Can_B_spots <- Can_B_spots %>% filter((Plasmablasts > 0.1 | `B-cells` >0.1) &
                                        `T-cells` > 0.1)
saveRDS(Can_B_spots, 'Can_B_spots.rds')
predf <- LR_Rowman %>% filter(Ligand.ApprovedSymbol %in% rownames(temp_spa)&
                                Receptor.ApprovedSymbol %in% rownames(temp_spa))

LRE <- LRP_enriched[LRP_enriched$Freq >= 5,]
rownames(LRE)
paired_lists <- Map(list, brObj_list, weight_list)

Spa_LRdf_list <- lapply(paired_lists, function(x){
  temp_spa <- x[[1]]
  tempw <- x[[2]]
  tempw <- as.data.frame(tempw)
  tempw <- normalize_weights(tempw)
  Can_B_spots <- tempw[, c('Plasmablasts', 'B-cells', 'T-cells')]
  Can_B_spots <- Can_B_spots %>% filter((Plasmablasts > 0.05 | `B-cells` >0.05) &
                                          `T-cells` > 0.05)
  
  spot_BC <- temp_spa[, rownames(Can_B_spots)]
  predf <- LR_Rowman %>% filter(Ligand.ApprovedSymbol %in% rownames(temp_spa)&
                                  Receptor.ApprovedSymbol %in% rownames(temp_spa))
  my_LRE <- intersect(rownames(LRE), predf$Pair.Name)
  DefaultAssay(spot_BC) <- 'SCT'
  predf <- LR_Rowman %>% filter(Ligand.ApprovedSymbol %in% rownames(temp_spa)&
                                  Receptor.ApprovedSymbol %in% rownames(temp_spa))
  
  my_LRE <- intersect(rownames(LRE), predf$Pair.Name)
  my_enrich <-sapply(my_LRE, function(x){
    LRSplit <- str_split(x,'_')
    mylig <- LRSplit[[1]][1]
    myrec <- LRSplit[[1]][2]
    ligcell <- spot_BC@assays$SCT@data[mylig,]
    ligmean <- mean(ligcell)
    reccell <- spot_BC@assays$SCT@data[myrec,]
    recmean <- mean(reccell)
    ligexp <-as.numeric(ligcell > ligmean)
    recexp <- as.numeric(reccell > recmean)
    LRexp <- ligexp + recexp
  })
  rownames(my_enrich) <- colnames(spot_BC)
  return(my_enrich)
})

names(Spa_LRdf_list) <- NULL

saveRDS(Spa_LRdf_list, 'Spa_LRdf_list.rds')

Spa_df_list <- lapply(Spa_LRdf_list, function(x){
  x <- x ==2
  LRP_enriched <- as.data.frame(colSums(x))
  colnames(LRP_enriched) <- 'Freq'
  LRP_enriched$percent <- LRP_enriched$Freq * 100/length(rownames(x))
  LRP_enriched <- LRP_enriched %>% rownames_to_column('LRP')
  return(LRP_enriched)
}) 

Spa_LRdf <- do.call(rbind, Spa_df_list)

total_spots <- sum(sapply(Spa_LRdf_list, function(x) length(rownames(x)))) #1173
Spa_LRdf <- Spa_LRdf %>% group_by(LRP) %>% mutate(total_freq = sum(Freq),
                                                  final_percent = sum(Freq)*100/total_spots)

Spa_LRdf_final <- Spa_LRdf %>% select(c(LRP, total_freq, final_percent))
Spa_LRdf_final <- unique(Spa_LRdf_final)

Spa_LRdf_final <-Spa_LRdf_final[order(Spa_LRdf_final$final_percent, decreasing = T),]

saveRDS(Spa_LRdf_final, 'Spa_LRdf_final.rds')

p1<-ggplot(Spa_LRdf_final[1:50,], aes(x=LRP, y=final_percent, fill = LRP)) + 
  geom_bar(stat = 'identity') + 
  labs(title = ' ',
       y= 'spot proportion (%)', x=NULL)+
  theme_bw() + coord_flip() + NoLegend()
png("/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/breast cancer/plots/Spatial_LRP.png", width = 15, height = 10.5, units = 'in', res = 600)
p1
dev.off()
