##########################################intestine data predictions###############################
###############Andrews et al singlet data
setwd("/mnt/8TB/users/shameed/shameed/Doublet predictions/new data/Newmans")
set.seed(101324)
int_sig <- GetSignature(int_singData, ident_col = int_singData$Cell_Type)
saveRDS(int_sig, '/mnt/8TB/users/shameed/shameed/Doublet predictions/new data/Newmans/int_sig.rds')
my_scores <- GetCellScores(seurat_obj = int_singData, signatures = int_sig, assay = 'RNA3', slot = 'data')
my_ass <- GetCellAssignments(score_data = my_scores)
my_scores_intSing <- my_scores
my_ass_intSing <- my_ass
int_singData <- AddMetaObject(int_singData, cell_class_df = my_ass)
colnames(int_singData@meta.data)
my_mult <- GetMultiplet(int_singData)
my_mult_filt <- FilterMultiplet(int_singData)
my_network_df <- GetNodeDF(mat = my_mult_filt$multSummaryFilt)
p1<-PlotNetwork(network_df = my_network_df, edge_width_factor = 10, node_text_size = 10,
                legend_text_size = 25, legend_title_size = 25, main_size =28 ) + theme(plot.background = element_rect(color = "black", size = 1, linetype = "dashed"))

p1 

#dir.create('plots')
png("/mnt/8TB/users/shameed/shameed/Doublet predictions/new data/Newmans/plots/intSing.png", width = 27, height = 13.5, units = 'in', res = 600)
p1
dev.off()

p1<- DimPlot(int_singData, reduction = 'umap', label=T, group.by = 'Cell_Type') + ggtitle('Intestinal Tissue')
png("/mnt/8TB/users/shameed/shameed/Doublet predictions/new data/Newmans/plots/umap_intSing.png", width = 8.5, height = 4.5, units = 'in', res = 600)
p1
dev.off()

###########################################Andrews et al multiplet data
set.seed(101324)
int_sig <- GetSignature(int_singData, ident_col = int_singData$Cell_Type)
my_scores <- GetCellScores(seurat_obj = int_multData, signatures = int_sig, assay = 'RNA3', slot = 'data')
my_ass <- GetCellAssignments(score_data = my_scores)
#int_multData@meta.data<- int_multData@meta.data %>% dplyr:: select(-c("count_ulm" , "celltype_ulm", "avg_pvalue", "avg_score"))
my_scores_intMult <- my_scores
my_ass_intMult <- my_ass
int_multData <- AddMetaObject(int_multData, cell_class_df = my_ass)
colnames(int_multData@meta.data)
my_mult <- GetMultiplet(int_multData)
my_mult_filt <- FilterMultiplet(int_multData,minFreq = 7)
my_network_df <- GetNodeDF(mat = my_mult_filt$multSummaryFilt)
p1<-PlotNetwork(network_df = my_network_df, edge_width_factor = 10, node_text_size = 8,
                legend_text_size = 22, legend_title_size = 22, main_size =25 )


png("/mnt/8TB/users/shameed/shameed/Doublet predictions/new data/Newmans/plots/intMult.png", width = 15, height = 10.5, units = 'in', res = 600)
p1
dev.off()

saveRDS(int_multData, 'int_multData.rds')
saveRDS(int_singData, 'int_singData.rds')

###################################Manco et al multiplet data 
setwd("/mnt/8TB/users/shameed/shameed/Doublet predictions/new data/clumpseq")
set.seed(101324)
int_sig <- GetSignature(int_singData, ident_col = int_singData$Cell_Type)

#####clumps/multiplets
set.seed(151024)
my_scores <- GetCellScores(seurat_obj = clumpsData, signatures = int_sig, assay = 'RNA3', slot = 'data')
my_ass <- GetCellAssignments(score_data = my_scores)
my_scores_intClump <- my_scores
my_ass_intClump <- my_ass
clumpsData <- AddMetaObject(clumpsData, cell_class_df = my_ass)
colnames(clumpsData@meta.data)
my_mult <- GetMultiplet(clumpsData)
my_mult_filt <- FilterMultiplet(clumpsData,minFreq = 10)
my_network_df <- GetNodeDF(mat = my_mult_filt$multSummaryFilt)
p1<-PlotNetwork(network_df = my_network_df, edge_width_factor = 10, node_text_size = 10,
                legend_text_size = 25, legend_title_size = 25, main_size =30 )

#dir.create('plots')
png("/mnt/8TB/users/shameed/shameed/Doublet predictions/new data/clumpseq/plots/Mult.png", width = 27.5, height = 14.5, units = 'in', res = 600)
p1
dev.off()

saveRDS(clumpsData, 'clumpsData.rds')

##################################lung data predictions###############################
setwd("/mnt/8TB/users/shameed/shameed/Doublet predictions/new data/Newmans")
set.seed(151024)
lung_sig <- GetSignature(lung_singData, ident_col = lung_singData$Cell_Type)
saveRDS(lung_sig, 'lung_sig.rds')
my_scores <- GetCellScores(seurat_obj = lung_singData, signatures = lung_sig, assay = 'RNA3', slot = 'data')
my_ass <- GetCellAssignments(score_data = my_scores)
my_scores_lungSing <- my_scores
my_ass_lungSing <- my_ass
lung_singData <- AddMetaObject(lung_singData, cell_class_df = my_ass)
colnames(lung_singData@meta.data)
my_mult <- GetMultiplet(lung_singData)
my_mult_filt <- FilterMultiplet(lung_singData, minFreq = 10)
my_network_df <- GetNodeDF(mat = my_mult_filt$multSummaryFilt)
p1<-PlotNetwork(network_df = my_network_df, edge_width_factor = 10, node_text_size = 10,
                legend_text_size = 24, legend_title_size = 24, main_size =27 ) +
  theme(plot.background = element_rect(color = "black", size = 1, linetype = "dashed"))

#dir.create('plots')
png("/mnt/8TB/users/shameed/shameed/Doublet predictions/new data/Newmans/plots/LungSing.png", width = 24, height = 12.5, units = 'in', res = 600)
p1
dev.off()

p1<- DimPlot(lung_singData, reduction = 'umap', label=T, group.by = 'Cell_Type', repel = T) + ggtitle('Lung Tissue')
png("/mnt/8TB/users/shameed/shameed/Doublet predictions/new data/Newmans/plots/umap_LungSing.png", width = 10.5, height = 5.5, units = 'in', res = 600)
p1
dev.off()

########################multiplet data
set.seed(1510242)
lung_sig <- GetSignature(lung_singData, ident_col = lung_singData$Cell_Type)
saveRDS(lung_sig, 'lung_sig.rds')
my_scores <- GetCellScores(seurat_obj = lung_multData, signatures = lung_sig, assay = 'RNA3', slot = 'data')
my_ass <- GetCellAssignments(score_data = my_scores)
my_scores_lungMult <- my_scores
my_ass_lungMult <- my_ass
lung_multData <- AddMetaObject(lung_multData, cell_class_df = my_ass)
colnames(lung_multData@meta.data)
my_mult <- GetMultiplet(lung_multData)
my_mult_filt <- FilterMultiplet(lung_multData)
my_network_df <- GetNodeDF(mat = my_mult_filt$multSummaryFilt)
p1<-PlotNetwork(network_df = my_network_df, edge_width_factor = 10, node_text_size = 10,
                legend_text_size = 25, legend_title_size = 25, main_size =30 )

#dir.create('plots')
png("/mnt/8TB/users/shameed/shameed/Doublet predictions/new data/Newmans/plots/LungMullt.png", width = 32, height = 12.5, units = 'in', res = 600)
p1
dev.off()

saveRDS(lung_singData, 'lung_singData.rds')
saveRDS(lung_multData, 'lung_multData.rds')

###################################ovarian cancer prediction#####################################
setwd("/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/ovarian cancer")
set.seed(1610242)
OC_sig <- GetSignature(ElenaData, ident_col = ElenaData$cell_type)
saveRDS(OC_sig, 'OC_sig.rds')
my_scores <- GetCellScores(seurat_obj = ElenaData, signatures = OC_sig, assay = 'RNA3', slot = 'data')
my_ass <- GetCellAssignments(score_data = my_scores, cut_off = 1)
my_scores_OC <- my_scores
my_ass_OC <- my_ass
#ElenaData@meta.data<- ElenaData@meta.data %>% dplyr:: select(-c("count_ulm" , "celltype_ulm", "avg_pvalue", "avg_score"))
ElenaData <- AddMetaObject(ElenaData, cell_class_df = my_ass)
colnames(ElenaData@meta.data)
# my_mult <- GetMultiplet(ElenaData)
my_mult_filt <- FilterMultiplet(ElenaData, minFreq = 5)
my_network_df <- GetNodeDF(mat = my_mult_filt$multSummaryFilt)
p1<-PlotNetwork(network_df = my_network_df, edge_width_factor = 11, node_text_size = 10,
                legend_text_size = 25, legend_title_size = 25, main_size =30 )

#dir.create('plots')
png("/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/ovarian cancer/plots/ulm_net.png", width = 25, height = 15.5, units = 'in', res = 600)
p1
dev.off()

saveRDS(ElenaData, 'ElenaData.rds')

################################breast cancer prediction##############################################
setwd("/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/breast cancer")

sunnyData$celltype_ref <- ifelse(sunnyData$celltype_major=='Plasmablasts' |
                                   sunnyData$celltype_major=='B-cells', 'B-Plasma',
                                 sunnyData$celltype_major)
table(sunnyData$celltype_ref)
set.seed(1610243)
BC_sig <- GetSignature(sunnyData, ident_col = sunnyData$celltype_ref)
saveRDS(BC_sig, 'BC_sig.rds')
my_scores <- GetCellScores(seurat_obj = sunnyData, signatures = BC_sig, assay = 'RNA3', slot = 'data')
my_ass <- GetCellAssignments(score_data = my_scores, cut_off = 1)
my_scores_BC <- my_scores
my_ass_BC <- my_ass
#sunnyData@meta.data<- sunnyData@meta.data %>% dplyr:: select(-c("count_ulm" , "celltype_ulm", "avg_pvalue", "avg_score"))
sunnyData <- AddMetaObject(sunnyData, cell_class_df = my_ass)
colnames(sunnyData@meta.data)
# my_mult <- GetMultiplet(sunnyData)
my_mult_filt <- FilterMultiplet(sunnyData)
my_network_df <- GetNodeDF(mat = my_mult_filt$multSummaryFilt)
p1<-PlotNetwork(network_df = my_network_df, edge_width_factor = 10, node_text_size = 8,
                legend_text_size = 22, legend_title_size = 22, main_size =25 )
p1
#dir.create('plots')
png("/mnt/8TB/users/shameed/shameed/Doublet predictions/spatial/breast cancer/plots/ulm_net.png", width = 24, height = 10.5, units = 'in', res = 600)
p1
dev.off()

saveRDS(sunnyData, 'sunnyData.rds')

