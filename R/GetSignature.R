#' generating gene signature
#' @description
#'  'GetSignature()' generates cell type-specific gene signatures from scRNAseq data
#' @param seurat_obj a prepossessed Seurat object storing the scRNAseq data
#' @param ident_col  a column in the Seurat object metadata which is a character vector storing cell names or labels. If not specified, the default ident of the Seurat object will be used.
#' @param n a numeric value specifying the number of genes to be used for each cell signature. Default is 100 genes per cell type.
#' @param p_val  a numeric value specifying the adjusted p-value cut-off
#' @returns a dataframe of cell type signatures.
#' @examples
#' data(int_singData)
#' int_sig <- GetSignature(seurat_obj = int_singData[,1:1000], ident_col = int_singData$Cell_Type)
#' head(int_sig)
#'
#' @importFrom Seurat Idents<- FindAllMarkers
#' @importFrom dplyr filter group_by slice_max select mutate
#' @export
#'
#'

GetSignature = function(seurat_obj, ident_col = NULL, n = 100, p_val = 0.05){
  if (!class(seurat_obj) %in% c('Seurat',  "SeuratObject") ) {
    stop('Input data is not a Seurat object, please provide a valid Seurat object')
  }
  if(!is.null(ident_col)){
    message(paste('using the specified seurat ident to generate signatures'))
    Idents(seurat_obj) <- ident_col
  }
  else{
    message('using the default seurat ident to generate signatures')
  }
  Markers <- FindAllMarkers(seurat_obj, only.pos = T)
  Markers <- Markers %>% filter(p_val_adj <0.05) %>% group_by(cluster) %>%
    slice_max(., order_by = avg_log2FC, n = n) %>% dplyr:: select(c(cluster, gene)) %>%
    mutate(mor= 1, cluster= as.character(cluster))
  colnames(Markers) <- c('source', 'target', 'mor')
  return(Markers)

}
