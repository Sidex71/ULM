#' adding cell assignments to seurat object
#' @description
#'  'AddMtaObject()' a function to add the predicted cell labels to the metadata of scRNAseq
#' @param seurat_obj a prepossessed Seurat object storing the scRNAseq data
#' @param cell_class_df a data frame of cell barcodes and cell type assignments, ideally obtained as the output from the GetCellAssignments() function.
#' @returns a new Seurat object with the updated metadata containing predicted cell labels in the "celltype_ulm" column
#' @examples
#' data(int_multData)
#' data(int_signature)
#' my_scores <- GetCellScores(seurat_obj = int_multData[,1:1000], signatures = int_signature, assay = 'RNA', slot = 'data')
#' my_ass <- GetCellAssignments(score_data = my_scores)
#' new_obj <- AddMetaObject(seurat_obj = int_multData[,1:1000], cell_class_df = my_ass)
#' head(new_obj$celltype_ulm, 20)
#'
#' @importFrom dplyr left_join
#' @importFrom tibble rownames_to_column column_to_rownames
#' @export
#'
#'
AddMetaObject <- function(seurat_obj, cell_class_df) {
  new_meta <- seurat_obj@meta.data %>% rownames_to_column('barcode_ulm') %>%
    left_join(cell_class_df, by= c('barcode_ulm'= 'barcode')) %>% column_to_rownames('barcode_ulm')
  seurat_obj@meta.data <- new_meta
  return(seurat_obj)
}
