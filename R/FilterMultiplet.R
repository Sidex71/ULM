#' filtering multiplets
#' @description
#'   'FilterMultiplet()' a function to filter multiplets of a predefined frequency
#' @param seurat_obj a Seurat object with the metadata containing predicted cell labels in the "celltype_ulm" column and the number of cells in the "count_ulm" column. Ideally the output from the AddMetObject() function or the multiplet Seurat object from the GetMultiplet() function
#' @param minCells a numeric value specifying the minimum number of cells. Default is 2 to include doublets and/or higher order multiplets
#' @param minFreq a numeric value specifying the minimum frequency of a multiplet type for it to be retained. Default is 10.
#' @returns  a filtered list containing a Seurat object of multiplets and a dataframe of multiplet distribution summary.
#' @examples
#'
#' data(int_multData)
#' data(int_signature)
#' my_scores <- GetCellScores(seurat_obj = int_multData[,1:1000], signatures = int_signature, assay = 'RNA', slot = 'data')
#' my_ass <- GetCellAssignments(score_data = my_scores)
#' new_obj <- AddMetaObject(seurat_obj = int_multData[,1:1000], cell_class_df = my_ass)
#' my_mult <- GetMultiplet(seurat_obj = new_obj)
#' my_mult_filt <- FilterMultiplet(seurat_obj = new_obj)
#' my_mult_filt
#'
#' @importFrom Seurat Idents<-
#' @export
#'
#'
#'

FilterMultiplet <- function(seurat_obj, minCells = 2, minFreq = 10){
  multObj <- subset(seurat_obj, subset = count_ulm >= minCells)
  multSummary <- as.data.frame(table(multObj$celltype_ulm))
  colnames(multSummary) <- c('multipletType', 'frequency')
  multSummary <- multSummary[multSummary$frequency >= minFreq,]
  Idents(multObj) <- multObj$celltype_ulm
  celltypes <- unique(multSummary$multipletType)
  multObj <- subset(multObj, idents = celltypes)
  return(list(multSummaryFilt = multSummary, multObjFilt = multObj))
}
