#' scoring cells for gene signatures
#' @description
#'   'GetCellScores()' scores each cell in the scRNAseq data for cell type-specific gene signatures
#'
#' @param seurat_obj a prepossessed Seurat object storing the scRNAseq data
#' @param signatures a data frame of cell type signatures with cell types in the source column, genes in the target column, and weights in the mor column. Ideally the output from the GetSignature() function.
#' @param assay a character specifying assay in the Seurat object containing count matrix.
#' @param slot a character specifying the slot to draw counts from. This can be 'counts' for raw counts, 'data' for normalized counts, or 'scaled' for scaled counts.
#' @returns a data frame of cell barcodes and gene signature scores
#' @examples
#' data(int_singData)
#' data(int_signature)
#' my_scores <- GetCellScores(seurat_obj = int_singData[,1:1000], signatures = int_signature, assay = 'RNA', slot = 'data')
#' head(my_scores)
#'
#' @importFrom Seurat GetAssayData
#' @importFrom decoupleR run_ulm
#' @export
#'
#'
GetCellScores <- function(seurat_obj, signatures,  assay = 'RNA', slot = 'data'){
  mat <- GetAssayData(seurat_obj, assay = assay, slot = slot )
  acts <- run_ulm(mat=mat, net= signatures, .source='source', .target='target',
                  .mor='mor', minsize = 5)
  colnames(acts)[2:3] <- c('celltype', 'barcode')
  acts <- acts[, c('barcode', 'celltype', 'score', 'p_value', 'statistic')]
  return(acts)
}
