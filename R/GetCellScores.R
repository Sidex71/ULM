#' scoring cells for gene signatures
#' @description
#'   'GetCellScores()' scores each cell in the scRNAseq data for cell type-specific gene signatures
#'
#' @param seurat_obj a preprocessed Seurat object storing the scRNAseq data
#' @param signatures a data frame of cell type signatures with cell types in the source column, genes in the target column, and weights in the mor column. Ideally the output from the GetSignature() function.
#' @param assay a character specifying assay in the Seurat object containing count matrix.
#' @param layer a character specifying the layer to draw counts from (Seurat v5). This can be "counts" for raw counts, "data" for normalized counts, or "scaled" for scaled counts. Default is "data"
#' @param layer a character specifying the layer to draw counts from (Seurat v3/v4). This can be "counts" for raw counts, "data" for normalized counts, or "scaled" for scaled counts. Default is NULL.
#' @returns a data frame of cell barcodes and gene signature scores
#' @examples
#' data(int_singData)
#' data(int_signature)
#' my_scores <- GetCellScores(seurat_obj = int_singData[,1:1000],
#'                            signatures = int_signature,
#'                            assay = 'RNA',
#'                            layer = 'data')
#' head(my_scores)
#'
#' @importFrom Seurat GetAssayData 
#' @importFrom decoupleR run_ulm
#' @export
#'
GetCellScores <- function(seurat_obj, signatures, assay = 'RNA', slot = NULL, layer = 'data') {
  # Check for invalid argument combinations
  if (is.null(slot) & is.null(layer)) {
    stop("You must specify either `slot` (Seurat v3/v4) or `layer` (Seurat v5).")
  }
  if (!is.null(slot) & !is.null(layer)) {
    stop("Please specify only one of `slot` or `layer`, not both.
         For example, you need to set `layer = NULL` if you are using slot, and vice versa")
  }
  
  # Get assay data depending on Seurat version
  if (!is.null(slot)) {
    mat <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  } 
  else {
    mat <- GetAssayData(seurat_obj, assay = assay, layer = layer)
  } 
  
  
  # Score cells using Univariate Linear Models
  acts <- run_ulm(mat = mat, net = signatures, .source = "source", .target = "target", .mor = "mor", minsize = 5)
  
  colnames(acts)[2:3] <- c("celltype", "barcode")
  acts <- acts[, c("barcode", "celltype", "score", "p_value", "statistic")]
  
  return(acts)
}