#' getting final cell assignments
#' @description
#'  'GetCellAssignments()' a function that assigns cell type labels to each barcode based on the signature scores
#' @param score_data a data frame of cell barcodes and gene signature scores, ideally obtained as the output from the GetCellScores() function.
#' @param p_val  a numeric value specifying the p value cut-off to filter significant signature scores (default: 0.05)
#' @param cut_off a numeric value specifying the cut-off for signature scores (default: 1)
#' @returns a data frame of cell barcodes and cell type assignments. A barcode may be assigned a single or multi cell type assignment depending on signature enrichment scores.
#' @examples
#' data(int_singData)
#' data(int_signature)
#' my_scores <- GetCellScores(seurat_obj = int_singData[,1:1000], signatures = int_signature, assay = 'RNA', slot = 'data')
#' my_ass <- GetCellAssignments(score_data = my_scores)
#' head(my_ass)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter mutate group_by select distinct
#' @importFrom stringr str_replace_all
#' @export
#'
#'
GetCellAssignments <- function(score_data, p_val = 0.05, cut_off = 1){
  acts_filt <- score_data %>% filter(p_value <= p_val & score > cut_off) %>%
    mutate(count_ulm = 1,
           celltype = str_replace_all(celltype, '_', ' '))
  cell_class <- acts_filt %>%
    group_by(barcode) %>%
    mutate(count_ulm= sum(count_ulm),
           celltype_ulm= paste(celltype, collapse = '_'),
           avg_pvalue = mean(p_value),
           avg_score = mean(score)) %>%
    dplyr::select(-c(celltype, p_value, score)) %>% distinct() %>% as.data.frame()

  return(cell_class)
}
