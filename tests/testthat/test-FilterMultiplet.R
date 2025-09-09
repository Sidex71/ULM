

test_that("FilterMultiplet filters multiplets correctly", {
  # load example data included in package
  data("int_multData", package = "ULM")
  data("int_signature", package = "ULM")
  library(ULM)
  # compute scores
  my_scores <- GetCellScores(
    seurat_obj = int_multData[, 1:200],   # keep smaller subset for speed
    signatures = int_signature,
    assay = "RNA",
    layer = "data"
  )
  
  # assign cells
  my_ass <- GetCellAssignments(score_data = my_scores)
  
  # add metadata
  new_obj <- AddMetaObject(
    seurat_obj = int_multData[, 1:200],
    cell_class_df = my_ass
  )
  

  # filter multiplets — expect a warning about removed cells
  expect_warning(
    my_mult_filt <- FilterMultiplet(seurat_obj = new_obj, minCells = 2, minFreq = 2),
    regexp = "Removing .* cells missing data"
  )
  # expectations
  expect_type(my_mult_filt, "list")
  expect_named(my_mult_filt, c("multSummaryFilt", "multObjFilt"))
  
  expect_true(all(c("multipletType", "frequency") %in% colnames(my_mult_filt$multSummaryFilt)))
  expect_s4_class(my_mult_filt$multObjFilt, "Seurat")
  
  # sanity check: filtered frequencies are >= minFreq
  expect_true(all(my_mult_filt$multSummaryFilt$frequency >= 2))
})
