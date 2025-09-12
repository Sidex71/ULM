test_that("GetNodeDF generates correct edges from multiplet data", {
  ######### Generate a dummy multiplet data frame
  mat <- data.frame(
    multilpetType = c(
      paste("A", "B", "C", sep = "_"),
      paste("A", "B", sep = "_"),
      paste("B", "C", "A", "E", sep = "_"),
      paste("D", "C", "E", sep = "_")
    ),
    frequency = rep(50, 4),
    stringsAsFactors = FALSE
  )
  
  ###########Run function to generate edges
  network_df <- GetNodeDF(mat)
  
  ### Basic checks on expected output
  expect_s3_class(network_df, "data.frame")                        ## output must be a data frame
  expect_true(all(c("Cell1", "Cell2", "n_cells") %in% colnames(network_df))) ## required columns
  expect_true(all(network_df$n_cells >= 50))                       ## all counts must be at least input frequency
  expect_gt(nrow(network_df), 0)                                   ## must return some edges
  expect_true(all(network_df$Cell1 != network_df$Cell2))           ## no self-loops
  
  ## Check that known edge exists
  expect_true(any(network_df$Cell1 == "A" & network_df$Cell2 == "B"))
})
