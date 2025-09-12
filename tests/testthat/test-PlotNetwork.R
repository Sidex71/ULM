test_that("PlotNetwork produces a ggraph plot", {
  ######### Generate a dummy multuplet data frame
  my_network_df <- data.frame(
    Cell1 = c("A", "B", "C", "D", "E", "F"),
    Cell2 = c("D", "A", "F", "C", "F", "B"),
    n_cells = c(20, 40, 60, 80, 100, 120)
  )
  
  ## Run function (suppress printed plot for testing purpose)
  net_plot <- suppressMessages(
    PlotNetwork(network_df = my_network_df)
  )
  
  # Basic checks
  expect_s3_class(net_plot, "ggplot")   ## returns a ggplot object
  # Title check
  expect_equal(net_plot$labels$title, "Network Plot")
})
