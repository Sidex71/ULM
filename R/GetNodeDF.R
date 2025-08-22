#' getting nodes and edges
#' @description
#'  'GetNodeDF()' a function to generate nodes and edges from multiplet data
#' @param mat a data frame of multiplets types and their frequencies, ideally the summary data frame obtained from the FilterMultiplet() or GetMultiplet() function.
#' @returns a data frame of cell-cell edges and nodes and the corresponding frequencies
#' @examples
#' mat <- data.frame('multilpetType' = c(paste('A', 'B', 'C', sep = '_'),
#'                                       paste('A', 'B', sep = '_'),
#'                                       paste('B', 'C', 'A', 'E', sep = '_'),
#'                                       paste('D', 'C', 'E', sep = '_')),
#'                   frequency = rep(50, 4))
#'
#' mat
#' GetNodeDF(mat)
#'
#' @importFrom stringr str_split
#' @importFrom utils combn
#' @importFrom dplyr mutate
#' @importFrom stats aggregate
#' @export
#'
#'

GetNodeDF <- function(mat){
  mat <- as.data.frame(mat)
  my_network <- apply(mat, 1, function (x){
    vec_split <- str_split(x[1], '_', simplify =F)
    vec_df <- as.data.frame(vec_split)
    colnames(vec_df)[1] <- 'cells'
    vec_comb <- t(combn(vec_df$cells, 2)) %>%
      as.data.frame() %>%
      mutate(n_cells = as.numeric(x[2])
      )
  })
  network_df <- do.call(rbind, my_network)
  colnames(network_df) <- c('Cell1', 'Cell2', 'n_cells')
  network_df <- aggregate(n_cells~Cell1+Cell2, data= network_df, FUN = sum)
  return(network_df)
}
