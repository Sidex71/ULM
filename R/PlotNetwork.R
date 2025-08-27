#' plotting network
#' @description
#'  'PlotNetwork()' a function to plot a customizable cell-cell interaction network graph
#'
#'
#' @param network_df a data frame containing the network data, ideally the output from the GetNodeDF() function containing the edges, nodes and weights.
#' @param node_size a numeric value specifying the size of the nodes in the plot (default: 20).
#' @param node_color a character defining the color of the nodes (default: "blue").
#' @param node_text_size a numeric value specifying the size of the text labels for nodes (default: 4).
#' @param edge_width_factor a numeric value specifying the scaling factor applied to edge weights to determine edge width (default: 50).
#' @param edge_color a character specifying the color of the edges (default: "red").
#' @param network_layout a character defining the layout algorithm to use for positioning nodes (default: "fr"). Common options include `"fr"` (Fruchterman–Reingold), `"kk"` (Kamada–Kawai), or Large Graph Layout `"lgl"`
#' @param legend_title a character specifying the title of the legend (default: "Scaled Counts").
#' @param legend_position a character specifying the position of the plot legend (default: "bottom"). Options include `"bottom"`, `"top"`, `"left"`, `"right"` .
#' @param min_edge_width a numeric value defining the minimum width of the edges (default: 0.5).
#' @param max_edge_width a numeric value defining the maximum width of the edges (default: 3).
#' @param main a character specifying the main title of the plot (default: "Network Plot").
#' @param main_size a numeric value specifying the text size of the main plot title (default: 15).
#' @param hjust a numeric value defining the horizontal justification of the main title (default: 0.5, centered).
#' @param legend_text_size a numeric value specifying the text size of the legend (default: 12).
#' @param legend_title_size a numeric value specifying the size of the legend title (default: 14).
#'
#' @return a `ggplot`  `ggraph` object representing the cell-cell interaction network visualization.
#' @examples
#' my_network_df <- data.frame(Cell1 =c('A', 'B', 'C', 'D', 'E', 'F'),
#'                             Cell2= c('D', 'A', 'F', 'C', 'F', 'B'),
#'                             n_cells = c(20, 40, 60, 80, 100, 120))
#'
#' PlotNetwork(network_df = my_network_df)
#'
#' @importFrom igraph graph_from_data_frame
#' @importFrom tidygraph as_tbl_graph activate
#' @importFrom dplyr mutate
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text scale_edge_width_continuous
#' @importFrom ggplot2 theme_void ggtitle theme element_text guides aes
#' @export
#'
PlotNetwork <- function(network_df, node_size = 20, node_color = 'blue',
                        node_text_size = 4, edge_width_factor = 50,
                        edge_color = 'red', network_layout = 'fr',
                        legend_title = 'Scaled Counts', legend_position = 'bottom',
                        min_edge_width = 0.5, max_edge_width = 3,
                        main = 'Network Plot', main_size = 15,
                        hjust = 0.5, legend_text_size = 12, legend_title_size = 14) {
  g <- graph_from_data_frame(network_df, directed = FALSE)

  tg <- as_tbl_graph(g)

  # Set node and edge aesthetics
  tg <- tg %>%
    activate(nodes) %>%
    mutate(node_size = node_size) %>%
    activate(edges) %>%
    mutate(edge_width = n_cells / edge_width_factor,
           edge_color = edge_color)

  # Create the network plot using ggraph
  set.seed(837100)
  net_plot <- ggraph(tg, layout = network_layout) +
    geom_edge_link(aes(width = edge_width), color = edge_color, show.legend = TRUE) +
    geom_node_point(aes(size = node_size), color = node_color, show.legend = FALSE) +
    geom_node_text(aes(label = name), vjust = 1.8, size = node_text_size, show.legend = FALSE) +
    scale_edge_width_continuous(range = c(min_edge_width, max_edge_width), name = legend_title) +
    theme_void() +
    ggtitle(main) +
    theme(
      legend.position = legend_position,
      legend.text = element_text(size = legend_text_size),
      legend.title = element_text(size = legend_title_size),
      plot.title = element_text(size = main_size, hjust = hjust)
    ) +
    guides(size = "none")
  print(net_plot)
}
