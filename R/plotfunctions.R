#'Visualizes RNAMagnet output as a signaling network
#'
#'@param rnamagnet An object of class \code{\link{rnamagnet}}, as produced by \code{\link{RNAMagnetSignaling}}
#'@param score One of \code{specificity} or \code{volume}, see detail below.
#'@param threshold Threshold interaction strength for including an edge in the network.
#'@param colors A named vector of colors to be used for plotting; names correspond to levels of \code{rnamagnet@class}
#'@param useLabels Plot the population labels on the graph?
#'@param ... Parameters passed to \code{\link[igraph]{plot.igraph}}, e.g. a layout can be stored as described in the example code below
#'@details If \code{score} is set to \code{specificity}, two cell types are connected if they are enriched for interactions.
#'@details If \code{score} is set to \code{volume}, two cell types are connected if there are many interactions between them.
#'@details To understand the difference, consider the following example of three cell types A, B, and C, where A expresses 50 different cytokines to which both B and C express cognate receptors. B expresses 10 cytokines, for which only C expresses cognate receptors. If plotting the number of interactions, A would be connected to B and C; if plotting their specificity, B would be connected to C.
#'@examples
#'\dontrun{
#'network <- PlotSignalingNetwork(rnamagnet_result)
#'layout <- layout_nicely(network)
#'PlotSignalingNetwork(rnamagnet_result, layout = layout)
#'}
#'@return Plots the network and invisibly returns an object of class \code{\link[igraph]{igraph}}
#'@export
PlotSignalingNetwork <- function(rnamagnet, score = "specificity", threshold = ifelse(score == "specificity", 0.01,10), colors = NULL, useLabels=T, ...) {

  pops <- colnames(rnamagnet@specificity)

  if (tolower(score) == "specificity"){
  mean_by_pop <- sapply(pops, function(id) {
    sapply(pops, function(id2) {
      mean(rnamagnet@specificity[rnamagnet@celltype == id2,id])
    })
  })
  factor <- 1000
  } else if (tolower(score) == "volume"){
    mean_by_pop <- sapply(pops, function(id) {
      sapply(pops, function(id2) {
        mean(rnamagnet@interaction[rnamagnet@celltype == id2,id])
      })
    })
    factor <- 0.1
  } else stop("Score must be set to one of specificity or volume")

  graphf <- reshape2::melt(t(mean_by_pop))

  colnames(graphf) <- c("sender", "receiver", "value")
  graphf$sender <- as.character(graphf$sender); graphf$receiver <- as.character(graphf$receiver);
  graphf$width <- graphf$value^1.5 *factor

  vertices <- data.frame(id = unique(graphf$sender), size = 5, stringsAsFactors = F, label = unique(graphf$sender))

if(!is.null(colors)){
  vertices$label.color <- colors[vertices$id]
  vertices$color <- colors[vertices$id]
  graphf$color <- colors[graphf$sender]
}

  plotgraph <- igraph::graph_from_data_frame(subset(graphf,value > threshold ), vertices = vertices)
  if (useLabels) plot(plotgraph,vertex.frame.color=NA, edge.arrow.mode = "-",  edge.curved=0.3, ...) else plot(plotgraph,vertex.frame.color=NA, edge.arrow.mode = "-", vertex.label = NA, edge.curved=0.3, ...)
  invisible(plotgraph)
}
