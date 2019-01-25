#'Visualizes RNAMagnet output as a signaling network
#'
#'@param rnamagnet An object of class \code{\link{rnamagnet}}, as produced by \code{\link{RNAMagnetSignaling}}
#'@param return Determines object to return; one of "summary" or "rnamagnet-class"
#'@param threshold Threshold interaction strength for including an edge in the network.
#'@param colors A named vector of colors to be used for plotting; names correspond to levels of \code{rnamagnet@class}
#'@param useLabels Plot the population labels on the graph?
#'@return Plots the network and invisibly returns an object of class \code{\link[igraph]{igraph}}
#'@export
PlotSignalingNetwork <- function(rnamagnet, threshold = 0.01, colors = NULL, useLabels=T) {

  pops <- colnames(rnamagnet@specificity)
  mean_by_pop <- sapply(pops, function(id) {
    sapply(pops, function(id2) {
      mean(rnamagnet@specificity[rnamagnet@celltype == id2,id])
    })
  })

  graphf <- reshape2::melt(t(mean_by_pop))

  colnames(graphf) <- c("sender", "receiver", "value")
  graphf$sender <- as.character(graphf$sender); graphf$receiver <- as.character(graphf$receiver);
  graphf$width <- graphf$value^1.5 *1000

  vertices <- data.frame(id = unique(graphf$sender), size = 5, stringsAsFactors = F, label = unique(graphf$sender))

if(!is.null(colors)){
  vertices$label.color <- colors[vertices$id]
  vertices$color <- colors[vertices$id]
  graphf$color <- colors[graphf$sender]
}

  plotgraph <- igraph::graph_from_data_frame(subset(graphf,value > 0.01 ), vertices = vertices)
  if (useLabels) plot(plotgraph,vertex.frame.color=NA, edge.arrow.mode = "-",  edge.curved=0.3) else plot(plotgraph,vertex.frame.color=NA, edge.arrow.mode = "-", vertex.label = NA, edge.curved=0.3)
  invisible(plotgraph)
}
