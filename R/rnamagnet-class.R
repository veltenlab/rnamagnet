
#'RNAMagnet class
#'
#'This class contains representation for the output of \code{\link{RNAMagnetBase}}.
#'
#'@slot anchors Mean gene expression matrix of anchor populations
#'@slot interaction Interaction score matrix between cells and anchor populations
#'@slot specificity Specificity of interaction scores relative to similar cells
#'@slot adhesiveness Total interaction potential per cell
#'@slot celltype Per-cell cell type annotation, copied over from \code{seurat@ident}
#'@slot mylr The output of \code{\link{getLigandsReceptors}} that this instance of RNAMagnet is based on
#'@slot params Named list of paramers used in \code{\link{RNAMagnetBase}}
#'@export
rnamagnet <- setClass(
  "rnamagnet",
  slots = c(
    anchors = "matrix",
    weight = "matrix",
    interaction = "matrix",
    specificity = "matrix",
    adhesiveness = "numeric",
    celltype = "factor",
    mylr = "data.frame",
    params = "list"

  )

)
