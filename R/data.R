#' 10x dataset from mouse bone marrow
#'
#' @format An object of class \code{\link[Seurat]{seurat}}, containing clustering and dimensionality reduction
"NicheData10x"

#' LCM dataset from mouse bone marrow
#'
#' @format Read count matrix, see \code{\link{NicheMetaDataLCM}} for per-sample metadata.
"NicheDataLCM"

#' Meta data for LCM dataset
#'
#' @format Data frame, column \code{id} corresponds to the column names of \code{\link{NicheDataLCM}}. Other colummns describe various parameters related to potential sources of batch effects (e.g. the microscopy slide, the day of sample collection and processing, the sequencing lane and the size of the area sampled) and biology (Basic sample class, presence of sinusoids and distance from the endosteum)
#' @describeIn NicheDataLCM
"NicheMetaDataLCM"

#' Markers for populations defined by 10x genomics
#'
#' @format Data frame of cell type specific markers - the output of running \code{FindMarkersAll{\link{NicheData10x}, method="roc"}}
#' @describeIn NicheData10x
"NicheMarkers10x"
