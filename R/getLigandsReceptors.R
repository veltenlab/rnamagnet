#'Retrieve a database of ligand-receptor pairs
#'
#'@param version Currently supports the following values: \itemize{
#'   \item{stable} points to \code{1.0.0}
#'   \item{latest} points to \code{2.2.0}
#'   \item{1.0.0} contains manual annotation for all genes expressed in bone marrow. This version was used for analysis is the Baccin et al paper. Information on heterodimeric receptors is only included for integrins.
#'   \item{2.0.0} contains manual annotation for all genes in the geneome, January 2019
#'   \item{2.1.0} integrates heterodimer information for non-integrin receptors from cellphoneDB, Janary 2019
#'   \item{2.2.0} updated to remove missing entries and several additions, April 2019
#'   \item{3.0.0} adds a translation of the complete cellphoneDB to mouse homologues.
#'   \item Alternatively, a data frame with the same column names as \code{ligandsReceptors_1.0.0} can be used.
#'}
#'@param cellularCompartment A character vector to select permitted localizations of the \strong{ligand}. Valid entries are \code{Membrane}, \code{ECM}, \code{Secreted} and \code{Both}, which refers to membrane-bound ligands that can also be secreted.
#'@param manualAnnotation A character vector to select permitted annotation status of the interaction. Valid entries are \code{Correct}, \code{Incorrect}, \code{Not Expressed} (i.e. entries not annotated in v1.0.0) and \code{Irrelevant} (i.e. entries we consider correct but not relevant at homeostasis, e.g. interactions involving activated components of the complement system)
#'@param ligandClass A character vector to select permitted classes of ligands. Defaults to all, i.e. \code{c("Other","Cytokine","Chemokine","GrowthFactor","Interleukin")}
#'@return Returns a data frame with 11 columns: \itemize{
#'\item{Pair.Name}{Unique identifier for the ligand-receptor pair}
#'\item{Ligand.Mouse}{MGI symbol of the ligand. }
#'\item{Receptor.Mouse}{MGI symbol of the receptor. In case multiple receptors are required for complex formation, multiple MGI symbols separated by \strong{&}, as in Itgal&Itgb2. }
#'\item{Source}{Ramilowski for entries taken from \emph{Ramilowski et al., Nature Communications 2015}; Baccin for entries added as part of this work}
#'\item{Reference}{Pubmed ID or similar}
#'\item{ManualAnnotation}{The literature evidence for all entries from Ramilowski et al with expression in bone where checked manually.}
#'\item{Ligand.CC}{Cellular Compartment of the ligand - one of Secreted, Membrane, Both, or ECM. Annotation is based on Gene Ontology, but was checked manually for entries with expression in bone}
#'\item{Ligand.GO}{Classification of the ligand by gene ontology into one of the following classes: GrowthFactor, Cytokine, Chemokine, Interleukin, or Other}
#'}
#'@export
getLigandsReceptors <- function(version = "latest", cellularCompartment = c("Membrane","ECM","Both","Secreted"), manualAnnotation = "Correct", ligandClass = c("Other","Cytokine","Chemokine","GrowthFactor","Interleukin")) {
  if (is.character(version)) out <- switch(version,
                latest = ligandsReceptors_2.2.0,
                stable = ligandsReceptors_1.0.0,
                "1.0.0" = ligandsReceptors_1.0.0,
                "2.0.0" = ligandsReceptors_2.0.0,
                "2.1.0" = ligandsReceptors_2.1.0,
                "2.2.0" = ligandsReceptors_2.2.0,
                "3.0.0" = ligandsReceptors_3.0.0) else if (is.data.frame(version)) out <- version else stop("Version needs to be a character string or a data frame.")
  if(is.null(out)) error("Version", version, "not supported. See documentation.")

  out <- subset(out, Ligand.CC %in% cellularCompartment & ManualAnnotation %in% manualAnnotation & Ligand.GO %in% ligandClass)
  out
}


#'Extracts receptor-ligand pairs mediating interaction between two cell types
#'
#'@param rnamagnet An object of class \code{\link{rnamagnet}}
#'@param pop_l Sending population to consider
#'@param pop_r Receiving population to consider
#'@param thresh Threshold score for ligand and receptor expression
#'@return A data frame conaining the score and the receptor-ligand pair
#'@export
getRNAMagnetGenes <- function(rnamagnet, pop_l, pop_r, thresh = 0.05) {
  rnamagnet@mylr$expression_ligand <- rnamagnet@anchors[rnamagnet@mylr$Ligand.Mouse, pop_l]
  rnamagnet@mylr$expression_receptor <- rnamagnet@anchors[rnamagnet@mylr$Receptor.Mouse, pop_r]
    out <- data.frame(score = kernel(rnamagnet@mylr$expression_ligand, k =rnamagnet@params[[".k"]], x0 = rnamagnet@params[[".x0"]]) * kernel(rnamagnet@mylr$expression_receptor, k =rnamagnet@params[[".k"]], x0 = rnamagnet@params[[".x0"]] ),
                      pair = rnamagnet@mylr$Pair.Name)
  subset(out[order(out$score, decreasing = T),], score > thresh)
}
