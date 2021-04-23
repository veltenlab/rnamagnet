#'Runs RNAMagnet for identifying localization of single cells to anchor populations
#'
#' RNAMagnet comes in two flavors: \code{RNAMagnetAnchors} and \code{\link{RNAMagnetSignaling}}. This function is meant to identify, for each single cell from a \code{\link[Seurat]{seurat}} object, the propensity to physically adhere to a set of anchor populations. For example, this function can identify if a single cell is more likely to bind to arteriolar or sinusoidal vessels.
#'@param seurat An object of class \code{\link[Seurat]{seurat}} containing a valid clustering and t-SNE information. For information on how to create such an object, see https://satijalab.org/seurat/get_started.html
#'@param anchors A character vector of anchor populations. Entries must be levels of \code{seurat@@ident}
#'@param return Determines object to return; one of "summary" or "rnamagnet-class"
#'@param neighborhood.distance See detail
#'@param neighborhood.gradient See detail
#'@param ... For explanation of all further parameters, see \code{\link{RNAMagnetBase}}.
#'@return Returns a data frame containing, for each cell, the propensity to physically interact with the various anchor populations as well as an overall adhesiveness score and prefered interaction partner. Alternatively, if \code{return} ist set to \code{rnamagnet-class}, an object of class \code{\link{rnamagnet}}.
#'@details Highly similar cell types can localize to different physical structures. For example, one type of pericyte may localize to sinusoids and another type may localize to arterioles. To increase the resolution in this scenario, \code{RNAMagnetAnchors} for each pair of single cells and anchor populations therefore computes a specificity score that describes how the single cell differs from \emph{similar} cells.
#'@details In our hands, the behavior of RNAMagnet is largely insensitive to the parameters \code{neighborhood.distance} and \code{neighborhood.gradient} that define what actually constitutes a similar cell. To explore how these parameters affect the weight each single cell carries in specificity score computation, see the example code.
#'@examples
#'\dontrun{
#'result <- RNAMagnetAnchors(NicheData10x, anchors = c("Sinusoidal ECs","Arteriolar ECs","Smooth muscle","Osteoblasts"))
#'qplot(x =NicheData10x@dr$tsne@cell.embeddings[,1], y=NicheData10x@dr$tsne@cell.embeddings[,2], \
#'      color = direction,size=I(0.75),alpha= adhesiveness,data=result) + \
#'      scale_color_brewer(name = "RNAMagnet\nLocation",palette= "Set1") + \
#'      scale_alpha_continuous(name = "RNAMagnet\nAdhesiveness")
#'
#'#To understand the effect of the neighborhood.distance and neighborhood.gradient parameters
#'#consider the following snippet:
#'myMagnet <- RNAMagnetAnchors(NicheData10x, return = "rnamagnet-class", \
#'            anchors = c("Sinusoidal ECs","Arteriolar ECs","Smooth muscle","Osteoblasts"))
#'use <- 1234 #select some cell of interest
#'kernel <- function(x,k=10, x0=0.5) 1/(1+exp(-k * (x-x0))) #defines weighing function
#'dbig <- as.matrix(1-cor(t(seurat0h@dr$pca@cell.embeddings[,1:15])))
#'plf <- data.frame(x = seurat@dr$tsne@cell.embeddings[,1],y = seurat@dr$tsne@cell.embeddings[,2], \
#'                  weight = 1-kernel(dbig[use,],k=neighborhood.gradient,x0=neighborhood.distance))
#'qplot(x = x, y= y, color = weight, data=plf) + \
#'     scale_color_gradientn(name = "Weight in local neighborhood", colours = c("#EEEEEE","#999999","blue","red")) + \
#'     geom_point(color = "black", shape = 17, size= 3, data=plf[use,])
#'}
#'@export
RNAMagnetAnchors <- function(seurat, anchors, return = "summary", neighborhood.distance = 0.7, neighborhood.gradient = 3, .k = 10, .x0 = 0.5, .minExpression = 0, .minMolecules = 1, .version = "1.0.0", .cellularCompartment = c("Membrane","ECM","Both"), .manualAnnotation = "Correct" ) {

  myMagnet <- RNAMagnetBase(seurat, anchors, neighborhood.distance,neighborhood.gradient, .k, .x0, .minExpression, .minMolecules, .version, .cellularCompartment, .manualAnnotation,TRUE)
  if (return=="rnamagnet-class") myMagnet else data.frame(direction = as.factor(colnames(myMagnet@specificity)[apply(myMagnet@specificity,1,which.max)]), adhesiveness = myMagnet@adhesiveness, myMagnet@specificity[,anchors])

}

#'Runs RNAMagnet for identifying signaling interactions between cells
#'
#' RNAMagnet comes in two flavors: \code{RNAMagnetSignaling} and \code{\link{RNAMagnetAnchors}}. This function is meant to identify, for each cell type from a \code{\link[Seurat]{seurat}} object, potential signaling interactions with other cell types.
#'@param seurat An object of class \code{\link[Seurat]{seurat}} containing a valid clustering and t-SNE information. For information on how to create such an object, see https://satijalab.org/seurat/get_started.html
#'@param ... For explanation of all further parameters, see \code{\link{RNAMagnetBase}}.
#'@return Returns an objects of class \code{\link{rnamagnet}}. \code{\link{PlotSignalingNetwork}} or \code{\link{getRNAMagnetGenes}} can be used for further analyses.
#'@export
RNAMagnetSignaling <- function(seurat, neighborhood.distance = NULL, neighborhood.gradient = NULL, .k = 10, .x0 = 0.5, .minExpression = 10, .minMolecules = 1, .version = "1.0.0", .cellularCompartment = c("Secreted","Both"), .manualAnnotation = "Correct" ) {

  RNAMagnetBase(seurat, anchors = NULL, neighborhood.distance,neighborhood.gradient, .k, .x0, .minExpression, .minMolecules, .version, .cellularCompartment, .manualAnnotation, FALSE)

}

#'Low-level function to run core RNA magnet steps
#'
#' Users are advised to use the top-level functions \code{\link{RNAMagnetAnchors}} and \code{\link{RNAMagnetSignaling}} which appropriately set default parameters and return user-friendly return values.
#' This is a low level function for development purposes.
#'@param seurat An object of class \code{\link[Seurat]{seurat}} containing a valid clustering and t-SNE information
#'
#'. For information on how to create such an object, see https://satijalab.org/seurat/get_started.html
#'@param anchors A character vector of anchor populations. Entries must be levels of the seurat identity vector. If \code{NULL}: All entries of the seurat identity vector are used as anchors.
#'@param neighborhood.distance See \code{\link{RNAMagnetAnchors}}
#'@param neighborhood.gradient See \code{\link{RNAMagnetAnchors}}
#'@param .k Fuzzification parameter, see detail. Recommended to leave at the default value.
#'@param .x0 Fuzzification parameter, see detail. Recommended to leave at the default value.
#'@param .minExpression Minimal expression level of genes to be included, specified as the number of cells in the dataset that express the gene.
#'@param .minMolecules Number of molecules per cell required to count the cells as positive.
#'@param .version The version of the underlying ligand-receptor database. See \code{\link{getLigandsReceptors}}.
#'@param .cellularCompartment Types of ligands to be included. For physical interactions, defaults to \code{ c("Membrane","ECM","Both")}.  See \code{\link{getLigandsReceptors}}.
#'@param .manualAnnotation Annotation status of ligands to be included. Default to \code{"Correct"}. See \code{\link{getLigandsReceptors}}.
#'@param .symmetric Assume that if A is a receptor for B, B is also a receptor for A
#'@details The algorithm takes the following steps: \enumerate{
#'\item Ligand-receptor pairs are selected based on the parameters \code{.version}, \code{.cellularCompartment} and \code{.manualAnnotation}. Choice of \code{.cellularCompartment} is crucial for determining the algorithm's behavior, e.g. if set to \code{c("Secreted","Both")}, paracrine signaling interactions involving soluble ligands are investigated.
#'\item Dropout values in the expression levels of ligands and receptors are imputed using \code{\link[Rmagic]{magic}}
#'\item Mean expression level of ligands and receptors is computed for all anchor populations
#'\item For each cell or anchor population, the expression of each ligand and receptor is encoded as a fuzzy logic variable
#'\item Fuzzy logic AND is used to compute probabilities for a given interaction to be active between a single cell and an anchor population
#'\item An interaction score is computed as the sum of interaction probabilities across all possible ligand-receptor pairs
#'\item Specificty scores are computed by comparing interaction scores to average interaction scores in a local neighborhood.
#'}
#'@details Add the methods section of the paper here!
#'@return Returns an object of class \code{\link{rnamagnet}}
#'@export
RNAMagnetBase <- function(seurat, anchors=NULL,neighborhood.distance=NULL, neighborhood.gradient =NULL, .k = 10, .x0 = 0.5, .minExpression,.minMolecules=1, .version = "1.0.0", .cellularCompartment, .manualAnnotation = "Correct", .symmetric = F) {
  cat("Setting everything up...\n")

  if (grepl("^3|^4", Biobase::package.version("Seurat"))) {
    seurat.ident <- Idents(seurat)
    seurat.cell.names <- colnames(seurat)
    seurat.pca <- Embeddings(seurat, reduction = "pca")
    seurat.raw.data <- GetAssayData(seurat, slot = "counts")

  } else {
    seurat.ident <- seurat@ident
    seurat.cell.names <- seurat@cell.names
    seurat.pca <- seurat@dr$pca@cell.embeddings
    seurat.raw.data <- seurat@raw.data
  }

  if (is.null(anchors)) anchors <- as.character(unique(seurat.ident))

  out <- new("rnamagnet", celltype = seurat.ident, params = list("neighborhood.distance"=neighborhood.distance, "neighborhood.gradient" =neighborhood.gradient, ".k" = .k, ".x0" = .x0, ".minExpression" = .minExpression, ".minMolecules" = .minMolecules, ".cellularCompartment" = .cellularCompartment, ".manualAnnotation" = .manualAnnotation, ".symmetric" = .symmetric))

  #compute cell-cell similarity
  similarity <- as.matrix(1-cor(t(seurat.pca[,1:15])))

  #prepare database
  ligrec <- getLigandsReceptors(.version, .cellularCompartment, .manualAnnotation)
  if (.symmetric) ligrec <- makeSymmetric(ligrec) #for physical interactions: i a binds b, b binds a.


  #prepare genes included into MAGIC
  filteredGenes <- rownames(seurat.raw.data)[apply(seurat.raw.data[,seurat.cell.names] >=.minMolecules ,1,sum) > .minExpression]

  genes <- unique(c(ligrec$Receptor.Mouse,
                    ligrec$Ligand.Mouse))

  genes <- genes[sapply(genes, function(x) {
    entries <- strsplit(x, "[&|]")
    if (grepl("&", x))
      all(entries[[1]] %in% filteredGenes)
    else any(entries[[1]] %in% filteredGenes)
  })]

  genes_formagic <- unlist(strsplit( genes,"[&|]"))


  formagic <- Matrix::t(seurat.raw.data[,seurat.cell.names])
  formagic <- Rmagic::library.size.normalize(formagic)
  formagic <- sqrt(formagic)

  #Run MAGIC
  cat("Now running MAGIC to impute dropout values...\n")
  mymagic <- Rmagic::magic(formagic, genes = genes_formagic, seed = 0xbeef)
  mymagic <- as.data.frame(mymagic)

  #handle expression levels of heterodimers
  resolvedRawData <- resolveData(t(mymagic), genes)
  resolvedRawData <- t(apply(resolvedRawData, 1, function(x) (x - min(x )) / (max(x) - min(x)) ))


  annotated_genes <- rownames(resolvedRawData)
  out@mylr <- subset(ligrec, Receptor.Mouse %in% annotated_genes &
                   Ligand.Mouse %in% annotated_genes)

  stepf<- Vectorize(function(x) if (x<0) 0 else x)

  cat("Now running RNAMagnet...\n")

  #compute mean gene expression level per population
  out@anchors <- do.call(cbind, lapply(anchors, function(id) {
    apply(resolvedRawData[,seurat.ident == id],1,mean)
  }));
  colnames(out@anchors) <- anchors

  #use fuzzy logic and operations to compute interaction score
  out@interaction <- sapply(anchors, function(pop_l) {
    out@mylr$expression_ligand <- out@anchors[out@mylr$Ligand.Mouse, pop_l]
    sapply(seurat.cell.names, function(cell_r) {
      out@mylr$expression_receptor <-resolvedRawData[out@mylr$Receptor.Mouse, cell_r]
      sum(kernel(out@mylr$expression_ligand, k =.k, x0 = .x0) * kernel(out@mylr$expression_receptor, k = .k, x0=.x0 )) #kernel performs fuzzification
    })
  })

  out@specificity <- t(sapply(rownames(out@interaction), function(cell) {
    x <- out@interaction[cell,]
    beta <- x/sum(x)
    if (!is.null(neighborhood.distance)) alpha <- apply(out@interaction * (1-kernel(similarity[cell,],neighborhood.gradient,x0=neighborhood.distance )),2,sum) else alpha <- apply(out@interaction,2,mean)
    alpha <- alpha / sum(alpha)
    beta- alpha
  }))
  rownames(out@specificity) <- rownames(out@interaction)

  out@adhesiveness <- apply(mymagic, 1,function(x) sum(kernel(x,x0 = .x0, k = .k)))
  return(out)

}


makeSymmetric <- function(ligrec) {
  toadd <- list()
  for (i in 1:nrow(ligrec)) {
    if (! any(ligrec$Ligand.Mouse[i] == ligrec$Receptor.Mouse & ligrec$Receptor.Mouse[i] == ligrec$Ligand.Mouse )) {
      toadd[[length(toadd)+1]] <- ligrec[i,]
      toadd[[length(toadd)]]$Receptor.Mouse <- ligrec[i,"Ligand.Mouse"]
      toadd[[length(toadd)]]$Ligand.Mouse <- ligrec[i,"Receptor.Mouse"]
    }
  }
  rbind(ligrec, do.call(rbind, toadd))
}

#function to handle heterodimeric complexes by applying fuzzy logic AND operations
resolveData <- function(rawdata, lr) {
  #as ligands/receptor can contain gene lists (complex components), apply fuzzy logic to gene expression values
  singleentries <- lr[!grepl("[&|]",lr)]
  use_normdata <- rawdata[singleentries,]

  doubleentries <- lr[grepl("[&|]",lr)]

  add_normdata <- do.call(rbind,lapply(doubleentries, function(x) {

    if (grepl("&",x)) {
      entries <- strsplit(x, "&",fixed = T)[[1]]
      apply(rawdata[entries,],2,min)
    }  else {
      entries <- strsplit(x, "|",fixed=T)
      entries <- entries[[1]][entries[[1]] %in% rownames(rawdata)]
      if(length(entries) > 1) apply(rawdata[entries,],2,max) else rawdata[entries,]
    }
  }))
  rownames(add_normdata) <- doubleentries

  use_normdata <- rbind(use_normdata, add_normdata)
  use_normdata[lr,]
}

kernel <- function(x,k=10, x0=0.5) 1/(1+exp(-k * (x-x0)))
