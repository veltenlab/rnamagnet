## ---- echo = T, message=F,warning=F--------------------------------------
require(RNAMagnet)
require(ggplot2)

ligrec <- getLigandsReceptors("1.0.0",cellularCompartment = c("Secreted","Both"),manualAnnotation = "Correct")
head(ligrec)

## ---- echo = F-----------------------------------------------------------
NicheDataColors <-
c(Erythroblasts = "#bc7c7c", Chondrocytes = "#a6c7f7", Osteoblasts = "#0061ff", 
`Fibro/Chondro p.` = "#70a5f9", `pro-B` = "#7b9696", `Arteriolar ECs` = "#b5a800", 
`B cell` = "#000000", `large pre-B.` = "#495959", `Sinusoidal ECs` = "#ffee00", 
Fibroblasts = "#70a5f9", `Endosteal fibro.` = "#264570", `Arteriolar fibro.` = "#567fba", 
`Stromal fibro.` = "#465f82", `small pre-B.` = "#323d3d", `Adipo-CAR` = "#ffb556", 
`Ng2+ MSCs` = "#ab51ff", Neutrophils = "#1f7700", `T cells` = "#915400", 
`NK cells` = "#846232", `Schwann cells` = "#ff00fa", `Osteo-CAR` = "#ff0000", 
`Dendritic cells` = "#44593c", Myofibroblasts = "#dddddd", Monocytes = "#8fff68", 
`Smooth muscle` = "#ff2068", `Ery prog.` = "#f9a7a7", `Mk prog.` = "#f9e0a7", 
`Ery/Mk prog.` = "#f9cda7", `Gran/Mono prog.` = "#e0f9a7", `Neutro prog.` = "#c6f9a7", 
`Mono prog.` = "#f4f9a7", LMPPs = "#a7f9e9", `Eo/Baso prog.` = "#a7b7f9", 
HSPCs = "#c6f9a7")

## ---- echo=TRUE, message =F, warning=F, fig.width=8, fig.height=5--------

qplot(x = NicheData10x@dr$tsne@cell.embeddings[,1], y = NicheData10x@dr$tsne@cell.embeddings[,2], color =NicheData10x@ident) + scale_color_manual(values = NicheDataColors) + theme_bw() + theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank())



## ---- echo=TRUE, message =F, warning=F, fig.width=8, fig.height=6--------
  result <- RNAMagnetSignaling(NicheData10x, .version = "1.0.0")

## ---- echo=TRUE, message =F, warning=F, fig.width=6, fig.height=4.5------
qplot(x =NicheData10x@dr$tsne@cell.embeddings[,1], y=NicheData10x@dr$tsne@cell.embeddings[,2], color = result@specificity[,"Adipo-CAR"] ,size=I(0.75)) + scale_color_gradientn(name = "Strength of signal\nfrom Adipo-CAR",colours = c("black","blue","red")) + theme_bw() + theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank())


## ---- echo=TRUE, message =F, warning=F, fig.width=8, fig.height=6--------
PlotSignalingNetwork(result, threshold = 0.01)

## ---- echo=TRUE, message =F, warning=F, fig.width=8, fig.height=6--------
getRNAMagnetGenes(result, "Osteo-CAR", "pro-B")


