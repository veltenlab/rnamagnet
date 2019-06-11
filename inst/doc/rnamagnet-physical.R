## ---- echo = T, message=F,warning=F--------------------------------------
require(RNAMagnet)
require(ggplot2)

ligrec <- getLigandsReceptors("1.0.0",cellularCompartment = c("ECM","Surface","Both"),manualAnnotation = "Correct")
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

## ---- echo=TRUE, message =F, warning=F, fig.width=8, fig.height=6--------

qplot(x = NicheData10x@dr$tsne@cell.embeddings[,1], y = NicheData10x@dr$tsne@cell.embeddings[,2], color =NicheData10x@ident) + scale_color_manual(values = NicheDataColors) + theme_bw() + theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank())



## ---- echo=TRUE, message =F, warning=F, fig.width=8, fig.height=6--------
  result <- RNAMagnetAnchors(NicheData10x, anchors = c("Sinusoidal ECs","Arteriolar ECs","Smooth muscle","Osteoblasts"), .version = "1.0.0")

## ---- echo=TRUE, message =F, warning=F, fig.width=8, fig.height=6--------
head(result)

## ---- echo=TRUE, message =F, warning=F, fig.width=6, fig.height=4.5------
qplot(x =NicheData10x@dr$tsne@cell.embeddings[,1], y=NicheData10x@dr$tsne@cell.embeddings[,2], color = direction,size=I(0.75),alpha= adhesiveness,data=result) + scale_color_brewer(name = "RNAMagnet\nLocation",palette= "Set1") + scale_alpha_continuous(name = "RNAMagnet\nAdhesiveness") + theme_bw() + theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank())


## ---- echo=FALSE, message =F, warning=F, fig.width=4, fig.height=4.5-----
require(plyr)
require(reshape2)
require(pheatmap)

result$id <- NicheData10x@ident

summarised <- ddply(result, c("id"), summarise,  RNAmagnet.score = table(direction) / length(direction), n = rep(length(direction),4), RNAmagnet.adhesiveness = rep(mean(adhesiveness),4), experiment = names(table(direction)))

castMagnet <- dcast(subset(summarised, RNAmagnet.adhesiveness > 35), id ~ experiment, value.var = "RNAmagnet.score")
rownames(castMagnet) <- castMagnet[,1]
castMagnet <- castMagnet[,-1]
castMagnet <- t(apply(castMagnet,1,function(x) (x-min(x)) / (max(x)-min(x))))

pheatmap(castMagnet, cluster_cols = F, annotation_legend = F, annotation_names_col = F, color = colorRampPalette(c("white","white","blue","red"))(100), fontsize=8, treeheight_row = 20)

