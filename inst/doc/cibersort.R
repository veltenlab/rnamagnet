## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

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
HSPC = "#c6f9a7")

## ---- echo =T, warning=F, message=F, fig.width=6,fig.height=4------------
require(RNAMagnet)
require(Seurat)
require(ggplot2)

n.pca <- prcomp(t(DESeq2::varianceStabilizingTransformation(as.matrix(NicheDataLCM)) ))
qplot(x = n.pca$x[,1], y = n.pca$x[,2], color = NicheMetaDataLCM$biological.class) + scale_color_discrete(name="Sample type") + xlab("PC1") + ylab("PC2")
outliers <- n.pca$x[,1] > 70 | n.pca$x[,2] < -40
remove <- colnames(NicheDataLCM)[outliers]


## ---- echo =T, warning=F, message=F, fig.width=4.5,fig.height=3.8--------
for (pop in c("Ery/Mk prog.","Neutro prog.","Mono prog.","Gran/Mono prog.","LMPPs","Mk prog.","Eo/Baso prog.","Ery prog.")) NicheData10x <- RenameIdent(NicheData10x, pop, "HSPC")

usegenes <- unique(NicheMarkers10x$gene[(NicheMarkers10x$myAUC > 0.8 |NicheMarkers10x$myAUC < 0.2) ])

mean_by_cluster <- do.call(cbind, lapply(unique(NicheData10x@ident), function(x) {
  apply(NicheData10x@raw.data[usegenes,NicheData10x@cell.names][,NicheData10x@ident == x], 1,mean )
}))
colnames(mean_by_cluster) <- unique(NicheData10x@ident)

## ---- echo =T, warning=F, message=F, fig.width=4.5,fig.height=3.8--------
#character vector that maps column names of NicheDataLCM to sample type
LCM_design <- NicheMetaDataLCM$biological.class
names(LCM_design) <- NicheMetaDataLCM$id

CIBER <- runCIBERSORT(NicheDataLCM, mean_by_cluster, LCM_design, mc.cores=3)
head(CIBER)

## ---- echo =T, warning=F, message=F, fig.width=8,fig.height=6------------

CIBER <- subset(CIBER,CellType %in% c("Adipo-CAR","Ng2+ MSCs","Osteoblasts", "Arteriolar fibro.", "Sinusoidal ECs", "Osteo-CAR","Chondrocytes","Endosteal fibro.", "Fibro/Chondro p.", "Stromal fibro.", "Arteriolar ECs", "Smooth muscle") &  !SampleID %in% remove)

labeler <- c("ARTERIES" = "Arteriolar", "ENDOSTEUM" = "Endosteal", "HIGH SINUSOIDS" = "Sinusoidal", "LOW SINUSOIDS" = "Non-vascular", "SUB-ENDOSTEUM" = "Sub-Endosteal", "Other" = "darkgrey")

ggplot(aes(x = SampleClass, y= Fraction,color = CellType),data=CIBER) + geom_point(stat="summary", fun.y=mean) + facet_wrap(~ CellType, scales="free_y") + 
  theme_bw(base_size=12) + theme(axis.text.x = element_text(angle=90, color="black"),axis.text.y = element_blank(), panel.grid = element_blank()) + 
  geom_errorbar(stat="summary", fun.ymin=function(x) mean(x)+sd(x)/sqrt(length(x)),fun.ymax=function(x) mean(x)-sd(x)/sqrt(length(x)), width=0.2) + 
  ylab("CIBERSORT estimate (a.u.)") + scale_color_manual(values = NicheDataColors, guide=F) + xlab("Niche") + scale_x_discrete(labels = labeler)



