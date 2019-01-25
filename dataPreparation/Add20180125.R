require(RNAMagnet)

lr <- ligandsReceptors_1.0.0

julia <- read.csv("Julia20180125.csv", stringsAsFactors = F)

#check that entries with no entry in Julia TRUE/FALSE also don't have a value in the other fields
noentry <- julia$Julia.TRUE.FALSE.NA == ""
noentry[is.na(noentry)] <- F

wrongFormat <- julia$Pair.Name[which(noentry & apply(julia[,12:ncol(julia)],1,function(x) any(x != "")))]

naentries <- julia$Pair.Name[(is.na(julia$Julia.TRUE.FALSE.NA) | grepl("NA",julia$Julia.TRUE.FALSE.NA)) & julia$Annotation == "Not Expressed"]


#abgleichen: entries mit annotation von Chiara und trotzdem eintrag von Julia
#Ist Location die gleiche zwischen Julia und der original DB?
#Julia should double check: Are entries with Location Ligand %in% c(ECM,Membrane) really mediating interaction?

julia <- julia[,-(12:ncol(julia))]
julia$Julia.TRUE.FALSE.NA[grepl("NA",julia$Julia.TRUE.FALSE.NA) | is.na(julia$Julia.TRUE.FALSE.NA)] <- ""
julia$Julia.TRUE.FALSE.NA <- sapply(julia$Julia.TRUE.FALSE.NA, function(x) switch(as.character(x),
                                                                                  "FALSE" = "Incorrect",
                                                                                  "Irrelevant" = "Irrelevant",
                                                                                  "TRUE" = "Correct",
                                                                                  "WAHR" = "Correct",
                                                                                  ""))
colnames(julia)[1:8] <- colnames(lr)

lr <- merge(lr, julia, by = c("Pair.Name","Ligand.Mouse","Receptor.Mouse","Source"),all.x = T, suffixes = c("",".julia") )

#changed location
lr$Ligand.CC.julia[!lr$Ligand.CC.julia %in% unique(lr$Ligand.CC)] <- lr$Ligand.CC[!lr$Ligand.CC.julia %in% unique(lr$Ligand.CC)]
chloc <- which(lr$Ligand.CC != lr$Ligand.CC.julia)
lr$Ligand.CC[chloc] <- lr$Ligand.CC.julia[chloc]

#entries that she re-annotated
reanno <- which(lr$ManualAnnotation != "Not Expressed" & lr$Julia.TRUE.FALSE.NA != "" & lr$Julia.TRUE.FALSE.NA != lr$ManualAnnotation)
lr$ManualAnnotation[reanno] <- lr$Julia.TRUE.FALSE.NA[reanno]

#entries that were not previously annotated

notanno <- which(lr$ManualAnnotation == "Not Expressed" & lr$Julia.TRUE.FALSE.NA != "")
lr$ManualAnnotation[notanno] <- lr$Julia.TRUE.FALSE.NA[notanno]
lr$Source[notanno] <- "Julia Schnell"
lr$Reference[notanno] <- paste(lr$Reference[notanno], lr$Julia.literature.reference[notanno], sep ="; ")

ligandsReceptors_2.0.0 <- lr[,colnames(ligandsReceptors_1.0.0)]

save(ligandsReceptors_2.0.0, file = "/g/steinmetz/velten/Scripts/rnamagnet/data/ligandsReceptors_2.0.0.rda")


