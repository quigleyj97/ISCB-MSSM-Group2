// Expects the dataset to be in the working directory
load("BTBR.data.clean.Rdata")

sirt1_islet <- islet.rz[, annot[grep("Sirt1", annot$gene1), 1]]
foxo_islet <- islet.rz[, annot[grep("Foxo1", annot$gene1), 1]]
