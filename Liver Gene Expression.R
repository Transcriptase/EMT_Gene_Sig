library("lumi")
library("lumiMouseIDMapping")
library("lumiMouseAll.db")
library("annotate")
library("mogene10stprobeset.db")

PATH = paste(getwd(), "/Raw Liver Data", sep = "")
file_list <- list.files(path = PATH)
file_list <- paste(PATH, file_list, sep = "/")

x.lumi <- lumiR.batch(file_list, lib.mapping = "lumiMouseIDMapping")
lumi.T <- lumiT(x.lumi, method = "log2")
dataMatrix <- exprs(lumi.T)
presentCount <- detectionCall(lumi.T)
selDataMatrix <- dataMatrix[presentCount > 0,]
probeList <- rownames(selDataMatrix)


oncogene <- c("M",
              "Normal",
              "Normal",
              "M",
              "MT",
              "MT",
              "MT",
              "MT",
              "MT",
              "MT",
              "MT",
              "MT",
              "MT",
              "MT",
              "MT",
              "MT",
              "MT",
              "MT",
              "M",
              "M")
site <- c("Pri",
          "Norm",
          "Norm",
          "Pri",
          "Pri",
          "Met",
          "Pri",
          "Met",
          "Pri",
          "Met",
          "Pri",
          "Met",
          "Pri",
          "Met",
          "Met",
          "Pri",
          "Met",
          "Met",
          "Pri",
          "Pri")
annot <- cbind(seq(from = 1, to = 20, by = 1),
               colnames(selDataMatrix),
               oncogene,
               site,
               c(rep("A", 10), rep("B", 10))
               )
colnames(annot) <- c("sample", "Id", "oncogene", "site", "batch")

batch.trts <- factor(annot[,"batch"])
site.trts