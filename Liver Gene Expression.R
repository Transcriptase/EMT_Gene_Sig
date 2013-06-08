library("lumi")
library("lumiMouseIDMapping")
library("lumiMouseAll.db")
library("annotate")
library("mogene10stprobeset.db")
library("limma")

WD = "E:/Russell/Documents/GitHub/EMT_Gene_Sig"
setwd(WD)

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
annot <- as.data.frame(annot)
annot$exp <- paste(annot$site, annot$oncogene, sep = "_")

batch.trts <- factor(annot$batch)
exp.trts <- factor(annot$exp)

batch_design <- model.matrix(~0 + exp.trts + batch.trts)
batch_fit <- lmFit(selDataMatrix, batch_design)
batch_corrected <- selDataMatrix -
    batch_fit$coefficients[,'batch.trtsB']%*%t(batch_design[,'batch.trtsB'])

MT_M.contrasts <- makeContrasts(exp.trtsPri_MT-exp.trtsPri_M,
                                levels=batch_design)

fit2 <- eBayes(contrasts.fit(batch_fit, MT_M.contrasts))

fit2$genes$Symbol <- getSYMBOL(fit2$genes$ID,'lumiMouseAll.db')
#fit2$genes$PMID <- getPMID(fit2$genes$ID,'lumiMouseAll.db')
fit2$genes$PMID <- NULL
fit2$genes$EG <- getEG(fit2$genes$ID,'lumiMouseAll.db')
fit2$genes$RefSeq <- nuID2RefSeqID(fit2$genes$ID, lib.mapping = "lumiMouseIDMapping")

M_MT_results <- topTable(fit2,number=length(which(topTable(fit2,adjust.method='BH',
                                                             number=length(fit2$genes$ID))$B<0)),
                           adjust.method='BH')