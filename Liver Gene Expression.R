library("lumi")
library("lumiMouseIDMapping")
library("lumiMouseAll.db")
library("annotate")
library("mogene10stprobeset.db")
library("limma")
library("ggplot2")
library("sva")

#home wd
#WD = "E:/Russell/Documents/GitHub/EMT_Gene_Sig"
#work wd
#WD = "C:/Users/rwill127/Documents/GitHub/EMT_Gene_Sig"
setwd(WD)

PATH = paste(getwd(), "/Raw Liver Data", sep = "")
file_list <- list.files(path = PATH)
file_list <- paste(PATH, file_list, sep = "/")

# read raw data
x.lumi <- lumiR.batch(file_list, lib.mapping = "lumiMouseIDMapping")

# assess quality of RNA in raw data
boxplot(log2(exprs(x.lumi)),las=2)

# filter bad quality samples

# play around with normalization!!!!
lumi.T <- lumiT(x.lumi, method = "log2") # may want to try vst, log2 just log transforms

# continue from here!!!
dataMatrix <- exprs(lumi.T)
presentCount <- detectionCall(lumi.T)
selDataMatrix <- dataMatrix[presentCount > 0,]
probeList <- rownames(selDataMatrix)

### evaluate based upon distribution of expression values
boxplot(selDataMatrix,las=2,names=annot$sample)

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

rawclust <- hclust(dist(t(selDataMatrix)))
plot(rawclust,
     labels = paste(annot$exp[which(annot$Id == colnames(selDataMatrix))],
                    annot$batch[which(annot$Id == colnames(selDataMatrix))],
                    sep = " "
     )
)

# clustering based upon correlation (less senstive to changes in scale)
plot(standard.pearson(selDataMatrix),
          labels = paste(annot$exp[which(annot$Id == colnames(selDataMatrix))],
                    annot$batch[which(annot$Id == colnames(selDataMatrix))],
                    sep = " "
     )
)

batch_design <- model.matrix(~0 + exp.trts + batch.trts)
batch_fit <- lmFit(selDataMatrix, batch_design)
batch_corrected <- selDataMatrix -
    batch_fit$coefficients[,'batch.trtsB']%*%t(batch_design[,'batch.trtsB'])

batchclust <- hclust(dist(t(batch_corrected)))
plot(batchclust,
     labels = paste(annot$exp[which(annot$Id == colnames(selDataMatrix))],
                    annot$batch[which(annot$Id == colnames(selDataMatrix))],
                    sep = " "
     )
)

comcorrected <- ComBat(selDataMatrix, annot$batch, annot$exp)
comclust <- hclust(dist(t(comcorrected)))
plot(comclust,
     labels = paste(annot$exp[which(annot$Id == colnames(selDataMatrix))],
                    annot$batch[which(annot$Id == colnames(selDataMatrix))],
                    sep = " "
     )
)


MT_M.contrasts <- makeContrasts(exp.trtsPri_MT-exp.trtsPri_M,
                                levels=batch_design)

fit2 <- eBayes(contrasts.fit(batch_fit, MT_M.contrasts))

fit2$genes$Symbol <- getSYMBOL(fit2$genes$ID,'lumiMouseAll.db')
#fit2$genes$PMID <- getPMID(fit2$genes$ID,'lumiMouseAll.db')
fit2$genes$PMID <- NULL
fit2$genes$EG <- getEG(fit2$genes$ID,'lumiMouseAll.db')
fit2$genes$RefSeq <- nuID2RefSeqID(fit2$genes$ID, lib.mapping = "lumiMouseIDMapping")

M_MT_results <- topTable(fit2,number=length(which(topTable(fit2,adjust.method='BH',
                                                             number=length(fit2$genes$ID))$adj.P.Val<1)),
                           adjust.method='BH')

raw_p_plot <- ggplot(M_MT_results, aes(x = P.Value))
raw_p_plot + stat_bin()

adj_p_plot <- ggplot(M_MT_results, aes(x = adj.P.Val))
adj_p_plot + stat_bin()

comfit <- lmFit(comcorrected, batch_design)
comfit2 <- eBayes(contrasts.fit(comfit, MT_M.contrasts))
comfit2$genes$Symbol <- getSYMBOL(comfit2$genes$ID,'lumiMouseAll.db')
#fit2$genes$PMID <- getPMID(fit2$genes$ID,'lumiMouseAll.db')
comfit2$genes$PMID <- NULL
comfit2$genes$EG <- getEG(comfit2$genes$ID,'lumiMouseAll.db')
comfit2$genes$RefSeq <- nuID2RefSeqID(comfit2$genes$ID, lib.mapping = "lumiMouseIDMapping")

comM_MT_results <- topTable(comfit2,number=length(which(topTable(comfit2,adjust.method='BH',
                                                              number=length(comfit2$genes$ID))$adj.P.Val<1)),
                            adjust.method='BH')

raw_p_plot <- ggplot(comM_MT_results, aes(x = P.Value))
raw_p_plot + stat_bin()

adj_p_plot <- ggplot(comM_MT_results, aes(x = adj.P.Val))
adj_p_plot + stat_bin()
