library("lumi")
library("lumiMouseIDMapping")
library("lumiMouseAll.db")
library("annotate")
library("mogene10stprobeset.db")
library("limma")
library("ggplot2")
library("sva")
library("simpleaffy")

#home wd
#WD = "E:/Russell/Documents/GitHub/EMT_Gene_Sig"
#work wd
WD = "C:/Users/rwill127/Documents/GitHub/EMT_Gene_Sig"
setwd(WD)

PATH = paste(getwd(), "/Raw Liver Data", sep = "")
file_list <- list.files(path = PATH)
file_list <- paste(PATH, file_list, sep = "/")

# read raw data
x.lumi <- lumiR.batch(file_list, lib.mapping = "lumiMouseIDMapping")

# assess quality of RNA in raw data
boxplot(log2(exprs(x.lumi)),las=2)

# filter bad quality samples
bad_samples <- c(15, 20)
good_samples <- c(1:20)[-bad_samples]
filtered.lumi <- x.lumi[, good_samples]
boxplot(log2(exprs(filtered.lumi)), las = 2)

# play around with normalization!!!!
lumi.T <- lumiT(filtered.lumi, method = "vst") # may want to try vst, log2 just log transforms

# continue from here!!!
dataMatrix <- exprs(lumi.T)
presentCount <- detectionCall(lumi.T)
selDataMatrix <- dataMatrix[presentCount > 0,]
probeList <- rownames(selDataMatrix)

### evaluate based upon distribution of expression values


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
annot <- annot[good_samples,]

boxplot(selDataMatrix,las=2,names=annot$sample)

pri_only <- which(annot$site == "Pri")
priDataMatrix <- selDataMatrix[,pri_only]
pri_annot <- annot[pri_only,]
pri_batch <- factor(pri_annot$batch)
pri_exp <- factor(pri_annot$exp)

batch.trts <- factor(annot$batch)
exp.trts <- factor(annot$exp)

rawclust <- hclust(dist(t(selDataMatrix)))
plot(rawclust,
     labels = paste(annot$exp,
                    annot$batch,
                    sep = " "
     )
)

# clustering based upon correlation (less senstive to changes in scale)
plot(standard.pearson(selDataMatrix),
          labels = paste(annot$exp,
                    annot$batch,
                    sep = " "
     )
)

batch_design <- model.matrix(~0 + exp.trts + batch.trts)
batch_fit <- lmFit(selDataMatrix, batch_design)
batch_corrected <- selDataMatrix -
    batch_fit$coefficients[,'batch.trtsB']%*%t(batch_design[,'batch.trtsB'])

batchclust <- standard.pearson(batch_corrected)
plot(batchclust,
     labels = paste(annot$exp,
                    annot$batch,
                    sep = " "
     )
)

comcorrected <- ComBat(selDataMatrix, annot$batch, annot$exp)
comclust <- standard.pearson(comcorrected)
plot(comclust,
     labels = paste(annot$exp,
                    annot$batch,
                    sep = " "
     )
)



design <- model.matrix(~0 + exp.trts)
fit <- lmFit(selDataMatrix, design)
MT_M.contrasts <- makeContrasts(exp.trtsPri_MT-exp.trtsPri_M,
                                levels=design)

fit2 <- eBayes(contrasts.fit(fit, MT_M.contrasts))

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

comfit <- lmFit(comcorrected, design)
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

priDataMatrix <- normalizeBetweenArrays()
pri_com <- ComBat(priDataMatrix, pri_annot$batch, pri_annot$exp)
pricom_clust <- standard.pearson(pri_com)
plot(pricom_clust, labels = pri_annot$oncogene)
pri_design <- model.matrix(~0 + pri_exp + pri_batch)
pri_fit <- lmFit(pri_com, pri_design)
priMT_M.contrasts <- makeContrasts(pri_expPri_MT-pri_expPri_M,
                                levels=pri_design)
pri_fit2 <- eBayes(contrasts.fit(pri_fit, priMT_M.contrasts))
priM_MT_results <- topTable(pri_fit2,number=length(which(topTable(pri_fit2,adjust.method='BH',
                                                                 number=length(pri_fit2$genes$ID))$adj.P.Val<1)),
                            adjust.method='BH')

raw_p_plot <- ggplot(priM_MT_results, aes(x = P.Value))
raw_p_plot + stat_bin()

adj_p_plot <- ggplot(priM_MT_results, aes(x = adj.P.Val))
adj_p_plot + stat_bin()

#ID Myc-only samples
myc_samples <- colnames(priDataMatrix)[which(pri_annot$oncogene == "M")]

#init vector to track twist expression
twist.trts <- rep(0, ncol(priDataMatrix))
names(twist.trts) <- colnames(priDataMatrix)

#get Twist values
twist_values <- priDataMatrix[which(getSYMBOL(rownames(priDataMatrix), "lumiMouseAll.db")== "Twist1"),]

twiDataMatrix <- priDataMatrix[, which(colnames(priDataMatrix) != "DB09_D4399_PriM")]
twi_annot <- pri_annot[1:8,]
twi_batch <- factor(twi_annot$batch)
twi_exp <- factor(twi_annot$exp)
twi_clust <- standard.pearson(twiDataMatrix)
plot(twi_clust)

twi_design <- model.matrix(~0 + twi_exp + twi_batch)
twiMT_M.contrasts <- makeContrasts(twi_expPri_MT-twi_expPri_M,
                                   levels=twi_design)

fit <- lmFit(twiDataMatrix, twi_design)
fit2 <- eBayes(contrasts.fit(fit, twiMT_M.contrasts))

fit2$genes$Symbol <- getSYMBOL(twi_fit2$genes$ID,'lumiMouseAll.db')
fit2$genes$RefSeq <- nuID2RefSeqID(twi_fit2$genes$ID, lib.mapping = "lumiMouseIDMapping")
results <- topTable(fit2,number=length(which(topTable(fit2,adjust.method='BH',
                                                      number=length(fit2$genes$ID))$adj.P.Val<.05)),
                    adjust.method='BH')
genes_only <- results[which(results$Symbol != "NA"),]
write.csv(genes_only, "Liver Differential Genes.csv")

twiDataMatrix <- normalizeBetweenArrays(twiDataMatrix, method  = "cyclicloess")
twi_com <- ComBat(twiDataMatrix, twi_annot$batch, twi_annot$exp)

twi_fit <- lmFit(twi_com, twi_design)

twi_fit2 <- eBayes(contrasts.fit(twi_fit, twiMT_M.contrasts))
twi_fit2$genes$Symbol <- getSYMBOL(twi_fit2$genes$ID,'lumiMouseAll.db')
#fit2$genes$PMID <- getPMID(fit2$genes$ID,'lumiMouseAll.db')
twi_fit2$genes$PMID <- NULL
twi_fit2$genes$EG <- getEG(twi_fit2$genes$ID,'lumiMouseAll.db')
twi_fit2$genes$RefSeq <- nuID2RefSeqID(twi_fit2$genes$ID, lib.mapping = "lumiMouseIDMapping")
twiM_MT_results <- topTable(twi_fit2,number=length(which(topTable(twi_fit2,adjust.method='BH',
                                                                  number=length(twi_fit2$genes$ID))$B>0)),
                            adjust.method='BH')
genes_only <- twiM_MT_results[which(twiM_MT_results$Symbol != "NA"),]
write.csv(genes_only, "Liver Differential Genes.csv")