#Process Lung mRNA data to prepare for sylamer analysis

#load libraries
library("lumi")
library("lumiMouseIDMapping")
library("lumiMouseAll.db")
library("annotate")
library("mogene10stprobeset.db")

#find files
CTR_PATH <- ("M:/MRSTran/Lab Notebook/Array/Lung 2009/Raw Data/CT and CTR 6-24-09/Copy of Ptran_for JH_012710/Data/Ptran_not Normalized_Sample Probe Profile.txt")
CR_PATH <- ("M:/MRSTran/Lab Notebook/Array/Lung 2009/Raw Data/CM, CR and CMR 6-12-09/NOT Norm sample Probe Profile-061209.txt")
OTHER_DATA_PATH <- ("M:/MRSTran/Lab Notebook/Array/Lung 2009/Raw Data/CT-LSL On, CT-LSL Off and LSL 12-7-09/Not Normalized_Sample Probe Profile_120709.txt")
LIST_PATH <- ("M:/MRSTran/Lab Notebook/Array/Lung 2009/Analysis/Single set GSEA - All 3-2010/TranP_annot.txt")

file_list <- c(CTR_PATH, CR_PATH, OTHER_DATA_PATH)
x.lumi <- lumiR.batch(file_list, lib.mapping = "lumiMouseIDMapping", sampleInfoFile = LIST_PATH)
lumi.T <- lumiT(x.lumi, method = "log2")
dataMatrix <- exprs(lumi.T)
presentCount <- detectionCall(lumi.T)
selDataMatrix <- dataMatrix[presentCount > 0,]
probeList <- rownames(selDataMatrix)

annot <- as.data.frame(read.table(LIST_PATH, header = TRUE))
rownames(annot) <- annot[,"ID"]
trimAnnot <- annot[which(rownames(annot)%in%colnames(selDataMatrix)),]
exp <- factor(trimAnnot[,"group"])
names(exp) <- trimAnnot[,"ID"]
design <- model.matrix(~0 + exp)
fit <- lmFit(selDataMatrix, design)
CR_CRT.contrasts <- makeContrasts(exprtTATwistRas-expras,
                                  levels=design)
fit2 <- eBayes(contrasts.fit(fit, CR_CRT.contrasts))

fit2$genes$Symbol <- getSYMBOL(fit2$genes$ID,'lumiMouseAll.db')
#fit2$genes$PMID <- getPMID(fit2$genes$ID,'lumiMouseAll.db')
fit2$genes$PMID <- NULL
fit2$genes$EG <- getEG(fit2$genes$ID,'lumiMouseAll.db')
fit2$genes$RefSeq <- nuID2RefSeqID(fit2$genes$ID, lib.mapping = "lumiMouseIDMapping")

gene_list <- cbind(fit2$genes$RefSeq, fit2$coefficients)
write.csv(gene_list, "CRT_CR_ordered_gene_list.csv")
CR_CRT_results <- topTable(fit2,number=length(which(topTable(fit2,adjust.method='BH',
                                                                      number=length(fit2$genes$ID))$adj.P.Val<0.05)),
                      adjust.method='BH')
write.csv(CR_CRT_results, "CRT_CR_results.csv")

MR_R.contrasts <- makeContrasts(expmycras-expras, levels = design)
fit3 <- eBayes(contrasts.fit(fit, MR_R.contrasts))
fit3$genes$Symbol <- getSYMBOL(fit3$genes$ID,'lumiMouseAll.db')
fit3$genes$RefSeq <- nuID2RefSeqID(fit3$genes$ID, lib.mapping = "lumiMouseIDMapping")

MR_R_gene_list <- cbind(fit3$genes$RefSeq, fit2$coefficients)
write.csv(gene_list, "CRT_CR_ordered_gene_list.csv")