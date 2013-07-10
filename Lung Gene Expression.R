library("lumi")
library("lumiMouseIDMapping")
library("lumiMouseAll.db")
library("annotate")
library("mogene10stprobeset.db")
library("limma")

WD <- "C:/Users/rwill127/Documents/GitHub/EMT_Gene_Sig"
#find files
CTR_PATH <- paste(WD, 
                  "/Raw Lung Data/Ptran_not Normalized_Sample Probe Profile.txt",
                  sep = "")
CR_PATH <- paste(WD,
                 "/Raw Lung Data/NOT Norm sample Probe Profile-061209.txt",
                 sep ="")
OTHER_DATA_PATH <- paste(WD, 
                         "/Raw Lung Data/Not Normalized_Sample Probe Profile_120709.txt",
                         sep = "")
LIST_PATH <- paste(getwd(),
                   "/Raw Lung Data/TranP_annot.txt", 
                   sep = "")

file_list <- c(CTR_PATH, CR_PATH, OTHER_DATA_PATH)
x.lumi <- lumiR.batch(file_list, lib.mapping = "lumiMouseIDMapping", sampleInfoFile = LIST_PATH)
lumi.T <- lumiT(x.lumi, method = "log2")
dataMatrix <- exprs(lumi.T)
presentCount <- detectionCall(lumi.T)
selDataMatrix <- dataMatrix[presentCount > 0,]
probeList <- rownames(selDataMatrix)

annot <- as.data.frame(read.table(LIST_PATH, header = TRUE))
rownames(annot) <- annot[,"ID"]
batch <- c(rep("A", 16), rep("B", 8), rep("C", 8))
trimAnnot <- annot[which(rownames(annot)%in%colnames(selDataMatrix)),]
trimAnnot$batch <- batch
batch.trts <- factor(trimAnnot[,"batch"])
names(batch) <- trimAnnot[,"ID"]
exp <- factor(trimAnnot[,"group"])
names(exp) <- trimAnnot[,"ID"]

mapping_info <- as.data.frame(nuID2RefSeqID(probeList,
                              lib.mapping='lumiMouseIDMapping',
                              returnAllInfo = TRUE))
hk.trts <- rep(0, ncol(selDataMatrix))
names(hk.trts) <- colnames(selDataMatrix)
hk.names <- c("Gapdh",
              "Ppp2r1a",
              "Gsk3a",
              "Ndufa1",
              "Actb",
              "Rpl3",
              "Ppia",
              "Col6a1",
              "Mapkapk2")

gene_values <- function(gene_sym){
  probes <- rownames(mapping_info[which(mapping_info$Symbol == gene_sym),])
  gene_values <- selDataMatrix[probes,]
  if (length(probes) > 1) gene_values <- apply(gene_values, 2, mean)
  med_cent <- gene_values - median(gene_values)
  return(med_cent)
}

hk.values <- sapply(hk.names, gene_values)
hk.trts <- apply(hk.values, 1, mean)

batch_design <- model.matrix(~0 + exp + hk.trts)
batch_fit <- lmFit(selDataMatrix, batch_design)
batch_corrected <- selDataMatrix -
  batch_fit$coefficients[,'hk.trts']%*%t(batch_design[,'hk.trts'])

CR_CRT.contrasts <- makeContrasts(exprtTATwistRas-expras,
                                  levels=batch_design)
fit2 <- eBayes(contrasts.fit(batch_fit, CR_CRT.contrasts))

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

full_results <- topTable(fit2, number = length(which(topTable(fit2,adjust.method='BH',
                                                              number=length(fit2$genes$ID))$P.Val<1)),
                         adjust.method='BH')

raw_p_plot <- ggplot(full_results, aes(x = P.Value))
raw_p_plot <- raw_p_plot + 
  geom_histogram(binwidth = .02,
                 aes(fill = ..count..)) +
  scale_fill_gradient("Count", 
                      low = "dark green",
                      high = "dark red")
raw_p_plot

adj_p_plot <- ggplot(full_results, aes(x = adj.P.Val))
adj_p_plot <- adj_p_plot + 
  geom_histogram(binwidth = .01,
                 aes(fill = ..count..)) +
  scale_fill_gradient("Count", 
                      low = "dark green",
                      high = "dark red") +
  geom_vline(xintercept = .05)
adj_p_plot

