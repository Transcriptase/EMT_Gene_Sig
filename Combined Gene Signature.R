library('GSA')
library('GSEABase')
library('limma')

#WD <- "G:/R Workspace"
WD <- "E:/Russell/Documents/GitHub/EMT_Gene_Sig"


setwd(WD)
liver <- read.csv("Liver LM v LMT.csv")
prostate <- read.csv("MycCaP WT Twist v Vector.csv")
lung <- read.csv("CRT_CR_results.csv")

e1 <- new.env()
e1$liver <- liver$Symbol
e1$prostate <- prostate$Symbol
e1$lung <- lung$Symbol
AllComps <- as.list(e1)
#counts <- venn(AllComps)

shared_genes <- intersect(e1$liver, intersect(e1$lung, e1$prostate))
#write.csv(shared_genes, "Liver Lung Prostate Shared Gene Sig.csv")

build_table <- function(gene_list){
    full_table<-data.frame(row.names = gene_list)
    full_table$p.lung <- lung$adj.P.Val[which(unique(lung$Symbol)%in%shared_genes)]
    full_table$fc.lung <- 2^lung$logFC[which(unique(lung$Symbol)%in%shared_genes)]
    #full_table$direction.lung
    full_table$p.liver <- liver$Corrected.p.value[which(unique(liver$Symbol)%in%shared_genes)]
    full_table$fc.liver <- liver$FCAbsolute[which(unique(liver$Symbol)%in%shared_genes)]
    full_table$direction.liver <- liver$regulation[which(unique(liver$Symbol)%in%shared_genes)]
    full_table$p.prostate <- prostate$adj.P.Val[which(unique(prostate$Symbol)%in%shared_genes)]
    full_table$fc.prostate <- 2^prostate$logFC[which(unique(lung$Symbol)%in%shared_genes)]
    #full_table$direction.prostate <- apply(prostate$logFC, 1, function(x) if (x>0) return ("up") else return ("down"))
    return(full_table)
}
info <- build_table(shared_genes)

probes2genes <- as.list(mogene10sttranscriptclusterSYMBOL)
all_genes <- unique(unlist(probes2genes))

gsea_fisher <- function(sig_genes, geneset){
    geneset <- unlist(geneset)
    sig_in_set <- length(intersect(sig_genes, geneset))
    sig_not_in_set <- length(sig_genes) - sig_in_set
    not_sig_in_set <- length(geneset) - sig_in_set
    not_sig_not_in_set <- length(all_genes) - not_sig_in_set
    data <- matrix(c(sig_in_set, not_sig_in_set, sig_not_in_set, not_sig_not_in_set), nrow = 2)
    test <- fisher.test(data,alternative='greater')
    return(test$p.value)
}

subset_gsea_test <- function(sig_genes){
    gsea_tests <- p.adjust(sapply(GSEASets$genesets, function(x){gsea_fisher(names(sig_genes), x)}), method = "BH")
    names(gsea_tests) <- GSEASets$geneset.names
    sig_genesets <-gsea_tests[which(gsea_tests < .05)]
    return(sig_genesets)
}
GSEASets <- GSA.read.gmt('c2.all.v3.0.symbols.gmt')
subset_gsea_test(toupper(shared_genes))