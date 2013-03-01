library("annotate")
library("mogene10sttranscriptcluster.db")

group_comparison <- function(contrast){
    contrasts <- makeContrasts(contrasts = contrast,
                               levels=batchDesign)
    contrast.fit <- eBayes(contrasts.fit(batchFit[probesSelect,],
                                         contrasts))
    contrast.fit$genes$Symbol <- getSYMBOL(
        contrast.fit$genes$ID,
        'mogene10sttranscriptcluster.db')
    results <- topTable(contrast.fit,
                        number=length(which(topTable(
                            contrast.fit,adjust.method='BH',                                                                         number=length(probesSelect))$B>0)),
                        adjust.method='BH')
    return(results)
}

full_gene_lookup <- function(gene_symbol){
    gene_probes <- genes2Probes[[gene_symbol]]
    gene_data <- batchFit[gene_probes,]
    gene_expression <- gene_data$coefficients[1:5]
    names(gene_expression) = c("AQA", "DQD", "F191G", "WT", "VEC")
    return(gene_expression)
}

pick_samples <- function(mutants){
    compSamples <- which(Annotation[colnames((expressionData)),'knock-in']%in%c('VEC','WT',mutants))
    return(compSamples)
}

gene_lookup <- function (gene, samples) {
    probes <- genes2Probes[[gene]]
    data <- batchCorrected[probes[1], samples]
    return(data)
}

make_gene_table <- function(filedata, samples){
    names_list <- as.vector(filedata)
    genedata <- sapply(names_list, gene_lookup, samples)
    return(genedata)
}

gsea_fisher <- function(sig_genes, geneset){
    geneset <- unlist(geneset)
    sig_in_set <- length(intersect(sig_genes, geneset))
    sig_not_in_set <- length(sig_genes) - sig_in_set
    not_sig_in_set <- length(geneset) - sig_in_set
    not_sig_not_in_set <- length(batchCorrected[,1]) - not_sig_in_set
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

prep_for_GSEA <- function(subset_names){
    #Attaches probe IDs to gene names and capitalizes names
    #so that human GSEA sets will be recognized
    subset_probes <- probesSelect[subset_names]
    names(subset_probes) <- toupper(names(subset_probes))
    return(subset_probes)
}

find_genes_of_interest <- function(geneset_name, comparison_list){
    gs_members <- (GSEASets$genesets[which(GSEASets$geneset.names
                                           == geneset_name)])
    genes_of_interest <- intersect(unlist(toupper(comparison_list)),
                                   unlist(gs_members))
    return(genes_of_interest)
}

comparison_gene_lookup <- function(gene_symbol, comparison){
    return(comparison[which(
        toupper(comparison$Symbol) == gene_symbol),])
}

gsoi_lookup <- function(geneset_name, comparison_list, comparison){
    table <- t(sapply(find_genes_of_interest(geneset_name, comparison_list),
                      comparison_gene_lookup, comparison))
}

gs_lookup_output <- function(gs_lookup, filename){
    for (i in 1:length(gs_lookup)){
        write.table(names(gs_lookup[i]), file = filename, sep = ",",
                    append = TRUE, col.names = FALSE)
        write.table(gs_lookup[[i]], file = filename, sep = ",",
                    append = TRUE, col.names = (i < 2))
    }
}

