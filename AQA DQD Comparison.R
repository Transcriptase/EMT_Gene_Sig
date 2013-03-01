DQDResults <- group_comparison("exp.trtsDQD - exp.trtsVEC")
AQAResults <- group_comparison("exp.trtsAQA - exp.trtsVEC")
DQDWTResults <- group_comparison("exp.trtsDQD - exp.trtsWT")

AQA_DQD_venn <- venn(list(WTResults$Symbol,
                          AQAResults$Symbol,
                          DQDResults$Symbol))

all_diff_genes <- union(union(WTResults$Symbol,
                        AQAResults$Symbol),
                        DQDResults$Symbol)

samples <- pick_samples(c("AQA", "DQD"))
gene_table <- make_gene_table(all_diff_genes, samples)
gene_matrix <-  t(data.matrix(gene_table))
colnames(gene_matrix) <- paste(Annotation[colnames(gene_matrix),'knock-in'],
                               Annotation[colnames(gene_matrix),'replicate'],
                               sep=".")

hm_lmat <- matrix(c(0, 0, 0, 0, 0,
                    0, 1, 1, 1, 2,
                    0, 0, 4, 0, 3), 3, 5, byrow = TRUE)
hm_lwid <- c(1, 6, 8, 6, 1)
hm_lhei <- c(1, 12, 3)

heatmap.2(gene_matrix, 
          col = redgreen (75),
          scale = "row",
          trace = "none",
          keysize = 0.5,
          margins = c(5,0),
          density.info = "none",
          dendrogram = "none",
          labRow = "",
          lmat = hm_lmat,
          lhei = hm_lhei,
          lwid = hm_lwid)

GSEASets <- GSA.read.gmt('c2.all.v3.0.symbols.gmt')
WT_not_AQA <- prep_for_GSEA(
    setdiff(WTResults$Symbol, AQAResults$Symbol))
DQD_diff_WT <- prep_for_GSEA(setdiff(DQDResults$Symbol, WTResults$Symbol))
WT_not_AQA_gs <- subset_gsea_test(WT_not_AQA)
DQD_diff_WT_gs <- subset_gsea_test(DQD_diff_WT)

AQA_gsoi <- c("PUIFFE_INVASION_INHIBITED_BY_ASCITES_DN",
              "MARKEY_RB1_ACUTE_LOF_UP",
              "PEREZ_TP53_TARGETS",
              "FRIDMAN_SENESCENCE_UP",
              "BROCKE_APOPTOSIS_REVERSED_BY_IL6",
              "MARTINEZ_RB1_TARGETS_UP",
              "MARTINEZ_TP53_TARGETS_UP",
              "MILI_PSEUDOPODIA_CHEMOTAXIS_DN",
              "MILI_PSEUDOPODIA_HAPTOTAXIS_DN",
              "DEBIASI_APOPTOSIS_BY_REOVIRUS_INFECTION_DN",
              "REACTOME_UNFOLDED_PROTEIN_RESPONSE")

AQA_gsoi_tables <- sapply(AQA_gsoi, gsoi_lookup, names(WT_not_AQA), WTResults)
gs_lookup_output(AQA_gsoi_tables, "WT not AQA Lookup.csv")
DQD_gsoi_tables <- gsoi_lookup(names(DQD_diff_WT_gs), names(DQD_diff_WT), DQDResults)
gs_lookup_output(DQD_gsoi_tables, "DQD not WT Lookup.csv")