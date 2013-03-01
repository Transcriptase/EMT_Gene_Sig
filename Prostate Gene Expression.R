library(oligo)
library('ClassDiscovery')
library('AnnotationDbi')
library('limma')
library('pd.mogene.1.0.st.v1')
library('gplots')
library('KEGG.db')
library('org.Mm.eg.db')
library('GO.db')
library('annotate')
library('mogene10sttranscriptcluster.db')
library('mogene10stprobeset.db')
library('GSA')
library('GSEABase')

#Get paths of CEL files
celFiles <- list.files(path=
    'G:/Project - PCa Twist1 Mutant Gene Expression/Raw CEL Files',
    pattern='cel',ignore.case=T,full.names=T)
#Use oligo library function to read CEL files into raw intensity data
rawExpressionData <- read.celfiles(celfile.path=celFiles)
#Use oligo library RMA fuction to do preprocessing
expressionData <- rma(rawExpressionData,target='core')
#Use AnnotationDbi function to map between gene symbols and probe numbers
#as a Bimap object
genes2Probes <- revmap(as.list(mogene10sttranscriptclusterSYMBOL))
#Use BioBase function to format the data as an expression set with only
#the probes representing genes.
expressionData <- exprs(expressionData)[unique(unlist(genes2Probes)),]
#Annotate the samples so that the files can be id'd by various aspects
Annotation <- cbind(seq(from=1,to=15,by=1),
                    rep(c('VEC','WT','AQA','DQD','F191G'),3),
                    rep(seq(1,3),each=5),
                    c('JF_A1_1_(MoGene-1_0-st-v1).CEL', 
                      'JF_A2_2_(MoGene-1_0-st-v1).CEL',
                      'JF_A3_3_(MoGene-1_0-st-v1).CEL', 
                      'JF_A4_04_(MoGene-1_0-st-v1).CEL', 
                      'JF_A5_5_(MoGene-1_0-st-v1).CEL', 
                      'JF_A6_6_(MoGene-1_0-st-v1).CEL', 
                      'JF_A7_7_(MoGene-1_0-st-v1).CEL', 
                      'JF_A8_08_(MoGene-1_0-st-v1).CEL',
                      'JF_B1_9_(MoGene-1_0-st-v1).CEL',
                      'JF_B2_10_(MoGene-1_0-st-v1).CEL', 
                      'JF_B3_11_(MoGene-1_0-st-v1).CEL',
                      'JF_B4_12_(MoGene-1_0-st-v1).CEL', 
                      'JF_B5_13_(MoGene-1_0-st-v1).CEL',
                      'JF_B6_14_(MoGene-1_0-st-v1).CEL',
                      'JF_B7_15_(MoGene-1_0-st-v1).CEL'),
                    c(rep('A',8),rep('B',7)))

colnames(Annotation) <- c('sample','knock-in','replicate','CelFile','batch')
row.names(Annotation) <- Annotation[,'CelFile']

#create vectors id'ing samples by batch and by knockin
batch.trts <- factor(Annotation[colnames((expressionData)),'batch'])
exp.trts <- factor(Annotation[colnames((expressionData)),'knock-in'])

#pull the filenames of the vector samples
VECSamples <- colnames(expressionData)[grep('VEC',Annotation[colnames(expressionData),'knock-in'])]

#initialize a vector to track Twist expression
twst.trts <- rep(0, ncol(expressionData))
names(twst.trts) <- colnames(expressionData)

#extract probe-level data (previous extraction was gene-level) and map it
expressionProbeData <- rma(rawExpressionData,target='probeset')
genes2Probesset <- revmap(as.list(mogene10stprobesetSYMBOL))

#pull the Twist expression levels
twistValues <- exprs(expressionProbeData)[genes2Probesset[['Twist1']],]
#select the Twist probe number with the greatest internal range?
twistValues <- twistValues[which.max(apply(twistValues,1,function(x){max(x)-min(x)})),]
#Sets the value of the Twist status vector for the VEC knock-ins equal to the mean-centered
#twist expression of the samples
twst.trts[VECSamples] <- twistValues[VECSamples] - mean(twistValues[VECSamples])
#sets the value of the non-VEC samples to be centered to the mean 
#twist expression for the remaining samples
twst.trts[setdiff(colnames(expressionData),VECSamples)] <- twistValues[setdiff(colnames(expressionData),
   VECSamples)] - mean(twistValues[setdiff(colnames(expressionData),
    VECSamples)])

#sets the model matrix for the batch correction by combining the knockin, batch, and twist expression factors
batchDesign <- model.matrix(~0+exp.trts+twst.trts+batch.trts)
# linear model fit 
batchFit <- lmFit(expressionData, batchDesign)
#subtract out batch and twist expression effects
batchCorrected <- expressionData -
    batchFit$coefficients[,c('batch.trtsB','twst.trts')]%*%t(batchDesign[,c('batch.trtsB','twst.trts')])
#set up contrast matrix for WT v VEC
compSample <- which(Annotation[colnames((expressionData)),'knock-in']%in%c('VEC','WT','F191G'))

exprsRange <- apply(batchCorrected[,compSample],1,function(x){max(x)-min(x)})
probesSelect <- sapply(genes2Probes,function(x){x[which.max(exprsRange[x])]})
WTVec.contrasts <- makeContrasts(exp.trtsWT-exp.trtsVEC,
                                 levels=batchDesign)
#get est. coeefs and stderrs, smoothed with bayesian technique
WTVec.contrast.fit <- eBayes(contrasts.fit(batchFit[probesSelect,],
                                           WTVec.contrasts))
#not sure i understand this data structure, but this annotates the results with gene symbols
WTVec.contrast.fit$genes$Symbol <- getSYMBOL(WTVec.contrast.fit$genes$ID,'mogene10sttranscriptcluster.db')

#returns the top genes, benjamini-hochberg corrected, B > 0 is the criteria how is this related to adjsuted p-value?
WTResults <- topTable(WTVec.contrast.fit,number=length(which(topTable(WTVec.contrast.fit,adjust.method='BH',
                                                                      number=length(probesSelect))$B>0)),
                      adjust.method='BH')

#repeat process with F191G v WT
F191GWT.contrasts <- makeContrasts(exp.trtsF191G-exp.trtsWT,
                                   levels=batchDesign)
F191GWT.contrast.fit <- eBayes(contrasts.fit(batchFit[probesSelect,],
                                             F191GWT.contrasts))
F191GWT.contrast.fit$genes$Symbol <- getSYMBOL(F191GWT.contrast.fit$genes$ID,'mogene10sttranscriptcluster.db')

F191GWTResults <- topTable(F191GWT.contrast.fit,number=length(which(topTable(F191GWT.contrast.fit,adjust.method='BH',
                                                                             number=length(probesSelect))$B>0)),
                           adjust.method='BH')
#repeat method for F191G v Vec
F191GVec.contrasts <- makeContrasts(exp.trtsF191G-exp.trtsVEC,
                                    levels=batchDesign)
F191GVec.contrast.fit <- eBayes(contrasts.fit(batchFit[probesSelect,],
                                              F191GVec.contrasts))
F191GVec.contrast.fit$genes$Symbol <- getSYMBOL(F191GVec.contrast.fit$genes$ID,'mogene10sttranscriptcluster.db')

F191GResults <- topTable(F191GVec.contrast.fit,number=length(which(topTable(F191GVec.contrast.fit,adjust.method='BH',
                                                                            number=length(probesSelect))$B>0)),
                         adjust.method='BH')
#Perform analysis of overlaps
e1 <- new.env()
e1$F191GvWT <- F191GWTResults$Symbol
e1$F191GvVec <- F191GResults$Symbol
e1$WTvVec <- WTResults$Symbol
AllComps <- as.list(e1)
counts <- venn(AllComps)
#works for now, can make prettier with VennDiagram package

#get list of WT-Vec-Only genes
WTVecOnly <- setdiff(WTResults$Symbol, union(F191GWTResults$Symbol, F191GResults$Symbol))
WT_not_F191G <- setdiff(WTResults$Symbol, intersect(F191GResults$Symbol, WTResults$Symbol))

#make heat map
pick_samples <- function(mutant){
    compSamples <- which(Annotation[colnames((expressionData)),'knock-in']%in%c('VEC','WT',mutant))
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

samples <- pick_samples("F191G")
gene_table <- make_gene_table(WT_not_F191G, samples)
gene_matrix <-  t(data.matrix(gene_table))
colnames(gene_matrix) <- paste(Annotation[colnames(gene_matrix),'knock-in'],
                               Annotation[colnames(gene_matrix),'replicate'],
                               sep=".")
hm_lmat <- matrix(c(0,0,0,0,0,0,1, 1, 1, 2, 0, 0, 4, 0, 3), 3, 5, byrow = TRUE)
hm_lwid <- c(1, 6, 8, 6, 1)
hm_lhei <- c(1, 12, 3)

heatmap.2(gene_matrix, 
          col = redgreen (75),
          scale = "row",
          trace = "none",
          keysize = 0.5,
          margins = c(5,0),
          density.info = "none",
          dendrogram = "column",
          labRow = "",
          #lmat = hm_lmat,
          #lhei = hm_lhei,
          #lwid = hm_lwid
          )

#GSEA/GO Analysis
GSEASets <- GSA.read.gmt('c2.all.v3.0.symbols.gmt')
WTVecOnly_probes <- probesSelect[WTVecOnly$Symbol]
AllWT_probes <- probesSelect[WTResults$Symbol]
F191G_and_WT_probes <- probesSelect[intersect(WTResults$Symbol, F191GResults$Symbol)]
WT_not_F191G_probes <- probesSelect[setdiff(WTResults$Symbol, intersect(WTResults$Symbol, F191GResults$Symbol))]
names(WTVecOnly_probes) <- toupper(names(WTVecOnly_probes))
names(AllWT_probes) <- toupper(names(AllWT_probes))
names(F191G_and_WT_probes) <-toupper(names(F191G_and_WT_probes))
names(WT_not_F191G_probes) <- toupper(names(WT_not_F191G_probes))

subsets_to_test <- list(WTVecOnly_probes, AllWT_probes, F191G_and_WT_probes,
                        WT_not_F191G_probes)

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

gsea_results <- sapply(subsets_to_test, subset_gsea_test)
gsea_overlap <- venn(list(names(gsea_results[[1]]),
                          names(gsea_results[[2]]),
                          names(gsea_results[[3]]),
                          names(gsea_results[[4]])))