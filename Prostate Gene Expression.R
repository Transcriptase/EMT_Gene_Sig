library(oligo)
library('AnnotationDbi')
library('limma')
library('pd.mogene.1.0.st.v1')
library('AnnotationDbi')
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
PATH = paste(getwd(), "/Raw Prostate Data", sep = "")
celFiles <- list.files(path= PATH,
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

compSample <- which(Annotation[colnames((expressionData)),'knock-in']%in%c('VEC','WT'))
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
prostate_results <- topTable(WTVec.contrast.fit,number=length(which(topTable(WTVec.contrast.fit,adjust.method='BH',
                                                                      number=length(probesSelect))$B>0)),
                      adjust.method='BH')
write.csv(prostate_results, "MycCap Twist v Vector Gene Expression.csv")