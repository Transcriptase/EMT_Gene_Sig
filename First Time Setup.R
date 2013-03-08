source("http://bioconductor.org/biocLite.R")

cran_libraries <- c("gplots",
                   "RUnit")


bioc_libraries <- c("oligo",
                    'ClassDiscovery',
                    'AnnotationDbi',
                    'limma',
                    'pd.mogene.1.0.st.v1',
                    'AnnotationDbi',
                    'KEGG.db',
                    'org.Mm.eg.db',
                    'GO.db',
                    'annotate',
                    'mogene10sttranscriptcluster.db',
                    'mogene10stprobeset.db',
                    'GSA',
                    'GSEABase',
                    "lumi",
                    "lumiMouseIDMapping",
                    "lumiMouseAll.db")

install.packages(cran_libraries)
biocLite(bioc_libraries)
