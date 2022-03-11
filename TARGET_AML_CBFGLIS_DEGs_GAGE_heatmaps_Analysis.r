#Jenny Smith
#original Author: Hamid Bolouri
#February 15, 2017 
#Purpose: Perform DE analysis of batch corrected, low pass, RNAseq data from 2016April_BCCA_illumina_data



# #installthe R package GAGE.
# source("https://bioconductor.org/biocLite.R")
# biocLite("gage")
# #install the R package GageData.
# source("https://bioconductor.org/biocLite.R")
# biocLite("gageData")
# #install R package Pathview. 
# source("https://bioconductor.org/biocLite.R")
# biocLite("pathview")
# source("https://bioconductor.org/biocLite.R")
# biocLite("edgeR")
# biocLite("limma")
# source("https://bioconductor.org/biocLite.R")
# biocLite("EnsDb.Hsapiens.v79")
# biocLite("Homo.sapiens")
# biocLite("AnnotationHub")
# biocLite("BSgenome.Hsapiens.UCSC.hg19")
# source("https://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
# install.packages("VennDiagram")


#load the Gage library 
library("gage")
#load the library GageData
library("gageData")
#load the library Pathview. 
library("pathview")
library(limma)
library(edgeR)
library(biomaRt)
library(VennDiagram)
library(Homo.sapiens)
library(AnnotationHub)
library(dplyr)
library(EnsDb.Hsapiens.v79)


getwd()

################ prepare low depth Datasets ############

#setwd("X:/fast/meshinchi_s/RNAseq_BCCA28Apr2016/gene_coverage_GSC-1367/")

setwd("~/RNA_seq_Analysis/2017.02.15_CBF-GLIS_DEG/ExpressionData")


#From Hamids previous analysis which used batch correction methods on the 
#raw read counts - which also, are not intergers, but fractional counts that were rounded. 
correctedExp <- get(load('combatAdjustedExp.RData'))

head(correctedExp[,1:5])
class(correctedExp)
dim(correctedExp) #486 Samples

#must have been filtered for read counts - 30,048 genes. 
dimnames(correctedExp) 

#any values less than zero (negative) must be made to be zero counts. 
correctedExp[which(correctedExp < 0, arr.ind=TRUE)] <- 0
dim(correctedExp)

#ROWSUMS HERE to remove very low count samples
correctedExp <- correctedExp[ rowSums(correctedExp) > 10, ] #30,048
dim(correctedExp)

############### Prepare High Depth Datasets ################
# 
setwd("~/RNA_seq_Analysis/2017.02.15_CBF-GLIS_DEG/ExpressionData/")

#read in the tsv file. 
Counts <- read.csv("rawCounts_TARGET_AML_Aug2014.csv")


head(Counts)
names <- Counts$ensemblSymbol
Counts <- Counts[,3:ncol(Counts)]

Dx_counts <- subset(Counts, select= grepl("T.+0[3,9]A.+", names(Counts)))

dim(Dx_counts) #58450

#Convert the target barcode into only the target USI
Dx_ID <- gsub("T.+([A-Z]+{6}).+", "\\1", colnames(Dx_counts))
Dx_ID <- gsub("\\.", "-", Dx_ID)
Dx_ID

colnames(Dx_counts) <- Dx_ID
rownames(Dx_counts) <- names

#Characterize the dataset
head(Dx_counts)
dim(Dx_counts) #58450

#ROWSUMS HERE to remove very low count samples
Dx_counts <- Dx_counts[rowSums(Dx_counts) > 10, ]
dim(Dx_counts) #42020


###################### Same for both datasets below this ####################

#log2 transformation of batch corrected read counts. This puts all expression values on the same scale. 
log2Low <- apply(correctedExp, 2, function(x) log2(x + 1))
log2High <- apply(Dx_counts, 2, function(x) log2(x + 1))


colnames(log2Low)
colnames(log2High)

####TMM Normalization method #######

lowCounts <- DGEList(correctedExp)

highCounts <- DGEList(Dx_counts)


dge_low <- calcNormFactors(lowCounts)

head(dge_low$counts[,1:5])


#wow very big range of normalization factors! 
range(dge_low$samples$norm.factors)
#0.1220644 5.4637443

dge_high <- calcNormFactors(highCounts)

dge_high


range(dge_high$samples$norm.factors)
#0.2228303 1.8429547

#prepare for DE using CPM to normalize
logCPM_low <- cpm(dge_low, log=TRUE, prior.count = 3)

logCPM_low


logCPM_high <- cpm(dge_high, log=TRUE, prior.count = 3)


logCPM_high[,1:10]


#read in the file with the CBF positive samples 
CBFA2T3_GLIS2_RegNo <- read.csv("~/reference_mapping-files/CBFA2T3_GLIS2_positive_samples.csv")[,1:4]
# CBFA2T3_GLIS2_RegNo <- read.csv("/Users/Shared/2017.02_MeshinchiLab/reference_mapping-files/CBFA2T3_GLIS2_positive_samples.csv")[,1:4]

CBFA2T3_GLIS2_RegNo


regNo_highDepth <- subset(CBFA2T3_GLIS2_RegNo, grepl("Hi", CBFA2T3_GLIS2_RegNo$X))
regNo_lowDepth <- subset(CBFA2T3_GLIS2_RegNo, grepl("Low", CBFA2T3_GLIS2_RegNo$X))

regNo_highDepth 
regNo_lowDepth

#subset the batch corrected data into CBF-GLIS positive and CBF-GLIS negative groups
posH <- as.character(regNo_highDepth$USI)
posH <- unique(posH)
posL <- as.character(regNo_lowDepth$USI)
posL <- unique(posL)

negH <- colnames(log2High[ , !(colnames(log2High) %in% regNo_highDepth$USI)])
negH <- negH[!grepl("BM", negH)]

negL <- colnames(log2Low[ , !(colnames(log2Low) %in% regNo_lowDepth$USI)])
negL <- negL[!grepl("BM", negL)]

# colnames(logCPM_high)
# negH <- colnames(logCPM_high[ , !(colnames(logCPM_high) %in% regNo_highDepth$USI)])
# negL <- colnames(logCPM_low[ , !(colnames(logCPM_low) %in% regNo_lowDepth$USI)])
# negL <- negL[!grepl("BM", negL)] 

posH
negH
posL
negL

length(posH)
length(negH)

length(posL)
length(negL)


#Check out the dimensions of the data subset that will be analyzed. 
# dim(logCPM_low[ ,c(posL, negL)])
# dim(logCPM_high[ ,c(posH, negH)])

################# Differential Expression Analysis ###################

calcDE<-function(expData, g1, g2) {
  library(limma)
  
  designMatrix <- matrix(0, nrow=dim(expData)[2] , ncol=2)
  colnames(designMatrix) <- c("g1","g2")
  rownames(designMatrix) <- colnames(expData)
  designMatrix[g1, 1] <- 1
  designMatrix[g2, 2] <- 1
  
  fit <- lmFit(expData,designMatrix)
  tmp <- paste("g1","g2",sep="-")
  
  cont.matrix <- makeContrasts(contrasts=tmp,
                               levels=designMatrix)
  fit2<-contrasts.fit(fit, cont.matrix)
  fit2<-eBayes(fit2)
  DE<-topTable(fit2,adjust.method="BH",sort.by="P",
               number=20000,p.value=0.05)
  return(DE)
}



#Call the function identify DEGs for low counts 
DE_low <- calcDE(log2Low[ ,c(posL, negL)], posL, negL)							#
DE_low <- DE_low[which(abs(DE_low$logFC) > 1), ]									#1830 cpm versus 1960 log2

dim(DE_low)

#Call the function identify DEGs for high counts 
DE_high <- calcDE(log2High[ ,c(posH, negH)], posH, negH)							#
DE_high <- DE_high[which(abs(DE_high$logFC) > 1), ]									#1211 cpm versus 1240 log2

dim(DE_high)

write.table(DE_low, file="", sep="\t",
            col.names=TRUE, row.names=TRUE, quote=FALSE)


write.table(DE_high, file="", sep="\t",
            col.names=TRUE, row.names=TRUE, quote=FALSE)

write.csv(DE_low, file="DE_lowDepth_TMMCPM_2017.02.24.csv")


write.csv(DE_high, file="DE_highDepth_TMMCPM_2017.02.24.csv")




###########KEGG Ontology################

#Read in the limma results.
DE_highDepth <- read.delim("DE_highDepth_TMMCPM_2017.02.24.csv", sep = ',', stringsAsFactors = FALSE)
DE_lowDepth <- read.delim("DE_lowDepth_TMMCPM_2017.02.24.csv", sep = ",", stringsAsFactors = FALSE)


head(DE_highDepth)
dim(DE_highDepth)


head(DE_lowDepth)
dim(DE_lowDepth)

range(DE_highDepth$AveExpr)


#conver the ensmeble IDs to the gene symbol
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v79)

genesL <- DE_lowDepth$X
  
genesH <- DE_highDepth$X


#print the keytypes
# keytypes(org.Hs.eg.db)
# 
# columns(org.Hs.eg.db)

keytypes(EnsDb.Hsapiens.v79)

columns(EnsDb.Hsapiens.v79)


#geneID = ensebmble gene identifier in ENSG# format. 
gene_listH <- select(EnsDb.Hsapiens.v79, key=genesH, 
                     columns = c("ENTREZID", "SYMBOL"),
                     keytype = "GENEID")

#geneID = ensebmble gene identifier in ENSG# format. 
gene_listL <- select(EnsDb.Hsapiens.v79, key=genes, 
                     columns = c("ENTREZID", "SYMBOL"),
                     keytype = "GENEID")

#geneID = ensebmble gene identifier in ENSG# format. #FUCKING AWFUL 
# gene_listL <- select(org.Hs.eg.db, key=genesL, 
#                      columns = c("GENENAME"),
#                      keytype = "GENENAME")

#geneID = ensebmble gene identifier in ENSG# format. IT HAS HORRIBLE RESULTSSSS
# gene_listH <- select(org.Hs.eg.db, key=genesH, 
#                      columns = c("SYMBOL"),
#                      keytype = "ENSEMBL")


#WTF WHY IS EVERYTHING NA?
gene_listL
dim(gene_listL)


#check out the genelist
head(gene_listH)
dim(gene_listH) #1164 of 1211 DEGs 
tail(gene_listH)

dim(gene_listH[which(duplicated(gene_listH$SYMBOL)), ])

dim(gene_listH[which(!is.na(gene_listH$SYMBOL)), ])

gene_listH[which(is.na(gene_listH$SYMBOL)), ]

uniqueGeneListH <- gene_listH[which(!duplicated(gene_listH$GENEID)), ]
dim(uniqueGeneListH) #1164


#check out the genelist
head(gene_listL)
dim(gene_listL)
tail(gene_listL)

#Look at the column names - will be used for merging the dataframes. 
colnames(DE_highDepth)
colnames(DE_lowDepth)

colnames(gene_listH)
colnames(gene_listL)

dim(gene_listH[which(is.na(gene_listH$SYMBOL)), ])



#merge the  data frame with the gene_list data frame using column names. 
DE_highDepth <- merge(DE_highDepth, uniqueGeneListH, by.x = "X",
             by.y = "GENEID", all.x = TRUE)



#merge the  data frame with the gene_list data frame using column names. 
# DE_lowDepth <- merge(DE_lowDepth, uniqueGeneListL, by.x = 0,
#                       by.y = "external_gene_name", all.x = TRUE)




#check that the merge worked. 
colnames(DE_highDepth)
dim(DE_highDepth) #1211
head(DE_highDepth)

#check that the merge worked. 
colnames(DE_lowDepth)
dim(DE_lowDepth)
head(DE_lowDepth)


#rearrange the columns
DE_highDepth <- DE_highDepth[, c(1,8,9,2:7)]
head(DE_highDepth)

length(which(is.na(DE_highDepth$SYMBOL)))




# DE_lowDepth <- DE_lowDepth[, c(1,8,2:7)]
# head(DE_lowDepth)

#NOTE: Same number of rows as input datasets which were the original logFCs. 
dim(DE_highDepth)
dim(DE_lowDepth)


#save the dataframe.
#remove results wihtout a matched official gene symbol. This is temporary.
#MUST manually find the remaining ~75 genes (47 genes using EnsDb.Hsapiens.v79)

DE_highDepth <- DE_highDepth[which(!is.na(DE_highDepth$SYMBOL)), ]

dim(DE_highDepth)

write.csv(DE_highDepth, "DE_highDepth_TMMCPM_withGeneSymbol_2017.02.24.csv", quote = FALSE, 
          row.names = FALSE, col.names = TRUE)

#save the dataframe.
# write.csv(DE_lowDepth, "DE_lowDepth_withGeneSymbol_2017.02.20.csv", quote = FALSE, 
#           row.names = FALSE, col.names = TRUE)


getwd()
setwd("H:/RNA_seq_Analysis/2017.02.15_CBF-GLIS_DEG/")


#Read in the files to be analyzed. 
# highDepth <- read.delim("DE_highDepth_withGeneSymbol_2017.02.20.csv", sep= ",", stringsAsFactors = FALSE)
# lowDepth <- read.delim("DE_lowDepth_withGeneSymbol_2017.02.20.csv", sep= ",", stringsAsFactors = FALSE)

highDepth <- DE_highDepth
lowDepth <- DE_lowDepth

#WHY ARE THERE SO MANY NAs??? TWO DIFFERNET METHODS HAVE ALMOST THE SAME OUT PUT WITH NAs??? 
# highDepth <- highDepth[which(!is.na(highDepth$entrezgene)), ]
# 
# lowDepth <- lowDepth[which(!is.na(lowDepth$entrezgene)), ]


head(highDepth)
head(lowDepth)

dim(highDepth)
dim(lowDepth)

nrow(highDepth[which(highDepth$logFC > 1), ]) #911 cpm
nrow(highDepth[which(highDepth$logFC < -1), ]) #300 cpm

nrow(lowDepth[which(lowDepth$logFC > 1), ]) #1286
nrow(lowDepth[which(lowDepth$logFC < -1), ]) #544


#Load the kegg human datasets 
data("kegg.gs")
data("go.sets.hs")
data("egSymb")


#create objects to hold the information from the datasets loaded above
kg.hsa=kegg.gsets("hsa")

#sig = signaling, met = metabolic pathways. 
kegg.gs=kg.hsa$kg.sets[kg.hsa$sigmet.idx]

#convert to gene symbols
kegg.gs.sym <- lapply(kegg.gs, eg2sym)

#save the data 
save(kegg.gs.sym, file="kegg.hsa.sigmet.sym.gsets.RData")

lapply(kegg.gs.sym[1:3],head)
length(kegg.gs.sym)


#Hedhod gene list in entrez gene identifiers. 
SHH_path <- kegg.gs.sym$`hsa04340 Hedgehog signaling pathway`
WNT_path <- kegg.gs.sym$`hsa04310 Wnt signaling pathway`	
TGFb_path <- kegg.gs.sym$`hsa04350 TGF-beta signaling pathway`

length(SHH_path)



#construct a vector of log2 fold change values and use the 
#names() function to name each component with a ko ID 
#(if it has one). This will allow us to do GSEA for all 
#of the KEGG pathways, based on log2 fold change as our 
#per- gene measure of differential expression.

#make a variable to hold only the fold_changes column values
highDepth_foldchanges <- highDepth$logFC
lowDepth_foldchanges <- lowDepth$logFC


head(highDepth_foldchanges)


#names function will set the names for each  fold-change value as the 
#entrez gene IDs  or gene symbols
names(highDepth_foldchanges) <- highDepth$SYMBOL
names(lowDepth_foldchanges) <- lowDepth$X

which(is.na(names(highDepth_foldchanges)))


highDepth_foldchanges[1:10]
lowDepth_foldchanges[1:10]

length(highDepth_foldchanges)
length(lowDepth_foldchanges)


#run gage on vector of log2 fold changes to test for 
#enrichment of genes with extreme values 
DE_highDepth_test <- gage(highDepth_foldchanges, gsets = kegg.gs.sym, 
                   same.dir = FALSE)

DE_highDepth_SHH <- gage(highDepth_foldchanges, gsets = list(SHH_path), 
                         same.dir = FALSE)

lapply(DE_highDepth_test, head)
DE_highDepth_test

DE_highDepth[which(DE_highDepth$SYMBOL %in% SHH_path), ]
DE_highDepth[which(DE_highDepth$SYMBOL %in% WNT_path), ]
DE_highDepth[which(DE_highDepth$SYMBOL %in% TGFb_path), ]



#enrichment of genes with extreme values 
DE_lowDepth_test <- gage(lowDepth_foldchanges, gsets = kegg.gs.sym, 
                          same.dir = FALSE)

DE_lowDepth_SHH <- gage(lowDepth_foldchanges, gsets = list(SHH_path), 
                         same.dir = FALSE)


lapply(DE_lowDepth_test, head, 30)


#Save the Results. 
KEGG_high <- subset(DE_highDepth_test$greater, !is.na(DE_highDepth_test$greater[,4]))

dim(KEGG_high)

KEGG_low <- subset(DE_lowDepth_test$greater, DE_lowDepth_test$greater[,3] < 0.15)

dim(KEGG_low)

rownames(DE_highDepth_test$greater)
rownames(KEGG_low)

#save the files 
write.csv(KEGG_high, file = "DE_CBFGLIS_highDepth_KEGGpaths_2017.02.20.csv")

write.csv(KEGG_low, file = "DE_CBFGLIS_lowDepth_KEGGpaths_2017.02.20.csv")


# pathways <- rownames(DE_highDepth_test$greater)[1:10]
# pathways <- rownames(KEGG_high)

pathways <- rownames(KEGG_low)

#start and stop are the character positions in the 
#string "KO00000" for KO ids
ids <- substr(pathways, start = 1, stop = 8)
ids

length(ids)


#Draw the pathway with a color scale that reflects log2 
#fold change for each gene. The pathview() function will 
#write 3 files to your working directory. 
#ko04610.pathview.png is the one you want.

setwd("H:/RNA_seq_Analysis/2017.02.15_CBF-GLIS_DEG/GAGE/")
getwd()


i <- 1
pathList_high <- NULL
for (path in ids){
  print(path)
  name <- paste(path, i, sep = "_")
  pathList_high <- append(pathList_high, name)
  tmp <- pathview(gene.data=highDepth_foldchanges,
           pathway.id=path, species="hsa", new.signature=FALSE,
           trans.fun = list(gene = NULL, cpd = NULL),
           low = list(gene = "green", cpd = "blue"),
           mid = list(gene = "yellow", cpd = "gray"),
           high = list(gene = "red", cpd = "yellow"),
           na.col = "transparent", 
           same.layer = FALSE)
  assign(name, tmp)
  i <- i + 1
}

pathList_high

hsa04974_1$plot.data.gene
hsa04974_1$plot.data.cpd
hsa04151_2$plot.data.gene



i <- 1
pathList_low <- NULL
for (path in ids){
  print(path)
  name <- paste(path, i, sep = "_")
  pathList_low <- append(pathList_low, name)
  tmp <- pathview(gene.data=lowDepth_foldchanges,
                  pathway.id=path, species="hsa", new.signature=FALSE,
                  trans.fun = list(gene = NULL, cpd = NULL),
                  low = list(gene = "green", cpd = "blue"),
                  mid = list(gene = "yellow", cpd = "gray"),
                  high = list(gene = "red", cpd = "yellow"), 
                  na.col = "transparent", 
                  same.layer = FALSE)
  assign(name, tmp)
  i <- i + 1
}

pathList_low

hsa04640_1

#original one gene at a time. TRY THIS WITH LAYERING 
# pathview(gene.data=highDepth_foldchanges,
#          pathway.id="hsa04974", species="hsa", new.signature=FALSE,
#          trans.fun = list(gene = NULL, cpd = NULL),
#          low = list(gene = "green", cpd = "blue"),
#          mid = list(gene = "yellow", cpd = "gray"),
#          high = list(gene = "red", cpd = "yellow"), na.col = "transparent")


#513 subset of genes for GAGE?



##########Comparison of High and Low Depth genes in common##############

head(highDepth)
head(lowDepth)

dim(highDepth)
dim(lowDepth)

commonDEGs <- merge(highDepth[,c(3:4,8)], lowDepth[,c(1,3,7)], by.x = "external_gene_name", by.y = "Row.names")
head(commonDEGs)

# commonDEGs <- commonDEGs[, c(1:2,4:5,8:9)]

colnames(commonDEGs) <- c("Gene_Symbol", "LogFC_highDepth", "adj.P.Val_highDepth", "logFC_lowDepth", "adj.P.Val_lowDepth")
head(commonDEGs)
dim(commonDEGs)


#save the file. 
write.csv(commonDEGs, file = "commonDEGs_CBFGLIS_2017.02.20.csv", quote = FALSE, 
          row.names = FALSE)

common <- read.csv("commonDEGs_CBFGLIS_2017.02.20.csv", header = TRUE)



head(common)
nrow(common)


numBothUp <- nrow(common[which(common$LogFC_highDepth > 1 & common$logFC_lowDepth > 1), ])
numBothdown <- nrow(common[which(common$LogFC_highDepth < 1 & common$logFC_lowDepth < 1), ])
opposite <- common[which(common$LogFC_highDepth > 1 & common$logFC_lowDepth < 1), ]


genesH <- highDepth$external_gene_name
genesL <- lowDepth$Row.names


venny <- venn.diagram(list(highDepth = genesH, lowDepth = genesL), "CommonDEGs_CBFGLIS_2017.02.20.tiff", fill = c("cornflowerblue", "darkorchid1"), 
                      print.mode = c("raw","percent"), force.unique = FALSE, cat.cex = 0.8)





#####Comparison of Previously published ###########
setwd("H:/RNA_seq_Analysis/2017.02.15_CBF-GLIS_DEG")


#pupblished RNAseq data
Gruber <- read.csv("CBFA2T3-GLIS2_Fusion_RNAseq_T.Gruber_2012.csv", header = TRUE)

head(Gruber)
head(highDepth)
head(lowDepth)

dim(highDepth)
dim(lowDepth)


colnames(highDepth)
colnames(lowDepth)
colnames(Gruber)

upRegTranscripts_fromGruberSupp <- c("CBFA2T3", "GATA2", "HOXA9", "MN1", "FLI1",
                                      "NIPBL", "HOXB9", "NUP98", "KDM5A", "GRB10", "SDK1")


#high depth 

genesH
length(genesH)

genesGuber <- Gruber$Gene
length(genesGuber)

highDepth_Gruber <- merge(highDepth, Gruber, by.x = "external_gene_name", by.y = "Gene")

head(highDepth_Gruber)
dim(highDepth_Gruber)


numBothUp <- nrow(highDepth_Gruber[which(highDepth_Gruber$logFC >1 & highDepth_Gruber$log2.fold.change. > 1), ])

numBothdown <- nrow(highDepth_Gruber[which(highDepth_Gruber$logFC <1 & highDepth_Gruber$log2.fold.change. < 1), ])

opposite <- nrow(highDepth_Gruber[which(highDepth_Gruber$logFC > 1 & highDepth_Gruber$log2.fold.change. < 1), ])


#save the file 
write.csv(highDepth_Gruber, file="highDepth_TARGET_withGruber_2017.02.20.csv")

#create venn diagram 
venny <- venn.diagram(list(highDepth = genesH, T.Gruber = genesGuber), "CommonDEGs_TARGET_withGruber_2017.02.20.tiff", fill = c("cornflowerblue", "darkorchid1"), 
                      print.mode = c("raw","percent"), force.unique = FALSE, cat.cex = 0.8)


#low depth
genesL
length(genesL)

genesGuber <- Gruber$Gene
length(genesGuber)


lowDepth_Gruber <- merge(lowDepth, Gruber, by.x = "Row.names", by.y = "Gene")


head(lowDepth_Gruber)
dim(lowDepth_Gruber)

numBothUp <- nrow(lowDepth_Gruber[which(lowDepth_Gruber$logFC >1 & lowDepth_Gruber$log2.fold.change. > 1), ])

numBothdown <- nrow(lowDepth_Gruber[which(lowDepth_Gruber$logFC <1 & lowDepth_Gruber$log2.fold.change. < 1), ])

opposite <- nrow(lowDepth_Gruber[which(lowDepth_Gruber$logFC > 1 & lowDepth_Gruber$log2.fold.change. < 1), ])


#save the file
write.csv(lowDepth_Gruber, file="lowDepth_TARGET_withGruber_2012.02.20.csv")


#create a venn Diagram. 
venny <- venn.diagram(list(lowDepth = genesL, T.Gruber = genesGuber), "CommonDEGs_lowDepth_TARGET_withGruber_2017.02.20.tiff", fill = c("cornflowerblue", "darkorchid1"), 
                      print.mode = c("raw","percent"), force.unique = FALSE, cat.cex = 0.8)

head(lowDepth_Gruber)
dim(lowDepth_Gruber)



########### Heatmap construction ##############


#install the packages gplots, RColorBrewer, and dendextend
# biocLite("gplots")
# biocLite("RColorBrewer")
# biocLite("dendextend")

#load the libraries
library("gplots")
library("RColorBrewer")
library("dendextend")


#### high Depth 

setwd("H:/RNA_seq_Analysis/2017.02.15_CBF-GLIS_DEG/ExpressionData/")

#read in the tsv file. 
Counts <- read.csv("Dx_rawcounts_FilteredLowCounts_TARGET_AML_Aug2014.csv")

head(Counts)
names <- Counts$ensemblSymbol
Counts <- Counts[,3:ncol(Counts)]

Dx_counts <- subset(Counts, select= grepl("T.+0[3,9]A.+", names(Counts)))


dim(Dx_counts) #58450



#Convert the target barcode into only the target USI
Dx_ID <- gsub("T.+([A-Z]+{6}).+", "\\1", colnames(Dx_counts))
Dx_ID <- gsub("\\.", "-", Dx_ID)
Dx_ID

colnames(Dx_counts) <- Dx_ID
rownames(Dx_counts) <- names


#ROWSUMS HERE to remove very low count samples
Dx_counts <- Dx_counts[rowSums(Dx_counts) > 10, ]
dim(Dx_counts) #42020 by 160 


#sort the count data by CBF-GLIS groups
Dx_counts <- Dx_counts[, c(posH,negH)]

#vector of phenotype information
posH
negH

pheno <- c(rep("pos", length(posH)), rep("neg", length(negH)))
length(pheno)


#Characterize the dataset
head(Dx_counts)
dim(Dx_counts) #42020

#CPM normalization of the read counts 
dge <- DGEList(counts = Dx_counts, group = pheno)

dge$samples

dim(dge$counts) #42,020


#From EDGER manual to remove very low count samples. 
#Did not follow this becuase the DE analysis above did not, and instead filtered for genes with low counts (less than 10 in across all samples)
# keep <- rowSums(cpm(dge) >= 1) >= 2
# 
# dge <- dge[keep, ,keep.lib.sizes=FALSE]
# 
# 
# dim(dge$counts) #20,503


#TMM Normalization of Counts 
#b y finding a set of scaling factors for the library sizes that minimize the 
#log-fold changes between the samples for most genes
dge <- calcNormFactors(dge)

dge$samples
dim(dge$counts)


#Create an object to hold the log2 CPM 
#prior count adds 2 to every count to avoid taking the log of zero
#using larger count (usually 1), bc this increase shrinkage
TMMCPM <- cpm(dge, normalized.lib.sizes = TRUE)


head(TMMCPM)
dim(TMMCPM)
rownames(TMMCPM)


#Use the ~500 genes in common to both sets to subset the genes analyzed. 

common <- read.csv("commonDEGs_CBFGLIS_2017.02.20.csv", header = TRUE)
head(common)

common <- common$Gene_Symbol

head(common)
length(common)

ID_Map <- read.csv("DE_highDepth_withGeneSymbol_2017.02.20.csv", header = TRUE)

dim(ID_Map)
head(ID_Map)

ID_Map <- ID_Map[,c(1,3)]

head(ID_Map)


#subset ID_MAP only those 513 genes in common to both low and high depth
ID_Map <- ID_Map[ID_Map$external_gene_name %in% common, ]


#subset the logCPM data frame
TMMCPM <- subset(TMMCPM, rownames(TMMCPM) %in% ID_Map$Row.names)

head(TMMCPM)
dim(TMMCPM)


names_TMMCPM <- NULL


#vector of official gene names. Used match() to ensure they are in the same order as 
#the CPM dataframe. 
names_TMMCPM <- ID_Map[match(rownames(TMMCPM), ID_Map$Row.names), 2]

head(names_TMMCPM)

#Add a count 0.01 to avoid taking the log of zero
TMMCPM <- TMMCPM + 0.01

head(TMMCPM)

#log2 transform
TMMCPM <- as.data.frame(log2(TMMCPM))

#transpose for row-wise scaling. Scale to center and range the CPM counts
TMMCPM.n <- scale(t(TMMCPM))

head(TMMCPM.n[1:10, 1:10])


#put back in original orientation
TMMCPM.tn <- t(TMMCPM.n)

head(TMMCPM.tn)

#Calculate Euclidean Distances between samples
highDepth.d1 <- dist(TMMCPM.n, method = "euclidean", diag = FALSE,
                     upper = FALSE)


#calculate euclidean distances between genes
highDepth.d2 <- dist(TMMCPM.tn, method = "euclidean", diag = FALSE,
                     upper = TRUE)

#cluster samples based on ward linkage clustering 
highDepth.c1 <- hclust(highDepth.d1, method = "ward.D2", members = NULL)

#Cluster genes based on ward linkage clustering
highDepth.c2 <- hclust(highDepth.d2, method = "ward.D2", members = NULL)


head(highDepth.c2)



#plot the dendrograms. 
par(mfrow=c(3,1),cex=0.5) 
plot(lowDepth.c1) 
plot(rotate(as.dendrogram(highDepth.c1), order = c(160:1)))
plot(highDepth.c2)

highDepth.c2$labels
highDepth.c1$labels


##try to color the dendrograms
#http://stackoverflow.com/questions/18802519/label-and-color-leaf-dendrogram-in-r
# install.packages('dendextend')
library(dendextend)

head(pheno)
length(pheno)

colorCodes <- c(pos="red",neg="dark blue")
dend <- as.dendrogram(highDepth.c1)
labels_colors(dend) <- colorCodes[pheno][order.dendrogram(dend)]

minimized <- cut(dend, h=125)

# mar=c(bottom,left,top, right)

par(mfrow=c(2,1), cex=0.3, mar=c(6, 6, 2, 2), pty="m")
plot(dend, hang= -1,
     xlab = " ", ylab=" ",  main=" ", 
     axes=TRUE, cex.axis=1.5,
     type = "rectangle", 
     horiz = FALSE)
par(cex=0.8, cex.main = 1, cex.lab = 0.85)
title(xlab="Patient Sample", ylab="Euclidean Distance",
      main="Unsupervised Clustering of AML Samples")

par(cex=0.85, mar=c(5, 6, 2, 2), pty="m")
plot(minimized$lower[[2]], hang= -1,
     xlab = " ", ylab=" ",  main=" ", 
     axes=TRUE, cex.axis=1,
     type = "rectangle", 
     horiz = FALSE)
par(cex=0.8, cex.main = 1, cex.lab=0.85)
title(xlab=" ", ylab="Euclidean Distance",
      main="")



#print the order of the samples in the dendrogram 
highDepth.c1$order

colorPal <- colorRampPalette(c("green", "yellow", "red"))(n=299)

#plot with heatmap.2
par(cex.main=1, cex=0.75, font=2, font.axis=1, lend=1) 
heatmap.2(TMMCPM.tn, 
          Colv=rotate(as.dendrogram(highDepth.c1), 
                      order = c(160:1)), 
          Rowv=as.dendrogram(highDepth.c2), 
          labRow="",
          labCol = "",
          ColSideColors = c(rep("red", length(posH)), rep("dark blue", length(negH))),
          density.info="none",
          trace="none",
          scale="none",
          col = colorPal,
          cexRow=0.3,
          margins=c(2,10), 
          lwid=c(.8,3), lhei=c(.8,3), srtCol=75, 
          adjCol=c(1,1),
          keysize=1.0, 
          key.par = list(cex=0.85),
          main="Unsupervised Clustering of AML Patient Samples")


########## low Depth 

# setwd("/Volumes/meshinchi_s/RNAseq_BCCA28Apr2016/gene_coverage_GSC-1367")
# getwd()

setwd("Desktop/2017.02_MeshinchiLab/2017.02.15_CBF-GLIS_DEG/")

correctedExp <- get(load('combatAdjustedExp.RData'))

head(correctedExp)
class(correctedExp)
dim(correctedExp)

#must have been filtered for read counts - 30,048 genes. 
dimnames(correctedExp) 

#any values less than zero (negative) must be made to be zero counts. 
correctedExp[which(correctedExp < 0, arr.ind=TRUE)] <- 0
dim(correctedExp)

#ROWSUMS HERE to remove very low count samples
correctedExp <- correctedExp[ rowSums(correctedExp) > 10, ] #30,048
dim(correctedExp)

#sort the count data by CBF-GLIS groups
correctedExp <- correctedExp[, c(posL,negL)]

head(correctedExp)
dim(correctedExp)
posH

#vector of phenotype information
posL
negL

pheno <- c(rep("pos", length(posL)), rep("neg", length(negL)))
length(pheno)

#CPM normalization of the read counts 
dge <- DGEList(counts = correctedExp, group = pheno)

dge$samples

dim(dge$counts) #30,048

#TMM Normalization of Counts 
#b y finding a set of scaling factors for the library sizes that minimize the 
#log-fold changes between the samples for most genes
dge <- calcNormFactors(dge)

dge$samples
dim(dge$counts)


#Create an object to hold the log2 CPM 
#prior count adds 2 to every count to avoid taking the log of zero
#using larger count (usually 1), bc this increase shrinkage
TMMCPM <- cpm(dge, normalized.lib.sizes = TRUE)


head(TMMCPM)
dim(TMMCPM)
rownames(TMMCPM)


#Use the ~500 genes in common to both sets to subset the genes analyzed. 

common <- read.csv("commonDEGs_CBFGLIS_2017.02.20.csv", header = TRUE)
head(common)

common <- common$Gene_Symbol

head(common)
length(common)

#subset the logCPM data frame
TMMCPM <- subset(TMMCPM, rownames(TMMCPM) %in% common)

head(TMMCPM[,1:10])

dim(TMMCPM)

names_TMMCPM <- NULL

#vector of official gene names.
names_TMMCPM <- rownames(TMMCPM)

head(names_TMMCPM)

#Add a count 0.01 to avoid taking the log of zero
TMMCPM <- TMMCPM + 0.01

head(TMMCPM)

#log2 transform
TMMCPM <- as.data.frame(log2(TMMCPM))

#transpose for row-wise scaling. Scale to center and range the CPM counts
TMMCPM.n <- scale(t(TMMCPM))

head(TMMCPM.n[1:10, 1:10])


#put back in original orientation
TMMCPM.tn <- t(TMMCPM.n)

head(TMMCPM.tn)

#Calculate Euclidean Distances between samples
lowDepth.d1 <- dist(TMMCPM.n, method = "euclidean", diag = FALSE,
                     upper = FALSE)


#calculate euclidean distances between genes
lowDepth.d2 <- dist(TMMCPM.tn, method = "euclidean", diag = FALSE,
                     upper = TRUE)

#cluster samples based on ward linkage clustering 
lowDepth.c1 <- hclust(lowDepth.d1, method = "ward.D2", members = NULL)

#Cluster genes based on ward linkage clustering
lowDepth.c2 <- hclust(lowDepth.d2, method = "ward.D2", members = NULL)

head(lowDepth.c2)

#plot the dendrograms. 
par(mfrow=c(3,1),cex=0.5) 
plot(lowDepth.c1) 
plot(rotate(as.dendrogram(lowDepth.c1), order = c(160:1)))
plot(lowDepth.c2)

lowDepth.c2$labels
lowDepth.c1$labels


##try to color the dendrograms
#http://stackoverflow.com/questions/18802519/label-and-color-leaf-dendrogram-in-r
install.packages('dendextend')
library(dendextend)

head(pheno)
length(pheno)

colorCodes <- c(pos="red",neg="dark blue")
dend <- as.dendrogram(lowDepth.c1)
labels_colors(dend) <- colorCodes[pheno][order.dendrogram(dend)]

minimized <- cut(dend, h=100)

# mar=c(bottom,left,top, right)
# pdf("lowDepth_dendrograms.pdf", 7, 6)
par(mfrow=c(2,1), cex=0.3, mar=c(6, 7, 8, 2), pty="m")  #mpg=c(3,4,0)
plot(dend,
     xlab = " ", ylab=" ",  main=" ", 
     axes=TRUE, cex.axis=1.5,
     type = "rectangle", 
     horiz = FALSE)
par(cex=0.8, cex.main = 1, cex.lab = 0.85)
# axis(1, line = 3)
title(xlab="Patient Sample", ylab="Euclidean Distance",
      main="Unsupervised Clustering of AML Samples", line = 4.5)
par(cex=0.85, mar=c(8, 6, 2, 2), pty="m")
plot(minimized$lower[[1]],
     xlab = " ", ylab=" ",  main=" ", 
     axes=TRUE, cex.axis=1,
     type = "rectangle", 
     horiz = FALSE)
par(cex=0.8, cex.main = 1, cex.lab=0.85)
title(xlab="Patient Sample", ylab="Euclidean Distance",
      main="", line = 4.5)
# dev.off()


#print the order of the samples in the dendrogram 
lowDepth.c1$order

colorPal <- colorRampPalette(c("green", "yellow", "red"))(n=299)

#plot with heatmap.2
# pdf("LowDepth_heatmap.pdf", height=10, width=10, bg="white")
png(filename = "lowDepth_heatmap.png", 800, 900)
par(cex.main=1, cex=0.75, font=2, font.axis=1, lend=1,  pty="m") #mar=c(5.1,4.1,4.1,2.1),
heatmap.2(TMMCPM.tn, 
          Colv=rotate(as.dendrogram(lowDepth.c1), 
                      order = c(465:1)), 
          Rowv=as.dendrogram(lowDepth.c2), 
          labRow="",
          labCol = "",
          ColSideColors = c(rep("red", length(posL)), rep("dark blue", length(negL))),
          density.info="none",
          trace="none",
          scale="none",
          col = colorPal,
          cexRow=0.3,
          margins=c(1,10), 
          lwid=c(.8,3), lhei=c(.8,3), srtCol=75, 
          adjCol=c(1,1),
          keysize=0.5, 
          key.par = list(cex=0.75),
          main="Unsupervised Clustering of AML Patient Samples")
dev.off()



######### GSEA #################

Jenny Smith 

#February 16th, 2017

#purpose: Convert TPM files into GSEA broad institute format for CBF-GLIS patients GSEA. 


setwd("H:/RNA_seq_Analysis/2017.02.15_CBF-GLIS_DEG/")


####### ID conversion #################
#read in the TPM only file. 
df1 <- read.csv(file = "TPM_TARGET_AML.csv", header = TRUE,
                sep=",")



#prepare the TPM dataframe
names <- df1$ensemblSymbol 
remLowCounts <- df1[,2:ncol(df1)] #keep gene expression in TPM
rownames(remLowCounts) <- names #Make the ensemble IDs to be rownames.



dim(remLowCounts)
head(remLowCounts[,1:10])


#remove any samples with a TPM of less than 1. Not stringent. Essentially only removes genes with zeros across all samples. 
remLowCounts <- remLowCounts[rowSums(remLowCounts) > 1, ]

dim(remLowCounts)
head(remLowCounts[,1:10])


#different, but correct ways to download a dataset to use. 
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#Ensembl IDs to convert
genes <- rownames(remLowCounts)


#use getBM to query the dataset/base(?) for the associated official gene symbol. 
gene_list <- getBM(attributes= c('ensembl_gene_id', 'external_gene_name'), 
                   filters = 'ensembl_gene_id', 
                   values = genes, 
                   mart = mart)


#check out the genelist
head(gene_list)
dim(gene_list)
tail(gene_list)

#remove an one to multiple ensemble ID mappings
uniqueGeneList <- gene_list[which(!duplicated(gene_list$ensembl_gene_id)), ]

dim(uniqueGeneList)

#Look at the column names - will be used for merging the dataframes. 
colnames(remLowCounts)
colnames(gene_list)

#merge the TPM data frame with the gene_list data frame using column names. 
df2 <- merge(remLowCounts,uniqueGeneList , by.x = 0,
             by.y = "ensembl_gene_id", all.x = TRUE)


#check that the merge worked. 
colnames(df2)
dim(df2)
head(df2)


#rearrange df2 
df3 <- df2[ ,c(1,209,2:208)]

head(df3[,1:10])


#save the dataframe.
write.csv(df3, "TPM_withGeneSymbol_FilteredLowCounts_TARGET_AML.csv", quote = FALSE, 
          row.names = FALSE)

###################### Convert to GSEA format #############

#read in the TPM only file. 
df1 <- read.csv(file = "TPM_withGeneSymbol_FilteredLowCounts_TARGET_AML.csv", header = TRUE,
                sep=",")

head(df1)
dim(df1)

df1 <- df1[,2:ncol(df1)]



#Prepare for GSEA analysis. Remove genes whose only official gene symbol was "NA"
df2 <- df2[which(!is.na(df2$external_gene_name)), ]

dim(df2)

#Subset for only diagnostic samples
Dx_GSEA <- subset(df2, select= grepl("T.+0[3,9]A.+", names(df2)))

head(Dx_GSEA)
dim(Dx_GSEA)

#Convert the target barcode into only the target USI
Dx_ID <- gsub("T.+([A-Z]+{6}).+", "\\1", colnames(Dx_GSEA))
Dx_ID <- gsub("\\.", "-", Dx_ID)

#Make column names to be diagnostic IDs (patient USI)
colnames(Dx_GSEA) <- Dx_ID

#Create a Name column for the gene symbols and a description column which can be NAs. 
#http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GCT:_Gene_Cluster_Text_file_format_.28.2A.gct.29
Dx_GSEA$NAME <- df2$external_gene_name
Dx_GSEA$Description <-  c(rep("na",51942))

head(Dx_GSEA)


Dx_GSEA <- Dx_GSEA[,c(161,162,1:160)]

head(Dx_GSEA[,1:10])

#Still must convert this to *.gct - a tab delimited file using excel
#must improve this. 
write.csv(Dx_GSEA, "TPM_GSEAfmt_TARGET_AML.csv")

##################GSEA####################


### High Depth #####
#save the dataframe.
df3 <- read.table("TPM_withGeneSymbol_TARGET_AML.csv", header = TRUE, sep=",")

head(df3)
dim(df3)

#remove unidentified genes
df3 <- df3[which(!is.na(df3$external_gene_name)), ]

dim(df3)
head(df3[,1:10])


#Ensure that the dataset is filtered for genes with low counts. 
#very low stringency right now. Only greater than 1 TPM for all samples
df4 <- df3[rowSums(df3[,3:ncol(df3)]) > 1, ]

head(df4)
dim(df4)

Dx_GSEA <- subset(df4, select= grepl("T.+0[3,9]A.+", names(df4)))

head(Dx_GSEA)
dim(Dx_GSEA)


#Convert the target barcode into only the target USI
Dx_ID <- gsub("T.+([A-Z]+{6}).+", "\\1", colnames(Dx_GSEA))
Dx_ID <- gsub("\\.", "-", Dx_ID)
Dx_ID

#update the column names
colnames(Dx_GSEA) <- Dx_ID


posH
negH

#sort the data frame by CBF-GLIS pos or negative
Dx_GSEA <- Dx_GSEA[, c(posH,negH)]

head(Dx_GSEA)

#Create a name column and desription column. 
#update the column names ot USI only
Dx_GSEA$NAME <- df4$external_gene_name
Dx_GSEA$Description <-  c(rep("na",34735))

dim(Dx_GSEA)

#Reorder the dataframe
Dx_GSEA <- Dx_GSEA[,c(161,162,1:160)]

head(Dx_GSEA[,1:10])

#create a phenotype list for GSEA program. 
cat("160 2 1\n# POS NEG\n", pheno, file = "GSEAfmt_CBFGLIS_phenotype.cls")

#save the what will become the *.gct file. 
write.table(Dx_GSEA, "TPM_GSEAfmt_FilteredLowCounts_TARGET_AML.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)


###low Depth

head(correctedExp)
dim(correctedExp)

#list of positive and negative BCFGLIS Samples
posL
negL

gsea_low <- NULL

dim(correctedExp[, c(posL, negL)])

names <- rownames(correctedExp)

gsea_low <- as.data.frame(correctedExp[, c(posL, negL)], row.names = FALSE)

head(gsea_low)
dim(gsea_low)

gsea_low$NAME <- names 
gsea_low$Description <-  c(rep("na",30048))

head(gsea_low)
dim(gsea_low)

#rearrange the dataframe
gsea_low <- gsea_low[, c(466, 467, 1:465)]

head(gsea_low[,1:20])
dim(gsea_low)

#vector of phenotypes 
pheno <- c(rep("pos", length(posL)), rep("neg", length(negL)))
length(pheno)

#create a phenotype list for GSEA program. 
cat("465 2 1\n# POS NEG\n", pheno, file = "batchCorrectedExpression_GSEAfmt_CBFGLIS_phenotype.cls")

getwd()

#save the what will become the *.gct file. 
write.table(gsea_low, "batchCorrectedExpression_GSEAfmt_TARGET_AML.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

