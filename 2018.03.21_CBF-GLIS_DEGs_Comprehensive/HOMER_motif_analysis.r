#For Motif Analysis 

#Use Modules for analysis
#ml Homer/v4.9-foss-2016b

setwd("~/RNA_seq_Analysis/2018.03.21_CBF-GLIS_DEGs_Comprehensive/")
options(stringsAsFactors = FALSE)
dir.create("HOMER")

library(dplyr)



IDmap <- read.csv("~/RNA_seq_Analysis/0000.00.02_Reference_GeneInfo/GeneSymbol_EnsemblID_Conversion_GRCh37.69_FromBCCA.csv")
dim(IDmap)
head(IDmap)


#DEGs

DEGs <- read.csv("DEGs/TARGET_AML_1031_CBFGLIS_vs_OtherAML_DEGs_updateAnno.csv") %>%
  left_join(., IDmap, by=c("gene"="geneSymbol"))


dim(DEGs) #3926   12
head(DEGs)


DEGs %>%
  select(gene_id, gene, logFC,adj.P.Val) %>%
  filter(logFC > 1) %>%
  write.table(., file="HOMER/TARGET_AML_CBFGLIS_vs_OtherAML_UpReg.txt", sep = "\t", quote = FALSE)
  
DEGs %>%
  select(gene_id, gene, logFC,adj.P.Val) %>%
  filter(logFC < -1) %>%
  write.table(., file="HOMER/TARGET_AML_CBFGLIS_vs_OtherAML_DnReg.txt", sep = "\t", quote = FALSE)


#Interactions DEGs

int.up <- read.csv("miRNA-mRNA_Interactions/anamiR/TARGET_AML_CBFGLIS_vs_OtherAML_UpRegDEGs_DnRegMIR_miRNA-mRNA_interactions_.csv")
head(int.up)
dim(int.up)

int.up %>%
  select(Ensembl, everything()) %>%
  write.table(., file="HOMER/TARGET_AML_CBFGLIS_vs_OtherAML_UpDEGs_Interactions.txt", sep = "\t", quote = FALSE)
  
int.dn <- read.csv("miRNA-mRNA_Interactions/anamiR/TARGET_AML_CBFGLISvsOtherAML_DnRegDEGs_UpRegMIR_miRNA-mRNA_interactions_.csv")
head(int.dn)
dim(int.dn)


int.dn %>%
  select(Ensembl, everything()) %>%
  write.table(., file="HOMER/TARGET_AML_CBFGLIS_vs_OtherAML_DnDEGs_Interactions.txt", sep = "\t", quote = FALSE)










