#Jenny Smith
#Dec 4, 2019

library(dplyr)
library(readr)
library(stringr)
library(magrittr)


#### Request

#18 Nov 2019
# Building networks with the ARACNe/VIPER approach requires to have a minimum of between 100 and 200 data points. So, if OK with you, maybe we can try to aim for a total of 150 to 200 patient RNAseq data points to be on the safe side.
# We also need some diversity in the profiles to be able to compute the networks, so we should include patients with different genetic alterations. 
# Here are the patients data that would be useful (the indicated number is the number found in your Nature Genetics or CCR papers…maybe they do not reflect the exact number of available RNAseq data but this is just to have an idea)
# -CBFA2T3-GLIS2 patients (both AMKL and other AML): 37 samples
# -RBM15-MKL1 patients: 10 samples
# -NUP98-KDM5A patients: 17 samples
# Then we would compared them to about 90 to 140 other pediatric AML that could include the following:
#   -RUNX1-RUNX1T1 patients: maybe 30-40 samples if possible
# -MLL/KMT2A fusions patients: maybe 30-40 samples if possible
# -40 to 60 other patients: that could include samples taken randomly or those you believe would be important to further explore functionally. 
#If we find anything that pops up on those, we would, of course, let you know along the way. 
#Of particular interest as it relates to CBFA2T3-GLIS2 are those with an ERG alteration.
# To see if there is some correlation between transcription factor activity and karyotype, 
#it would be great to get the karyotypes of these patients if possible.



#### October 19, 2020
# Dear Jenny, Dear Soheil,
# 
# I hope that things are going well for you and your families during this troubled period (and I really hope to be able to meet in person again soon).
# 
# In the meantime, here are some results that appeared during the analyses done by Elie (cc’ed). We are quite impacted by the COVID situation in our work organization so could not do everything we wanted yet.
# In brief, two groups of CBFA2T3-GLIS2 patients (we call it ETO2-GLIS2 here) appeared upon clustering. 
# I thought it would be nice to see whether the two groups are correlated with:
#   -different phenotypes of the blasts ? ...as reported by Masetti/Locatelli, Blood 2013
# -different blast % ? In fact, the xCell analysis strongly suggests that it may in fact just be a matter of leukemic cell purity and cross-contamination with other blood cells.
# -or other parameters (e.g. sex, age at diagnosis, karyotype) ?
#   -maybe most importantly, different survival ?
#   
#   We have attached the excel file with the classification of the ETO2-GLIS2 patients in these two groups.
# Would you be willing to share these additional information to see if there are any correlations?
#   
#   Elie and I would be very happy to further discuss if you need additional information.


#### Counts/TPMS
 

TPMs <- read_csv(file.path(HOME,
                           "0000.00.03_Expression_Matrices/TARGET_AML_0531_1031_Ribodepleted_RNAseq_dupGenesRemoved_TPM.csv"))
head(TPMs[,1:5])
dim(TPMs[,])

cts <- read_csv(file.path(HOME,
                          "0000.00.03_Expression_Matrices/TARGET_AML_0531_1031_Ribodepleted_RNAseq_dupGenesRemoved_FractionalCounts.csv"))
head(cts[,1:5])
dim(cts)


##### Subset CDEs
USIs <- str_split_fixed(colnames(TPMs),pattern = "\\.", n=5)[,3]
CDEs <- readRDS(file.path(CDE,"Merged/TARGET_AML_0531_1031_merged_CDEs_9.4.19.RDS"))

CDEs2 <- CDEs %>% 
  filter(!is.na(USI)) %>% 
  filter(USI %in% USIs) %>%
  mutate(RBM15.MKL1=case_when(
    grepl("RBM15-MKL1", Primary.Fusion.CNV) ~ "Yes", 
    ScreenedForFusion == "No" ~ "Unknown", 
    TRUE ~ "No"),
    KMT2A=case_when(
      grepl("KMT2A", Primary.Fusion.CNV) | grepl("KMT2A", Additional.Fusions.CNV) ~ "Yes", 
      ScreenedForFusion == "No" ~ "Unknown", 
      TRUE ~ "No")) %>%
  mutate(Groups=case_when(
    CBFA2T3.GLIS2 == "Yes" ~ "CBFA2T3.GLIS2",
    NUP98.KDM5A== "Yes" ~ "KDM5A",
    RBM15.MKL1== "Yes" ~ "MKL1",
    KMT2A== "Yes" ~"KMT2A",
    grepl("RUNX1-RUNX1T1",Primary.Fusion.CNV) ~"t.8.21",
    grepl("ERG",Primary.Fusion.CNV) ~ "ERG",
    TRUE ~ "Others")) %>%
  
  group_by(Groups) %>%
  mutate(N=n()) %>%
  sample_n(size=ifelse(unique(N)<40, unique(N),40),replace=FALSE) %>%
  ungroup() %>%
  
  select(USI,Protocol,Age.in.years,ISCN, WBC..x10.3.MicroLiter..levels,
         Bone.marrow.leukemic.blast.percentage....,Peripheral.blasts....,
         FAB_or_WHO.Classification,
         CBFA2T3.GLIS2,NUP98.KDM5A,RBM15.MKL1,KMT2A,
         Primary.Fusion.CNV,Additional.Fusions.CNV,
         NPM.mutation., CEBPA.mutation.,WT1.mutation.,
         FLT3.ITD.positive.,FLT3.ITD.allelic.ratio) %>%
  as.data.frame() %>%
  set_rownames(.$USI)

head(CDEs2[,1:5])
dim(CDEs2) 

table(CDEs2$CBFA2T3.GLIS2) #38
table(CDEs2$NUP98.KDM5A) #30
table(CDEs2$RBM15.MKL1) #17
table(CDEs2$KMT2A) # 40
table(CDEs2$Protocol)

# CBFA2T3.GLIS2 ERG         KDM5A         KMT2A 
# 38             19            30            40 
# MKL1        Others        t.8.21 
# 10            40            40 

##Subset Expn Data for 1031 only
subset_data <- function(CDE,ExpnData,NBM=FALSE){
  #https://stackoverflow.com/questions/30604107/r-conditional-evaluation-when-using-the-pipe-operator
  if(NBM){
    rm.samples <- "Kas|MV4|MPN[0-9]|Sort|repl"
  }else{
    rm.samples <- "Kas|MV4|MPN[0-9]|Sort|repl|BM[0-9]|RO[0-9]"
  }
  print(c("Selected NBM?", NBM))
  
  rm.samples <- grep(rm.samples,colnames(ExpnData), ignore.case = TRUE)
  if(length(rm.samples)>0){ExpnData <- ExpnData[,-rm.samples]}
  if(!any(grepl("[A-Za-z]", rownames(ExpnData)))){colnames(ExpnData)[1] <- "Gene_Name"}
  
  if(any(grepl("\\.", colnames(ExpnData)))){
    USI <- str_split_fixed(colnames(ExpnData),pattern = "\\.", n=5)[,3]
  }else{
    USI <- colnames(ExpnData)
  }
  
  colNames <- c("Gene_Name",intersect(USI, CDE$USI))%>% 
    {if(NBM) c(., grep("BM[0-9]|RO[0-9]",colnames(ExpnData), value=TRUE)) else .} %>%
    paste(.,  collapse = "|")
  
  CDE <- CDE[grep(colNames, CDE$USI, value = TRUE),]
  ExpnData <- ExpnData[,grep(colNames, colnames(ExpnData), value = T)]
  
  return(list(ExpnData, CDE))
}


cts.subset <-  subset_data(CDE = CDEs2, ExpnData = cts)
dim(cts.subset[[1]])  #51573   218
dim(cts.subset[[2]]) #217  19
head(cts.subset[[1]][,1:5])
table(grepl("rep",colnames(cts.subset[[1]]))) 

table(cts.subset[[2]]$CBFA2T3.GLIS2) #38 positive
table(cts.subset[[2]]$NUP98.KDM5A) #30
table(cts.subset[[2]]$RBM15.MKL1) #10
table(cts.subset[[2]]$KMT2A, useNA = 'ifany') #251
table(cts.subset[[2]]$Protocol)


tpms.subset <- subset_data(CDE = CDEs2,ExpnData = TPMs, NBM=FALSE)
dim(tpms.subset[[1]]) 
dim(tpms.subset[[2]])
head(tpms.subset[[1]][,1:5])
table(grepl("rep",colnames(tpms.subset[[1]])))

table(tpms.subset[[2]]$CBFA2T3.GLIS2) #38 positive
table(tpms.subset[[2]]$NUP98.KDM5A) #30
table(tpms.subset[[2]]$RBM15.MKL1) #10
table(tpms.subset[[2]]$Protocol)



### Save the Data 

# write.csv(tpms.subset[[1]],
#           file.path(SCRATCH,"jlsmith3/TARGET_AML_Ribodepleted_RNAseq_dupGenesRemoved_TPMs.csv"),
#           row.names = FALSE)
# write.csv(cts.subset[[1]],
#           file.path(SCRATCH,"jlsmith3/TARGET_AML_Ribodepleted_RNAseq_dupGenesRemoved_FractionalCounts.csv"),
#           row.names = FALSE)
# 
# write.csv(cts.subset[[2]],
#           file.path(SCRATCH,"jlsmith3/TARGET_AML_Clinical_Data.csv"),
#           row.names = FALSE)

table(cts.subset[[2]]$Primary.Fusion.CNV) %>% 
  as.data.frame() %>% 
  arrange(desc(Freq)) %>% 
  head()

#            Var1 Freq
# 1          None  323
# 2 RUNX1-RUNX1T1  146
# 3    CBFB-MYH11  104
# 4   KMT2A-MLLT3   89
# 5  KMT2A-MLLT10   54
# 6    NUP98-NSD1   54


### October 20, 2020

setwd(file.path(PROJHOME,"2017.02.15_CBF-GLIS_DEG/2018.03.21_CBF-GLIS_DEGs_Comprehensive"))

CDEs <- readRDS(file.path(CDE,"Merged/TARGET_AML_0531_1031_merged_CDEs_9.18.20.RDS"))
CDEs <- CDEs %>% 
  filter(USI != "Unknown", !is.na(USI))

head(CDEs[,1:5])
dim(CDEs)

manifest <- read.csv(file.path(TARGET, "SequencingDataMatrix/TARGET_AML_Ribodepleted_Manifest_10.08.20.csv"),
                     row.names = 1) 

dim(manifest)
head(manifest)

toFill <- xlsx::read.xlsx("Thomas_Mercher/metadataSeattle_withEGclusters.xls",
                          sheetIndex = 1)

head(toFill)
dim(toFill)

Filled <- toFill %>% 
  select(Sample=NA.,everything(),-Sex) %>%
  left_join(., select(manifest,Sample,USI,Batch,Time_point,Tissue),
            by=c("Sample")) %>%
  left_join(.,select(CDEs, USI, Age.in.years,Sex,
                     ISCN,
                     Blast.percent..by.flow.,
                     Bone.marrow.leukemic.blast.percentage....,
                     Peripheral.blasts....,
                     CR.status.at.end.of.course.1,
                     CR.status.at.end.of.course.2,
                     OS.event.ID,OS.time..days.,
                     EFS.event.type.ID,EFS.time..days.),
            by="USI") %>%
  select_if(~!all(is.na(.))) 

head(Filled)
dim(Filled)
# 
# write.csv(Filled,"Thomas_Mercher/TARGET_AML_metadataSeattle_withEGclusters.csv",
#           row.names = FALSE)


