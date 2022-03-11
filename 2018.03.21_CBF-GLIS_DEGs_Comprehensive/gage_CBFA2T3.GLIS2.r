library(stringr)
library(magrittr)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(gage)
library(gageData)
library(methods)

setwd('~/RNA_seq_Analysis/2018.03.21_CBF-GLIS_DEGs_Comprehensive/')
source("~/scripts/RNAseq_Analysis/DifferentialExpn_PathwayAnalysis/GAGE_GSEA_Function.r")


CBFGLIS <- readRDS("DEGs/TARGET_AML_AllCohorts_CBFGLISvsOtherAML_DEGs_10.23.18.RDS")
filename <- "TARGET_AML_hd.1031_CBFGLIS_vs_OtherAML"
C2.KEGG <- readRDS("~/RNA_seq_Analysis/0000.00.01_GSEA_geneSets_gmt/c2.cp.kegg.v6.0.symbols.RDS")


print("starting1")

GSA <- lapply(CBFGLIS, gage_from_pipeline, type="expn",geneset=C2.KEGG)
  save(GSA,file=paste0(filename, "_expn_C2.KEGG_11.23.18.RData"))
  rm(GSA)
  gc()

print("done1")

print("starting2")

GSA.FC <- lapply(CBFGLIS, gage_from_pipeline, type="FC",geneset=C2.KEGG)
  save(GSA.FC,file=paste0(filename, "_FC_C2.KEGG_11.23.18.RData"))
  rm(GSA.FC)
  gc()


rm(C2.KEGG)
print("done2")

print("starting3")


C2.All <- readRDS("~/RNA_seq_Analysis/0000.00.01_GSEA_geneSets_gmt/c2.all.v6.0.symbols.RDS")

GSA.C2.All <- lapply(CBFGLIS, gage_from_pipeline, type="expn",geneset=C2.All)
save(GSA.C2.All,file=paste0(filename, "_expn_C2.All_11.23.18.RData"))
rm(GSA.C2.All)
gc()

rm(GSA.C2.All)
print("done3")


print("starting4")

GSA.KEGG <- lapply(CBFGLIS, 
                   gage_from_pipeline,
                   type="expn",
                   geneset=NULL, 
                   include.Disease=TRUE)
  save(GSA.KEGG,file=paste0(filename, "_expn_HSA.KEGG_w.Diseases_11.23.18.RData"))
  rm(GSA.KEGG)
  gc()

print("done4")

C5 <- readRDS("~/RNA_seq_Analysis/0000.00.01_GSEA_geneSets_gmt/c5_list_SetSize_GE.50_LE.300_v6.1.symbols.RDS")

print("starting5")
GSA.GO.BioProcess <- lapply(CBFGLIS, gage_from_pipeline, type="expn", geneset=C5[["c5.bp"]])
    save(GSA.GO.BioProcess, file=paste0(filename, "_expn_C5.BioProcess_SetSize50to300_11.23.18.RData"))
    rm(GSA.GO.BioProcess)
    gc()
print("done5")



print("starting6")
GSA.GO.CellComp <- lapply(CBFGLIS, gage_from_pipeline, type="expn", geneset=C5[["c5.cc"]])
    save(GSA.GO.CellComp, file=paste0(filename, "_expn_C5.CellComp_SetSize50to300_11.23.18.RData"))
    rm(GSA.GO.CellComp)
    gc()
print("done6")



print("starting7")
GSA.GO.MolFunc <- lapply(CBFGLIS, gage_from_pipeline, type="expn", geneset=C5[["c5.mf"]])
    save(GSA.GO.MolFunc, file=paste0(filename, "_expn_C5.MolFunc_SetSize50to300_11.23.18.RData"))
    rm(GSA.GO.MolFunc)
    gc()
print("done7")


print("starting8")
GSA.GO.BioProcess.FC <- lapply(CBFGLIS, gage_from_pipeline, type="FC", geneset=C5[["c5.bp"]])
  save(GSA.GO.BioProcess.FC, file=paste0(filename, "_FC_C5.BioProcess_SetSize50to300_11.23.18.RData"))
  rm(GSA.GO.BioProcess.FC)
  gc()
print("done8")



print("starting9")
GSA.GO.CellComp.FC <- lapply(CBFGLIS, gage_from_pipeline, type="FC", geneset=C5[["c5.cc"]])
  save(GSA.GO.CellComp.FC, file=paste0(filename, "_FC_C5.CellComp_SetSize50to300_11.23.18.RData"))
  rm(GSA.GO.CellComp.FC)
  gc()
print("done9")



print("starting10")
GSA.GO.MolFunc.FC <- lapply(CBFGLIS, gage_from_pipeline, type="FC", geneset=C5[["c5.mf"]])
  save(GSA.GO.MolFunc.FC, file=paste0(filename, "_FC_C5.MolFunc_SetSize50to300_11.23.18.RData"))
  rm(GSA.GO.MolFunc.FC)
  gc()
print("done10")


C3 <- readRDS("~/RNA_seq_Analysis/0000.00.01_GSEA_geneSets_gmt/c3.tft.v6.2.symbols.RDS")

print("starting11")
GSA.TFB <- lapply(CBFGLIS, gage_from_pipeline,
                  type="expn",
                  geneset=C3)
  save(GSA.TFB, file=paste0(filename, "_expn_C3_TFbindingSites.RData"))
  rm(GSA.TFB)
  gc()
print("done11")


print("starting12")
GSA.TFB.FC <- lapply(CBFGLIS, gage_from_pipeline,
                     type="FC",
                     geneset=C3)
  save(GSA.TFB.FC, file=paste0(filename, "_FC_C3_TFbindingSites.RData"))
  rm(GSA.TFB.FC)
  gc()
print("done12")


print("starting13")
hall <- readRDS("~/RNA_seq_Analysis/0000.00.01_GSEA_geneSets_gmt/hallmark_MdbSig/h.all.v6.2.symbols.RDS")

GSA.hallmark <- lapply(CBFGLIS,
                       gage_from_pipeline,
                       type="expn",
                       geneset=hall)
    save(GSA.hallmark, file=paste0(filename, "_expn_MSigDB_Hallmark_11.23.18.RData"))
    rm(GSA.hallmark)
    gc()

print("done13")

    