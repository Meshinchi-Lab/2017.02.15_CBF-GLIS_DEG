#Jenny Smith 
#12/27/18

#Purpose: Create a pipeline for detecting anti-correlated miRNA-mRNAs and query databases for experimentally validated or known interactions. 
setwd('~/RNA_seq_Analysis/2018.03.21_CBF-GLIS_DEGs_Comprehensive/')

corr.miRNA.mRNA <- function(miRNA.Expn, gene.Expn, phenovector,ref){
  library(ggplot2)
  library(psych)
  # miRNA.Expn is a named numeric vector of miRNA log2 RPMs or TMM normalized counts
  #gene.Expn is a named numeric vector of mRNA log2 normalized expn (CPM, TPM, RPKM, etc)
  #phenovector is named characer vector of the groups (eg pos, neg)
  #ref is a charcter vector length 1, which has the group from phenovector that was "control group"
  
  g1 <- names(which(phenovector != ref))
  g2 <- names(which(phenovector == ref))
  
  #Subset for the group of interest and ensure same order in expression sets
  RPM <- t(miRNA.Expn[,g1])
  GE <- t(gene.Expn[,g1])
  
  #combine the expression sets
  expn.g1 <- RPM %>%
    as.data.frame() %>%
    rownames_to_column("USI") %>%
    inner_join(., 
               rownames_to_column(as.data.frame(GE),"USI"), 
               by="USI")
  
  expn.g2 <- t(miRNA.Expn[,names(phenovector)]) %>%
    as.data.frame() %>%
    rownames_to_column("USI") %>%
    inner_join(., 
               rownames_to_column(as.data.frame(t(gene.Expn[,names(phenovector)])),"USI"), 
               by="USI") %>%
    mutate(Group=phenovector) %>%
    dplyr::select(USI,Group, everything())
  
  #correlation of miRNA to mRNA in group 1, CBF-GLIS in this case. 
  corr <- corr.test(RPM,GE, method="spearman",
                    adjust="BH",ci=FALSE)
  
  #Format the results 
  pearson <- corr$r %>%
    as.data.frame() %>%
    rownames_to_column("MIMAT") %>%
    gather(Gene,SpearmanRho, -MIMAT) 
  
  pvals <- corr$p %>%
    as.data.frame() %>%
    rownames_to_column("MIMAT") %>%
    gather(Gene,Adj.P.val, -MIMAT)
  
  res <- merge(pearson,pvals, by=c("MIMAT","Gene")) %>%
    separate(MIMAT,into = c("Mir","MIMAT"), sep="\\.")
  
  list <- list(expn.g1,expn.g2, res, corr)
  names(list) <- c("Merged_Expression","Merged_Expression_Ref", "Results", "correlation_list")
  
  return(list)
  
}

queryAnamir <- function(pairs, cutoff=2){
  db <- RMySQL::dbConnect(RMySQL::MySQL(), user = "visitor",
                          password = "visitor",
                          dbname = "visitor",
                          host = "anamir.cgm.ntu.edu.tw")
  
  interaction <- list()
  
  for (i in seq_len(nrow(pairs))){
    # print(i)
    mirna <- pairs[i,1]
    gene <- pairs[i,2]
    corr <- pairs[i,3]
    
    query <- paste0("SELECT `miRNA_21`, `Gene_symbol`, `Ensembl`,
                    `Gene_ID`, `DIANA_microT_CDS`, `EIMMo`, `Microcosm`,
                    `miRDB`, `miRanda`, `PITA`, `rna22`, `Targetscan`,
                    `Sum`, `miRecords`, `miRTarBase`,
                    `Validate` FROM `all_hsa` WHERE miRNA_21 like '",
                    mirna, "' AND gene_symbol like '", gene, "' ;")
    
    tmp <- DBI::dbGetQuery(db, query)
    
    #If statement to avoid adding empty query results to the interactions list
    if (nrow(tmp) == 0 && cutoff > 0) {
      next 
    } else {
      # print(colnames(tmp)) #16 columns, #13 is sum and #16 is validate
      tmp <- c(tmp, corr, 0)
      interaction[[i]] <- tmp
    }
  }
  
  # interaction <- data.frame(do.call(rbind, interaction))
  interaction <- do.call(rbind, interaction)
  
  # disconnect db
  cons <- RMySQL::dbListConnections(RMySQL::MySQL())
  for (con in cons) RMySQL::dbDisconnect(con)
  
  # add column de novo
  colnames(interaction)[ncol(interaction)] <- c("de novo")
  
  if (cutoff > 1 && nrow(interaction) > 0) {
    del_row <- c()
    for (i in seq_len(nrow(interaction))) {
      if (interaction[i, 13] < cutoff && interaction[i, 16] %in% "FALSE") {
        del_row <- c(del_row, i)
      }
    }
    interaction <- interaction[-(del_row), ]
  }
  
  return(interaction)
}


reformat.Anamir <- function(anamir.res,UpGenePairs,DEGs, DEMirs){
  require(dplyr)
  
  anamir.res <- anamir.res %>%
    apply(., 2, as.character) %>% #all columns are "list" class, so need this line
    as.data.frame(.) %>%
    
    mutate(Gene_symbol=toupper(Gene_symbol)) %>%
    merge(x = ., y = dplyr::select(UpGenePairs, 
                            miRNA_21=Mir.V.21,
                            Gene_symbol=Gene,
                            MIMAT,
                            Adj.P.val), 
          all.x=TRUE, 
          by=c("miRNA_21", "Gene_symbol")) %>%
    dplyr::select(miRNA_21,MIMAT,Gene_symbol,Ensembl, Gene_ID, Correlation,Adj.P.val, everything()) %>%
    
    left_join(., dplyr::select(DEGs,
                               Gene_symbol=gene, 
                               logFC_Gene=logFC,  
                               adj.P.Val_Gene=adj.P.Val), 
              by="Gene_symbol") %>%
    
    left_join(., dplyr::select(DEMirs.1031,
                               MIMAT, 
                               logFC_miR=logFC,  
                               adj.P.Val_miR=adj.P.Val),
              by="MIMAT") %>%
    
    mutate_all(funs(ifelse(grepl("^[0-9]|^\\-[0-9]", .), as.numeric(.), .))) %>%
    filter(Validate == TRUE | Sum >= 2) %>%
    
    dplyr::select(miRNA_21, MIMAT,logFC_miR,adj.P.Val_miR,
                  Gene_symbol,logFC_Gene ,adj.P.Val_Gene,Ensembl,
                  Gene_ID,Correlation,Adj.P.val_Corr=Adj.P.val, everything())

}



reformat_ints <- function(multimir.res, gene.miR.corrs, DEGs, DEMirs, type="summary"){
  
  if(type=="summary"){
    reformatted <- multimir.res@summary %>%
      filter(all.sum >= 2 | validated.sum > 0) #only have
    
  }else if (type=="data"){
    reformatted <- Up.Gene.Ints@data
  }
  
  reformatted <- reformatted %>%
    inner_join(., dplyr::select(gene.miR.corrs,
                                Mir.V.21,
                                Gene, 
                                SpearmanRho, 
                                adj.P.Val_Corr=Adj.P.val),
               by=c("mature_mirna_id"="Mir.V.21",
                    "target_symbol"="Gene")) %>%
    left_join(., dplyr::select(DEGs,
                               gene,
                               logFC_Gene=logFC,
                               adj.P.Val_Gene=adj.P.Val,
                               contains("Compartment")),
              by=c("target_symbol"= "gene")) %>%
    left_join(., dplyr::select(DEMirs,
                               miR,
                               logFC_miRNA=logFC,
                               adj.P.Val_miRNA=adj.P.Val,
                               gene),
              by=c("mature_mirna_id"= "miR")) %>%
    dplyr::select(mature_mirna_acc, mature_mirna_id,logFC_miRNA,adj.P.Val_miRNA,
                  target_symbol,logFC_Gene,adj.P.Val_Gene,
                  SpearmanRho, adj.P.Val_Corr,
                  everything())
  
  return(reformatted)
  
}


extract_core_genes <- function(gage_from_pipeline.res,direction="up", ret.genesets=FALSE){
  
  gage.res.sig <- gage_from_pipeline.res[sapply(gage_from_pipeline.res,length) == 5]
  
  
  if(direction=="up"){
    idx <- "essSets.Up"
    idx.num <- 4
  }else if (direction=="down"){
    idx <- "essSets.Dn"
    idx.num <- 5
  }
  
  #Remove extraneous gene-sets from the MSigDB C2.All set
  c2.all <- grep("OtherAML_expn_C2.All", names(gage.res.sig))
  i <- sapply(names(gage.res.sig[[c2.all]][[idx]][["coreGeneSets"]]), 
              function(x) grep("^PID_|^REACTOME|^BIOCARTA", x)) 
  i <- unlist(lapply(i, function(x) length(x) > 0))
  
  gage.res.sig[[c2.all]][[idx]][["coreGeneSets"]] <- gage.res.sig[[c2.all]][[idx]]$coreGeneSets[i]
  
  #select the correct index for the pathway core genesets
  core.genes <- lapply(gage.res.sig, `[[`, idx.num) 
  core.genes <- lapply(core.genes, `[[`, 7) 
  
  if(ret.genesets){
    return(core.genes)
  }
  
  #Define an empty list to append to 
  paths.core.genes <- list()
  #index position for the new lis
  n <- 1
  
  #for loop to un-nest the core.genes, so that they are not in a depth == 2 list
  for (i in 1:length(core.genes)){
    
    gs <- core.genes[[i]]
    
    if(i == 1){
      stop=length(gs)
    }else{
      stop <- n + c(length(gs)-1)
    }
    
    items <- n:stop
    paths.core.genes[items] <- gs
    names(paths.core.genes)[items] <- names(gs)
    
    n <- stop+1
  }
  
  return(paths.core.genes)
  
}



find.Pathway.Ints <- function(pathway.genes, Interaction.df, pathway.name=NULL){
  
  Ints <- unlist(lapply(pathway.genes, function(x)
    grep(paste0("^",x,"$"), Interaction.df$Gene_symbol,value = TRUE, ignore.case = TRUE)))
  
  pairs <- Interaction.df %>%
    filter(Gene_symbol %in% Ints) %>%
    mutate(Pair=paste(miRNA_21, Gene_symbol, sep=", ")) %>%
    dplyr::select(Pair)
  
  if(!is.null(pathway.name)){
    colnames(pairs) <- pathway.name
  }
  
  return(pairs)
  
}




miRNA.mRNA_Interactsions <- function(DEGs.res, DEMiRs.res,mir.IDmap=NULL,
                                     gage.res=NULL,
                                     FC.cutoff=1, corr=0, corr.FDR=0.05){
  library(tibble)
  library(dplyr)
  library(anamiR, lib.loc="/home/jlsmith3/R/x86_64-pc-linux-gnu-library/3.5")
  
  #Define Functions
  filter.FC <- function(res, cut.off){
    
    df <- res$DE$DE %>%
      rownames_to_column("gene") %>%
      mutate(gene=toupper(as.character(gene)))
    
    if(cut.off > 0){
      df <- filter(df,logFC > cut.off)}
    else{
      df <- filter(df,logFC < cut.off)
    }
  }
  
  
  reFormat.anaMir <- function(filt.cor.res){
    cor.mat <- filt.cor.res %>%
      dplyr::select(miRNA=Mir.V.21, Gene, Correlation=SpearmanRho) %>%
      as.matrix()
    return(cor.mat)
  }
  
  reformatEnrichment <- function(mat){
    df <- as.data.frame(mat, stringsAsFactors=FALSE) %>%
      set_colnames(gsub(" |\\-", "_", colnames(.))) %>%
      mutate(Raw_P_Value=as.numeric(as.character(Raw_P_Value)),
             Empirical_P_Value=as.numeric(as.character(Empirical_P_Value))) %>%
      mutate(FDR=p.adjust(Empirical_P_Value, method="BH")) %>%
      filter(FDR < 0.05) %>%
      arrange(FDR)
    
    return(df)
  }
  
  #Read in the ID map for miRNAs if necessary
  if(is.null(mir.IDmap)){
    ID.map <- read.csv("/fh/fast/meshinchi_s/workingDir/TARGET/Reference_Data/miRBase_v21/hsa_gff3_IDMap.csv")
  }
  
  
  #Define genes/mirs as up or down regulated
  Up.DE <- filter(DEGs.res, cut.off=FC.cutoff)
  Dn.DE <- filter(DEGs.res, cut.off=FC.cutoff*-1)
  Up.DEMirs <- filter(DEMiRs.res, cut.off=FC.cutoff)
  Dn.DEMirs <- filter(DEMiRs.res, cut.off=FC.cutoff*-1)
  
  #Calculate correlations
  corrs.UpGenes <- corr.miRNA.mRNA(miRNA.Expn = DEMiRs.res$DE$Voom$E[Dn.DEMirs[["gene"]],],
                                   gene.Expn = DEGs.res$DE$Voom$E[Up.DE[["gene"]],],
                                   phenovector = DEGs.res$phenovector,
                                   ref="GroupB")
  
  corrs.DnGenes <- corr.miRNA.mRNA(miRNA.Expn = DEMiRs.res$DE$Voom$E[Up.DEMirs[["gene"]],],
                                   gene.Expn = DEGs.res$DE$Voom$E[Dn.DE[["gene"]],],
                                   phenovector = DEGs.res$phenovector,
                                   ref="GroupB")
  
  #Filter for Significance
  UpGenePairs <- corrs.UpGenes$Results %>%
    left_join(., dplyr::select(ID.map, Mir.V.21=miR, everything()),
              by=c("MIMAT"="MIMAT.ID")) %>%
    mutate(Gene_symbol=toupper(Gene_symbol)) %>%
    filter(Adj.P.val < corr.FDR & SpearmanRho < corr) %>%
    arrange(Adj.P.val)
  
  DnGenePairs <- corrs.DnGenes$Results %>%
    left_join(., dplyr::select(ID.map, Mir.V.21=miR, everything()),
              by=c("MIMAT"="MIMAT.ID")) %>%
    mutate(Gene_symbol=toupper(Gene_symbol)) %>%
    filter(Adj.P.val < corr.FDR & SpearmanRho < corr) %>%
    arrange(Adj.P.val)
  
  
  #Reformat Dataframes for R Package
  up.pairs.fmt <- reFormat.anaMir(UpGenePairs)
  dn.pairs.fmt <- reFormat.anaMir(DnGenePairs)
  
  #Query Databases
  Interaction.Up <- queryAnamir(up.pairs.fmt)
  Interaction.Dn <- queryAnamir(dn.pairs.fmt)
  
  #Clean up the Results and Merge in the FC, FDR, and other Info. 
  Int.Up.Clean <- reformat.Anamir(Interaction.Up,UpGenePairs,
                                  DEGs = Up.DE, DEMirs = Dn.DEMirs)
  
  Int.Dn.Clean <- reformat.Anamir(Interaction.Dn,DnGenePairs,
                                  DEGs = Dn.DE, DEMirs = Up.DEMirs)

  #Gene Set Enrichment
  up.paths <- enrichment(data_support = Interaction.Up, org = "hsa")
  dn.paths <- enrichment(data_support = Interaction.Dn, org="hsa")
  
  #Reformat the GSEA Results
  up.paths.fmt <-  reformatEnrichment(up.paths)
  dn.paths.fmt <- reformatEnrichment(dn.paths)
  
  #create a list with the results
  res <- list("UpGenes"=c("spearman"=corrs.UpGenes,"AntiCorr_Pairs"=UpGenePairs,"Interactions"=Int.Up.Clean, "GSA"=up.paths.fmt),
              "DnGenes"=c("spearman"=corrs.DnGenes,"AntiCorr_Pairs"=DnGenePairs,"Interactions"=Int.Dn.Clean,"GSA"=dn.paths.fmt))
  
  
  if(!is.null(gage.res)){
    up.core.sets <- extract_core_genes(gage_from_pipeline.res = gage.res, direction = "up")
    dn.core.sets <- extract_core_genes(gage_from_pipeline.res = gage.res, direction = "down")
    
    UpPath.Ints<- lapply(up.core.sets, 
                         find.Pathway.Ints, 
                         Interaction.df=Int.Up.Clean)
    
    UpPath.Ints <- UpPath.Ints[sapply(UpPath.Ints, function(x) nrow(x) >0 )]
    DnPath.Ints <- lapply(dn.core.sets, 
                          find.Pathway.Ints,
                          Interaction.df=Int.Dn.Clean)
    
    DnPath.Ints <- DnPath.Ints[sapply(DnPath.Ints, function(x) nrow(x) >0 )]
    
    res[["UpPath.Ints"]] <- UpPath.Ints
    res[["DnPath.Ints"]] <- DnPath.Ints
  }
  
  return(res)
  
}









