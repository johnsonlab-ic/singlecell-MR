

check_colnames=function(df,colnames){
    
    df_colnames=colnames(df)
    input_colnames=colnames
    
    #check if all input colnames exist in df colnames
    missing=setdiff(input_colnames,df_colnames)
    
    if(length(missing)==0){
        message("All column names present.")
    } else {
        stop(paste0("Columns missing from GWAS: ",missing,". Please check formatting."))
    }
    
}



select_regions=function(gwas,
                        pval,
                        window,
                        parallel=FALSE,
                        local=TRUE,
                        plink_bin=genetics.binaRies::get_plink_binary(),
                        path_to_binaries="/rds/general/user/ah3918/projects/single_cell_eqtl/live/MRC/REFERENCE_DATASETS/clumping_ref/EUR/EUR"){


    error_catch_clump=function(gwas,clump_kb,
        clump_r2,clump_p,plink_bin,bfile){

        tryCatch({test<-ieugwasr::ld_clump_local(gwas,clump_kb=clump_kb,clump_r2=clump_r2,clump_p=clump_p,
        plink_bin=plink_bin,bfile=bfile);},error=function(e){test<<-NA});
        return(test)

    }
    message(paste0("local: ",local))
    options(warn=-1)

    # window is each side of the index snp. So "window" x 2 is the full window considered.
    window=window
    message(paste0(Sys.time(),": Splitting into chromosomes.."))
    gwas_chr_list<-split(gwas,gwas$chr)

    gwas_chr_list=lapply(gwas_chr_list,function(x){

    b=x[x$pval<pval,]
    if(nrow(b)==0){
        return(NULL)
    }else{
        return(x)
    }
    })

  gwas_chr_list=Filter(Negate(is.null), gwas_chr_list)

  kb_window=window/1e3
  message(paste0(Sys.time(),": Selecting lead index snps per region.."))

  if(local==TRUE){

      gwas_clump<-suppressMessages(lapply(gwas_chr_list,error_catch_clump,clump_kb=kb_window,clump_r2=0,clump_p=pval,plink_bin=plink_bin,
      bfile=path_to_binaries))
  } else{
      gwas_clump<-suppressMessages(lapply(gwas_chr_list,ieugwasr::ld_clump,clump_kb=kb_window,clump_r2=0,clump_p=pval))

  }

  gwas_clump<-Filter(function(a) any(!is.na(a)),gwas_clump)
  gwas_chr_list<-gwas_chr_list[names(gwas_clump)]
  message(paste0(Sys.time(),": Populating regions.."))
  populate_clumps=function(clumped_gwas_chr,gwas_chr,window=window,pval=pval){

    signif_snps<-clumped_gwas_chr
    x<-gwas_chr

    signif_snps$pos<-as.numeric(signif_snps$pos)
    signif_snps$start<-signif_snps$pos-window
    signif_snps$end<-signif_snps$pos+window

    start<-signif_snps[1,]$start
    end<-signif_snps[1,]$end
    df<-subset(x,pos>start & pos<end)
    df$n_snps_region<-rep(nrow(df),nrow(df))
    df$signif_snp_region<-rep(signif_snps[1,]$rsid,nrow(df))

    if(nrow(signif_snps)>1){
      for(i in 2:nrow(signif_snps)){
        start<-signif_snps[i,]$start
        end<-signif_snps[i,]$end
        df_2<-subset(x,pos>start & pos<end)
        df_2$n_snps_region<-rep(nrow(df_2),nrow(df_2))
        df_2$signif_snp_region<-rep(signif_snps[i,]$rsid,nrow(df_2))

        df<-rbind(df,df_2)
      }
      return(df)
    } else {
      return(df)
    }
  }


  window=window/2
  test=list()
  for(i in 1:length(gwas_chr_list)){
    test[[i]]<-populate_clumps(clumped_gwas=gwas_clump[[i]],gwas_chr=gwas_chr_list[[i]],window=window)
    }

  test<-do.call(rbind,test)
  return(test)


}



preprocess_mateqtlouts=function(mateqtlouts,GWAS){

    message("Keeping SNPs appearing in clumped GWAS..")

    mateqtlouts<-lapply(mateqtlouts,function(x){x<-x[which(x$SNP %in% GWAS$rsid),]
        return(x)})

    
    
    mateqtlouts<-lapply(mateqtlouts,function(x){
        row.names(x)<-paste0(x$SNP,"_",x$gene)
        return(x)
    })

    message("Calculating standard errors")
    mateqtlouts<-lapply(mateqtlouts,function(x){x<-x[!x$t.stat==0,]
        x<-x[!x$t.stat==-Inf,]
        x <- x[!is.infinite(x$t.stat), ]
        x$std_error=x$beta/x$t.stat
        return(x)})

    return(mateqtlouts)
}



### define cellCOLOC_object classes. simplifies downstream processing

# Helper function to create tibbles
createTibble <- function(regionData, SNP, gene, value) {
  DF <- tidyr::spread(regionData[, c(SNP, gene, value)], key = gene, value = value)
  rownames(DF) <- DF$SNP
  return(tibble::as_tibble(DF))
}

createCellTypeObject <- function(mateqtl, gwasList, genesToKeepList, maf_df) {
  regionsList <- list()

  for (regionName in names(gwasList)) {
    gwasSNPs <- gwasList[[regionName]]$rsid
    eQTLsnps <- unique(mateqtl$SNP)
    commonSNPs <- intersect(gwasSNPs, eQTLsnps)
    regionData <- mateqtl[mateqtl$SNP %in% commonSNPs, ]

    if (!is.null(genesToKeepList[[regionName]])) {
      genesToKeep <- genesToKeepList[[regionName]]
      regionData <- regionData[regionData$gene %in% genesToKeep, ]
    }

    if (nrow(regionData) == 0) next

    maf_info = maf_df[match(commonSNPs, maf_df$snp), ]

    regionsList[[regionName]] <- new("regionData", 
                                     beta = createTibble(regionData, "SNP", "gene", "beta"), 
                                     std_error = createTibble(regionData, "SNP", "gene", "std_error"),
                                     pvalue = createTibble(regionData, "SNP", "gene", "p.value"),
                                     FDR = createTibble(regionData, "SNP", "gene", "FDR"),
                                     maf_info = tibble::as_tibble(maf_info))
  }

  return(new("cellType", regions = regionsList))
}

createCellCOLOCobject <- function(mateqtl, gwasList, genesToKeepList, indiv_numbers, maf_df,gwas_type,gene_locs,snp_locs) {
  cellTypesList <- list()
    
  if(!gwas_type %in% c("cc", "quant")) {
    stop("gwasInfo should be either 'cc' or 'quant'")
  }

  for (cellType in names(mateqtl)) {
    cellObj = createCellTypeObject(mateqtl[[cellType]], gwasList, genesToKeepList,maf_df=maf_df)
    cellTypesList[[cellType]] <- cellObj

    n_indivs <- as.numeric(indiv_numbers[rownames(indiv_numbers) == cellType, "n_indivs_pseudobulk"])
    cellTypesList[[cellType]]@n_indivs <- n_indivs
  }

  gwas_N=as.numeric(gwasList[[1]]$N_total[1])
  if(gwas_type=="cc"){
    gwas_S=as.numeric(gwasList[[1]]$case.prop[1])
  }else{
    gwas_S=0
  }

  gwasList <- lapply(gwasList, function(gwasData) {
    alleles_df <- as.data.frame(gwasData[, c("rsid", "A1", "A2")])
      
    if ("MAF" %in% colnames(gwasData)) {
        maf_tibble = as_tibble(gwasData[, c("rsid", "MAF")])
    } else {
        warning("warning: MAF info missing from GWAS. Taking MAF from eQTL study")
        commonSNPs <- intersect(gwasData$rsid, maf_df$snp)
        maf_tibble = maf_df[match(commonSNPs, maf_df$snp), c("snp", "maf")]
        colnames(maf_tibble) <- c("rsid", "MAF")
        maf_tibble=as_tibble(maf_tibble)
    }
    new("gwasData", 
        beta = createTibble(gwasData, "rsid", "b"),
        std_error = createTibble(gwasData, "rsid", "se"),
        pvalue = createTibble(gwasData, "rsid", "pval"),
        maf=as_tibble(maf_tibble),
        alleles=as_tibble(alleles_df)
    )
  })

  return(new("cellCOLOCObj", cellTypes = cellTypesList, gwas = gwasList,gwasInfo=gwas_type,gwas_N=gwas_N,gwas_S=gwas_S,snp_locs=snp_locs,gene_locs=gene_locs))
}




##define COLOC functions 

run_coloc_specific = function(cellCOLOC_obj, cellTypeName, regionName, geneName){
    
    
  # Extracting the specified region and cell type
  # cellTypeName = "Microglia"
  # regionName = "rs6733839"
  # geneName = "BIN1"

    # Accessing the regionData for the specified cellType and region
    regionData <- cellCOLOC_obj@cellTypes[[cellTypeName]]@regions[[regionName]]

    # Generating eqtl_trait list
    eqtl_trait = list(
        beta = regionData@beta[[geneName]],  
        varbeta = regionData@std_error[[geneName]]^2,  
        snp = regionData@beta$SNP,
        type = "quant",
        N = cellCOLOC_obj@cellTypes[[cellTypeName]]@n_indivs,
        MAF = regionData@maf_info$maf,  
        sdY=1
    )
    
    

    # Accessing the gwasData for the specified region
    gwasData <- cellCOLOC_obj@gwas[[regionName]]

    # Generating gwas_trait list using the provided column names
    
    commonSNPs = eqtl_trait$snp[eqtl_trait$snp %in% gwasData@beta$rsid]
    idx = which(gwasData@beta$rsid %in% commonSNPs)
    maf_filt=gwasData@maf
    maf_filt=maf_filt[match(commonSNPs,maf_filt$rsid),]
    
    type=cellCOLOC_obj@gwasInfo

    if(type=="cc"){
      gwas_trait = list(
            beta = as.numeric(gwasData@beta$b[idx]),  # Subset and convert to numeric vector
            varbeta = (as.numeric(gwasData@std_error$se[idx]))^2,  # Subset, convert and square
            snp = gwasData@beta$rsid[idx],
            type = cellCOLOC_obj@gwasInfo, # Assuming the type is stored here
            N = cellCOLOC_obj@gwas_N,  # GWAS data doesn't provide this, you might fill it in later
            MAF = maf_filt$MAF, # Adjusted MAF column name
            s=cellCOLOC_obj@gwas_S
        )
    }else if (type=="quant"){
            gwas_trait = list(
            beta = as.numeric(gwasData@beta$b[idx]),  # Subset and convert to numeric vector
            varbeta = (as.numeric(gwasData@std_error$se[idx]))^2,  # Subset, convert and square
            snp = gwasData@beta$rsid[idx],
            type = cellCOLOC_obj@gwasInfo, # Assuming the type is stored here
            N = cellCOLOC_obj@gwas_N,  # GWAS data doesn't provide this, you might fill it in later
            MAF = maf_filt$MAF, # Adjusted MAF column name
            sdY=1
        )

    }    ### TODO: Run colocalization (to be provided)

    # Returning the traits for now, will add more returns as function evolves
    capture.output(coloc_result <-suppressWarnings(suppressMessages(coloc::coloc.abf(eqtl_trait, gwas_trait))),file=nullfile())
    resvec=coloc_result$summary
    lead_snp_data = coloc_result$results
    lead_snp_data = lead_snp_data[order(lead_snp_data$SNP.PP.H4, decreasing=T), ]
    
    lead_snp = lead_snp_data$snp[1]
    snph4 = lead_snp_data$SNP.PP.H4[1]
    
    # Extract eQTL FDR, eQTL p-value, and GWAS p-value for the lead SNP
    eqtl_fdr = regionData@FDR[regionData@beta$SNP == lead_snp, geneName]
    eqtl_pval = regionData@pvalue[regionData@beta$SNP == lead_snp, geneName]
    gwas_pval = gwasData@pvalue$pval[gwasData@beta$rsid == lead_snp]

    # Append these values to resvec
    resvec = c(regionName, cellTypeName, geneName, lead_snp, snph4, eqtl_pval, eqtl_fdr, gwas_pval, resvec)
    names(resvec) = c("region", "celltype", "gene", "lead_snp","SNP.PP.H4", "eQTL_pval", "eQTL_FDR", "GWAS_pval", "nsnps","PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4")
    resvec=unlist(resvec)

    return(resvec)
}


run_coloc_all=function(cellCOLOC_obj){
    df_all_results <- data.frame()

    cellTypes <- names(cellCOLOC_obj@cellTypes)
    for (cellType in cellTypes) {
        regions <- names(cellCOLOC_obj@cellTypes[[cellType]]@regions)
        for (region in regions) {
            genes_in_region <- colnames(cellCOLOC_obj@cellTypes[[cellType]]@regions[[region]]@beta)
            genes_in_region <- setdiff(genes_in_region, "SNP") # remove the SNP column
            for (gene in genes_in_region) {
                resvec <- run_coloc_specific(cellCOLOC_obj, cellType, region, gene)

                # Convert the resvec vector into a dataframe
                df_resvec <- as.data.frame(t(resvec))
                colnames(df_resvec) <- names(resvec)

                # Bind the result to the main dataframe
                df_all_results <- rbind(df_all_results, df_resvec)
            }
        }
    }


    df_all_results<- df_all_results %>% mutate(across(c("SNP.PP.H4","eQTL_FDR","eQTL_pval","GWAS_pval","PP.H0","PP.H1","PP.H2","PP.H3","PP.H4"),as.numeric))
    numeric_cols <- sapply(df_all_results, is.numeric)
    df_all_results[numeric_cols] <- lapply(df_all_results[numeric_cols], signif, 3) # 3 is the number of significant digits
    return(df_all_results)

}


##define MR functions

run_cellMR=function(cellCOLOC_obj,
use_coloc_lead_snp=TRUE,
coloc_res,pph4_cutoff=0.7,
eqtl_FDR_cutoff=0.05,r2_cutoff=0.01,
plink_bin=genetics.binaRies::get_plink_binary(),
path_to_binaries="/path_to_binaries/EUR/EUR"){
    
    coloc_res=coloc_res %>% filter(PP.H4>pph4_cutoff)
    if (nrow(coloc_res) == 0) {
      message("No rows left after the pph4 cutoff")
        return(NA)
    }
    
    MR_results=list()
    message(paste0(Sys.time(),": Preparing MR inputs. Pruning SNPs and harmonizing alleles.."))
    if(use_coloc_lead_snp==TRUE){
      message("Using coloc lead SNP as index SNP to prune around")
    }
    for(i in 1:nrow(coloc_res)){
        
        hit=coloc_res[i,]
        region=hit$region
        celltype=hit$celltype
        gene=hit$gene
        lead_snp=hit$lead_snp
        

        df=cellCOLOC_obj@cellTypes[[celltype]]@regions[[region]]@FDR[,c("SNP",gene)]
        df=df %>% filter((!!sym(gene)) <eqtl_FDR_cutoff) %>%
          select(SNP, !!sym(gene)) %>%
          rename_with(~ c("rsid", "pval")) %>% as.data.frame()
        
        if (nrow(df) == 0) {
        # warning("No SNPs left after filtering for eQTL FDR in iteration ", i)
            MR_results[[i]] = c(celltype,gene,NA,NA,NA)
            next
          }

        if(use_coloc_lead_snp==TRUE){

          ##this ensure that SNPS are pruned around the lead SNP by setting the p-value to zero.
          ## Occasionnally, the lead SNP proposed by COLOC is not the strongest eQTL in the region.
          if(use_coloc_lead_snp) {
            df <- df %>%
          mutate(pval = ifelse(rsid == lead_snp, 0, pval))
            }
            
        }

        ##do the LD pruning
        suppressMessages(capture.output(snps_tokeep<-ieugwasr::ld_clump_local(df,clump_r2=r2_cutoff,clump_p=eqtl_FDR_cutoff,plink_bin=plink_bin,
        bfile=path_to_binaries,clump_kb=1e6),file=nullfile()))

        trait_beta <- cellCOLOC_obj@cellTypes[[celltype]]@regions[[region]]@beta[, c("SNP",gene)] %>% filter(SNP %in% snps_tokeep$rsid)
        trait_se <- cellCOLOC_obj@cellTypes[[celltype]]@regions[[region]]@std_error[, c("SNP",gene)] %>% filter(SNP %in% snps_tokeep$rsid)
        trait_alleles<-cellCOLOC_obj@cellTypes[[celltype]]@regions[[region]]@maf_info[, c("snp","ref","alt")] %>% filter(snp %in% snps_tokeep$rsid)


        gwas_beta<-cellCOLOC_obj@gwas[[region]]@beta %>% filter(rsid %in% snps_tokeep$rsid)
        gwas_se <- cellCOLOC_obj@gwas[[region]]@std_error %>% filter(rsid %in% snps_tokeep$rsid)
        gwas_alleles=cellCOLOC_obj@gwas[[region]]@alleles %>% filter(rsid %in% snps_tokeep$rsid)


        ### MungeSumstats considers A2 to be the effect allele and not A1.
        ### In the below code, A1 = effect allele (as by convention).
        # https://neurogenomics.github.io/MungeSumstats/articles/MungeSumstats.html
        ### Changing so that "gwas_A1" is actually gwas_alleles$A2 .
        #Additionally, the SeqArray genotype data is in the form of 0,1,2 where 0 is alternate allele dosage.
        #In MatrixEQTL additive model, the ref allele is therefore the effect allele.

        mr_df=data.frame(SNP=trait_beta$SNP,eqtl_beta=trait_beta[,2],
                         eqtl_se=trait_se[,2],
                         eqtl_A1=trait_alleles$ref,
                         eqtl_A2=trait_alleles$alt,
                         gwas_beta=gwas_beta$b,
                         gwas_se=gwas_se$se,
                         gwas_A1=gwas_alleles$A2,
                        gwas_A2=gwas_alleles$A1)

        colnames(mr_df)=c("SNP","eqtl_beta","eqtl_se","eqtl_A1","eqtl_A2","gwas_beta","gwas_se","gwas_A1","gwas_A2")

        #harmonize alleles
        mr_df=mr_df %>%
          mutate(
            eqtl_beta = ifelse(eqtl_A1 != gwas_A1 & eqtl_A1 == gwas_A2, -eqtl_beta, eqtl_beta)
          )
        MRInputObject<-MendelianRandomization::mr_input(bx=mr_df$eqtl_beta,
                                                        bxse=mr_df$eqtl_se,by=mr_df$gwas_beta,byse=mr_df$gwas_se)
    
      res=MendelianRandomization::mr_ivw(MRInputObject)
      res=data.frame(Method="IVW",Estimate=res@Estimate,`P-value`=res@Pvalue,check.names=F,std_error=res@StdError)

        mr_results=res %>% filter(Method == "IVW") %>%
          mutate(
            IVW_beta = signif(Estimate, digits = 4),
            IVW_pval = signif(`P-value`, digits = 4),
            IVW_se=std_error
          ) %>%
          select(IVW_beta, IVW_pval,IVW_se)

        mr_results$IVs <- paste(unique(mr_df$SNP), collapse = ",")
        mr_results$gene=gene
        mr_results$celltype=celltype
        mr_results=mr_results %>% select(celltype,gene,IVs,IVW_beta,IVW_pval,IVW_se)
        
        MR_results[[i]]=mr_results
    }
    MR_results=as.data.frame(do.call(rbind,MR_results))
    return(MR_results)
    
    
    
}


### S4 classes for cellCOLOC_obj

setClass(
  "regionData",
  slots = list(
    beta = "tbl_df",
    std_error = "tbl_df",
    pvalue = "tbl_df",
    FDR = "tbl_df",
    maf_info = "tbl_df"
  )
)



# Define the cellType class
setClass(
  "cellType",
  representation(
    regions = "list",
      n_indivs = "numeric"# A list of regionData objects
  )
)

setClass(
  "cellCOLOCObj",
  slots = list(
    cellTypes = "list",
      gwas = "list",
      gwasInfo = "character",
      gwas_N="numeric",
      gwas_S="ANY",
      snp_locs="ANY",
      gene_locs="ANY"
  )
)

setClass(
  "gwasData",
  slots = list(
    beta = "ANY",
    std_error = "ANY",
      pvalue="ANY",
      maf="ANY",
       alleles="ANY"
  )
)
