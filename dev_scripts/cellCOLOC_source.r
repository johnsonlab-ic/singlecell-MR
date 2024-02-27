
cellCOLOC=function(GWAS,
  mateqtlouts,
  gene_loc_file=list.files(pattern="gene_locations","../../MatrixEQTL_IO/",full.names=TRUE)[1],
  maf_file="../../MatrixEQTL_IO/MAF_mat.csv",
  snp_locs_file="../../MatrixEQTL_IO/snp_chromlocations.csv",
  indiv_numbers_file="../../MatrixEQTL_IO/n_indivs_pseudobulk.txt",
  intersected_mateqtlouts_file="GWAS_intersected_mateqtlouts.rds",
  preprocess_mateqtlouts=TRUE,
  preprocess_GWAS=TRUE,
  processed_GWAS_file="GWAS_clumped.rds",
  GWASsignif=5e-8,
  GWAS_window=1e6,
  GWAS_type=c("quant","cc"),
  GWAS_name="AD",
  cellMR=TRUE,
  pph4_cutoff=0.7,
  use_coloc_lead_snp=TRUE){

  
  start_time<-Sys.time()
  message(paste0(Sys.time(),": Starting cellCOLOC pipeline."))
  
  if(preprocess_mateqtlouts==TRUE){
    if(preprocess_GWAS==TRUE){
      message(paste0(Sys.time(),": preprocess_GWAS=TRUE. Reading in unprocessed GWAS. Selecting regions for colocalisation based on supplied parameters: \n GWAS window=",GWAS_window,
      "\n GWAS signif=",GWASsignif))
      message("This step is server-based and can therefore take a long time (>2h).\nTime will also depend on the number of loci passing thresholds.")
      message(paste0("You can supply your own clumped GWAS file, as long as it is saved in outdir/'run_name'/ as '",processed_GWAS_file,"' and you specify preprocess_GWAS=FALSE.\nSee tutorial for formatting. "))



      GWAS<-as.data.frame(suppressWarnings(data.table::fread(GWAS)))


      message(paste0(Sys.time(),": Checking column names.."))
      
      if(GWAS_type=="cc"){
        check_colnames(GWAS,c("rsid","pval","b","pos","chr","se","A2","A1","case.prop","N_total"))
      }else{
        check_colnames(GWAS,c("rsid","pval","b","pos","chr","se","A2","A1","N_total"))
      }
      if(any(GWAS$pval<GWASsignif)==FALSE){
          stop(paste0("No associations below specified p-value: ",GWASsignif,". Change 'GWASsignif' if you want to proceed."))
      }

      ##Filter out zero values
      GWAS=GWAS %>% dplyr::filter(!b %in% 0 & !se %in% 0 & !pval %in% 0)
      GWAS=GWAS %>% dplyr::arrange(rsid, pval) %>% dplyr::distinct(rsid, .keep_all = TRUE)

      GWAS<-select_regions(GWAS,window=GWAS_window,pval=GWASsignif)
      saveRDS(GWAS,processed_GWAS_file)
      message(paste0(Sys.time(),": Processed GWAS saved under '",processed_GWAS_file,"'."))


    }
      message(paste0(Sys.time(),": preprocess_GWAS=FALSE, Reading in pre-processed GWAS..."))
      GWAS<-readRDS(processed_GWAS_file)

      message(paste0(Sys.time(),": preprocess_mateqtlouts=TRUE, reading in and filtering mateqtlouts.."))
      mateqtlouts<-readRDS(mateqtlouts)
      mateqtlouts<-preprocess_mateqtlouts(mateqtlouts,GWAS)

    saveRDS(mateqtlouts,"GWAS_intersected_mateqtlouts.rds")

  } else {
    message(paste0(Sys.time(),": preprocess=FALSE. Reading in preprocessed GWAS and mateqtlouts."))
    mateqtlouts<-readRDS(intersected_mateqtlouts_file)
    GWAS<-readRDS(processed_GWAS_file)
  }
  if(length(grep("ch",GWAS$chr))==0){
    GWAS$chr<-paste0("chr",GWAS$chr)
  }
  processed_gwas<-split(GWAS,GWAS$signif_snp_region)

  nregions<-length(processed_gwas)
  message(paste0(Sys.time(),": ",nregions," regions were selected."))
  

  # -============================================================================================================
  # -=-=-=-=-=-=-=- STEP 2. Generating cellCOLOC objects and run COLOC (+ MR if specified) -=-=-=-=-=-=-===
  # -============================================================================================================

  date=Sys.Date()
  message(paste0(Sys.time(),": Prepping coloc data.."))

  # first, get the genes to keep per region
  gene_locs<-read.table(gene_loc_file)
  genestokeep=lapply(processed_gwas,get_genes_per_region,gene_locs=gene_locs)
  indiv_numbers=read.table(indiv_numbers_file)
  maf_df=as.data.frame(data.table::fread(maf_file))
  snp_locs=as.data.frame(data.table::fread(snp_locs_file))

  # create the cellCOLOC object
  cellCOLOC_obj=createCellCOLOCobject(mateqtlouts, processed_gwas, 
  genestokeep,indiv_numbers = indiv_numbers,
  maf_df=maf_df,gwas_type=GWAS_type,snp_locs=snp_locs,gene_locs=gene_locs)

  ##filter out regions with less than N snps
  min_snps=100

  cellCOLOC_obj@gwas <- remove_small_regions(cellCOLOC_obj@gwas,min_snps=min_snps)
  cellTypes_names <- names(cellCOLOC_obj@cellTypes)
  for (cellType in cellTypes_names) {
    cellCOLOC_obj@cellTypes[[cellType]]@regions <- remove_small_regions(cellCOLOC_obj@cellTypes[[cellType]]@regions,min_snps=min_snps)
  }
  saveRDS(cellCOLOC_obj,paste0(GWAS_name,"_cellCOLOC_obj.rds"))


  ### Run COLOC analysis

  coloc_df=run_coloc_all(cellCOLOC_obj)
  coloc_df=coloc_df[order(coloc_df$PP.H4,decreasing=T),]
  write.table(coloc_df,paste0(GWAS_name,"_COLOC_results.txt"))
  message(paste0(Sys.time(),": cellCOLOC run complete."))

  if(cellMR==TRUE){
    mr_results=run_cellMR(cellCOLOC_obj,use_coloc_lead_snp=use_coloc_lead_snp,coloc_df,pph4_cutoff=pph4_cutoff,eqtl_FDR_cutoff=0.05,r2_cutoff=0.01)
    if(length(mr_results)>1){
        write.table(mr_results,paste0(GWAS_name,"_MR_results.txt"))
    }
    message(paste0(Sys.time(),": cellMR run complete."))
  }

}


get_genes_per_region=function(processed_gwas_region,gene_locs){

  gene_locs<-split(gene_locs,gene_locs$chr)
  gene_locs<-gene_locs[which(names(gene_locs)==unique(processed_gwas_region$chr))]
  gene_locs<-do.call(rbind,gene_locs)
  rownames(gene_locs)<-rep(1:nrow(gene_locs))

  startpos<-min(processed_gwas_region$pos)
  endpos<-max(processed_gwas_region$pos)

  #this is to see whether the end of the gene is within region
  gene_locs_1<-subset(gene_locs,right>startpos & left<endpos)

  #this is to see whether the start of the gene is within region
  gene_locs<-subset(gene_locs,left>startpos & left<endpos)

  #combine the two
  gene_locs_1<-setdiff(gene_locs_1,gene_locs)
  gene_locs<-rbind(gene_locs,gene_locs_1)

  return(gene_locs$geneid)

}

preprocess_mateqtlouts=function(mateqtlouts,GWAS){

    message("Keeping SNPs appearing in clumped GWAS..")

    mateqtlouts<-lapply(mateqtlouts,function(x){x<-x[which(x$SNP %in% GWAS$rsid),]
        return(x)})

    
    
    mateqtlouts<-lapply(mateqtlouts,function(x){
        row.names(x)<-paste0(x$SNP,"_",x$gene)
        return(x)
    })


    # hyprcoloc can't handle zero values. These are filtered out
    message("Calculating standard errors")
    mateqtlouts<-lapply(mateqtlouts,function(x){x<-x[!x$t.stat==0,]
        x<-x[!x$t.stat==-Inf,]
        x <- x[!is.infinite(x$t.stat), ]
        x$std_error=x$beta/x$t.stat
        return(x)})

  #   rowsdiff=rowsbefore-nrow(mateqtlouts[[1]])
  #   message(paste0(Sys.time(),": Intersection done. ",rowsdiff, " cis-EQTLs removed from analysis."))
    return(mateqtlouts)
}


# Redefine the regionData class
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


createCellTypeObject <- function(mateqtl, gwasList, genesToKeepList, maf_df) {
  regionsList <- list()

  for (regionName in names(gwasList)) {
    gwasSNPs <- gwasList[[regionName]]$rsid
    eQTLsnps <- unique(mateqtl$SNP)
    
    # Take the intersection of GWAS and eQTL SNPs
    commonSNPs <- intersect(gwasSNPs, eQTLsnps)

    # Filter the mateqtl data based on the commonSNPs
    regionData <- mateqtl[mateqtl$SNP %in% commonSNPs, ]

    # Check if genes are provided for this region
    if (is.null(genesToKeepList[[regionName]])) next

    # Further filter the regionData by genes in the genesToKeepList
    genesToKeep <- genesToKeepList[[regionName]]
    regionData <- regionData[regionData$gene %in% genesToKeep, ]

    # Check if after filtering, the regionData is empty
    if (nrow(regionData) == 0) next

    # Create beta, std_error, pval, and FDR tibbles with genes as columns and SNPs as rows
    betaDF <- tidyr::spread(regionData[, c("SNP", "gene", "beta")], key = "gene", value = "beta")
    stdErrorDF <- tidyr::spread(regionData[, c("SNP", "gene", "std_error")], key = "gene", value = "std_error")
    pvalDF <- tidyr::spread(regionData[, c("SNP", "gene", "p.value")], key = "gene", value = "p.value")
    FDR_DF <- tidyr::spread(regionData[, c("SNP", "gene", "FDR")], key = "gene", value = "FDR")
    maf_info = maf_df[match(commonSNPs, maf_df$snp), ]

    # Ensure SNPs are the row names
    rownames(betaDF) <- betaDF$SNP
    rownames(stdErrorDF) <- stdErrorDF$SNP

    regionsList[[regionName]] <- new("regionData", 
                                     beta = tibble::as_tibble(betaDF), 
                                     std_error = tibble::as_tibble(stdErrorDF),
                                     pvalue = tibble::as_tibble(pvalDF),
                                     FDR = tibble::as_tibble(FDR_DF),
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

    # if (all(c("A1", "A2") %in% colnames(gwasData))) {
        alleles_df <- as.data.frame(gwasData[, c("rsid", "A1", "A2")])
      # } else {
      #   warning("warning: A1 and A2 allele info missing from GWAS")
      #   alleles_df <- data.frame(rsid = character(0), A1 = character(0), A2 = character(0))  # Empty data frame
      # }
      
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
          beta = as_tibble(gwasData[, c("rsid", "b")]),
          std_error = as_tibble(gwasData[, c("rsid", "se")]),
          pvalue = as_tibble(gwasData[, c("rsid", "pval")]),
          maf=as_tibble(maf_tibble),
          alleles=as_tibble(alleles_df)
      )
    })
  

  return(new("cellCOLOCObj", cellTypes = cellTypesList, gwas = gwasList,gwasInfo=gwas_type,gwas_N=gwas_N,gwas_S=gwas_S,snp_locs=snp_locs,gene_locs=gene_locs))
}


setMethod(
  "show",
  "cellCOLOCObj",
  function(object) {
    cat("A cellCOLOC data object containing colocalisation input data for multiple cell types.\n")

    # Number and names of cell types
    cellTypeNames <- names(object@cellTypes)
    cat(length(cellTypeNames), "cell types:", toString(cellTypeNames), "\n")

    allRegionsCounts <- sapply(object@cellTypes, function(cellType) length(cellType@regions))
    cat("Total regions:", allRegionsCounts[1], "\n")

    # Information about genes and SNPs in regions
    allGenesCounts <- sapply(object@cellTypes, function(cellType) {
      sapply(cellType@regions, function(region) {
        ncol(region@beta)
      })
    })

    allSNPsCounts <- sapply(object@cellTypes, function(cellType) {
      sapply(cellType@regions, function(region) {
        nrow(region@beta)
      })
    })
    cat("------------------------------ \n")
    cat("GWAS Type:", object@gwasInfo, "\n")
    cat("Max number of genes in a region:", max(unlist(allGenesCounts)), "\n")
    cat("Min number of genes in a region:", min(unlist(allGenesCounts)), "\n")
    cat("Max number of SNPs in a region:", max(unlist(allSNPsCounts)), "\n")
    cat("Min number of SNPs in a region:", min(unlist(allSNPsCounts)), "\n")
  }
)
remove_small_regions <- function(regions,min_snps) {
  regions_to_keep <- sapply(regions, function(region) {
    beta_rows <- nrow(region@beta)
    beta_rows >= min_snps
  })
  regions[regions_to_keep]
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

run_coloc_specific = function(cellCOLOC_obj, cellTypeName, regionName, geneName){
    
    
    # Extracting the specified region and cell type
  # cellTypeName = "Microglia"
  # regionName = "rs6733839"
  # geneName = "BIN1"

    # Accessing the regionData for the specified cellType and region
    regionData <- cellCOLOC_obj@cellTypes[[cellTypeName]]@regions[[regionName]]

    # Generating eqtl_trait list
    eqtl_trait = list(
        beta = regionData@beta[[geneName]],  # All values for the specified gene
        varbeta = regionData@std_error[[geneName]]^2,  # Squaring all values for the specified gene
        snp = regionData@beta$SNP,
        type = "quant",
        N = cellCOLOC_obj@cellTypes[[cellTypeName]]@n_indivs,
        MAF = regionData@maf_info$maf,  # Adjusted MAF column name
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

    # this is the local version
    error_catch_clump=function(gwas,clump_kb,
      clump_r2,clump_p,plink_bin,bfile){

      tryCatch({test<-ieugwasr::ld_clump_local(gwas,clump_kb=clump_kb,clump_r2=clump_r2,clump_p=clump_p,
      plink_bin=plink_bin,bfile=bfile);},error=function(e){test<<-NA});
      return(test)

    }
    message(paste0("local: ",local))
    options(warn=-1)
  #   error_catch_clump=function(gwas,clump_kb,
  #     clump_r2,clump_p,plink_bin,bfile){

  #     tryCatch({test<-suppressMessages(ld_clump_local(gwas,clump_kb=clump_kb,clump_r2=clump_r2,clump_p=clump_p,
  #     plink_bin=plink_bin,bfile=bfile));},error=function(e){test<<-NA});
  #     return(test)

  #   }
  # library(ieugwasr)
  # library(genetics.binaRies)
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
    # if(parallel==TRUE){
    #   gwas_clump<-mclapply(gwas_chr_list,error_catch_clump,ld_clump_local,clump_kb=kb_window,clump_r2=0,clump_p=pval,plink_bin=plink_bin,
    #   bfile=path_to_binaries,mc.cores=parallel.cores)
    # }else{
      gwas_clump<-suppressMessages(lapply(gwas_chr_list,error_catch_clump,clump_kb=kb_window,clump_r2=0,clump_p=pval,plink_bin=plink_bin,
      bfile=path_to_binaries))

    # }
  } else{

    # if(parallel==TRUE){
    #   gwas_clump<-mclapply(gwas_chr_list,ld_clump,clump_kb=kb_window,clump_r2=0,clump_p=pval,mc.cores=parallel.cores)
    # }else{
      gwas_clump<-suppressMessages(lapply(gwas_chr_list,ieugwasr::ld_clump,clump_kb=kb_window,clump_r2=0,clump_p=pval))

    # }

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


run_cellMR=function(cellCOLOC_obj,
use_coloc_lead_snp=FALSE,
coloc_res,pph4_cutoff=0.7,
eqtl_FDR_cutoff=0.05,r2_cutoff=0.01,
plink_bin=genetics.binaRies::get_plink_binary(),
path_to_binaries="/rds/general/user/ah3918/projects/single_cell_eqtl/live/MRC/REFERENCE_DATASETS/clumping_ref/EUR/EUR"){
    
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

        suppressMessages(capture.output(snps_tokeep<-ieugwasr::ld_clump_local(df,clump_r2=r2_cutoff,clump_p=eqtl_FDR_cutoff,plink_bin=plink_bin,
      bfile=path_to_binaries,clump_kb=1e6),file=nullfile()))

        trait_beta <- cellCOLOC_obj@cellTypes[[celltype]]@regions[[region]]@beta[, c("SNP",gene)] %>% filter(SNP %in% snps_tokeep$rsid)
        trait_se <- cellCOLOC_obj@cellTypes[[celltype]]@regions[[region]]@std_error[, c("SNP",gene)] %>% filter(SNP %in% snps_tokeep$rsid)
        trait_alleles<-cellCOLOC_obj@cellTypes[[celltype]]@regions[[region]]@maf_info[, c("snp","ref","alt")] %>% filter(snp %in% snps_tokeep$rsid)


        gwas_beta<-cellCOLOC_obj@gwas[[region]]@beta %>% filter(rsid %in% snps_tokeep$rsid)
        gwas_se <- cellCOLOC_obj@gwas[[region]]@std_error %>% filter(rsid %in% snps_tokeep$rsid)
        gwas_alleles=cellCOLOC_obj@gwas[[region]]@alleles %>% filter(rsid %in% snps_tokeep$rsid)


        #### 28th Dec 2023 UPDATE - MungeSumstats considers A2 to be the effect allele and not A1.
        ### In the below code, A1 = effect allele (as by convention).
        # https://neurogenomics.github.io/MungeSumstats/articles/MungeSumstats.html
        ### Changing so that "gwas_A1" is actually gwas_alleles$A2 .

        mr_df=data.frame(SNP=trait_beta$SNP,eqtl_beta=trait_beta[,2],
                         eqtl_se=trait_se[,2],
                         eqtl_A1=trait_alleles$ref,
                         eqtl_A2=trait_alleles$alt,
                         gwas_beta=gwas_beta$b,
                         gwas_se=gwas_se$se,
                         gwas_A1=gwas_alleles$A2,
                        gwas_A2=gwas_alleles$A1)

        colnames(mr_df)=c("SNP","eqtl_beta","eqtl_se","eqtl_A1","eqtl_A2","gwas_beta","gwas_se","gwas_A1","gwas_A2")

        mr_df=mr_df %>%
          mutate(
            eqtl_beta = ifelse(eqtl_A1 != gwas_A1 & eqtl_A1 == gwas_A2, -eqtl_beta, eqtl_beta)
          )
        MRInputObject<-MendelianRandomization::mr_input(bx=mr_df$eqtl_beta,
                                                        bxse=mr_df$eqtl_se,by=mr_df$gwas_beta,byse=mr_df$gwas_se)
        #         res<-MendelianRandomization::mr_allmethods(MRInputObject, method = "ivw")@Values
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



######## Plotting functions

create_regional_association_plot <- function(object, gene, celltype, 
                                             highlight_snps = NULL, ylim = NULL,
                                             fill_color = "#B9D3FF", 
                                             plink_bin_path="/rds/general/user/ah3918/projects/single_cell_eqtl/live/MRC/TMP/plink", 
                                             path_to_binaries="/rds/general/user/ah3918/projects/single_cell_eqtl/live/MRC/REFERENCE_DATASETS/clumping_ref/EUR/EUR") {
  library(ggplot2)
  library(data.table)
  library(ieugwasr)
  library(ggrepel)
  
  cellTypeObj <- object@cellTypes[[celltype]]
    regionName <- NULL
  for (name in names(cellTypeObj@regions)) {
    if (gene %in% colnames(cellTypeObj@regions[[name]]@beta)) {
      regionName <- name
      break
    }
  }
  
  if (is.null(regionName)) {
    stop(paste("No region found containing gene:", gene))
  }

  regionDataObj <- cellTypeObj@regions[[regionName]]
  # Assuming gene names are used as keys for regions list
  snp_locs <- object@snp_locs
  gene_locs <- object@gene_locs

  trait <- paste0(celltype, ".", gene)

  # Prepare data for plotting
  eqtl_pvals=regionDataObj@pvalue[,c("SNP",gene)]
  gwas_pvals=object@gwas[[regionName]]@pvalue

  common_snps=intersect(eqtl_pvals$SNP,gwas_pvals$rsid)
  eqtl_pvals=eqtl_pvals[match(common_snps,eqtl_pvals$SNP),]
  gwas_pvals=gwas_pvals[match(common_snps,gwas_pvals$rsid),]

  df=data.frame(eqtl_pvals[,2],gwas_pvals[,2])
  colnames(df)=c(trait,"GWAS")
  rownames(df)=eqtl_pvals$SNP

  # Identify the lead SNP in the eQTL and GWAS data
  lead_snp_eqtl <- rownames(df)[which.min(df[, trait])]
  lead_snp_gwas <- rownames(df)[which.min(df[, "GWAS"])]

  # Compute LD matrix
  ldmat <- ieugwasr::ld_matrix_local(rownames(df),with_alleles=F,bfile=path_to_binaries,plink_bin=plink_bin_path)
  missing_snps <- setdiff(rownames(df), rownames(ldmat))
  ldmat <- rbind(ldmat, matrix(0, nrow = length(missing_snps), ncol = ncol(ldmat), dimnames = list(missing_snps, colnames(ldmat))))
  ldmat <- cbind(ldmat, matrix(0, nrow = nrow(ldmat), ncol = length(missing_snps), dimnames = list(rownames(ldmat), missing_snps)))
  ldmat=ldmat^2


  df$pos <- snp_locs[match(rownames(df), snp_locs$annot),]$pos
  chromosome <- unique(snp_locs[match(rownames(df), snp_locs$annot),]$chrom)

  # Transform p-values to -log10 scale
  df[,"GWAS"] <- -log10(df[,"GWAS"])
  df[,trait] <- -log10(df[,trait])
    # Define R2 vector
  r2 <- ldmat[lead_snp_gwas, ]
      
      ##
      if(!is.null(highlight_snps)){
          r2=ldmat[highlight_snps[1], ]
      }


    # Add r2 values to df
    df$r2 <- r2[rownames(df)]
    
    # Define R2 categories
    df$r2_category <- "miss"
    df$r2_category[df$r2>=0 & df$r2<0.2 & !is.na(df$r2)] <- "0.0-0.2"
    df$r2_category[df$r2>=0.2 & df$r2<0.4 & !is.na(df$r2)] <- "0.2-0.4"
    df$r2_category[df$r2>=0.4 & df$r2<0.6 & !is.na(df$r2)] <- "0.4-0.6"
    df$r2_category[df$r2>=0.6 & df$r2<0.8 & !is.na(df$r2)] <- "0.6-0.8"
    df$r2_category[df$r2>=0.8 & df$r2<=1 & !is.na(df$r2)] <- "0.8-1.0" 
    df$r2_category <- factor(df$r2_category, levels=c("miss", "0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"))
    
    # Identify maximum p-value for setting y-axis limit
  #   ymax <- max(max(df[,"GWAS"]), max(df[,trait])) * 1.1
      ymax <- max(max(df[,"GWAS"])) * 1.1
      ymax_eqtl=max(df[,trait])* 1.1
      
      if(!is.null(ylim)){
  #         ymax=ylim
          ymax_eqtl=ylim
      }
      
      #### colour palette
      
      colors <- c("#FF0000", "#FFFF00", "#00FF00", "#00FFFF", "#0000FF")
      gradient <- colorRampPalette(colors)(6)
      snpcolours <- gradient
      
  #     colors <- c("#DCDCDC", "#0000FF", "#FF0000")

  #     # Create a color ramp palette and generate the gradient of length 6
  #     gradient <- colorRampPalette(colors)(6)
  #     snpcolours=gradient
  #     snpcolours=c("#1D3A60", "#3C4D8A", "#5770B2", "#BD83B8", "#F5D7DB", "#F1916D")
      

      ### add highlighted snps
      
  highlight_data <- NULL
    if (!is.null(highlight_snps)) {
      highlight_data <- subset(df, rownames(df) %in% highlight_snps)
    }else{
        highlight_snps=c(unique(lead_snp_eqtl,lead_snp_gwas))
        
        highlight_data <- subset(df, rownames(df) %in% highlight_snps)
    }

    # eQTL plot
    marker_plot_eqtl <- ggplot(data = df, mapping = aes(x = pos, y = !!as.symbol(trait))) +
      geom_point(aes(fill = r2),stroke=0.35, pch = 21, size = 3.5/.pt) +
      scale_fill_gradient2(low="#49409a",mid="#f9eb5d",high="#ff5234",midpoint=0.5)
      
    if (is.null(highlight_data)) {
      marker_plot_eqtl <- marker_plot_eqtl +
        geom_point(data = subset(df, rownames(df) == lead_snp_eqtl), aes(fill = r2), pch = 23, size = 5/.pt, fill = "red") 
    } else {
      marker_plot_eqtl <- marker_plot_eqtl +
        geom_text_repel(data = highlight_data, aes(label = rownames(highlight_data)),max.overlaps=Inf,min.segment.length=0,size=5/.pt) +
        geom_point(data = highlight_data, aes(fill =r2), pch = 23, size = 5/.pt, fill = "red") +
        geom_segment(data = highlight_data, 
                    aes(xend = pos, yend = GWAS, color = rownames(highlight_data)),                                        
                    linetype = 2, color = "grey50")
    }
    trait=gsub("[.]"," / ",trait)
    marker_plot_eqtl <- marker_plot_eqtl +
      scale_y_continuous(limits = c(0, ymax_eqtl)) +
      labs(y = "eQTL\n-log10(p)",title=trait, x = NULL)+
    theme_bw() +
      theme(
       legend.position = c(0.92, 0.75),
      legend.background = element_rect(fill="transparent"), # to remove the legend background
          legend.text = element_text(size=5),
          legend.key.width = unit(0.15, "cm"),
          legend.key.height = unit(0.15, "cm"),
          legend.text.align = -1,
          legend.title=element_text(size=6),
          axis.title.y = element_text(vjust = 2.25, size = 7, face = "bold"),
            panel.grid.major = element_blank(),
            axis.text.y = element_text(size=5,colour="black"),
            axis.text.x = element_blank(),  # No x-axis labels
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=7, face = "bold"),
            panel.border = element_rect(colour = "black", fill = NA, size = 1))  # Full square border

  
    # GWAS plot
    marker_plot_gwas <- ggplot(data = df, mapping = aes(x = pos, y = GWAS)) +
      geom_point(aes(fill = r2),stroke=0.35, pch = 21, size = 3.5/.pt) + 
      scale_fill_gradient2(low="#49409a",mid="#f9eb5d",high="#ff5234",midpoint=0.5)
    
    if (is.null(highlight_data)) {
      marker_plot_gwas <- marker_plot_gwas +
        geom_point(data = subset(df, rownames(df) == lead_snp_gwas), aes(fill =r2), pch = 23, size = 5/.pt, fill = "red")  # different shape and color
    } else {
      marker_plot_gwas <- marker_plot_gwas + 
        geom_text_repel(data = highlight_data, aes(label = rownames(highlight_data)),max.overlaps=Inf,min.segment.length=0,size=5/.pt) +
        geom_point(data = highlight_data, aes(fill =r2), pch = 23, size = 5/.pt, fill = "red") +
        geom_segment(data = highlight_data, 
                    aes(xend = pos, yend = GWAS, color = rownames(highlight_data)),                                        
                    linetype = 2, color = "grey50")
    }

    marker_plot_gwas <- marker_plot_gwas +
      scale_y_continuous(limits = c(0, ymax)) +
      labs(y = "GWAS\n-log10(p)", x = NULL)+
      theme_bw() +
      theme(legend.position = c(0.92, 0.75),
      legend.background = element_rect(fill="transparent"), # to remove the legend background
          legend.text = element_text(size=5),
          legend.key.width = unit(0.15, "cm"),
          legend.key.height = unit(0.15, "cm"),
          legend.text.align = 1,
          legend.title=element_text(size=6),
            axis.title.y = element_text(vjust = 2.25, size = 7, face = "bold"),
            axis.text.y = element_text(size=5,colour="black"),
            panel.grid.major = element_blank(),
            axis.text.x = element_blank(),  # No x-axis labels
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, size = 1))  # Full square border
                                                

                                                                  
    # Combine the plots
      
    gene_plot=create_gene_plot(gene_df=gene_locs,
                              chromosome=chromosome,
                              min_pos=min(df$pos),
                              max_pos=max(df$pos),
                              target_gene=gene)
      

  combined_plot <- (marker_plot_eqtl /
                    marker_plot_gwas /
                    gene_plot)

  return(combined_plot)
}
create_gene_plot <- function(gene_df, chromosome, min_pos, max_pos, target_gene) {
  # Filter genes within the SNP range and on the correct chromosome
  gene_df <- gene_df[gene_df$chr == chromosome & gene_df$left >= min_pos & gene_df$right <= max_pos,]
    
    

  gene_plot <- ggplot() +
    geom_segment(data = gene_df, 
                 aes(x = left, xend = right, y = geneid, yend = geneid, 
                     color = ifelse(geneid == target_gene, "red", "black")), linewidth = 1) +
    geom_text(data = gene_df, aes(x = (left + right)/2, y = geneid, label = geneid), vjust = -1, size = 5/.pt) +
    scale_color_identity() + scale_x_continuous(limits=c(min_pos,max_pos))+
    labs(y = NULL, x = paste0("Position on ", chromosome)) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(vjust = -2.25, size = 7, face = "bold"),
          axis.text.x = element_text(size=5,colour="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1)) # Full square border


    gene_plot <- gene_plot + scale_y_discrete(expand = expansion(mult = c(0.1, 0.3)))
    

    
  return(gene_plot)
}