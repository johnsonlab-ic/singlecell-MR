#' and results summarization.
#'
#' @param outdir output directory. Make sure you use the same output directory as cellQTL.
#' @param run_name String to name the run. We recommend you name this similar to the GWAS you are using.
#' @param GWAS Path to input GWAS summary statistics file, read in with "data.table::fread"
#' @param mateqtlouts Path to input mateqtlouts.rds to use. Look for it in your latest cellQTL run.
#' @param gene_loc_file Path to gene location file - by default will take the first gene_loc file in MatrixEQTL_IO.
#' @param snp_loc_file Path to snp location file - by default will take snp_chromlocations.csv file in MatrixEQTL_IO
#' @param indiv_numbers_file Path to indiv numbers file - by default will take n_indivs_pseudobulk.txt file in MatrixEQTL_IO. Needed for COLOC.
#' @param preprocess_GWAS Default is TRUE. Takes the file specified in "GWAS" and performs region selection based on GWASsignif and GWAS_window parameters.
#' @param processed_GWAS_file Path to alternate clumped GWAS.
#' @param GWASsignif Numeric value. Default is 5e-8, which means the function will create regions centered around SNPs below this p-value.
#' @param GWAS_window Numeric value. Indicates window size (make sure you match this with your cellQTL run)
#' @param GWAS_type Either "cc" for case/control or "quant" for continuous trait. Needed for COLOC.
#' @param coloc.p12_prior Value for the coloc.p12_prior. This is the default COLOC value - we don't recommend changing this unless you know what you're doing!.
#' @export

cellCOLOC=function(GWAS,
  mateqtlouts,
  gene_loc_file=list.files(pattern="gene_locations","../../MatrixEQTL_IO/",full.names=TRUE)[1],
  maf_file="../../MatrixEQTL_IO/MAF_mat.csv",
  snp_loc_file="../../MatrixEQTL_IO/snp_chromlocations.csv",
  indiv_numbers_file="../../MatrixEQTL_IO/n_indivs_pseudobulk.txt",
  intersected_mateqtlouts_file="GWAS_intersected_mateqtlouts.rds",
  preprocess_mateqtlouts=TRUE,
  preprocess_GWAS=TRUE,
  processed_GWAS_file="GWAS_clumped.rds",
  GWASsignif=5e-8,
  GWAS_window=1e6,
  hyprcoloc=FALSE,
  coloc=TRUE,
  GWAS_type=c("quant","cc"),
  coloc.p12_prior=1e-5,
  hcoloc_prior.1=1e-4,
  hcoloc_prior.c=0.02,
  hcoloc_prior.12=NULL,
  run_name="test"){

  
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




  #make sure GWAS is a data frame to check chr format
  if(length(grep("ch",GWAS$chr))==0){
    GWAS$chr<-paste0("chr",GWAS$chr)
  }
  processed_gwas<-split(GWAS,GWAS$signif_snp_region)

  nregions<-length(processed_gwas)
  message(paste0(Sys.time(),": ",nregions," regions were selected."))
    
    # -============================================================================================================
  # -=-=-=-=-=-=-=- STEP 2. Generating the gene-specific SE and Beta matrices -=-=-=-=-=-=-===
  # -============================================================================================================
  
  gene_locs<-read.table(gene_loc_file)
  date=Sys.Date()
  message(paste0(Sys.time(),": Prepping coloc data.."))
  coloc_data<-prep_coloc_data(mateqtlouts,processed_gwas,gene_locs)
  coloc_data<-add_directionality(coloc_data)
  saveRDS(coloc_data,paste0(run_name,"_coloc_data.rds"))
  

  if(length(coloc_data$betas)!=length(processed_gwas)){
      regions<-names(coloc_data$betas)
      processed_gwas<-processed_gwas[regions]
      regions_lost=nregions-length(processed_gwas)
      message(paste0(Sys.time(),": ",regions_lost," regions were lost because no genes from gene_loc file were found within."))
  }
    
   
  message(paste0(Sys.time(),": Reading in all auxiliary files..."))
  
    indiv_numbers<-read.table(indiv_numbers_file)
    snplocs<-as.data.frame(data.table::fread(snp_loc_file))
    maf_df<-as.data.frame(data.table::fread(maf_file))
    


  if(coloc==TRUE){
    message(paste0(Sys.time(),": Running COLOC.."))
    coloc_df<-coloc_wrap(coloc_data,
      snplocs=snplocs,
      maf_df=maf_df,
      processed_gwas=processed_gwas,
      GWAS_type=GWAS_type,
      indiv_numbers=indiv_numbers,
      coloc.p12_prior=coloc.p12_prior)


  
    write.table(coloc_df,paste0(date,"_",run_name,"_COLOC_results.txt"))
  }
      
      
  if(hyprcoloc==TRUE){
    # message("Running hyprcoloc..")
    # hyprcoloc_df<-run_hyprcoloc(coloc_data,processed_gwas)



    # hyprcoloc_df<-hyprcoloc_df[grep("GWAS",hyprcoloc_df$traits),]
    # if(is.data.frame(hyprcoloc_df)==TRUE){
    #   hyprcoloc_df$index_snp<-do.call(rbind,strsplit(rownames(hyprcoloc_df),"\\."))[,1]
    #   hyprcoloc_df<-add_pvalues_hyprcoloc(hyprcoloc_df,coloc_data)
    #   write.table(hyprcoloc_df,paste0(date,"_",run_name,"_hyprcoloc_results.txt"))
    # }S
    
    message(paste0(Sys.time(),": Running hyprcoloc by gene.."))
    hyprcoloc_df_bygene=run_hyprcoloc_bygene(coloc_data,
    processed_gwas,prior.1 = hcoloc_prior.1,
    prior.c = hcoloc_prior.c,
    prior.12 = hcoloc_prior.12)

    message(paste0(Sys.time(),": Running hyprcoloc by celltype.."))
    hyprcoloc_df_bycelltype=run_hyprcoloc_bycelltype(coloc_data,
    processed_gwas,prior.1 = hcoloc_prior.1,
    prior.c = hcoloc_prior.c,
    prior.12 = hcoloc_prior.12)

    if(is.data.frame(hyprcoloc_df_bygene)==TRUE){
      write.table(hyprcoloc_df_bygene,paste0(date,"_",run_name,"_hyprcoloc_results_bygene.txt"))
    }

    if(is.data.frame(hyprcoloc_df_bycelltype)==TRUE){
      write.table(hyprcoloc_df_bycelltype,paste0(date,"_",run_name,"_hyprcoloc_results_bycelltype.txt"))
    }
  }
  
   end_time<-Sys.time()
  runtime=end_time-start_time
  message(paste0(Sys.time(),": Analysis complete. Runtime: ",runtime))
  








}



