
#test comment
#' @export

cellMR=function(regions,
                traits,
                run_name="test_preprocess",
                use_coloc_results=TRUE,
                allele_df_file="../../MatrixEQTL_IO/MAF_mat.csv",
                coloc_data_file=paste0("../../cellCOLOC/",run_name,"/coloc_data.rds"),
                coloc_df_file="2022-06-23_AD_2022_COLOC_results.txt",
                processed_GWAS_file=paste0("../../cellCOLOC/",run_name,"/GWAS_clumped.rds"),
                eqtl_FDR_cutoff=0.2,
                coloc_cutoff=0.5,
                prune_snps_r_cutoff=0.01,
                prune_local=FALSE,
               save_MR_input=FALSE
               ){
  start_time<-Sys.time()
  date=Sys.Date()
  # cellMRdir=paste0(outdir,"/cellMR/")
  # mr_outdir=paste0(outdir,"/cellMR/",run_name)
  # dir.create(cellMRdir,showWarnings=FALSE)
  # dir.create(mr_outdir,showWarnings=FALSE)
  # setwd(mr_outdir)
    
    
   
   
    
    if(file.exists(coloc_data_file)==TRUE && file.exists(processed_GWAS_file)==TRUE){
        coloc_data<-readRDS(coloc_data_file)
        processed_gwas<-readRDS(processed_GWAS_file)
        processed_gwas<-split(processed_gwas,processed_gwas$signif_snp_region)
    }else{
        message("No coloc input data found. Did you run cellCOLOC?")
    }
    
    if(use_coloc_results==TRUE){
        message("use_coloc_results=TRUE. Reading in latest COLOC results..")
        # coloc_df_list<-list.files(paste0("../../cellCOLOC/",run_name,"/"),pattern="COLOC_results.txt",full.names=T)
        # coloc_df_file<-dplyr::last(coloc_df_list)
        coloc_df<-read.table(coloc_df_file)
        coloc_df<-coloc_df[coloc_df$PP.H4.abf>coloc_cutoff,]
        regions<-coloc_df$region
        # traits<-rownames(coloc_df)
        traits=coloc_df$trait
        if(nrow(coloc_df)==0){
          stop("No colocs above PP.H4.abf threshold")
        }
    }

    #filter by regions so we don't harmonize everything
    coloc_data=lapply(coloc_data,function(x){
      x=x[regions]
      return(x)
    })

    processed_gwas=processed_gwas[regions]


    
    message("Harmonizing alleles..")
    allele_df<-data.table::fread(allele_df_file)
    coloc_data<-harmonize_alleles_wrap(coloc_data,processed_gwas,allele_df)
      
    
    message("Pruning eQTLs based on LD..")
    MR_regions<-vector("list",length(regions))
    names(MR_regions)<-regions
    
    
    message(paste0("Filtering out eQTLs with FDR below ",eqtl_FDR_cutoff))
    
    
    for(i in 1:length(regions)){
        MR_regions[[i]]<-prune_snps(coloc_data=coloc_data,
          region=regions[i],
          r_cutoff=prune_snps_r_cutoff,local=prune_local,
        eqtl_FDR_cutoff=eqtl_FDR_cutoff,trait=traits[i])

    }
    traits<-traits[as.numeric(which(unlist(lapply(MR_regions,is.null))==FALSE))]
    MR_regions<-Filter(Negate(is.null), MR_regions)


  if(length(MR_regions)==0){
    message("No regions remain after eQTL filtering.")
  }else{
    message("Running MR...")
    MR_res_list<-list()
    for(i in 1:length(MR_regions)){
      MR_res_list[[i]]<-run_MR(MR_regions[[i]],
      trait=traits[i],
      processed_gwas=processed_gwas,save_input=save_MR_input)
    }
    names(MR_res_list)<-paste0(names(MR_regions),"_",traits)
    df<-do.call(rbind,MR_res_list)
    df$gene<-as.data.frame(do.call(rbind,strsplit(rownames(df),"[.]")))$V2
    df$celltype<-as.data.frame(do.call(rbind,strsplit(rownames(df),"_")))$V2
    write.table(df,paste0(date,"_",run_name,"_MR_res.txt"))

  }

  
  end_time<-Sys.time()
  runtime=end_time-start_time
  message(paste0("Analysis complete. Runtime: ",runtime))


    
    
}