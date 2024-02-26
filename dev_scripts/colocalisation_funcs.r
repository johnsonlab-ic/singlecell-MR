

#' This function harmonizes alleles to the GWAS of reference - if alleles don't match, flip the beta of the eQTL
#'
#' @export


harmonize_alleles=function(gwas_region,beta_df,direction_df,alleles_df){

  if(length(grep("Allele1",colnames(gwas_region)))==1 & length(grep("Allele2",colnames(gwas_region)))==1){
      names(gwas_region)[which(names(gwas_region)=="Allele1")]<-"A1"
      names(gwas_region)[which(names(gwas_region)=="Allele2")]<-"A2"
    }

    if(stringr::str_detect(gwas_region$A1,"[[:upper:]]")[1]==FALSE & stringr::str_detect(gwas_region$A2,"[[:upper:]]")[1]==FALSE){
      gwas_region$A1<-toupper(gwas_region$A1)
      gwas_region$A2<-toupper(gwas_region$A2)
    }

    gwas_region<-gwas_region[match(rownames(beta_df),gwas_region$rsid),]
    alleles_df<-alleles_df[match(rownames(beta_df),alleles_df$snp),]
    # direction_df<-direction_df[,trait]


    #in MatrixEQTL, the "ref" allele is the effect allele.
    beta_df_A1<-alleles_df$ref
    beta_df_A2<-alleles_df$alt

    #store GWAS values so they dont get flipped
    beta_df_GWAS<-beta_df$GWAS
    beta_df$GWAS<-NULL




    #now check if alleles match. If not, the whole row (and its betas) are flipped
    # by *(-1)

    counts<-0
    status<-vector()
    for(i in 1:nrow(beta_df)){
      if(beta_df_A1[i]!=gwas_region$A1[i]&&beta_df_A1[i]==gwas_region$A2[i]){

        beta_df[i,]<-as.numeric(beta_df[i,])*(-1)

        tmp_directiondf<-direction_df[i,]

        for(b in 1:ncol(tmp_directiondf)){
          if(tmp_directiondf[1,b]=="-"){
          tmp_directiondf[1,b]<-"+"
            }else if(tmp_directiondf[1,b]=="+"){
          tmp_directiondf[1,b]<-"-"
            }
        }
        direction_df[i,]<-tmp_directiondf
      }else{
          status[i]<-"Same"
          }
        }

     beta_df$GWAS<-beta_df_GWAS

    # message(paste0("A total of ",counts," beta effects were flipped"))

    res_list<-list(beta_df,direction_df)
    names(res_list)<-c("betas","directionality")
    return(res_list)




}

#' @export

harmonize_alleles_wrap=function(coloc_data,processed_gwas,alleles_df){

  betas<-coloc_data$betas
  directionality<-coloc_data$directionality

    for(i in 1:length(betas)){

      res<-harmonize_alleles(gwas_region=processed_gwas[[i]],
      beta_df=betas[[i]],
      direction_df=directionality[[i]],
      alleles_df=alleles_df)


      betas[[i]]<-res$betas
      directionality[[i]]<-res$directionality
    }
    coloc_data$betas<-betas
    coloc_data$directionality<-directionality
    return(coloc_data)
}


#' @export
prune_snps=function(coloc_data,
  region,
  r_cutoff=0.01,
  eqtl_FDR_cutoff=1,
  trait,
  local=FALSE,
  plink_bin="/rds/general/user/ah3918/projects/single_cell_eqtl/live/MRC/TMP/plink",
  path_to_binaries="/rds/general/user/ah3918/projects/single_cell_eqtl/live/MRC/REFERENCE_DATASETS/clumping_ref/EUR/EUR"){

      biglist_beta<-coloc_data$beta[region][[1]]
      biglist_se<-coloc_data$std.errors[region][[1]]
      biglist_pval<-coloc_data$pvalues[region][[1]]
      biglist_FDR<-coloc_data$FDR[region][[1]]
      biglist_directionality<-coloc_data$directionality[region][[1]]


      fdrs<-biglist_FDR[,trait]
      fdrs<-fdrs[which(fdrs<eqtl_FDR_cutoff)]


      if(length(fdrs)==0){
        return(NULL)
      }else{



      df<-data.frame(rsid=rownames(biglist_beta),pval=biglist_FDR[,trait])

      
      suppressMessages(capture.output(snps_tokeep<-ieugwasr::ld_clump(df,clump_r2=r_cutoff,clump_p=eqtl_FDR_cutoff,
                   clump_kb=1e6),file=nullfile()))
      if(local==TRUE){
          suppressMessages(capture.output(snps_tokeep<-ieugwasr::ld_clump_local(df,clump_r2=r_cutoff,clump_p=eqtl_FDR_cutoff,
                   clump_kb=1e6,bfile=path_to_binaries,plink_bin=plink_bin),file=nullfile()))
      }

      tmp<-list(biglist_beta,biglist_se,biglist_pval,biglist_FDR,biglist_directionality)
      names(tmp)<-names(coloc_data)
      tmp<-lapply(tmp,function(x){
        x<-x[match(snps_tokeep$rsid,rownames(x)),]
        return(x)
      })

      test<-tmp$pvalues
      test<-test[test$GWAS<1,]
      snpstokeep<-rownames(test)

      tmp<-lapply(tmp,function(x){
      x<-x[snpstokeep,]
      return(x)})


      return(tmp)
    }
}

#' @export
filter_eqtls_MR=function(pruned_region,trait,eqtl_pval_cutoff){


  # message(paste0("Filtering out eQTLs with p-value below ",eqtl_pval_cutoff))

  pvals<-pruned_region$FDR
  snpstokeep<-rownames(pvals)
  pvals<-pvals[,trait]
  snpstokeep<-snpstokeep[which(pvals<eqtl_pval_cutoff)]

  pruned_region$betas<-pruned_region$betas[snpstokeep,]
  pruned_region$std.errors<-pruned_region$std.errors[snpstokeep,]
  pruned_region$FDR<-pruned_region$FDR[snpstokeep,]
  pruned_region$pvalues<-pruned_region$pvalues[snpstokeep,]
  pruned_region$directionality<-pruned_region$directionality[snpstokeep,]


  if(nrow(pruned_region$pvalues)==0){
      return(NULL)
  } else{
  return(pruned_region)
    }
}
#' @export
run_MR=function(pruned_region,trait,processed_gwas,save_input=FALSE){

  # nsnps<-nrow(pruned_region$betas)
  #  message(paste0(nsnps, " snps selected for MR."))
   gwas_beta<-as.data.frame(pruned_region$betas)$GWAS
   gwas_se<-as.data.frame(pruned_region$std.errors)$GWAS

   betas<-as.data.frame(pruned_region$betas)[,trait]
   ses<-as.data.frame(pruned_region$std.errors)[,trait]
   pval<-as.data.frame(pruned_region$pvalues)[,trait]
   FDR<-as.data.frame(pruned_region$FDR)[,trait]
   directionality<-as.data.frame(pruned_region$directionality)[,trait]

   MRInputObject<-MendelianRandomization::mr_input(bx=betas,bxse=ses,by=gwas_beta,byse=gwas_se)
   if(save_input==TRUE){
    dir.create("MR_inputs",showWarnings=FALSE)
    saveRDS(MRInputObject,paste0("MR_inputs/",trait,"_mr_input.rds"))
   }

   res<-MendelianRandomization::mr_allmethods(MRInputObject, method = "ivw")@Values
   IVW_beta<-res$Estimate[1]
   res<-as.data.frame(t(res))
   colnames(res)<-res["Method",]
   res<-res["P-value",]
  #  if(length(betas)==1){
  #    if(betas>0){
  #        eqtl_direction<-"+"
  #        }else{
  #        eqtl_direction<-"-"
  #        }
  #    if(gwas_beta>0){
  #        gwas_direction<-"+"
  #        }else{
  #        gwas_direction<-"-"
  #        }
  #    }else if(length(betas)>1){
  #    if(all(betas>0)){
  #        eqtl_direction<-"+"
  #        }else if(all(betas<0)){
  #        eqtl_direction<-"-"
  #        }else{
  #        eqtl_direction<-"NA"
  #        }

  #        if(all(gwas_beta>0)){
  #        gwas_direction<-"+"
  #        }else if(all(gwas_beta<0)){
  #        gwas_direction<-"-"
  #        }else{
  #        gwas_direction<-"NA"
  #        }
  #            }
  #  #
  #  direction<-paste0("eQTL: ",directionality[1]," -> GWAS: ",gwas_direction)


   rownames(res)<-trait
  #  res$direction<-direction
   res$IVW_beta<-IVW_beta
   res$IVs<-paste(rownames(pruned_region$betas),collapse=", ")
   res$eQTL_pval<-min(pval)
   res$eQTL_FDR<-min(FDR)

   #add gwas alleles
   if(class(processed_gwas)=="list"){
     processed_gwas<-do.call(rbind,processed_gwas)
   }

   a1<-processed_gwas[match(res$IVs,processed_gwas$rsid),]$A1
   a2<-processed_gwas[match(res$IVs,processed_gwas$rsid),]$A2

   res$effect_allele<-a1
   res$ref_allele<-a2

   return(res)
}

###########################################################
#################### COLOCALISATION #######################
###########################################################
# make sure pvalue column is named "pval" and variant name is "rsid"

##to get directionality, make sure this is BEFORE harmonization of alelles;
#this is because not all alleles will be flipped

#' @export
get_directionality=function(coloc_df,beta_df){

    directionality_vector<-vector()
    for(i in 1:nrow(coloc_df)){

        trait=rownames(coloc_df[i,])
        snp=coloc_df$lead_snp[i]
        snp_index<-which(rownames(beta_df)==snp)

        eqtl_beta<-beta_df[snp_index,trait]
        gwas_beta<-beta_df[snp_index,"GWAS"]

        if(eqtl_beta>0){
            gene="gene:---"
        }else{
            gene="gene:+++"
        }

        if(gwas_beta>0){
            risk="risk:+++"
        }else{
            risk="risk:---"
        }

        directionality=paste0(gene," -> ",risk)
        directionality_vector[i]<-directionality


    }
    return(directionality_vector)

}

#' @export
coloc_function=function(coloc.p12_prior,biglist_se_df,
  biglist_beta_df,
  biglist_pval_df,
  biglist_FDR_df,
  processed_gwas,
  GWAS_type=c("quant","cc"),
  maf_df,
  region,
  indiv_numbers,
  snplocs){


    if(class(processed_gwas)=="list"){
      processed_gwas<-do.call(rbind,processed_gwas)
    }

    se_df<-as.data.frame(biglist_se_df)
    beta_df<-as.data.frame(biglist_beta_df)
    pval_df<-as.data.frame(biglist_pval_df)
    FDR_df<-as.data.frame(biglist_FDR_df)

    se_df$MAF<-maf_df$maf[match(rownames(se_df),maf_df$snp)]
    se_df$MAF<-as.numeric(se_df$MAF)


    gwas<-processed_gwas
    GWAS_betas<-beta_df$GWAS
    GWAS_se<-se_df$GWAS

    #𝑣𝑎𝑟(𝛽)  is the standard error of 𝛽, squared - COLOC tutorial https://chr1swallace.github.io/coloc/articles/a02_data.html
    
    GWAS_varbeta<-GWAS_se^2
    GWAS_pval<-pval_df$GWAS
    if(length(grep("MAF",colnames(gwas))==1)){
      GWAS_MAF<-gwas$MAF[match(rownames(beta_df),gwas$rsid)]
    } else{
      GWAS_MAF<-se_df$MAF
    }

    GWAS_SNP<-as.character(row.names(beta_df))
    GWAS_position<-snplocs$position[match(GWAS_SNP,snplocs$annot)]

    coloc_trait2<-list(beta=GWAS_betas,snp=GWAS_SNP,
    position=GWAS_position,varbeta=GWAS_varbeta,
    type=GWAS_type,MAF=GWAS_MAF,
    pvalues=GWAS_pval)

    if(GWAS_type=="quant"){
      if(any(colnames(processed_gwas)=="N_total")){
        coloc_trait2[["N"]]<-round(mean(processed_gwas$N_total))
        # message(paste0("GWAS study of ",processed_gwas$N[1]," individuals"))
      }else{
        coloc_trait2[["N"]]<-50000
      }
    }else if(GWAS_type=="cc"){
      if(any(colnames(processed_gwas)=="N_total")){
        coloc_trait2[["N"]]<-round(mean(processed_gwas$N_total))
        # message(paste0("GWAS study of ",processed_gwas$N[1]," individuals"))
      }else{
        coloc_trait2[["N"]]<-50000
      }

      coloc_trait2[["s"]]<-mean(processed_gwas$case.prop)
    }






    #
    # if(GWAS_type=="cc"){
    #   coloc_trait2$type="cc"
    #   coloc_trait2$s<-mean(processed_gwas$CASES/processed_gwas$CONTROLS)
    # }
    beta_df$GWAS<-NULL


    indiv_numbers<-indiv_numbers

    results<-list()
    snp_list<-list()
    for(i in 1:ncol(beta_df)){

      ####!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      #this part is PROJECT SPECIFIC. its a few lines to properly add the number of
      trait<-colnames(beta_df[i])
      n_indiv_index<-which(!is.na(pmatch(rownames(indiv_numbers),trait)))
      n_indiv<-as.numeric(indiv_numbers[n_indiv_index,])
      # message(paste0(n_indiv, " individuals for the estimation"))

      astro_betas<-beta_df[,i]


      ## eQTL is always quant
      astro_varbeta<-se_df[,i]^2
      astro_pval<-beta_df[,i]
      astro_MAF<-se_df$MAF
      astro_SNP<-as.character(row.names(se_df))
      astro_position<-GWAS_position
      coloc_trait1<-list(beta=astro_betas,snp=astro_SNP,varbeta=astro_varbeta,
        position=astro_position,
        type="quant",MAF=astro_MAF,
        pvalues=astro_pval,N=n_indiv)
        options(warn=-1)
        saveRDS(coloc_trait1,"coloc_trait1.rds")
        saveRDS(coloc_trait2,"coloc_trait2.rds")
        capture.output(res<-suppressMessages(coloc::coloc.abf(coloc_trait1,coloc_trait2)),file=nullfile())
        results[[i]]<-res$summary
        snp_list[[i]]<-res$results

    }


    results<-do.call(rbind,results)
    rownames(results)<-colnames(beta_df)

    #add lead snp
    snp_list<-lapply(snp_list,function(x){x<-x[order(x$SNP.PP.H4,decreasing=T),]
                                      x<-x[1,c("snp","SNP.PP.H4")]
                                      return(x)})
    snp_list<-do.call(rbind,snp_list)
    results<-cbind(results,snp_list)


    #add eQTLs
    eqtl_pvals<-list()
    for(i in 1:nrow(results)){
        eqtl_pvals[[i]]<-pval_df[results$snp[i],rownames(results[i,])]
        }

        eqtl_pvals<-unlist(eqtl_pvals)
        results$eQTL_pval<-eqtl_pvals
    FDR_pvals<-list()
    for(i in 1:nrow(results)){
        FDR_pvals[[i]]<-FDR_df[results$snp[i],rownames(results[i,])]
        }
    FDR_pvals<-unlist(FDR_pvals)
    results$eQTL_FDR<-FDR_pvals

    #add GWAS pvalue
    results$GWAS_pval<-processed_gwas$pval[match(results$snp,processed_gwas$rsid)]

    names(results)[names(results)=="snp"]<-"lead_snp"
    results$region<-rep(region,nrow(results))

    beta_df$GWAS<-GWAS_betas
    results$directionality<-get_directionality(results,beta_df)
    return(results)








}
#' @export
add_directionality=function(coloc_data){

    get_dir<-function(column){
    sapply(column,function(x){

        #beta>0 is negative expression because matrixeqtl is flipped (ref allele additive)
     if(x>0){
        return("+")
        }else{
        return("-")
        }
    })
    }
    betas<-coloc_data$betas
    directionality<-lapply(betas,function(x){
        x<-x %>% dplyr::mutate_all(dplyr::funs(get_dir))
        return(x)
        })

    coloc_data$directionality<-directionality
    return(coloc_data)
}

#' @export

coloc_wrap=function(coloc_data,snplocs,maf_df,processed_gwas,coloc.p12_prior,GWAS_type=c("quant","cc"),indiv_numbers){

  maf_df<-maf_df
  snplocs<-snplocs
  gwas<-processed_gwas
  biglist_beta<-coloc_data$beta
  biglist_se<-coloc_data$std.errors
  biglist_pval<-coloc_data$pvalues
  biglist_FDR<-coloc_data$FDR

  full_coloc_res<-list()
  for(i in 1:length(biglist_se)){
      full_coloc_res[[i]]<-coloc_function(coloc.p12_prior=coloc.p12_prior,biglist_se_df=biglist_se[[i]],
             biglist_beta_df=biglist_beta[[i]],
             biglist_pval_df=biglist_pval[[i]],
              biglist_FDR_df=biglist_FDR[[i]],
              region=names(biglist_se[i]),
             processed_gwas=processed_gwas,
             GWAS_type=GWAS_type,
             maf_df=maf_df,snplocs=snplocs,
           indiv_numbers=indiv_numbers)


           ###split genes

           df=full_coloc_res[[i]]
           splitdf=as.data.frame(do.call(rbind,strsplit(rownames(df),"[.]")))
           df$gene=splitdf$V2
           df$celltype=splitdf$V1
           df$trait=paste0(df$celltype,".",df$gene)
           full_coloc_res[[i]]=df


  }

  finaldf<-do.call(rbind,full_coloc_res)
  finaldf<-finaldf[order(finaldf$PP.H4.abf,decreasing=T),]
  return(finaldf)


}
#' @export
select_regions=function(gwas,
                        pval,
                        window,
                        parallel=FALSE,
                        local=FALSE,
                        plink_bin=genetics.binaRies::get_plink_binary(),
                        path_to_binaries="/rds/general/user/ah3918/projects/single_cell_eqtl/live/MRC/TMP/plink"){

    #this is the local version
  #   error_catch_clump=function(gwas,clump_kb,
  #     clump_r2,clump_p,plink_bin,bfile){

  #     tryCatch({test<-ld_clump_local(gwas,clump_kb=clump_kb,clump_r2=clump_r2,clump_p=clump_p,
  #     plink_bin=plink_bin,bfile=bfile);},error=function(e){test<<-NA});
  #     return(test)

  #   }

    options(warn=-1)
  error_catch_clump=function(gwas,clump_kb,
    clump_r2,clump_p,plink_bin,bfile){

    tryCatch({test<-suppressMessages(ld_clump(gwas,clump_kb=clump_kb,clump_r2=clump_r2,clump_p=clump_p,
    plink_bin=plink_bin,bfile=bfile));},error=function(e){test<<-NA});
    return(test)

  }
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
    if(parallel==TRUE){
      gwas_clump<-mclapply(gwas_chr_list,error_catch_clump,ld_clump_local,clump_kb=kb_window,clump_r2=0,clump_p=pval,plink_bin=plink_bin,
      bfile=path_to_binaries,mc.cores=parallel.cores)
    }else{
      gwas_clump<-suppressMessages(lapply(gwas_chr_list,error_catch_clump,clump_kb=kb_window,clump_r2=0,clump_p=pval,plink_bin=plink_bin,
      bfile=path_to_binaries))

    }
  } else{

    if(parallel==TRUE){
      gwas_clump<-mclapply(gwas_chr_list,ld_clump,clump_kb=kb_window,clump_r2=0,clump_p=pval,mc.cores=parallel.cores)
    }else{
      gwas_clump<-suppressMessages(lapply(gwas_chr_list,ieugwasr::ld_clump,clump_kb=kb_window,clump_r2=0,clump_p=pval))

    }

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
#' @export
createdf=function(x,uniquegene){

  generows<-which(x$gene==as.character(uniquegene))

  df<-data.frame(SNP=as.character(x$SNP[generows]),
  beta=x$beta[generows],
  se=x$std_error[generows],
  GWAS_se=x$GWAS_se[generows],
  GWAS_beta=x$GWAS_beta[generows])
  return(df)

}

#' @export
appendGWAS_beta=function(x,GWAS){

  x$GWAS<-GWAS$b[match(rownames(x),GWAS$rsid)]
  return(x)
}

#' @export
add_pvalues_hyprcoloc=function(hyprcoloc_df,coloc_data){

  final<-hyprcoloc_df

  # add p-values from trait 1
  traits<-final$traits
  traits<-strsplit(traits,",")
  traits<-lapply(traits,function(x){
  x<-x[1]
  return(x)})


  pval<-list()
  FDR<-list()
  candidate_snps<-final$candidate_snp
  index_snps<-final$index_snp

  biglist_pval<-coloc_data$pvalues
  biglist_FDR<-coloc_data$FDR
  for(i in 1:length(index_snps)){

    tmpdf<-biglist_pval[[index_snps[i]]]
    pval[i]<-tmpdf[candidate_snps[i],traits[[i]]]
    tmpdf<-biglist_FDR[[index_snps[i]]]
    FDR[i]<-tmpdf[candidate_snps[i],traits[[i]]]
  }

  pval<-unlist(pval)
  FDR<-unlist(FDR)

  final$eqtl_pval_first_trait<-pval
  final$FDR<-FDR
  return(final)


}

#' @export
appendGWAS_se=function(x,GWAS){

  x$GWAS<-GWAS$se[match(rownames(x),GWAS$rsid)]
  return(x)
}
#' @export
appendGWAS_pval=function(x,GWAS){

  x$GWAS<-GWAS$pval[match(rownames(x),GWAS$rsid)]
  return(x)
}
#' @export
append_gwas_info=function(hyprcoloc_res,GWAS){

    x<-hyprcoloc_res
     x$index_snp<-GWAS$signif_snp_region[match(x$candidate_snp,GWAS$rsid)]
     x$n_snps_region_gwas<-GWAS$n_snps_region[match(x$candidate_snp,GWAS$rsid)]
    x$GWAS_pval<-GWAS$pval[match(x$candidate_snp,GWAS$rsid)]
    return(x)
}
#' @export
preprocess_mateqtlouts=function(mateqtlouts,GWAS){
    rowsbefore<-nrow(mateqtlouts[[1]])
    message(paste0(Sys.time(),": Beginning with pre-processing. Assuming all cell types share same number of cis-EQTLs, starting with ",rowsbefore, " cis-EQTLs"))
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
    x$std_error=x$beta/x$t.stat
    return(x)})

      message(paste0(Sys.time(),": Filtering outputs to only contain associations present in all cell-types"))
  commonrows<-Reduce(intersect,lapply(mateqtlouts,row.names))

  mateqtlouts<-lapply(mateqtlouts,function(x){x<-x[commonrows,]
    return(x)})
  rowsdiff=rowsbefore-nrow(mateqtlouts[[1]])




  message(paste0(Sys.time(),": Intersection done. ",rowsdiff, " cis-EQTLs removed from analysis."))
  return(mateqtlouts)
}#' @export
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
#' @export
prep_coloc_data=function(mateqtlouts,processed_gwas,gene_locs){

  # library(dplyr)
  # library(tidyr)

  if(class(processed_gwas)=="data.frame"){
    processed_gwas<-split(processed_gwas,processed_gwas$signif_snp_region)
  }
  biglist<-list()

  for(i in 1:length(processed_gwas)){

    # message(paste0("Selecting eQTLs for region ",i))
    tmp_gwas<-processed_gwas[[i]]
    genestokeep<-get_genes_per_region(tmp_gwas,gene_locs)

    biglist[[i]]<-lapply(mateqtlouts,function(x){
      x<-x[x$gene %in% genestokeep,]
      x<-x[x$SNP %in% tmp_gwas$rsid,]
      return(x)
    })
  }
  names(biglist)<-names(processed_gwas)


  biglist_se<-list()
  for(i in 1:length(biglist)){
    biglist_se[[i]]<-lapply(biglist[[i]],function(x){
          x<-x[c("SNP","gene","std_error")]
          x<-x%>% tidyr::pivot_wider(names_from = gene, values_from = std_error)
          x<-as.data.frame(x)
          rownames(x)<-x$SNP
          x$SNP<-NULL
          return(x)})
        }

  biglist_beta<-list()
  for(i in 1:length(biglist)){
    biglist_beta[[i]]<-lapply(biglist[[i]],function(x){
      x<-x[c("SNP","gene","beta")]
      x<-x%>% tidyr::pivot_wider(names_from = gene, values_from = beta)
      x<-as.data.frame(x)
      rownames(x)<-x$SNP
      x$SNP<-NULL
      return(x)})
    }

  biglist_pval<-list()
  for(i in 1:length(biglist)){
        biglist_pval[[i]]<-lapply(biglist[[i]],function(x){
        x<-x[c("SNP","gene","p.value")]
        x<-x%>% tidyr::pivot_wider(names_from = gene, values_from = p.value)
        x<-as.data.frame(x)
        rownames(x)<-x$SNP
        x$SNP<-NULL
        return(x)})
    }

  biglist_FDR<-list()
  for(i in 1:length(biglist)){
      biglist_FDR[[i]]<-lapply(biglist[[i]],function(x){
        x<-x[c("SNP","gene","FDR")]
        x<-x%>% tidyr::pivot_wider(names_from = gene, values_from = FDR)
        x<-as.data.frame(x)
        rownames(x)<-x$SNP
        x$SNP<-NULL
        return(x)})
    }


  for(i in 1:length(biglist_se)){
      biglist_se[[i]]<-do.call(cbind,biglist_se[[i]])
      biglist_beta[[i]]<-do.call(cbind,biglist_beta[[i]])
      biglist_pval[[i]]<-do.call(cbind,biglist_pval[[i]])
      biglist_FDR[[i]]<-do.call(cbind,biglist_FDR[[i]])
    }

  biglist_se<-lapply(biglist_se,function(x){
      if(length(unique(colnames(x)))==1){
          names<-paste0(names(mateqtlouts),".",colnames(x))
          colnames(x)<-names
          }
      return(x)
      })
  biglist_beta<-lapply(biglist_beta,function(x){
          if(length(unique(colnames(x)))==1){
              names<-paste0(names(mateqtlouts),".",colnames(x))
              colnames(x)<-names
              }
          return(x)
          })
  biglist_pval<-lapply(biglist_pval,function(x){
              if(length(unique(colnames(x)))==1){
                  names<-paste0(names(mateqtlouts),".",colnames(x))
                  colnames(x)<-names
                  }
              return(x)
              })
  biglist_FDR<-lapply(biglist_FDR,function(x){
                if(length(unique(colnames(x)))==1){
                    names<-paste0(names(mateqtlouts),".",colnames(x))
                    colnames(x)<-names
                    }
                return(x)
                })


  names(biglist_se)<-names(processed_gwas)
  names(biglist_beta)<-names(processed_gwas)
  names(biglist_pval)<-names(processed_gwas)
  names(biglist_FDR)<-names(processed_gwas)

  biglist_pval<-biglist_pval[sapply(biglist_pval, nrow)>0]
  biglist_FDR<-biglist_FDR[sapply(biglist_FDR, nrow)>0]
  biglist_se<-biglist_se[sapply(biglist_se, nrow)>0]
  biglist_beta<-biglist_beta[sapply(biglist_beta, nrow)>0]

  if(class(processed_gwas)=="list"){
    processed_gwas<-do.call(rbind,processed_gwas)
  }

  biglist_beta<-lapply(biglist_beta,appendGWAS_beta,GWAS=processed_gwas)
  biglist_se<-lapply(biglist_se,appendGWAS_se,GWAS=processed_gwas)
  biglist_pval<-lapply(biglist_pval,appendGWAS_pval,GWAS=processed_gwas)

  final_list<-list(biglist_beta,biglist_se,biglist_pval,biglist_FDR)
  names(final_list)<-c("betas","std.errors","pvalues","FDR")
  return(final_list)
}
#' @export
run_hyprcoloc=function(coloc_data,processed_gwas){

  biglist_beta<-coloc_data$beta
  biglist_se<-coloc_data$std.errors
  biglist_pval<-coloc_data$pvalues
  biglist_FDR<-coloc_data$FDR

  biglist_beta<-lapply(biglist_beta,as.matrix)
  biglist_se<-lapply(biglist_se,as.matrix)


  traits_list<-list()
  snp_id<-list()

  for(i in 1:length(biglist_beta)){
      traits_list[[i]]<-colnames(biglist_beta[[i]])
      snp_id[[i]]<-rownames(biglist_beta[[i]])
      }
  hyprcoloc_results<-list()
  for(i in 1:length(biglist_beta)){
        traits_list[[i]]<-colnames(biglist_beta[[i]])
        snp_id[[i]]<-rownames(biglist_beta[[i]])
        hyprcoloc_results[[i]]<-hyprcoloc::hyprcoloc(biglist_beta[[i]],
          biglist_se[[i]], trait.names=traits_list[[i]], snp.id=snp_id[[i]],snpscores=TRUE)
      }

  names(hyprcoloc_results)<-names(biglist_beta)
  saveRDS(hyprcoloc_results,"fullhyprcoloc_res.rds")

  final<-lapply(hyprcoloc_results,function(x){
    x<-as.data.frame(x$results)
    x<-x[order(x$posterior_prob,decreasing=T),]
    x<-x[which(x$posterior_prob>0),]
    return(x)})

  final<-final[sapply(final, nrow)>0]
  if(length(final)==0){
    message("No colocalisations found.")

  }else{

    final<-do.call(rbind,final)

    if(class(processed_gwas)=="list"){
      processed_gwas<-do.call(rbind,processed_gwas)
    }
    final<-append_gwas_info(final,processed_gwas)
    final<-final[c("traits","regional_prob","posterior_prob","posterior_explained_by_snp",
                 "candidate_snp","GWAS_pval","n_snps_region_gwas","index_snp")]
    final<-final[order(final$GWAS_pval),]

    #add n_snps
    index_snps<-final$index_snp
    n_snps_mateqtls<-vector()
    for(i in 1:length(index_snps)){
      n_snps_mateqtls[i]<-nrow(biglist_se[[index_snps[i]]])
    }
    final$n_snps_mateqtls<-n_snps_mateqtls
    #
    #
    return(final)
  }

}

#' @export
run_hyprcoloc_bycelltype=function(coloc_data,processed_gwas, prior.1 = 1e-04,
  prior.c = 0.02,
  prior.12 = NULL){
    


    biglist_beta<-coloc_data$betas
      biglist_se<-coloc_data$std.errors
      biglist_pval<-coloc_data$pvalues
      biglist_FDR<-coloc_data$FDR
    full_res=list()
    for(region in 1:length(biglist_beta)){


        beta_df=biglist_beta[[region]]
        se_df=biglist_se[[region]]
        pval_df=biglist_pval[[region]]
        FDR_df=biglist_FDR[[region]]

        #store GWAS values and remove to extract colnames
        gwas_beta=beta_df$GWAS
        beta_df$GWAS=NULL
        genes=colnames(beta_df)
        genes=as.data.frame(do.call(rbind,strsplit(genes,"[.]")))$V1


        genes=unique(genes)
        genes
        traits_list<-list()
        snp_id<-list()
        hyprcoloc_results<-list()

        for(i in 1:length(genes)){

            tmp_beta_df=beta_df[,grep(genes[i],colnames(beta_df)),drop=FALSE]
            #reattach GWAS beta
            tmp_beta_df$GWAS=gwas_beta

            ##using grep again, but attaching GWAS col thats already there (ncol(se_df))
            tmp_se_df=se_df[,c(grep(genes[i],colnames(se_df)),ncol(se_df))]

            #as matrix
            tmp_beta_df=as.matrix(tmp_beta_df)
            tmp_se_df=as.matrix(tmp_se_df)

            traits_list[[i]]<-colnames(tmp_beta_df)
            snps<-rownames(tmp_beta_df)
             hyprcoloc_results[[i]]<-hyprcoloc::hyprcoloc(tmp_beta_df,
                  tmp_se_df, trait.names=traits_list[[i]], snp.id=snps,snpscores=TRUE, prior.1 = prior.1,
  prior.c = prior.c,
  prior.12 = prior.12)
        }

        # hyprcoloc_results<-list()
        # for(i in 1:length(biglist_beta)){
        #         traits_list[[i]]<-colnames(biglist_beta[[i]])
        #         snp_id[[i]]<-rownames(biglist_beta[[i]])
        #         hyprcoloc_results[[i]]<-hyprcoloc(biglist_beta[[i]],
        #           biglist_se[[i]], trait.names=traits_list[[i]], snp.id=snp_id[[i]],snpscores=TRUE)
        # }

        final<-lapply(hyprcoloc_results,function(x){
            x<-as.data.frame(x$results)
            x<-x[order(x$posterior_prob,decreasing=T),]
            x<-x[which(x$posterior_prob>0),]
            return(x)})

          final<-final[sapply(final, nrow)>0]
          if(length(final)==0){
        

          }else{

            final<-do.call(rbind,final)

            if(class(processed_gwas)=="list"){
              processed_gwas<-do.call(rbind,processed_gwas)
            }
            final<-append_gwas_info(final,processed_gwas)
            final<-final[c("traits","regional_prob","posterior_prob","posterior_explained_by_snp",
                         "candidate_snp","GWAS_pval","n_snps_region_gwas","index_snp")]
            final<-final[order(final$GWAS_pval),]

            #add n_snps
            index_snps<-final$index_snp
            n_snps_mateqtls<-vector()
              n_traits<-vector()
            
  #             final$n_snps_mateqtls<-n_snps_mateqtls
              final$n_traits<-length(traits_list[[1]])
        }

        full_res[[region]]<-final
        
        }

          if(length(full_res)==0){
            message("No colocalisations found.")

          }else{
            full_res<-do.call(rbind,full_res)
        }
    full_res<-full_res[grep("GWAS",full_res$traits),]
        return(full_res)
}
#' @export
run_hyprcoloc_bygene=function(coloc_data,processed_gwas, prior.1 = 1e-04,
  prior.c = 0.02,
  prior.12 = NULL){


    


    biglist_beta<-coloc_data$betas
      biglist_se<-coloc_data$std.errors
      biglist_pval<-coloc_data$pvalues
      biglist_FDR<-coloc_data$FDR
    full_res=list()
    for(region in 1:length(biglist_beta)){


        beta_df=biglist_beta[[region]]
        se_df=biglist_se[[region]]
        pval_df=biglist_pval[[region]]
        FDR_df=biglist_FDR[[region]]

        #store GWAS values and remove to extract colnames
        gwas_beta=beta_df$GWAS
        beta_df$GWAS=NULL
        genes=colnames(beta_df)
        genes=as.data.frame(do.call(rbind,strsplit(genes,"[.]")))$V2


        genes=unique(genes)
        genes
        traits_list<-list()
        snp_id<-list()
        hyprcoloc_results<-list()

        for(i in 1:length(genes)){

            tmp_beta_df=beta_df[,grep(genes[i],colnames(beta_df)),drop=FALSE]
            #reattach GWAS beta
            tmp_beta_df$GWAS=gwas_beta

            ##using grep again, but attaching GWAS col thats already there (ncol(se_df))
            tmp_se_df=se_df[,c(grep(genes[i],colnames(se_df)),ncol(se_df))]

            #as matrix
            tmp_beta_df=as.matrix(tmp_beta_df)
            tmp_se_df=as.matrix(tmp_se_df)

            traits_list[[i]]<-colnames(tmp_beta_df)
            snps<-rownames(tmp_beta_df)
             hyprcoloc_results[[i]]<-hyprcoloc::hyprcoloc(tmp_beta_df,
                  tmp_se_df, trait.names=traits_list[[i]],snp.id=snps,snpscores=TRUE, prior.1 = prior.1,
                  prior.c = prior.c,
                  prior.12 = prior.12)
        }

        # hyprcoloc_results<-list()
        # for(i in 1:length(biglist_beta)){
        #         traits_list[[i]]<-colnames(biglist_beta[[i]])
        #         snp_id[[i]]<-rownames(biglist_beta[[i]])
        #         hyprcoloc_results[[i]]<-hyprcoloc(biglist_beta[[i]],
        #           biglist_se[[i]], trait.names=traits_list[[i]], snp.id=snp_id[[i]],snpscores=TRUE)
        # }

        final<-lapply(hyprcoloc_results,function(x){
            x<-as.data.frame(x$results)
            x<-x[order(x$posterior_prob,decreasing=T),]
            x<-x[which(x$posterior_prob>0),]
            return(x)})

          final<-final[sapply(final, nrow)>0]
          if(length(final)==0){
        

          }else{

            final<-do.call(rbind,final)

            if(class(processed_gwas)=="list"){
              processed_gwas<-do.call(rbind,processed_gwas)
            }
            final<-append_gwas_info(final,processed_gwas)
            final<-final[c("traits","regional_prob","posterior_prob","posterior_explained_by_snp",
                         "candidate_snp","GWAS_pval","n_snps_region_gwas","index_snp")]
            final<-final[order(final$GWAS_pval),]

            #add n_snps
            index_snps<-final$index_snp
        
    #             final$n_snps_mateqtls<-n_snps_mateqtls
              final$n_traits<-length(traits_list[[1]])
        }

        full_res[[region]]<-final
        
        }

          if(length(full_res)==0){
            message("No colocalisations found.")

          }else{
            full_res<-do.call(rbind,full_res)
        }
    full_res<-full_res[grep("GWAS",full_res$traits),]
        return(full_res)
}


liftover_hg19_hg38=function(GWAS){
  # library(gdsfmt)
  # library(SNPRelate)
  # library(SeqArray)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
  library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
  library(org.Hs.eg.db)

  # message("Reading in GWAS..")
  # gwas_name<-GWAS
  #
  #
  # GWAS<-read.table(GWAS,header=TRUE,stringsAsFactors=F)
  #
  # # if(length(grep("beta",colnames(GWAS)))==0 && length(grep("OR",colnames(GWAS)))==1){
  # #   GWAS$beta<-log(GWAS$OR)
  # # }

  if(length(grep("chr",GWAS$chr))==0){
    GWAS$chr<-paste0("chr",GWAS$chr)
  }

  chromlocs<-cbind(GWAS$chr,GWAS$pos)
  colnames(chromlocs)<-c("Chr","Pos")
  chromlocs<-as.data.frame(chromlocs)





  # message("Creating GRanges object..")
  # chromlocs<-as.character(GWAS$SNP)
  # chromlocs<-strsplit(chromlocs,":")
  # chromlocs<-do.call(rbind.data.frame,chromlocs)

    chronly<-strsplit(as.character(chromlocs$Chr),"r")
    chronly<-do.call(rbind.data.frame,chronly)

   final<-data.frame(chrom=chronly[,2],position=chromlocs$Pos)
   final<-final[complete.cases(final),]

  final$paste<-paste0(final$chrom,":",final$position)
  GWAS_gr<-makeGRangesFromDataFrame(final,
  seqnames.field="chrom",
  start.field="position",
  end.field="position")

  seqlevelsStyle(GWAS_gr)="UCSC"

  message("Performing liftover..")

  chain<-import.chain("/mnt/data/ALEX/DUMP/Tools/hg19ToHg38.over.chain")
  # chain<-import.chain("/rds/general/user/ah3918/home/Tools/hg19ToHg38.over.chain")

  hg38GWAS_gr<-liftOver(GWAS_gr,chain)
  hg38GWAS_gr<-unlist(hg38GWAS_gr)

  seqlevelsStyle(GWAS_gr)="NCBI"
  seqlevelsStyle(hg38GWAS_gr)="NCBI"

  hg19snp<-SNPlocs.Hsapiens.dbSNP144.GRCh37
  hg19rsids<-snpsByOverlaps(hg19snp,GWAS_gr)

  hg38snp<-SNPlocs.Hsapiens.dbSNP151.GRCh38
  hg38rsids<-snpsByOverlaps(hg38snp,hg38GWAS_gr)
  #
  hg19_rsids_df<-data.frame(rsid=as.data.frame(hg19rsids)$RefSNP_id,
                         chr_pos=paste0("chr",as.data.frame(hg19rsids)$seqnames,":",as.data.frame(hg19rsids)$pos))
  #
  hg38_rsids_df<-data.frame(rsid=as.data.frame(hg38rsids)$RefSNP_id,
                         chr_pos=paste0("chr",as.data.frame(hg38rsids)$seqnames,":",as.data.frame(hg38rsids)$pos))
  #
  GWAS_df<-as.data.frame(GWAS_gr)
  hg38GWAS_df<-as.data.frame(hg38GWAS_gr)

  GWAS_df<-GWAS_df[rownames(GWAS_df) %in% rownames(hg38GWAS_df),]
  #
  conversiontable<-data.frame(hg19_pos=paste0("chr",GWAS_df$seqnames,":",GWAS_df$start),
                                       hg38_pos=paste0("chr",hg38GWAS_df$seqnames,":",hg38GWAS_df$start))


  conversiontable$hg38_chr_only<-hg38GWAS_df$seqnames
  conversiontable$hg38_pos_only<-hg38GWAS_df$start

  test<-hg38_rsids_df[match(conversiontable$hg38_pos,hg38_rsids_df$chr_pos),]
  conversiontable$rsid_hg38<-test$rsid
  conversiontable$hg19_pos<-as.character(conversiontable$hg19_pos)
  conversiontable$hg38_pos<-as.character(conversiontable$hg38_pos)
  conversiontable$rsid_hg38<-as.character(conversiontable$rsid_hg38)

  n_snps_before<-nrow(GWAS)
  message(paste0("Started with ",n_snps_before," variants."))

  GWAS<-GWAS[match(conversiontable$rsid_hg38,GWAS$rsid),]
  GWAS$rsid_hg38<-conversiontable$rsid_hg38
  GWAS$hg38_pos<-conversiontable$hg38_pos
  GWAS$hg38_chr_only<-conversiontable$hg38_chr_only
  GWAS$hg38_pos_only<-conversiontable$hg38_pos_only

  GWAS<-GWAS[complete.cases(GWAS),]

  n_snps_after<-nrow(GWAS)

  snps_lost<-n_snps_before - n_snps_after

  GWAS$pos<-GWAS$hg38_pos_only
  message(paste0(snps_lost," variants were lost from conversion."))

  return(GWAS)
}

#' @export
prep_munge_out=function(gwas_path,gwas_out,add_n=TRUE,n_total=100,n_cases=NULL){
    
    message("Reading in GWAS...")
    gwas=data.table::fread(gwas_path)
    
    if(add_n==TRUE){
       gwas$N_total=n_total
        
        if(is.numeric(n_cases)){
            gwas$N_cases=n_cases
            gwas$case.prop=gwas$N_cases/gwas$N_total
            }
        
    }
  
    message("Adding columns...")
    names(gwas)[which(names(gwas)=="P")]="pval"
    names(gwas)[which(names(gwas)=="BETA")]="b"
    names(gwas)[which(names(gwas)=="SE")]="se"
    names(gwas)[which(names(gwas)=="CHR")]="chr"
    names(gwas)[which(names(gwas)=="BP")]="pos"
    names(gwas)[which(names(gwas)=="SNP")]="rsid"

    
    message("Saving file...")
    
    data.table::fwrite(gwas,gwas_out)
    message(paste0("GWAS prepped. Saved under ",gwas_out))
}

#' @export
run_hyprcoloc_wrap=function(GWAS,
  mateqtlouts,
  gene_loc_file,
  maf_file,
  snp_loc_file,
  allele_df_file,
  indiv_numbers_file,
  preprocess=FALSE,
  intersected_mateqtlouts_file="GWAS_intersected_mateqtlouts.rds",
  preprocess_mateqtlouts=FALSE,
  processed_mateqtlouts_file="filtered_mateqtlouts.rds",
  preprocess_GWAS=FALSE,
  processed_GWAS_file="GWAS_processed_hyprcoloc.rds",
  GWASsignif=5e-8,
  GWAS_window=1e6,
  run_name=NULL,
  harmonize_alleles=FALSE,
  hyprcoloc=FALSE,
  coloc=TRUE,
  GWAS_type=c("quant","cc"),
  coloc.p12_prior=1e-5,
  MR=FALSE,
  save_MR_input=FALSE,
  coloc_cutoff=0.5,
  prune_snps_r_cutoff=0.01,
  # gwas_pval_cutoff=1,
  eqtl_FDR_cutoff=0.2,
  # eqtl_pval_cutoff=5e-8,
  outdir="/."){


  # The GWAS input must contain the following column names: pval, rsid, pos, chr, beta, se
  # The mateqtlouts input is the path to a list object (*.rds), where each element (named after cell type/experiment)
  # is a matrixeqtl output. This must contain the FULL list of MatrixEQTLs (all associations without thresholding FDR/p-value).
  # The function assumes NO FILTERING has been done.
  # In addition, the functions assume that the GWAS and mateqtlouts are on the same genome build and positions are shared.


  suppressMessages(library(SNPlocs.Hsapiens.dbSNP151.GRCh38))
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(GenomicFeatures))
  suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
  suppressMessages(library(SNPlocs.Hsapiens.dbSNP151.GRCh38))
  suppressMessages(library(org.Hs.eg.db))
  suppressMessages(library(dplyr))
  suppressMessages(library(tidyr))
  # library(hyprcoloc)
  suppressMessages(library(coloc))
  suppressMessages(library(stringr))
  suppressMessages(library(ieugwasr))
  suppressMessages(library(parallel))
  # library(MendelianRandomization)

  # -============================================================================================================
  # -=-=-=-=-=-=-=- STEP 1. Reading in data and filtering -=-=-=-=-=-=-===
  # -============================================================================================================

  start_time<-Sys.time()
  setwd(outdir)
  if(preprocess==TRUE){

      message("preprocess_mateqtlouts=TRUE. Reading in and filtering mateqtlouts..")
      mateqtlouts<-readRDS(mateqtlouts)
      mateqtlouts<-preprocess_mateqtlouts(mateqtlouts)
      # filter meqtl matrices based on overlap with GWAS
      # message("Saving pre-processed list under 'mateqtlouts.rds' in current dir. You can re-start the analysis from here by specifying preprocess=FALSE.")

    if(preprocess_GWAS==TRUE){
      message(paste0("preprocess_GWAS=TRUE. Reading in unprocessed GWAS. Selecting regions for colocalisation based on supplied parameters: \n GWAS window=",GWAS_window,
      "\n GWAS signif=",GWASsignif))
      GWAS<-read.table(GWAS,head=T)
      GWAS<-select_regions(GWAS,window=GWAS_window,pval=GWASsignif)
      saveRDS(GWAS,"GWAS_processed_hyprcoloc.rds")
      message("Processed GWAS saved under 'GWAS_processed_hyprcoloc.rds'.")


    }else{
      message("Reading in pre-processed GWAS...")
      GWAS<-readRDS(processed_GWAS_file)
    }

    mateqtlouts<-lapply(mateqtlouts,function(x){x<-x[which(x$SNP %in% GWAS$rsid),]
      return(x)})

    # second step of intersection
    commonrows<-Reduce(intersect,lapply(mateqtlouts,row.names))
    mateqtlouts<-lapply(mateqtlouts,function(x){
      x<-x[commonrows,]
      return(x)
    })

    saveRDS(mateqtlouts,"GWAS_intersected_mateqtlouts.rds")

  } else {
    message("preprocess=FALSE. Reading in preprocessed GWAS and mateqtlouts.")
    mateqtlouts<-readRDS(intersected_mateqtlouts_file)
    GWAS<-readRDS(processed_GWAS_file)
  }




  #make sure GWAS is a data frame to check chr format
  if(class(GWAS)=="list"){
    GWAS<-do.call(rbind,GWAS)
  }

  if(length(grep("chr",GWAS$chr))==0){
    GWAS$chr<-paste0("chr",GWAS$chr)
  }
  processed_gwas<-split(GWAS,GWAS$signif_snp_region)

  nregions<-length(processed_gwas)
  message(paste0(nregions," regions were selected."))


  # -============================================================================================================
  # -=-=-=-=-=-=-=- STEP 2. Generating the gene-specific SE and Beta matrices -=-=-=-=-=-=-===
  # -============================================================================================================

  gene_locs<-read.table(gene_loc_file)
  date=Sys.Date()
  message("Prepping coloc data..")
  coloc_data<-prep_coloc_data(mateqtlouts,processed_gwas,gene_locs)
  coloc_data<-add_directionality(coloc_data)

  if(length(coloc_data$betas)!=length(processed_gwas)){
      regions<-names(coloc_data$betas)
      processed_gwas<-processed_gwas[regions]
      regions_lost=nregions-length(processed_gwas)
      message(paste0(regions_lost," regions were lost because no genes from gene_loc file were found within."))
  }

  saveRDS(coloc_data,"coloc_inputs.rds")


  if(hyprcoloc==TRUE){
    # message("Running hyprcoloc..")
    # hyprcoloc_df<-run_hyprcoloc(coloc_data,processed_gwas)



    # hyprcoloc_df<-hyprcoloc_df[grep("GWAS",hyprcoloc_df$traits),]
    # if(is.data.frame(hyprcoloc_df)==TRUE){
    #   hyprcoloc_df$index_snp<-do.call(rbind,strsplit(rownames(hyprcoloc_df),"\\."))[,1]
    #   hyprcoloc_df<-add_pvalues_hyprcoloc(hyprcoloc_df,coloc_data)
    #   write.table(hyprcoloc_df,paste0(date,"_",run_name,"_hyprcoloc_results.txt"))
    # }
    message("Running hyprcoloc by gene..")
    hyprcoloc_df_bygene=run_hyprcoloc_bygene(coloc_data,processed_gwas)

    message("Running hyprcoloc by celltype..")
    hyprcoloc_df_bycelltype=run_hyprcoloc_bycelltype(coloc_data,processed_gwas)

    if(is.data.frame(hyprcoloc_df_bygene)==TRUE){
      write.table(hyprcoloc_df_bygene,paste0(date,"_",run_name,"_hyprcoloc_results_bygene.txt"))
    }

    if(is.data.frame(hyprcoloc_df_bycelltype)==TRUE){
      write.table(hyprcoloc_df_bycelltype,paste0(date,"_",run_name,"_hyprcoloc_results_bycelltype.txt"))
    }
  }

  if(coloc==TRUE){

    message("Running COLOC..")
    indiv_numbers<-read.table(indiv_numbers_file)
    snplocs<-as.data.frame(data.table::fread(snp_loc_file))
    maf_df<-as.data.frame(data.table::fread(maf_file))
    coloc_df<-coloc_wrap(coloc_data,
      snplocs=snplocs,
      maf_df=maf_df,
      processed_gwas=processed_gwas,
      GWAS_type=GWAS_type,
      indiv_numbers=indiv_numbers,
      coloc.p12_prior=coloc.p12_prior)
    write.table(coloc_df,paste0(date,"_",run_name,"_COLOC_results.txt"))
  }

  if(MR==TRUE){
    message("Prepping MR..")
    library(MendelianRandomization)

    #Step 1: Determine regions of interest
    coloc_df<-coloc_df[coloc_df$PP.H4.abf>coloc_cutoff,]
    if(nrow(coloc_df)==0){
      message("No colocs pass threshold.")
    } else {

      regions<-coloc_df$region
      traits<-rownames(coloc_df)

      message(paste0(length(regions)," regions were selected for MR based on COLOC PP.H4 threshold."))
      #step 2: Harmonize alleles
      coloc_data$betas<-coloc_data$betas[regions]
      coloc_data$std.errors<-coloc_data$std.errors[regions]
      coloc_data$directionality<-coloc_data$directionality[regions]
      processed_gwas<-processed_gwas[regions]

      message("Harmonizing alleles..")
      allele_df<-read.table(allele_df_file)
      coloc_data<-harmonize_alleles_wrap(coloc_data,processed_gwas,allele_df)


      #Step 3: Prune snps
      message("Pruning eQTLs based on LD..")
      MR_regions<-vector("list",length(regions))
        names(MR_regions)<-regions
        message(paste0("Filtering out eQTLs with FDR below ",eqtl_FDR_cutoff))
      for(i in 1:length(regions)){

        # MR_regions[[i]]<-filter_eqtls_MR(pruned_region=MR_regions[[i]],
        #                                   trait=traits[i],
        #                                   eqtl_pval_cutoff=eqtl_pval_cutoff)

        MR_regions[[i]]<-prune_snps(coloc_data=coloc_data,
          region=regions[i],
          r_cutoff=prune_snps_r_cutoff,
        eqtl_FDR_cutoff=eqtl_FDR_cutoff,trait=traits[i])

      }
     traits<-traits[as.numeric(which(unlist(lapply(MR_regions,is.null))==FALSE))]
      MR_regions<-Filter(Negate(is.null), MR_regions)


      if(length(MR_regions)==0){
        message("No regions remain after eQTL filtering.")
      } else{
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
        saveRDS(df,paste0(date,"_",run_name,"_MR_res.rds"))

      }
    }
  }
  end_time<-Sys.time()
  runtime=end_time-start_time
  message(paste0("Analysis complete. Runtime: ",runtime))
}
