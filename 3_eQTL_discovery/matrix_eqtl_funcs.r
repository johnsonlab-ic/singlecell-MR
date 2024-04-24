

###main eqtl function 


calculate_ciseqtl=function(exp_mat,
  geno_mat,
  pcs,
  geno_pcs,
  exp_loc,
  geno_loc,
  name,
  cisDist=1e6,
  covadj=FALSE,
  covadj_pc=FALSE,
  geno_pcs_adj=TRUE,
  standardize=FALSE,
  cov_file,
  specify_genes=FALSE,
  specify_gene_list,
  covs_to_include,
  pvOutputThreshold=2e-5,
  filter_trans_FDR=FALSE,
  pvOutputThreshold_cis=5e-2,
  save_results=TRUE){

  exp_locs<-exp_loc
  geno_loc<-geno_loc

  common_names<-intersect(colnames(exp_mat),colnames(geno_mat))
  geno_mat<-geno_mat[common_names]
  exp_mat<-exp_mat[common_names]

  # this is in case the number of genes differ
  commongenes<-intersect(exp_locs$geneid,row.names(exp_mat))
  exp_mat<-exp_mat[commongenes,]

  # from MatrixEQTL manual: "The order of genes in gene locations
  # does not have to match the order of gene expression"


  cvrt=MatrixEQTL::SlicedData$new();
  if(covadj==TRUE){
    if(cov_file==""){
      covs<-as.data.frame(matrix(ncol=ncol(exp_mat),nrow=0))
      colnames(covs)<-exp_mat
    }else{
      covs<-read.table(cov_file)
      covs<-covs[covs_to_include,]
      covs<-covs[colnames(exp_mat)]
    }

    if(covadj_pc==TRUE){
      covs<-rbind(covs,pcs)
    }

    if(geno_pcs_adj==TRUE){
      covs=rbind(covs,geno_pcs)
    }


    message("Covs used: ",paste0(rownames(covs),sep=", "))
    cvrt$CreateFromMatrix(as.matrix(covs))
    saveRDS(covs,paste0(name,"covs_used.rds"))
  }

  snps=MatrixEQTL::SlicedData$new();
  snps$CreateFromMatrix(as.matrix(geno_mat))
  
  if(specify_genes==TRUE){
      exp_mat<-exp_mat[rownames(exp_mat) %in% specify_gene_list,]
  }

  if(standardize==TRUE){
    scaled<-scale(t(exp_mat),scale=T,center=F)
    exp_mat<-as.data.frame(t(scaled))
  }
  gene=MatrixEQTL::SlicedData$new();
  gene$CreateFromMatrix(as.matrix(exp_mat))
  n_indivs<-ncol(exp_mat)

  me<-suppressMessages(MatrixEQTL::Matrix_eQTL_main(snps,
                    gene,
                    cvrt = cvrt,
                    pvOutputThreshold = pvOutputThreshold,
                    useModel = MatrixEQTL::modelLINEAR,
                    errorCovariance = numeric(),
                    verbose = FALSE,
                    output_file_name=NULL,
                    output_file_name.cis =NULL,
                    pvOutputThreshold.cis =pvOutputThreshold_cis,
                    snpspos = geno_loc,
                    genepos = exp_locs,
                    cisDist = cisDist,
                    pvalue.hist = FALSE,
                    min.pv.by.genesnp = FALSE,
                    noFDRsaveMemory = FALSE))

  if(save_results==TRUE){


    if(pvOutputThreshold_cis>0){
      tmp<-me$cis$eqtls
      names(tmp)[which(names(tmp)=="statistic")]<-"t.stat"
      names(tmp)[which(names(tmp)=="pvalue")]<-"p.value"
      names(tmp)[which(names(tmp)=="snps")]<-"SNP"
      saveRDS(tmp,paste0(name,"_cis_MatrixEQTLout.rds"))
    }


    if(pvOutputThreshold>0){
      if(pvOutputThreshold_cis>0){
        tmp<-me$trans$eqtls
        names(tmp)[which(names(tmp)=="statistic")]<-"t.stat"
        names(tmp)[which(names(tmp)=="pvalue")]<-"p.value"
        names(tmp)[which(names(tmp)=="snps")]<-"SNP"
      if(filter_trans_FDR==TRUE){
          tmp<-tmp[tmp$FDR<0.2,]
        }
      saveRDS(tmp,paste0(name,"_trans_MatrixEQTLout_0.2FDR.rds"))
        }else{

          tmp<-me$all$eqtls
          if(filter_trans_FDR==TRUE){
              tmp<-tmp[tmp$FDR<0.2,]
            }

          saveRDS(tmp,paste0(name,"_trans_MatrixEQTLout_0.2FDR.rds"))
              }
      }
  }
    message(paste0("MatrixEQTL calculated for ",name, "."))
}



## PC optimisation eqtl
optimize_eqtl=function(exp_mat,
  geno_mat,
  pcs,
  geno_pcs,
  exp_loc,
  geno_loc,
  name,
  cisDist=1e6,
  filter_chr=TRUE,
  geno_pcs_adj=TRUE,
  chr="chr1",
  cov_file,
  covs_to_include,
  pvOutputThreshold=2e-5){

    
  exp_mat<-exp_mat
  geno_mat<-geno_mat
  exp_locs<-exp_loc
  geno_loc<-geno_loc

  common_names<-intersect(colnames(exp_mat),colnames(geno_mat))
  geno_mat<-geno_mat[common_names]
  exp_mat<-exp_mat[common_names]

  # this is in case the number of genes differ
  commongenes<-intersect(exp_locs$geneid,row.names(exp_mat))
  exp_mat<-exp_mat[commongenes,]


  ### get Max PCs - set to 100, less if exp_mat is small
  if(ncol(exp_mat)<100){
      
      
      down_signif <- function(x, digits = 0) {
          m <- 10^(ceiling(log(x, 10)) - digits)
          (x %/% m)*m
        }
      
      max_pcs=down_signif(ncol(exp_mat),1)/10
  }else{
      max_pcs=10
  }

  gene=MatrixEQTL::SlicedData$new();
  gene$CreateFromMatrix(as.matrix(exp_mat))
  n_indivs<-ncol(exp_mat)
  
  if(cov_file==""){
    covs<-as.data.frame(matrix(ncol=ncol(exp_mat),nrow=0))
    colnames(covs)<-exp_mat
  }else{
    covs<-read.table(cov_file)
    covs<-covs[covs_to_include,]
    covs<-covs[colnames(exp_mat)]
  }

  covs<-as.matrix(covs)
  ##add in pcs as provided by user
  covs<-rbind(covs,pcs)

  ##add Geno PCs
  if(geno_pcs_adj==TRUE){
      covs=rbind(covs,geno_pcs)
  }

  #filter to only contain chr1 to speed up optimisation
  if(filter_chr==TRUE){
      geno_loc=geno_loc[geno_loc$chrom %in% chr,]
      exp_loc=exp_locs[exp_locs$chr %in% chr,]
      geno_mat=geno_mat[rownames(geno_mat) %in% geno_loc$annot,]
    }
  
  snps=MatrixEQTL::SlicedData$new();
  snps$CreateFromMatrix(as.matrix(geno_mat))

  nsnps=nrow(geno_mat)
  ngenes=nrow(exp_mat)


  res_vector=vector()
  time_start=Sys.time()
  n_eqtls_results=vector()
  n_egenes_results=vector()
  pcs_vector=vector()
  
  
  for(i in 1:max_pcs){

        n_pcs=i*10
        pcs_vector=c(pcs_vector,n_pcs)

        covvec=c(covs_to_include,paste0(rep("PC."),1:(n_pcs)))
        tmp_covs=covs[rownames(covs) %in% covvec,]

        message(paste0(Sys.time()," Testing MatrixEQTL for ",n_pcs," PCs"))
        message(paste0(nsnps," SNPs and ",ngenes," genes used."))
        message("Covs used: ",paste0(rownames(tmp_covs),sep=", "))

        cvrt=MatrixEQTL::SlicedData$new();
        cvrt$CreateFromMatrix(as.matrix(tmp_covs))

  
        
        me<-suppressMessages(MatrixEQTL::Matrix_eQTL_main(snps,
                            gene,
                            cvrt=cvrt,
                            pvOutputThreshold = 0,
                            useModel = MatrixEQTL::modelLINEAR,
                            errorCovariance = numeric(),
                            verbose = FALSE,
                            output_file_name=NULL,
                            output_file_name.cis =NULL,
                            pvOutputThreshold.cis =5e-2,
                            snpspos = geno_loc,
                            genepos = exp_loc,
                            cisDist = cisDist,
                            pvalue.hist = FALSE,
                            min.pv.by.genesnp = FALSE,
                            noFDRsaveMemory = FALSE))

        eqtl_out=me$cis$eqtls
        n_eqtls=nrow(eqtl_out[eqtl_out$FDR<0.05,])
        n_egenes=length(unique(eqtl_out[eqtl_out$FDR<0.05,]$gene))
        message(paste0(n_egenes," eGenes discovered at < 5% FDR."))

        n_eqtls_results=c(n_eqtls_results,n_eqtls)
        n_egenes_results=c(n_egenes_results,n_egenes)
    }
  
  res_df=data.frame(n_eqtls_results,n_egenes_results,pcs=pcs_vector)
  res_df=res_df[order(res_df$n_egenes_results,decreasing=T),]
  res_df=res_df[1,]
  return(res_df)


}



## Function to filter the genes in the expression matrix. Takes either a minimum percentage of individuals or a minimum number of individuals
filter_pseudobulk = function(exp_mat, minimum_percentage=NULL, minimum_indivs=NULL){
  
  if(!is.null(minimum_percentage) && is.null(minimum_indivs)){

    message(paste0("filter_pseudobulk=TRUE. Keeping genes expressed in minimum ",minimum_percentage, " % individuals."))
    total_individuals = ncol(exp_mat)
    minimum_indivs = ceiling(total_individuals * (minimum_percentage / 100))
    exp_mat <- exp_mat[rowSums(exp_mat > 0) >= minimum_indivs, ]

  }else if(!is.null(minimum_indivs) && is.null(minimum_percentage)){
    exp_mat <- exp_mat[rowSums(exp_mat > 0) >= minimum_indivs, ]
     message(paste0("filter_pseudobulk=TRUE. Keeping genes expressed in minimum ",minimum_indivs, " number of individuals."))

  }else{
    stop("You must supply either a minimum_percentage or minimum_indivs value")
  }
  
  return(exp_mat)
}

## Filter genotype matrix. Retains snps with at least 2 individuals in each of the 2 genotypic categories.
filter_genotype_matrix=function(genotypemat,filter_type=c("2.2","2.3")){

  tmp<-genotypemat

  message(paste0("Filtering genotype matrix. Starting with ",nrow(genotypemat)," snps."))

  tmp$counter_0<-rowSums(tmp[1:ncol(tmp)]==0)
  tmp$counter_1<-rowSums(tmp[1:ncol(tmp)]==1)
  tmp$counter_2<-rowSums(tmp[1:ncol(tmp)]==2)

  if(filter_type=="2.2"){

    message("Retaining snps with at least 2 individuals in each of the 2 genotypic categories.")

    option1=tmp[tmp$counter_0>=2 & tmp$counter_1>=2,]
    option2=tmp[tmp$counter_0>=2 & tmp$counter_2>=2,]
    option3=tmp[tmp$counter_1>=2 & tmp$counter_2>=2,]

    intersected_rows<-unique(c(rownames(option1),rownames(option2),rownames(option3)))

  } else if (filter_type=="2.3"){

    message("Retaining snps with at least 2 individuals in each of the 3 genotypic categories.")
    tmp<-tmp[tmp$counter_0>=2 & tmp$counter_1>=2 & tmp$counter_2>=2,]
    intersected_rows<-row.names(tmp)
  }
  genotypemat<-genotypemat[intersected_rows,]
  message(paste0("Filtering complete. ",nrow(genotypemat), " snps retained."))
  return(genotypemat)
}


#function to get residuals (for now, specifically adds Random effects to Sample_Source and Diagnosis, and fixed effects to all other covariates)
get_residuals=function(exp_mat,covs_to_include,cov_file){

  covmat=read.table(cov_file)
  message(paste0("Correcting expression matrix for known covariates.\nUser-defined covariates: '",paste(covs_to_include,collapse=", "),"'"))  
  rownames(covmat)=covmat$Individual_ID
  ##
  if(length(grep(covs_to_include[1],rownames(covmat))>0)){
    covmat=as.data.frame(t(covmat))
  }

  ###relevel according to largest
  for (col_name in covs_to_include) {
    if (is.character(covmat[[col_name]])) {
      covmat[[col_name]] <- as.factor(covmat[[col_name]])
    }
    
    if (is.factor(covmat[[col_name]])) {
      # Determine the level with the highest frequency
      ref_level <- names(sort(table(covmat[[col_name]]), decreasing = TRUE))[1]
      
      # Relevel the factor with the most frequent level as the reference
      covmat[[col_name]] <- relevel(covmat[[col_name]], ref = ref_level)
    }
  }

  ### keep common names
  samples=colnames(exp_mat)
  common_samples=intersect(samples,rownames(covmat))

  covmat=covmat[common_samples,]
  exp_mat=exp_mat[,common_samples]



  lmmodel = paste0("gene ~ ", paste(covs_to_include, collapse = " + "))
  exp_mat=t(apply(exp_mat,1,function(x){
      
    lm_mat=data.frame(gene=x)
    lm_mat=cbind(lm_mat,covmat)
      # fit=lmerTest::lmer(lmmodel,lm_mat)
      fit=lm(lmmodel,lm_mat)
      resid(fit)
    
      
  }))
  return(exp_mat)

}


#function to get genotype pcs
get_geno_pcs=function(geno_mat){

  pcs=prcomp(geno_mat,scale=T)
  pcs=pcs$rotation
  pcs=pcs[,1:100]
  colnames(pcs)=paste0("geno.",colnames(pcs))
  return(t(pcs))



}
