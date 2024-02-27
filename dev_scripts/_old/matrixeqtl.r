
#' This function gathers the information from a MatrixEQTL data frame,
#' and outputs a table with some descriptive statistics (number of eGenes, eQTLs)
#' @export


get_info=function(x){
  tmp<-vector()


    n_snp_genepairs<-nrow(x)
    n_eGenes<-length(unique(x$gene))
    n_eSNPs<-length(unique(x$SNP))

    # same again for FDR
    fdr_0_2<-x[x$FDR<0.2,]

    # n_snp_genepairs_fdr_0_2<-nrow(fdr_0_2)
    n_eGenes_fdr_0_2<-length(unique(fdr_0_2$gene))
    n_eSNPS_fdr_0_2<-length(unique(fdr_0_2$SNP))

    fdr_0_1<-x[x$FDR<0.1,]

    # n_snp_genepairs_fdr_0_0_5<-nrow(fdr_0_0_5)
    n_eGenes_fdr_0_1<-length(unique(fdr_0_1$gene))
    n_eSNPS_fdr_0_1<-length(unique(fdr_0_1$SNP))

    fdr_0_0_5<-x[x$FDR<0.05,]

    # n_snp_genepairs_fdr_0_0_5<-nrow(fdr_0_0_5)
    n_eGenes_fdr_0_0_5<-length(unique(fdr_0_0_5$gene))
    n_eSNPS_fdr_0_0_5<-length(unique(fdr_0_0_5$SNP))

    tmp<-c(n_eSNPs,
      n_eGenes,
      n_eSNPS_fdr_0_2,
      n_eGenes_fdr_0_2,
      n_eSNPS_fdr_0_1,
      n_eGenes_fdr_0_1,
      n_eSNPS_fdr_0_0_5,
        n_eGenes_fdr_0_0_5)

        names(tmp)<-c("n_eSNPs",
        "n_eGenes",
        "n_eSNPs_fdr_0_2","n_eGenes_fdr_0_2",
        "n_eSNPS_fdr_0_1","n_eGenes_fdr_0_1",
        "n_eSNPs_fdr_0_0_5","n_eGenes_fdr_0_0_5")

        return(tmp)


}
#' This function gathers the information from a MatrixEQTL data frame,
#' and outputs a table with some descriptive statistics (number of eGenes, eQTLs)
#' @export
get_unique=function(x){

    #x is a mateqtlouts list

    genes<-lapply(x,function(x){
        return(unique(x$gene))
    })
    genes<-length(unique(unlist(genes)))
    snps<-lapply(x,function(x){
        return(unique(x$SNP))
    })
    snps<-length(unique(unlist(snps)))

    #20%fdr------

    x<-lapply(x,function(a){
        a<-a[a$FDR<0.2,]
    })
    genes_0_2_FDR<-lapply(x,function(x){
        return(unique(x$gene))
    })
    genes_0_2_FDR<-length(unique(unlist(genes_0_2_FDR)))
    snps_0_2_FDR<-lapply(x,function(x){
        return(unique(x$SNP))
    })
    snps_0_2_FDR<-length(unique(unlist(snps_0_2_FDR)))

    #10%fdr------


    x<-lapply(x,function(a){
        a<-a[a$FDR<0.1,]
    })
    genes_0_1_FDR<-lapply(x,function(x){
        return(unique(x$gene))
    })
    genes_0_1_FDR<-length(unique(unlist(genes_0_1_FDR)))
    snps_0_1_FDR<-lapply(x,function(x){
        return(unique(x$SNP))
    })
    snps_0_1_FDR<-length(unique(unlist(snps_0_1_FDR)))




    #5%fdr------

    x<-lapply(x,function(a){
        a<-a[a$FDR<0.05,]
    })
    genes_0_0_5_FDR<-lapply(x,function(x){
        return(unique(x$gene))
    })
    genes_0_0_5_FDR<-length(unique(unlist(genes_0_0_5_FDR)))
     snps_0_0_5_FDR<-lapply(x,function(x){
        return(unique(x$SNP))
    })
    snps_0_0_5_FDR<-length(unique(unlist(snps_0_0_5_FDR)))

    finalvec<-c(snps,genes,snps_0_2_FDR,genes_0_2_FDR,snps_0_1_FDR,genes_0_1_FDR,
      snps_0_0_5_FDR,genes_0_0_5_FDR)
    return(finalvec)


}

#' @export
get_unique_nonoverlapping=function(x){


  get_nonoverlapping=function(mateqtlouts,type=c("gene","snp")){
    
    if(type=="gene"){
        
        genes=lapply(mateqtlouts,function(x){
            return(unique(x$gene))
        })
        all_genes <- unlist(genes)

        # Find the unique genes for each element
        unique_genes <- lapply(gene_list, function(gene_vector) {
          unique_genes <- gene_vector[!gene_vector %in% all_genes[duplicated(all_genes)]]
          return(unique_genes)
        })
        
        return(sum(unlist(lapply(unique_genes,length))))

    } else if (type=="snp"){
        
        genes=lapply(mateqtlouts,function(x){
            return(unique(x$SNP))
        })
        all_genes <- unlist(genes)

        # Find the unique genes for each element
        unique_genes <- lapply(gene_list, function(gene_vector) {
          unique_genes <- gene_vector[!gene_vector %in% all_genes[duplicated(all_genes)]]
          return(unique_genes)
        })
        
        return(sum(unlist(lapply(unique_genes,length))))
        
    }
  }

    genes=get_nonoverlapping(mateqtlouts,type="gene")
    snps=get_nonoverlapping(mateqtlouts,type="snp")

    #20%fdr------
    mateqtlouts=lapply(mateqtlouts,function(x){x[x$FDR<0.2,]})
    genes_0_2_FDR=get_nonoverlapping(mateqtlouts,type="gene")
    snps_0_2_FDR=get_nonoverlapping(mateqtlouts,type="snp")
        
    #10%fdr------
    mateqtlouts=lapply(mateqtlouts,function(x){x[x$FDR<0.1,]})
    genes_0_1_FDR=get_nonoverlapping(mateqtlouts,type="gene")
    snps_0_1_FDR=get_nonoverlapping(mateqtlouts,type="snp")
    
    #5%fdr------
    mateqtlouts=lapply(mateqtlouts,function(x){x[x$FDR<0.05,]})
    genes_0_0_5_FDR=get_nonoverlapping(mateqtlouts,type="gene")
    snps_0_0_5_FDR=get_nonoverlapping(mateqtlouts,type="snp")

    #x is a mateqtlouts list)

    finalvec<-c(snps,genes,snps_0_2_FDR,genes_0_2_FDR,snps_0_1_FDR,genes_0_1_FDR,
      snps_0_0_5_FDR,genes_0_0_5_FDR)
    return(finalvec)


}

#' This function gathers the information from a MatrixEQTL data frame,
#' and outputs a table with some descriptive statistics (number of eGenes, eQTLs)
#' @export
summarize_outputs_cis=function(outdir){

  setwd(outdir)
  cellnames<-gsub("_cis_MatrixEQTLout.rds","",list.files(pattern="_cis_MatrixEQTLout"))
  message("Reading in MatrixEQTL outputs..")
  mateqtlouts<-lapply(list.files(pattern="_cis_MatrixEQTLout"),readRDS)
  names(mateqtlouts)<-cellnames

  message("Saving MatrixEQTL outputs as list object..")
  saveRDS(mateqtlouts,"mateqtlouts.rds")

  res<-lapply(mateqtlouts,get_info)
  res<-as.data.frame(do.call(cbind,res))
  colnames(res)<-cellnames

  unique<-get_unique(mateqtlouts)
  true_unique=get_unique_nonoverlapping(mateqtlouts)
  res$total<-rowSums(res)
  res$unique<-unique
  res$true_unique=true_unique

  mateqtlouts<-lapply(mateqtlouts,function(x){
    x<-x[x$FDR<0.2,]
    return(x)
  })
  saveRDS(mateqtlouts,"mateqtlouts_0.2FDR.rds")

  write.table(res,"summary_cis_eQTL_numbers.txt")


}
#' This function gathers the information from a MatrixEQTL data frame,
#' and outputs a table with some descriptive statistics (number of eGenes, eQTLs)
#' @export
summarize_outputs_trans=function(outdir){

  setwd(outdir)
  cellnames<-gsub("_trans_MatrixEQTLout.rds","",list.files(pattern="_trans_MatrixEQTLout"))
  message("Reading in MatrixEQTL outputs..")
  mateqtlouts<-lapply(list.files(pattern="_trans_MatrixEQTLout"),readRDS)
  names(mateqtlouts)<-cellnames

  message("Saving MatrixEQTL outputs as list object..")
  saveRDS(mateqtlouts,"mateqtlouts.rds")

  res<-lapply(mateqtlouts,get_info)
  res<-as.data.frame(do.call(cbind,res))
  colnames(res)<-cellnames

  unique<-get_unique(mateqtlouts)
  res$total<-rowSums(res)
  res$unique<-unique

  mateqtlouts<-lapply(mateqtlouts,function(x){
    x<-x[x$FDR<0.2,]
    return(x)
  })
  saveRDS(mateqtlouts,"mateqtlouts_trans_0.2FDR.rds")

  write.table(res,"summary_trans_eQTL_numbers.txt")


}

#' Main eQTL calculation function. Calculates eQTLs using the MatrixEQTL package based on 4 main input files:
#' genotype matrix, snp locations, expression matrix, gene locations. Outputs a cis and
#' if the option is selected, a trans-eQTL file with snp, gene, p-value, t-stat, and FDR.
#'
#' @param exp_mat The expression matrix in dataframe format. Rownames are genes, colnames are individuals. Output from agg_cpm_pseudobulk().
#' @param geno_mat The genotype matrix in dataframe format. Values are 0,1,2 of the ref allele (2 is homozygous ref). Output of get_genotype_matrix().
#' @param exp_loc A data frame of gene locations. Output from get_gene_locations().
#' @param name The prefix for the output. In the cellQTL main, this will be the cell type.
#' @param cisDist The distance from snp to gene for a "cis" association. This distance is set to either side of the gene as per MatrixEQTL description. Default is 1e6.
#' @param covadj If TRUE, include a covariate matrix in the MatrixEQTL calculation.
#' @param covadj_pc if TRUE, perform PC analysis on the expression matrix and include it to covariate matrix.
#' @param covadj_pcs_n if covadj_pc is TRUE, the number of PCs to adjust expression matrix for.
#' @param cov_file if covadj is TRUE, path to cov_file to read in. Rownames are covariates, colnames are sample IDs. If empty, the function will generate an matrix of 0 rows and N-individual columns.
#' @param covs_to_include is a character vector of covs to include.
#' @importClassesFrom MatrixEQTL SlicedData
#'
#' @export
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
  covadj_pcs_n=5,
  geno_pcs_adj=TRUE,
  standardize=FALSE,
  cov_file,
  specify_genes=FALSE,
  specify_gene_list,
  covs_to_include,
  pvOutputThreshold=2e-5,
  filter_trans_FDR=FALSE,
  pvOutputThreshold_cis=5e-2,
  save_results=TRUE,
  treeqtl_threshold=0.1,
  treeqtl_only=FALSE){


  exp_mat<-exp_mat
  geno_mat<-geno_mat
  exp_locs<-exp_loc
  geno_loc<-geno_loc

  common_names<-intersect(colnames(exp_mat),colnames(geno_mat))
  geno_mat<-geno_mat[common_names]
  exp_mat<-exp_mat[common_names]


  # colnames(geno_mat)<-colnames(exp_mat)


  # this is in case the number of genes differ
  commongenes<-intersect(exp_locs$geneid,row.names(exp_mat))
  exp_mat<-exp_mat[commongenes,]

  # from MatrixEQTL manual: "The order of genes in gene locations
  # does not have to match the order of gene expression"

  # library(MatrixEQTL)

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

      # pcs<-prcomp(exp_mat,scale=T,center=T)
      # pcs<-pcs$rotation
      # pcs<-pcs[,1:covadj_pcs_n]
      # pcs<-t(pcs)
      # pcs<-pcs[,colnames(exp_mat)]
      # rownames(pcs)<-paste0(rep("PC."),1:covadj_pcs_n)
      # covs<-as.matrix(covs)
      covs<-rbind(covs,pcs)


    }

    if(geno_pcs_adj==TRUE){

      # pcs=get_geno_pcs(geno_mat,pcs_n=geno_pcs_n)
      covs=rbind(covs,geno_pcs)
      

    }


    message("Covs used: ",paste0(rownames(covs),sep=", "))
    cvrt$CreateFromMatrix(as.matrix(covs))
    saveRDS(covs,paste0(name,"covs_used.rds"))
  }

  snps=MatrixEQTL::SlicedData$new();
  snps$CreateFromMatrix(as.matrix(geno_mat))
  
  #remove geno_mat from mem to clear space
  # rm(geno_mat)
  # rm(geno_mat_tmp)

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
  

  # mem=pryr::mem_used()
  # mem=mem/1e9
  # message(paste0(mem," Gb memory used."))
  
  # if(treeqtl_only==FALSE){
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
    # mem=pryr::mem_used()
    # mem=mem/1e9
    # message(paste0(mem," Gb memory used."))

    message(paste0("MatrixEQTL calculated for ",name, "."))
}
  


#' @export
get_residuals=function(exp_mat,covs_to_include,cov_file,ref_levels=NULL){

  covmat=read.table(cov_file)
  message(paste0("Correcting expression matrix for known covariates.\nUser-defined covariates: '",paste(covs_to_include,collapse=", "),"'"))  

  ##
  if(length(grep(covs_to_include[1],rownames(covmat))>0)){
    covmat=as.data.frame(t(covmat))
  }

  if(is.null(ref_levels)==FALSE){

    message(paste0("Ref levels set. Re-ordering data frame to set the reference factor levels as:\n"))
    print(data.frame(covs_to_include,ref_levels))

    for(i in 1:length(covs_to_include)){

      cov=covs_to_include[i]
      ref=ref_levels[i]


      cov1=covmat[covmat[,cov] %in% ref,]
      cov2=covmat[!covmat[,cov] %in% ref,]
      covmat=rbind(cov1,cov2)
      covmat[,cov]=factor(covmat[,cov],levels=unique(covmat[,cov]))
      covmat[,cov]=relevel(covmat[,cov],ref)

    }

    
  }

  samples=colnames(exp_mat)
  common_samples=intersect(samples,rownames(covmat))

  covmat=covmat[common_samples,]
  exp_mat=exp_mat[,common_samples]

  
  lmmodel=paste0("gene~",paste(covs_to_include,collapse="+"),collapse="")
  message(paste0("Obtaining residuals using lm(): "),lmmodel)



  exp_mat=t(apply(exp_mat,1,function(x){
      
    lm_mat=data.frame(gene=x)
    lm_mat=cbind(lm_mat,covmat)
      lmmodel=paste0("gene~",paste(covs_to_include,collapse="+"),collapse="")
      fit=lm(lmmodel,lm_mat)
      resid(fit)
    
      
  }))
  return(exp_mat)



}



#' @export
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
  geno_pcs_n=3,
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


  # colnames(geno_mat)<-colnames(exp_mat)


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


  ### add PCs
  # pcs<-prcomp(exp_mat,scale=T,center=T)
  # pcs<-pcs$rotation
  # pcs<-pcs[,1:(max_pcs*10)]
  # pcs<-t(pcs)
  # pcs<-pcs[,colnames(exp_mat)]
  # rownames(pcs)<-paste0(rep("PC."),1:(max_pcs*10))
  covs<-as.matrix(covs)
  covs<-rbind(covs,pcs)


  ##add Geno PCs
  if(geno_pcs_adj==TRUE){

      # pcs=get_geno_pcs(geno_mat,pcs_n=geno_pcs_n)
      # covs=rbind(covs,pcs)
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

         #stop looping once n_eGenes has peaked
        # if(i>(max_pcs/2)){
        #   if(n_egenes<n_egenes_results[i-1] && n_egenes<n_egenes_results[i-2]){
        #     message(paste0("Optimization peaked at ",(i-2)*10," PCs.")
        #     )
        #     break
        #   }
        # }
    }
  
  res_df=data.frame(n_eqtls_results,n_egenes_results,pcs=pcs_vector)
  res_df=res_df[order(res_df$n_egenes_results,decreasing=T),]
  res_df=res_df[1,]
  return(res_df)


}


#' @export
prep_covariates=function(covs,covadj=TRUE,
  cov_file="",
  geno_PC=TRUE,
  geno_pcs_n=3,
  covadj_pc=TRUE,
  covadj_pcs_n=10){


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

      pcs<-prcomp(exp_mat,scale=T,center=T)
      pcs<-pcs$rotation
      pcs<-pcs[,1:covadj_pcs_n]
      pcs<-t(pcs)
      pcs<-pcs[,colnames(exp_mat)]
      rownames(pcs)<-paste0(rep("PC."),1:covadj_pcs_n)
      covs<-as.matrix(covs)
      covs<-rbind(covs,pcs)


    }

    return(covs)
    
  }

}

#' @export
get_geno_pcs=function(geno_mat){

  pcs=prcomp(geno_mat,scale=T)
  pcs=pcs$rotation
  pcs=pcs[,1:100]
  colnames(pcs)=paste0("geno.",colnames(pcs))
  return(t(pcs))



}

