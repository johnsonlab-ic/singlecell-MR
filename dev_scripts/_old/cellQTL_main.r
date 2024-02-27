#' Main function to generate eQTLs from single-cell and genotype data, end to end. Includes pseudobulking, genotype matrix generation, filtering
#' and results summarization.
#'
#' @param step1 default is TRUE. This can be switched off after the first run. Takes a seurat object, splits it by cell type and pseudobulks. Make sure your seurat object meta data has "CellType" indicating the cell type classification for each cell you want to run the analysis on.
#' @param min.cells Minimum number of cells per individual in a cell type for a pseudobulk value to be included. Default is 100. You are encouraged to explore the stability of your pseudobulk values depending on this threshold.
#' @param summarize_outputs default is TRUE. Summarizes the eQTLs generated into a summary table.
#' @param input_seurat Only needed if Step1=TRUE. Make sure the cells are classified by cell type and that in the metadata (seuratobj[[]]) the column name for that classificaiton is "CellType"
#' @param input_seurat_list if TRUE, takes a list of seurat input files and generates pseudobulk items for each object before binding them together. Useful for large datasets.
#' @param celltype_colname column name for the celltype classification in the Seurat object Metadata. Default is "CellType"
#' @param run_id Name or ID for the eQTL run. Useful if you are generating multiple runs and testing out different parameters.
#' @param step2 default is TRUE. This takes as input a gds file or vcf to generate genotype matrix, snp locations and MAF files.
#' @param preprocess_vcf Default is FALSE. If you are providing vcfs, set to TRUE.
#' @param vcfs VCF files if you need to convert to gds.
#' @param gds_inputfile A seqArray gds inputfile. Supply this if preprocess_vcf=FALSE and you already have a gds object.
#' @param autosomalonly Default is TRUE. Only considers chromosomes 1:22.
#' @param maf Minimum minor allele frequency (MAF) for a SNP to be included in the genotype matrix.
#' @param filter_genotype_matrix Filter snps by genotype dosage as defined in the geno_filter_type.
#' @param geno_filter_type A string dictating how the snps should be filtered by genotype. "2.2" means that for a snp to be included, there needs to be at least 2 snps in 2 out of the 3 genotypic categories (0,1,2). "2.3" means at least 2 snps in all 3 categories.
#' @param geno_loc_file Default is "snp_chromlocations.csv" as generated in step 1. You can provide your own snp locations (as long as it has the same format) if you are looking at specific snps or regions.
#' @param filter_pseudobulk if TRUE, filters the expression matrices based on minimum individuals expressing the gene.
#' @param filter_pseudobulk_thresh minimum individuals expressing the gene for a gene to be expressed. For example, if you set '3', a gene will be included as long as 3 out of all individuals have a value other than 0.
#' @param specify_genes default is FALSE. If TRUE, supply a character vector of genes you are interested in.
#' @param specify_genes_list Character vector of genes to include in analysis.
#' @param specify_celltypes default is FALSE. If TRUE, runs the analysis on certain cell types.
#' @param specify_celltypes_names Character vector of cell types to run the analysis on. Make sure they match the prefix of the pseudobulk files.
#' @param pval_thresh_cis maximum p-value for a cis-eQTL association to be included in output. Set to 1 if you want all possible associations, and 0 if you don't want any cis-asociations.
#' @param pval_thresh_trans same as cis. Be careful about running both cis and trans in the same analysis - this will add a lot of time to the run.
#' @param covadj_matrixeqtl if TRUE, adjust for covariates in the MatrixEQTL calculation.
#' @param covs Character vector of covariates to include.
#' @param cov_file Path to covariate matrix to read in. Rows are covs, column names are individuals.
#' @param covadj_pc_matrixeqtl default is TRUE. Adjusts expression matrices for principal components.
#' @param covadj_pcs_n_matrixeqtl Number of PCs to include in correction if above is TRUE. It is recommended you perform multiple runs to optimize PC usage.
#' @export

### test comment for github

cellQTL=function(
  step1=TRUE,
  summarize_outputs=TRUE,
  input_seurat,
  input_seurat_list=FALSE,
  input_seurat_list_files,
  celltype_colname="CellType",
  indiv_colname="Sample_ID",
  outputdir,
  run_id=NULL,
  method=c("Merged","Indiv"),
  pseudobulk_method=c("agg_cpm"),
  standardize=TRUE,
  step2=TRUE,
  gds_inputfile,
  covs=c(""),
  autosomalonly=TRUE,
  maf=0.05,
  filter_genotype_matrix=FALSE,
  geno_filter_type=c("2.2","2.3"),
  geno_loc_file="snp_chromlocations.csv",
  geno_mat_file="genotype_012mat.csv",
  preprocess_vcf=FALSE,
  filter_pseudobulk=FALSE,
  filter_pseudobulk_thresh=3,
  min.cells=50,
  specify_genes=FALSE,
  specify_genes_list=NULL,
  specify_samples=FALSE,
  specify_samples_vector=NULL,
  get_exp_residuals=TRUE,
  ref_levels=NULL,
  pval_thresh_cis=1,
  pval_thresh_trans=0,
  filter_trans_FDR=TRUE,
  covadj_matrixeqtl=FALSE,
  covadj_pc_matrixeqtl=FALSE,
  covadj_pcs_n_matrixeqtl=5,
  geno_pc_adj=TRUE,
  geno_pcs_n=3,
  cov_file="",
  cisDist=1e6,
  optimize_pcs=FALSE,
  filter_chr=FALSE,
  use_optimization_res=TRUE,
  specify_celltypes=FALSE,
  specify_celltypes_names=NULL){


  pseudobulk_method<-match.arg(pseudobulk_method)
  geno_filter_type<-match.arg(geno_filter_type)  


  message(paste0(Sys.time(),": Starting cellQTL run. Parameters used:\n covadj_matrixeqtl=",covadj_matrixeqtl,
  "\n covs= ",paste(covs,collapse=" "),
  "\n filter genotype matrix=",filter_genotype_matrix,
  "\n filter genotype matrix threshold=",geno_filter_type,
  "\n filter_pseudobulk=",filter_pseudobulk,
  "\n filter_pseudobulk_thresh=",filter_pseudobulk_thresh,
  "\n specify_genes=",specify_genes,
  "\n specify_celltypes=",specify_celltypes))

  message(paste0("Additional parameters include;",
  "\noptimize_pcs=",optimize_pcs,
  "\nget_exp_residuals=",get_exp_residuals,
  "\ngeno_pc_adj=",geno_pc_adj,
  "\ngeno_pcs_n=",geno_pcs_n,
  "\nstandardize=",standardize))

  start_time<-Sys.time()
  #---------------------------------------------------
  #--------------------------------------------------
  # STEP 1:Pseudobulking
  #--------------------------------------------------
  #---------------------------------------------------
  if(step1==TRUE){
    message(paste0(Sys.time(),": Step 1: Processing of expression matrices.."))
    dir.create(paste0(outputdir,"/MatrixEQTL_IO"),showWarnings=FALSE)
      matrixeqtldir<-paste0(outputdir,"/MatrixEQTL_IO/")
      setwd(matrixeqtldir)
    if(input_seurat_list==FALSE){
    
      message("Reading in seurat input file..")
      seuratobj<-readRDS(input_seurat)
      celltypelist<-Seurat::SplitObject(seuratobj,split.by=celltype_colname)

      message("Pseudobulking by cell type.. ")


      if(pseudobulk_method=="SCTransform"){
        message("Using SCTransform normalised data and calculating mean..")
        pseudobulk_sct_cell(celltypelist)


      }else if(pseudobulk_method=="agg_cpm"){
        message("Summing counts and normalising by cpm..")
        agg_count_list<-pseudobulk_counts(celltypelist,min.cells=min.cells,indiv_col=indiv_colname)
        agg_count_list=normalize_pseudobulk(agg_count_list)
        agg_count_list=Filter(Negate(is.null),agg_count_list)
        for(i in 1:length(agg_count_list)){
          write.table(agg_count_list[[i]],paste0(names(agg_count_list[i]),"_pseudobulk.csv"))
        }

      }


    } else {
      
      message("First, making sure that all Seurat objects have the same genes..")
      so=readRDS(input_seurat_list_files[[1]])
      Seurat::DefaultAssay(so)="RNA"
      commongenes=rownames(so)
      for(i in 2:length(input_seurat_list_files)){
        so=readRDS(input_seurat_list_files[[i]])
        Seurat::DefaultAssay(so)="RNA"
        genes=rownames(so)
        commongenes=intersect(genes,commongenes)
      }

      
      n_seurat_objs=length(input_seurat_list_files)
      message(paste0(n_seurat_objs," seurat objects were provided. Reading in Seurat obj 1 .."))
      seuratobj<-readRDS(input_seurat_list_files[[1]])
      seuratobj=seuratobj[commongenes,]
      
      celltypelist<-Seurat::SplitObject(seuratobj,split.by=celltype_colname)
      agg_count_list_full<-pseudobulk_counts(celltypelist,min.cells=min.cells,indiv_col=indiv_colname)

      for(i in 2:length(input_seurat_list_files)){
        message(paste0(" Reading in Seurat obj ",i,".."))
        seuratobj<-readRDS(input_seurat_list_files[[i]])
        seuratobj=seuratobj[commongenes,]
        celltypelist<-Seurat::SplitObject(seuratobj,split.by=celltype_colname)
        agg_count_list<-pseudobulk_counts(celltypelist,min.cells=min.cells,indiv_col=indiv_colname)
        agg_count_list<-agg_count_list[names(agg_count_list) %in% names(agg_count_list_full)]
        agg_count_list<-agg_count_list[names(agg_count_list)]
        agg_count_list_full<-Map(cbind,agg_count_list_full,agg_count_list)
      }
      nullvec=unlist(lapply(agg_count_list_full,is.null))
      dropped_celltypes=names(nullvec[which(nullvec==TRUE)])
      agg_count_list_full=Filter(Negate(is.null),agg_count_list_full)


      agg_count_list_full=normalize_pseudobulk(agg_count_list_full)


      if(length(dropped_celltypes)>0){
        message(paste0(dropped_celltypes," celltype was dropped due to no individuals passing min.cells criteria."))
      }


      for(i in 1:length(agg_count_list_full)){
          write.table(agg_count_list_full[[i]],paste0(names(agg_count_list[i]),"_pseudobulk.csv"))
        }
    }

   

    matrixeqtldir<-paste0(outputdir,"/MatrixEQTL_IO/")
    setwd(matrixeqtldir)

    cellnames<-gsub("_pseudobulk.csv","",list.files(pattern="_pseudobulk.csv"))
    tmp<-lapply(list.files(pattern="_pseudobulk.csv"),read.table,check.names=FALSE)
    names(tmp)<-cellnames

    message("Getting gene locations for each pseudobulk matrix..")
    message(paste0(Sys.time(),": Step 1 completed. Pseudobulk matrices written in MatrixEQTL_IO"))


    for(i in 1:length(tmp)){
      genelocs<-get_gene_locations(tmp[[i]])
      write.table(genelocs,paste0(names(tmp[i]),"_gene_locations.csv"))
    }

    
  } else {
    message("Step 1 skipped. Using pre-computed pseudobulk matrices")
  }


  #---------------------------------------------------
  #--------------------------------------------------
  # STEP 2: Matrix of n x p (n=sample size, p=SNPs)
  #--------------------------------------------------
  #---------------------------------------------------


  if(step2==TRUE){
    matrixeqtldir<-paste0(outputdir,"/MatrixEQTL_IO/")
    setwd(matrixeqtldir)

    message(paste0(Sys.time(),": Step 2: Creating genotype matrix .."))
    get_genotype_matrix(vcfs=vcfs,
      inputfile=gds_inputfile,
      preprocess=preprocess_vcf,
      outdir=outputdir,
      autosomalonly=autosomalonly,
      minmaf=maf)

      
    
  } else {
    message("Step 2 skipped. Using pre-computed genotype matrix")
  }


  #---------------------------------------------------
  #---------------------------------------------------
  # STEP 3: Calculate MatrixEQTL
  #---------------------------------------------------
  #---------------------------------------------------

  message(paste0(Sys.time(),": Step 3: Calculating eQTLs.."))
  matrixeqtldir<-paste0(outputdir,"/MatrixEQTL_IO/")
  setwd(matrixeqtldir)

  exp_list<-lapply(list.files(pattern="_pseudobulk.csv"),read.table,check.names=FALSE)
  namelist<-gsub("_pseudobulk.csv","",list.files(pattern="_pseudobulk.csv"))
  names(exp_list)<-namelist
  exp_list<-lapply(exp_list,function(x){
      colnames(x)<-gsub("/",".",colnames(x))
      colnames(x)<-gsub("-",".",colnames(x))
      return(x)

  })

  n_indivs_pseudobulk<-unlist(lapply(exp_list,ncol))
  n_indivs_pseudobulk<-as.data.frame(n_indivs_pseudobulk)
  write.table(n_indivs_pseudobulk,"n_indivs_pseudobulk.txt")

  # message("Matching genotype and expression matrix column names..")

  geno_mat<-as.data.frame(data.table::fread(geno_mat_file))
  rownames(geno_mat)<-geno_mat$snp
  geno_mat$snp<-NULL
  colnames(geno_mat)<-gsub("/",".",colnames(geno_mat))
  colnames(geno_mat)<-gsub("-",".",colnames(geno_mat))

  #remove NAs
  geno_loc<-as.data.frame(data.table::fread(geno_loc_file))
  geno_loc<-geno_loc[,c("annot","chrom","position")]
  #fix geno_loc based on removed NAs
  row.names(geno_loc)<-geno_loc$annot
  geno_mat<-geno_mat[rownames(geno_loc),]
  geno_mat<-geno_mat[complete.cases(geno_mat),]
  geno_loc<-geno_loc[rownames(geno_mat),]



  row.names(geno_loc)<-rep(1:nrow(geno_loc))
  message("Making sure geno_loc and geno_mat match..")
  geno_mat<-geno_mat[rownames(geno_mat) %in% geno_loc$annot,]


  #project specific - rename sample
  if(any(colnames(geno_mat)=="S03.009.B")==TRUE){
    names(geno_mat)[names(geno_mat)=="S03.009.B"]<-"S03.009"
  }



  if(filter_genotype_matrix==TRUE){
      geno_mat<-filter_genotype_matrix(geno_mat,filter_type=geno_filter_type)
      row.names(geno_loc)<-geno_loc$annot
      geno_loc<-geno_loc[rownames(geno_mat),]
      row.names(geno_loc)<-rep(1:nrow(geno_loc))
  }


  if(filter_pseudobulk==TRUE){
    message(paste0("filter_pseudobulk=TRUE. Keeping genes with minimum ",filter_pseudobulk_thresh, " individuals."))
    exp_list<-lapply(exp_list,filter_pseudobulk,minimum_indivs=filter_pseudobulk_thresh)
  }

  ####!! 19th December 2022, the below Standardized function moved into MatrixEQTL wrapper

  # if(standardize==TRUE){
  #   message("Standardize=TRUE. Scaling and centering matrices.")
  #   exp_list<-lapply(exp_list,function(x){
  #     x<-t(x)
  #     x<-scale(x,center=FALSE,scale=T)
  #     return(as.data.frame(t(x)))
  #   })
  # }


  # if(specify_genes==TRUE){
  #   exp_list<-lapply(exp_list,function(x){
  #     return(x[rownames(x) %in% specify_genes_list,])
  #   })
  # }

  # each pseudobulk matrix has its own gene locations file
  exp_loc<-lapply(list.files(pattern="gene_locations.csv"),read.table)
  names(exp_loc)<-namelist

if(specify_celltypes==TRUE){
    message("specify_celltypes=TRUE. Specific cell types selected.")
    message(paste0("Performing analysis for ",specify_celltypes_names))

    namelist<-specify_celltypes_names
    exp_list<-exp_list[namelist]
    exp_loc<-exp_loc[namelist]
  }

if(specify_samples==TRUE){

  message("Specify_samples=TRUE. Using a subset of the data.")
  specify_samples_vector=gsub("/",".",specify_samples_vector)
  specify_samples_vector=gsub("-",".",specify_samples_vector)

  exp_list=lapply(exp_list,function(x){
    names=intersect(colnames(x),specify_samples_vector)
    x=x[,names]
    return(x)
  })



}


  #create output dir
meqtl_outdir=paste0(matrixeqtldir,Sys.Date(),"_",run_id,"_MatrixEQTL_output")
dir.create(meqtl_outdir,showWarnings=FALSE)
setwd(meqtl_outdir)


### get PCs - this must be done before residuals
pcs_list=list()
for(i in 1:length(exp_list)){
  exp_mat=exp_list[[i]]
  pcs<-prcomp(exp_mat,scale=T,center=T)
  pcs<-pcs$rotation
  pcs<-pcs[,1:(max_pcs*10)]
  pcs<-t(pcs)
  pcs<-pcs[,colnames(exp_mat)]
  rownames(pcs)<-paste0(rep("PC."),1:(max_pcs*10))
  pcs_list[[i]]=pcs


}

#### get residuals

if(get_exp_residuals==TRUE){
 
 if(cov_file==""){
  message("You must supply a covariate file to obtain corrected expression residuals.")
 }else{
  message(paste0("Correcting expression matrix for known covariates.\nUser-defined covariates: '",paste(covs,collapse=", "),"'"))  

   if(is.null(ref_levels)==FALSE){
    message(paste0("Ref levels set. Re-ordering data frame to set the reference factor levels as:\n"))
    print(data.frame(covs,ref_levels))
   }

  lmmodel=paste0("gene~",paste(covs,collapse="+"),collapse="")
  message(paste0("Obtaining residuals using lm(): "),lmmodel)

  exp_list=suppressMessages(lapply(exp_list,get_residuals,
  covs_to_include=covs,
  cov_file=cov_file,ref_levels=ref_levels))

  for(i in 1:length(exp_list)){
    fn=paste0(names(exp_list[i]),"_residuals_pseudobulk.csv")
    
    write.table(exp_list[[i]],fn)
  }
  exp_list=lapply(exp_list,as.data.frame)

 }
 
  
}


###PC optimization 

if(optimize_pcs==TRUE){
  message("Optimize=TRUE. Testing optimal PCs for MatrixEQTL.")
    pcs_df=data.frame()


 for(i in 1:length(exp_list)){
      message(paste0("Optimizing PCs for ",namelist[i]))
      
    
      exp_mat<-exp_list[[i]]
      pcs=pcs_list[[i]]

    nsamples=length(intersect(colnames(exp_mat),colnames(geno_mat)))
    message(paste0(nsamples, " individuals will be used."))
      if(get_exp_residuals==TRUE){
        cov_file=""
      }
      best_pcs=optimize_eqtl(exp_mat=exp_mat,
        exp_loc=exp_loc[[i]],
        geno_mat=geno_mat,
        pcs=pcs,
        geno_loc=geno_loc,
        filter_chr=filter_chr,
        cisDist=cisDist,
        name=namelist[i],
        cov_file=cov_file,
        geno_pcs_adj=geno_pc_adj,
        geno_pcs_n=geno_pcs_n,
        covs_to_include=covs)

        pcs_df=rbind(pcs_df,best_pcs)
    }
    rownames(pcs_df)=namelist
    write.table(pcs_df,"PC_optimisation_results.txt")

}


  #check whether user has run the optimisation steps.
  if(file.exists("PC_optimisation_results.txt") && use_optimization_res==TRUE){
    message("Found PC optimisation file. Using results to choose number of PCs.")
    pcs_df=read.table("PC_optimisation_results.txt")
    
  }

  for(i in 1:length(exp_list)){
    message(paste0(Sys.time(),": Calculating eQTLs for ",namelist[i],".."))

    exp_mat<-exp_list[[i]]
    nsnps=nrow(geno_mat)
    ngenes=nrow(exp_mat)
    nsamples=length(intersect(colnames(exp_mat),colnames(geno_mat)))
    message(paste0(nsamples, " individuals will be used."))
    message(paste0(nsnps," SNPs and ",ngenes," genes will be included."))

    #check if PC optimisation was read in. two steps here to avoid reading in table every time
    if(exists("pcs_df")){
      covadj_pcs_n_matrixeqtl=pcs_df[namelist[i],]$pcs[1]
    }
 
    if(get_exp_residuals==TRUE){
      cov_file=""
    }
    calculate_ciseqtl(exp_mat=exp_mat,
    exp_loc=exp_loc[[i]],
    geno_mat=geno_mat,
    geno_loc=geno_loc,
    cisDist=cisDist,
    name=namelist[i],
    specify_genes=specify_genes,
    specify_gene_list=specify_gene_list,
    pvOutputThreshold=pval_thresh_trans,
    pvOutputThreshold_cis=pval_thresh_cis,
    covadj=covadj_matrixeqtl,
    covadj_pc=covadj_pc_matrixeqtl,
    geno_pcs_adj=geno_pc_adj,
    geno_pcs_n=geno_pcs_n,
    standardize=standardize,
    covadj_pcs_n=covadj_pcs_n_matrixeqtl,
    cov_file=cov_file,
    covs_to_include=covs)
  }

  ##clear memory
  rm(geno_mat)
  rm(exp_list)
  rm(geno_loc)
  gc()

  if(summarize_outputs==TRUE){
    message("summarize_outputs=TRUE. Summarizing cis-eQTL outputs..")
    summarize_outputs_cis(meqtl_outdir)
    if(pval_thresh_trans>0){
      summarize_outputs_trans(meqtl_outdir)
    }
  }

  end_time<-Sys.time()

  duration=end_time - start_time

  message("Analysis complete!")
  message(paste0("Runtime: ",duration))

}
