



suppressMessages(library(argparse))
library(SCGsuite)


###########################################################
################ INPUT PARAMETERS #########################
###########################################################
parser <- ArgumentParser()


##### input files
###########################################################


parser$add_argument("--exp_mat", help = "Expression matrix path")
parser$add_argument("--geno_mat", help = "Genotype dosage matrix (0,1,2) path")
parser$add_argument("--exp_loc", help = "Gene locations path")
parser$add_argument("--geno_loc", help = "SNP locations path")
parser$add_argument("--cov_file", help = "Covariate file")
parser$add_argument("--run_name", help = "Run name")
parser$add_argument("--covs_to_include",nargs="+",help = "")


##### processing params
###########################################################

parser$add_argument("--filter_genotype_matrix", 
help = "Filter genotype matrix",action="store_true")
parser$add_argument("--filter_pseudobulk_matrix", 
help = "Filter pseudobulk matrix",action="store_true")
parser$add_argument("--exp_mat_thresh_percent", 
help = "Filter pseudobulk matrix min indivs")
parser$add_argument("--get_residuals", 
help = "Gets residuals, correcting for covariates",action="store_true")
parser$add_argument("--specify_samples_file", 
help = "File with samples to include. Will be read with ReadLines().")
parser$add_argument("--trans_eqtls", 
help = "Whether or not to calculate trans-eQTLs 'true' or 'false' .")



args <- parser$parse_args(commandArgs(TRUE))



exp_mat=args$exp_mat
name=gsub(".*/([^/]+)_pseudobulk.csv", "\\1", exp_mat)
name=gsub("_pseudobulk.csv","",name)

geno_mat=args$geno_mat
exp_loc=args$exp_loc
geno_loc=args$geno_loc
run_id=args$run_name
cov_file=args$cov_file
if(length(cov_file)==0){
  cov_file=""
  cov_file_included=FALSE
}else{
  cov_file_included=TRUE
}


covs=args$covs_to_include
print(covs)


##### processing params
###########################################################

trans_eqtls=args$trans_eqtls
if(trans_eqtls=="true"){
  trans_eqtls=TRUE
  message("trans_eqtls set to 'true'. Will calculate trans-eQTls at nominal (p<0.05) significance.")
}else{
  trans_eqtls=FALSE
}

# filter_genotype_matrix=TRUE
geno_filter_type="2.2"
# filter_pseudobulk=TRUE
filter_pseudobulk_thresh=as.numeric(args$exp_mat_thresh_percent)
# specify_genes=FALSE
# specify_genes_list=NULL
# specify_samples=FALSE
if (!is.null(args$specify_samples_file)) {
  # Read the file and store the samples in specify_samples_vector
  specify_samples_vector <- readLines(args$specify_samples_file)
}
get_exp_residuals=TRUE
ref_levels=NULL


###########################################################
###### Helper funcs
###########################################################
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



  ###run the model, extract residuals
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

##### matrixEQTL params
###########################################################


optimize_pcs=TRUE
filter_chr=FALSE
use_optimization_res=TRUE
cisDist=1e6

standardize=TRUE
pval_thresh_cis=1
if(trans_eqtls==TRUE){
  pval_thresh_trans=5e-7
  pval_thresh_cis=0
}else{
  pval_thresh_trans=0
  pval_thresh_cis=1
}

covadj_matrixeqtl=TRUE
covadj_pc_matrixeqtl=TRUE
covadj_pcs_n_matrixeqtl=5
geno_pc_adj=TRUE
geno_pcs_n=3



###########################################################
################ Read in files + reformat #################
###########################################################



exp_mat=read.table(exp_mat)
exp_loc=read.table(exp_loc)
geno_loc=as.data.frame(data.table::fread(geno_loc))
geno_mat=as.data.frame(data.table::fread(geno_mat))


rownames(geno_mat)<-geno_mat$snp
geno_mat$snp<-NULL
colnames(geno_mat)<-gsub("/",".",colnames(geno_mat))
colnames(geno_mat)<-gsub("-",".",colnames(geno_mat))


geno_loc<-geno_loc[,c("annot","chrom","position")]
row.names(geno_loc)<-geno_loc$annot
geno_mat<-geno_mat[rownames(geno_loc),]
geno_mat<-geno_mat[complete.cases(geno_mat),]
geno_loc<-geno_loc[rownames(geno_mat),]



row.names(geno_loc)<-rep(1:nrow(geno_loc))
message("Making sure geno_loc and geno_mat match..")
geno_mat<-geno_mat[rownames(geno_mat) %in% geno_loc$annot,]


message("Keeping common samples in all matrices..")
common_names=intersect(colnames(exp_mat),colnames(geno_mat))

if(cov_file_included==TRUE){
  covmat=read.table(cov_file)
  covmat=covmat[complete.cases(covmat),]
  # if(length(grep(covs[1],rownames(covmat))>0)){
  #   covmat=as.data.frame(t(covmat))
  # }
  common_names=intersect(common_names,covmat$Individual_ID)
}

exp_mat=exp_mat[,common_names]
geno_mat=geno_mat[,common_names]



###########################################################
###########################################################
#########        post-process / QC        #################
###########################################################
###########################################################



####################################
### specify samples if needed
####################################

if(!is.null(args$specify_samples_file)){

  message("Specify_samples=TRUE. Using a subset of the data.")
  specify_samples_vector=gsub("/",".",specify_samples_vector)
  specify_samples_vector=gsub("-",".",specify_samples_vector)
  
  common_names=intersect(colnames(exp_mat),specify_samples_vector)
  exp_mat=exp_mat[,common_names]
  geno_mat=geno_mat[,common_names]

}

message(paste0(length(common_names)," individuals retained. Double-check that this number is accurate."))

####################################
### filter expression matrix
####################################

if(args$filter_pseudobulk_matrix==TRUE){
    exp_mat=filter_pseudobulk(exp_mat,minimum_percentage=filter_pseudobulk_thresh)
}
####################################
### filter genotype matrix
####################################

if(args$filter_genotype_matrix==TRUE){
    geno_mat<-filter_genotype_matrix(geno_mat,filter_type=geno_filter_type)
    row.names(geno_loc)<-geno_loc$annot
    geno_loc<-geno_loc[rownames(geno_mat),]
    row.names(geno_loc)<-rep(1:nrow(geno_loc))
}

####################################
### STANDARDIZE EXPRESSION MATRIX
####################################
scaled<-scale(t(exp_mat),scale=T,center=F)
exp_mat<-as.data.frame(t(scaled))


####################################
### GET PCS ###
####################################

if(ncol(exp_mat)<100){
    down_signif <- function(x, digits = 0) {
        m <- 10^(ceiling(log(x, 10)) - digits)
        (x %/% m)*m
      }
    
    max_pcs=down_signif(ncol(exp_mat),1)/10
}else{
    max_pcs=10
}
message("Getting PCs...")
pcs<-prcomp(exp_mat,scale=F,center=F)
pcs<-pcs$rotation
pcs<-pcs[,1:(max_pcs*10)]
pcs<-t(pcs)
pcs<-pcs[,colnames(exp_mat)]
rownames(pcs)<-paste0(rep("PC."),1:(max_pcs*10))
write.table(pcs,paste0(name,"_full_expression_PCs.txt"))

####################################
### GET GENO PCS ###
####################################


geno_pcs=get_geno_pcs(geno_mat)
saveRDS(geno_pcs,paste0(name,"_geno_pcs.rds"))
geno_pcs=geno_pcs[1:geno_pcs_n,]



### GET RESIDUALS (2-STEP approach, if TRUE) ###
####################################

if(get_exp_residuals==TRUE){
 
 if(cov_file==""){
  message("You must supply a covariate file to obtain corrected expression residuals.")
 }else{
  message(paste0("Correcting expression matrix for known covariates.\nUser-defined covariates: '",paste(covs,collapse=", "),"'"))  

  #  if(is.null(ref_levels)==FALSE){
  #   message(paste0("Ref levels set. Re-ordering data frame to set the reference factor levels as:\n"))
  #   print(data.frame(covs,ref_levels))
  #  }

  lmmodel=paste0("gene~",paste(covs,collapse="+"),collapse="")
  message(paste0("Obtaining residuals using lm(): "),lmmodel)

  exp_mat=suppressMessages(get_residuals(exp_mat,
  covs_to_include=covs,
  cov_file=cov_file))
  exp_mat=as.data.frame(exp_mat)
  write.table(exp_mat,paste0(name,"_residuals_pseudobulk.csv"))
  }
}

####################################
### OPTIMIZE N_PCS #################
####################################

if(optimize_pcs==TRUE && covadj_pc_matrixeqtl==TRUE){

    message("Optimize PCs = TRUE. Running Optimisation..")

    if(get_exp_residuals==TRUE){
        message("Expression matrix already corrected for covariates. Setting cov_file as null.")
        cov_file=""
        # geno_pc_adj=FALSE
      }

    best_pcs=optimize_eqtl(exp_mat=exp_mat,
        exp_loc=exp_loc,
        geno_mat=geno_mat,
        pcs=pcs,
        geno_pcs=geno_pcs,
        geno_loc=geno_loc,
        filter_chr=filter_chr,
        cisDist=cisDist,
        name=name,
        cov_file=cov_file,
        geno_pcs_adj=geno_pc_adj,
        geno_pcs_n=geno_pcs_n,
        covs_to_include=covs)

    best_pcs=best_pcs$pcs[1]


}

message(paste0("Optimization complete. ",best_pcs," PCs will be used."))

###########################################################
#########        MatrixEQTL                #################
###########################################################

nsnps=nrow(geno_mat)
ngenes=nrow(exp_mat)
nsamples=length(intersect(colnames(exp_mat),colnames(geno_mat)))
message(paste0(nsamples, " individuals will be used."))
message(paste0(nsnps," SNPs and ",ngenes," genes will be included."))

if(optimize_pcs==TRUE){
  pcs=pcs[1:best_pcs,]
}


calculate_ciseqtl(exp_mat=exp_mat,
exp_loc=exp_loc,
geno_mat=geno_mat,
geno_loc=geno_loc,
cisDist=1e6,
name=name,
specify_genes=FALSE,
specify_gene_list=specify_gene_list,
pvOutputThreshold=pval_thresh_trans,
pvOutputThreshold_cis=1,
covadj=covadj_matrixeqtl,
covadj_pc=covadj_pc_matrixeqtl,
pcs=pcs,
geno_pcs=geno_pcs,
geno_pcs_adj=geno_pc_adj,
standardize=TRUE,
cov_file=cov_file)

