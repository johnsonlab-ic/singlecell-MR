

#load libraries
library(doParallel)
library(lme4)
library(Matrix)
library(lmerTest)
library(afex)
library(qvalue)


#define inputs

celltype="Excitatory_Neurons"
geno_mat="genotype_012mat.csv"



#these are the lognormalized counts (pre-correction for covariates)
exp_mat=paste0(celltype,"_pseudobulk.csv")
exp_mat=read.table(exp_mat)



# covariates used in matrixEQTL. These include;
# - expression PCs calculated on residuals (after correction for clinical covariates)
# - genotype PCs calculated on genotype data 


covs_used=paste0(celltype,"covs_used.rds")
covs_used=readRDS(covs_used)


##read in genotype data
geno_mat=as.data.frame(data.table::fread(args$geno_mat))
rownames(geno_mat)=geno_mat$snp
geno_mat$snp=NULL
colnames(geno_mat)<-gsub("/",".",colnames(geno_mat))
colnames(geno_mat)<-gsub("-",".",colnames(geno_mat))

##read in MatrixEQTL outputs, only testing significant associations 
meqtl_outs=readRDS("mateqtlouts_0.2FDR.rds")
