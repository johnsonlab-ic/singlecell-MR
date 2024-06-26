

#load libraries
library(argparse)
library(doParallel)
library(lme4)
library(Matrix)
library(lmerTest)
library(afex)
library(qvalue)
library(dplyr)

###########################################################
################ INPUT PARAMETERS #########################
###########################################################
parser <- ArgumentParser()


##### input files
###########################################################

parser$add_argument("--celltype", help = "celltype to test")
parser$add_argument("--geno_mat", help = "genotype dosage matrix path")
parser$add_argument("--exp_mat", help = "path to bulk expression matrix")
parser$add_argument("--cov_file", help = "Covariate file")
parser$add_argument("--covs_used", help = "Covs used file in matrixeqtl")
parser$add_argument("--meqtl_outs",help = "meqtl 20% FDR outs")
parser$add_argument("--num_cores", help = "number of cores for the parallelisation")


args <- parser$parse_args(commandArgs(TRUE))


#### read in data
###########################################################
celltype=args$celltype
num_cores=args$num_cores

##geno mat
geno_mat=as.data.frame(data.table::fread(args$geno_mat))
rownames(geno_mat)=geno_mat$snp
geno_mat$snp=NULL
colnames(geno_mat)<-gsub("/",".",colnames(geno_mat))
colnames(geno_mat)<-gsub("-",".",colnames(geno_mat))

exp_mat=read.table(args$exp_mat)

## cov_file
covs=read.table(args$cov_file)
covs=covs[complete.cases(covs),]

#these are the covariates used in matrixEQTL. contains all expression PCs + genotype PCs.
covs_used=readRDS(args$covs_used)

#meqtl results 20% FDR
meqtl_outs=readRDS(args$meqtl_outs)


#### filter/reformat
###########################################################

common_names=intersect(colnames(exp_mat),covs$Individual_ID)
common_names=intersect(common_names,colnames(geno_mat))
common_names=intersect(common_names,colnames(covs_used))

covs=covs %>% filter(Individual_ID %in% common_names)
exp_mat=exp_mat[,common_names]
geno_mat=geno_mat[,common_names]
covs=covs[match(common_names,covs$Individual_ID),]
covs_used=covs_used[,common_names]

meqtl_outs=meqtl_outs[[celltype]]
meqtl_outs=filter(meqtl_outs,FDR<0.05)
meqtl_outs=meqtl_outs[!duplicated(meqtl_outs$gene),]


input_list <- split(meqtl_outs, seq(nrow(meqtl_outs)))
input_list <- lapply(input_list, function(x) c(x$SNP, x$gene))
start_time=Sys.time()


#### Run models
###########################################################


res=mclapply(input_list,function(x){
tryCatch({
    gene=x[2]
    snp=x[1]

    message(snp)
    message(gene)

    E=as.numeric(exp_mat[gene,])
    G=as.numeric(geno_mat[snp,])
    SAMPLE_SOURCE=covs$Sample_Source
    AGE=covs$Age
    SEX=covs$Sex
    PMI=covs$PMI
    D=covs$Diagnosis

    data=data.frame(E,G,SAMPLE_SOURCE,AGE,SEX,PMI,D)
    #adding in expression PCs (calculated on RESIDUALS prior to MatrixEQTL) + genotype PCs
    data=cbind(data,t(covs_used))


    ### relevel covs to largest group
    Sample_Source_idx <- which(unique(data$SAMPLE_SOURCE) == "Imperial")
    Diagnosis_idx <- which(unique(data$D) == "Control")
    data$SEX<- factor(data$SEX, level = c(unique(data$SEX)[1], unique(data$SEX)[-1]))
    data$SAMPLE_SOURCE <- factor(data$SAMPLE_SOURCE, level = c(unique(data$SAMPLE_SOURCE)[Sample_Source_idx], unique(data$SAMPLE_SOURCE)[-Sample_Source_idx]))
    data$D <- factor(data$D, level = c(unique(data$D)[Diagnosis_idx], unique(data$D)[-Diagnosis_idx]))

    ########################
    ### M0, no genotype   ##
    ########################

    m0=paste("E~1+(1+D|SAMPLE_SOURCE)+D+AGE+SEX+PMI+",paste0(colnames(data)[8:ncol(data)],collapse="+"),sep="")
    m0=suppressMessages(suppressWarnings(lmer(m0,data=data)))

    ########################
    ### M1 with genotype ##
    ########################


    m1=paste("E~1+G+(1+D|SAMPLE_SOURCE)+D+AGE+SEX+PMI+",paste0(colnames(data)[8:ncol(data)],collapse="+"),sep="")
    m1=suppressMessages(suppressWarnings(lmer(m1,data=data)))

    ### LRT test
    anova=suppressMessages(anova(m0,m1,test = "LRT"))
    anova_pval=anova["m1",8]
    if(anova_pval<0.05){
        m1vsm0="PASS"
        message(m1vsm0)
        m1pval=coef(summary(m1))["G",5]
    }else{
        m1vsm0="FAIL"
        message(m1vsm0)
        m1pval=NA
    }

    ########################
    ### M2 with interaction ##
    ########################
    m2=paste("E~1+(1+D|SAMPLE_SOURCE)+D*G+AGE+SEX+PMI+",paste0(colnames(data)[8:ncol(data)],collapse="+"),sep="")
    m2=suppressMessages(suppressWarnings(lmer(m2,data=data)))

     ### LRT test
    anova=suppressMessages(anova(m1,m2,test = "LRT"))
    anova_pval=anova["m2",8]

    if(anova_pval<0.05){
        m2vsm1="INTERACTION"
    }else{
        m2vsm1="NO_INTERACTION"
    }
        message(m2vsm1)
        m2pval=coef(summary(m2))["G",5]
        if ("DPD:G" %in% rownames(coef(summary(m2)))){
            m2pval_PD=coef(summary(m2))["DPD:G",5]
        }else{
            m2pval_PD=NA

        }

        if ("DAD:G" %in% rownames(coef(summary(m2)))){
            m2pval_AD=coef(summary(m2))["DAD:G",5]
        }else{
            m2pval_AD=NA

        }

         if ("DMS:G" %in% rownames(coef(summary(m2)))){
            m2pval_MS=coef(summary(m2))["DMS:G",5]
        }else{
            m2pval_MS=NA

        }
   

    resvec=c(gene=gene,snp=snp,m1vsm0=m1vsm0,
             m1pval=m1pval,
             m2vsm1=m2vsm1,
            m2pval_G=m2pval,
             anova_pval=anova_pval,
             m2pval_PD=m2pval_PD,
             m2pval_AD=m2pval_AD,
             m2pval_MS=m2pval_MS
            )
     return(resvec)

       }, error = function(e) {
    message(paste("Error with input:", gene," ",snp))
    message("Error message:", conditionMessage(e))
    c(gene,snp,rep(NA,8))
  })
},mc.cores=num_cores)

resdf=as.data.frame(do.call(rbind,res))
resdf$matrixeqtl_pval=meqtl_outs$p.value
resdf$matrixeqtl_FDR=meqtl_outs$FDR
resdf$m1pval=as.numeric(resdf$m1pval)
resdf$m2pval_G=as.numeric(resdf$m2pval_G)
resdf$m2pval_AD=as.numeric(resdf$m2pval_AD)
resdf$m2pval_PD=as.numeric(resdf$m2pval_PD)
resdf$m2pval_MS=as.numeric(resdf$m2pval_MS)
resdf$anova_pval=as.numeric(resdf$anova_pval)

#calculate q-values
resdf$anova_qval=qvalue(resdf$anova_pval)$qvalues
resdf=resdf %>% mutate(m2vsm1 = ifelse(anova_qval < 0.05, "INTERACTION", "NO_INTERACTION"))

numeric_cols <- sapply(resdf, is.numeric)
resdf[numeric_cols] <- lapply(resdf[numeric_cols], signif, 3) # 3 is the number of significant digits


saveRDS(resdf,paste0(celltype,"_lmer_results.rds"))