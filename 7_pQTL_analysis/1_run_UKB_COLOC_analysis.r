
#### the following script assumes you have downloaded all of the UKB-PPP data
# https://metabolomips.org/ukbbpgwas/
# https://www.synapse.org/#!Synapse:syn51364943/files/


library(dplyr)
library(tidyverse)
library(data.table)
library(openxlsx)
source("../4_colocalisation/cellCOLOC_source.r")
#######################################
#### 1. READ IN MR RESULTS / CLEAN ####
#######################################

#declare outpu/input directories
mr_controldir="/path_to_MR_results/"
coloc_dir="/path_to_coloc_objs/"
meqtldir="/path_to_process_expression_outputs/"
pqtl_analysis_dir="/path_to_pQTL_folders//"
lmedir="/path_to_LME_output_dir/"
plink_bin=genetics.binaRies::get_plink_binary()
path_to_binaries="/path_to_binaries/EUR/EUR"



##read in MR results
mrresults=list.files(mr_controldir,pattern="MR_results",full.names=T)
gwas_names=list.files(mr_controldir,pattern="MR_results")
gwas_names=gsub("_MR_results.txt","",gwas_names)
mrresults=lapply(mrresults,read.table)
names(mrresults)=gwas_names

mrresults <- lapply(names(mrresults), function(df_name) {
  df <- mrresults[[df_name]]
  df <- df[complete.cases(df),]
  if(nrow(df) < 1) {
    return(NULL)
  } else {
    df$GWAS <- df_name
    return(df)
  }
})

mrresults <- Filter(Negate(is.null), mrresults)
mrdf=as.data.frame(do.call(rbind,mrresults))


##The following section filters out certain genes; we decided to exclude MAPT + HLA genes, as well as genes overlapping with our interaction QTLs.
#genes identified to be age-related and overlapping disease-interacting eQTLs
disease_genes =readRDS(paste0(lmedir,"/disease_overlapping_age_genes.rds"))
mrdf=mrdf %>% filter(IVW_pval<0.05) %>% filter(!(celltype %in% names(disease_genes) & gene %in% unlist(disease_genes)))


gene_locs=read.table(paste0(meqtldir,"/Astrocytes_gene_locations.csv"))

#chr6:28,510,120-33,480,577 https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC
hla=gene_locs %>% filter(chr == "chr6", left > 28510120-5e5, right < 433480577+5e5)
#chr17:45,309,498-46,903,773 https://www.ncbi.nlm.nih.gov/grc/human/regions/MAPT 
mapt=gene_locs %>% filter(chr == "chr17", left > 45309498-5e5, right < 46903773+5e5)

mrdf$gene=sapply(mrdf$gene,function(x) ifelse(x %in% hla$geneid,"HLA-region",x))
mrdf$gene=sapply(mrdf$gene,function(x) ifelse(x %in% mapt$geneid,"MAPT-region",x))
mrdf <- mrdf %>%
  group_by(GWAS, celltype, gene) %>%
  slice_max(IVW_pval) %>%
  ungroup()

mrdf <- mrdf %>% filter(!(gene %in% c("HLA-region","MAPT-region")))


##add in SNP info

snplocs=data.table::fread(paste0(meqtldir,"/snp_chromlocations.csv"),data.table = F)
alleles=data.table::fread(paste0(meqtldir,"/MAF_mat.csv"),data.table = F)
snp_info = snplocs %>% left_join(alleles, by = c("annot" = "snp")) %>% 
mutate(chrom=gsub("chr","",chrom)) %>%
mutate(SNP=paste0(chrom,":",position))

### keep only first IV. This is only to grab the chromosome for pQTL MR.
mrdf <- mrdf %>%
  mutate(IVs = str_split(IVs, pattern = ",", simplify = TRUE)[, 1]) 

#add it to the MR data frame
mrdf=mrdf %>% left_join(snp_info %>% dplyr::select(SNP,ref,alt,annot,chrom),by=c("IVs"="annot"))

genes=unique(mrdf$gene)
writeLines(genes,paste0(pqtl_analysis_dir,"genes.txt"))

## now keep the top cell-type per gene (to avoid duplicating calculations)
mrdf=mrdf %>%
  arrange(GWAS, gene, desc(IVW_pval), desc(IVW_beta)) %>%
  group_by(GWAS, gene) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  mutate(gene_trait=paste0(gene,"_",GWAS))



#######################################
#### 2. test UKB-PPP for  COLOC ###
#######################################


analysis_dir=pqtl_analysis_dir

COLOC_list=list()
for(hit in 1:nrow(mrdf)){



gene=mrdf$gene[hit]
if(gene=="HLA-region" | gene=="MAPT-region"){
    COLOC_list[[hit]]=list(eqtl_pqtl=NULL,pqtl_gwas=NULL,eqtl_gwas=NULL)
    next
}
message(paste("Processing gene ",gene,"and trait ",trait))

coloc_objs_dir=paste0(pipelinedir,"/COLOCMR_OUTS/controls_only/")
cellCOLOC_objects=list.files(coloc_objs_dir,pattern="cellCOLOC",full.names = TRUE)
cellCOLOC_file=paste0(coloc_objs_dir,"/",trait,"_coloc_ready_cellCOLOC_obj.rds")
cellCOLOC_file=readRDS(cellCOLOC_file)

object=cellCOLOC_file
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


# Extract betas, standard errors, p-values, and SNP names
extract_data <- function(data, columns) {
  data[, columns] %>% pull(gene)
}

betas <- extract_data(regionDataObj@beta, c("SNP", gene))
std_errors <- extract_data(regionDataObj@std_error, c("SNP", gene))
eqtl_pvals <- extract_data(regionDataObj@pvalue, c("SNP", gene))
maf_input=regionDataObj@maf_info$maf
effect_allele=regionDataObj@maf_info$ref
non_effect_allele=regionDataObj@maf_info$alt
snp_names <- regionDataObj@beta[,"SNP"]

# Filter snp_locs to only include snps in snp_names
snp_locs <- snp_locs %>% 
  filter(annot %in% snp_names$SNP) %>% 
  .[match(snp_names$SNP, .$annot), ]

# Add in pos:name 
snp_locs <- snp_locs %>% 
  mutate(chrom = gsub("chr", "", chrom),
         SNP = paste0(chrom, ":", position))

# Create eQTL coloc input data frame
eqtl_coloc_input <- data.frame(
  SNP = snp_names$SNP,
  betas = betas,
  std_errors = std_errors,
  eqtl_pvals = eqtl_pvals,
  maf=maf_input,
  CHRPOS=snp_locs$SNP,
  pos=snp_locs$pos,
  rsid=snp_locs$annot,
  eqtl_effect_allele=effect_allele,
  eqtl_ref_allele=non_effect_allele
)

gene=paste0(gene,"_") 
pqtl_dir=grep(gene,pqtl_dirs)
pqtl_dir=pqtl_dirs[pqtl_dir]

# Check if pqtl_dir is not empty
if (length(pqtl_dir) > 0) {
  # Now grab the right chromosomal file
  pqtl_files=list.files(paste0(analysis_dir,pqtl_dir),full.names = TRUE)
  pqtl_file=grep(chrom,pqtl_files)
  pqtl_file=pqtl_files[pqtl_file]

  # Check if pqtl_file is not empty
  if (length(pqtl_file) > 0) {
    
    if(gene=="LMOD1_"){
        pqtl_file=pqtl_file[grep("Neurology",pqtl_file)]
    }
    # Read in the file
    pqtl_summary=data.table::fread(pqtl_file)

  } 
}else{
    message("Skipping..")
    COLOC_list[[hit]]=list(eqtl_pqtl=NULL,pqtl_gwas=NULL,eqtl_gwas=NULL)
    next
}

pqtl_summary=pqtl_summary %>% 
mutate(SNP=paste0(CHROM,":",GENPOS)) %>% 
mutate(pvalue = 10^(-LOG10P)) %>% 
mutate(FDR=p.adjust(pvalue,method="fdr"))
pqtl_summary <- dplyr::rename(pqtl_summary, CHRPOS = SNP,pqtl_effect_allele=ALLELE1,pqtl_ref_allele=ALLELE0)

combined_data <- inner_join(eqtl_coloc_input, pqtl_summary, by = "CHRPOS") %>%
  dplyr::select(CHRPOS, betas, std_errors, eqtl_pvals, BETA, SE, pvalue, maf,pos,rsid,
  eqtl_effect_allele,
  eqtl_ref_allele,
  pqtl_effect_allele,
  pqtl_ref_allele)



# Rename the columns
combined_data <- dplyr::rename(combined_data, 
                        SNP = CHRPOS,
                        pQTL_beta = BETA,
                        pQTL_se = SE,
                        pQTL_pval = pvalue)

combined_data=combined_data %>% filter(!duplicated(rsid))



## now GWAS input 
gwas_pvals=object@gwas[[regionName]]@pvalue
gwas_betas=object@gwas[[regionName]]@beta
gwas_se=object@gwas[[regionName]]@std_error
gwas_alleles=object@gwas[[regionName]]@alleles
gwas_input=gwas_pvals %>% 
left_join(gwas_betas,by="rsid") %>% 
left_join(gwas_se,by="rsid") %>% 
left_join(gwas_alleles,by="rsid") %>% dplyr::rename(GWAS_pval=pval,GWAS_beta=b,GWAS_se=se,GWAS_ref_allele=A1,GWAS_effect_allele=A2)

combined_data=combined_data %>% left_join(gwas_input,by="rsid")


library(coloc)

eqtl_coloc_input=list(beta=combined_data[,"betas"]*1,
varbeta=combined_data[,"std_errors"]^2,
snp=combined_data[,"rsid"],
type = "quant",
N = object@cellTypes[[celltype]]@n_indivs,
MAF = combined_data$maf,  # Adjusted MAF column name
sdY=1
)

pqtl_input=list(beta=combined_data[,"pQTL_beta"],
varbeta=combined_data[,"pQTL_se"]^2,
snp=combined_data[,"rsid"],
type = "quant",
N=45000,
MAF = combined_data$maf,  # Adjusted MAF column name
sdY=1
)

type=object@gwasInfo

if(type=="cc"){
  gwas_input = list(
    beta = combined_data %>% pull(GWAS_beta),  # Subset and convert to numeric vector
    varbeta = (combined_data %>% pull(GWAS_se))^2,  # Subset, convert and square
    snp = combined_data %>% pull(rsid),
    type = type, # Assuming the type is stored here
    N = object@gwas_N,  # GWAS data doesn't provide this, you might fill it in later
    MAF = combined_data %>% pull(maf), # Adjusted MAF column name
    s=object@gwas_S
)

}else if (type=="quant"){
  gwas_input = list(
      beta = combined_data %>% pull(GWAS_beta),  # Subset and convert to numeric vector
      varbeta = (combined_data %>% pull(GWAS_se))^2,  # Subset, convert and square
      snp = combined_data %>% pull(rsid),
      type = type, # Assuming the type is stored here
      N = object@gwas_N,  # GWAS data doesn't provide this, you might fill it in later
      MAF = combined_data %>% pull(maf), # Adjusted MAF column name
      sdY=1
  )

}   

res1=coloc.abf(eqtl_coloc_input, pqtl_input,p12=1e-2)
res2=coloc.abf(pqtl_input, gwas_input,p12=1e-2)
res3=coloc.abf(eqtl_coloc_input, gwas_input,p12=1e-2)


COLOC_list[[hit]]=list(eqtl_pqtl=res1,pqtl_gwas=res2,eqtl_gwas=res3,
input=combined_data)


}
names(COLOC_list)=mrdf$gene_trait
COLOC_list_res = lapply(names(COLOC_list), function(n) {
  x = COLOC_list[[n]]
  
  if(is.null(x$eqtl_pqtl)){
    return(NULL)
  } else {
    eqtl_pqtl = c(x$eqtl_pqtl$summary,n)
    pqtl_gwas = c(x$pqtl_gwas$summary,n)
    eqtl_gwas = c(x$eqtl_gwas$summary,n)
    df = rbind(eqtl_pqtl, pqtl_gwas, eqtl_gwas)
    
    
    return(df)
  }
})

COLOC_res = as.data.frame(do.call(rbind, COLOC_list_res)) %>%
  mutate_at(vars("PP.H0.abf":"PP.H4.abf"), as.numeric) %>% 
  mutate_at(vars("PP.H0.abf":"PP.H4.abf"), round, digits = 3) %>%
  dplyr::rename(gene_trait=V7) %>% 
  mutate(test = rownames(.)) %>%
 mutate(test = sub("\\..*", "", test))

saveRDS(COLOC_res,paste0(pqtl_analysis_dir,"/COLOC_res_1e-2.rds"))