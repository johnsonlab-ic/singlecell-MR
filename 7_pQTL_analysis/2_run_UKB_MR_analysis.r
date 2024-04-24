

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
#### 2. test UKB-PPP for  MR ###
#######################################


MR_list=list()

gene_locs=read.table(paste0(meqtldir,"/Astrocytes_gene_locations.csv"))
pqtl_dirs=list.files(pqtl_analysis_dir)

for(i in 1:nrow(mrdf)){

    gene = paste0(mrdf$gene[i], "_")
    chrom = paste0("chr",mrdf$chrom[i],"_")
    trait = gsub(" ","_",mrdf$GWAS[i])
    celltype = mrdf$celltype[i]
    pqtl_dir=NULL
    pqtl_file=NULL

    message(paste("Processing gene ",gene,"and trait ",trait))
    pqtl_dir=grep(gene,pqtl_dirs)
    pqtl_dir=pqtl_dirs[pqtl_dir]

  if(length(pqtl_dir)==0){
    MR_list[[i]]="No pQTL file"
    next
  }else{
    pqtl_files=list.files(paste0(pqtl_analysis_dir,"/",pqtl_dir),full.names = TRUE)
    pqtl_file=grep(chrom,pqtl_files)
    pqtl_file=pqtl_files[pqtl_file]
    if(gene=="LMOD1_"){
        pqtl_file=pqtl_file[grep("Neurology",pqtl_file)]
    }
    pqtl_summary=data.table::fread(pqtl_file)
    pqtl_summary=pqtl_summary %>% 
      arrange(desc(LOG10P)) %>% 
      mutate(SNP=paste0(CHROM,":",GENPOS)) %>% 
      mutate(pvalue = 10^(-LOG10P)) %>% 
      mutate(FDR=p.adjust(pvalue,method="fdr")) %>% filter(A1FREQ>0.05)

      if(nrow(pqtl_summary[pqtl_summary$FDR<0.05])==0){
        MR_list[[i]]="No significant pQTLs"
        next
    }else{

    gene=gsub("_","",gene)
    #we obtain the COLOC objects to center our analysis around the GWAS index SNP, and not the gene. 
    #This ensures the analysis is directly comparable to our own MR results.
    coloc_objs_dir=coloc_dir
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

 region_snps=object@gwas[[regionName]]@beta$rsid
 snp_locs_tmp=object@snp_locs %>% filter(annot %in% region_snps)
region_snps=paste0(gsub("chr","",snp_locs_tmp$chrom),":",snp_locs_tmp$position)

pqtl_summary_cis=pqtl_summary %>% filter(SNP %in% region_snps)

  # get GWAS
  GWAS_DIR="/rds/general/user/ah3918/projects/roche/live/ALEX/ANALYSIS/GWAS_STUDIES/COLOC_INPUTS/"
  gwaslist=list.files(full.names = T,pattern = "coloc_ready",paste0(GWAS_DIR,"/CC/"))
  gwaslist=c(gwaslist,list.files(full.names = T,pattern = "coloc_ready",paste0(GWAS_DIR,"/QUANT/")))


  gwas_file=gwaslist[grep(trait,gwaslist)]
  gwas=data.table::fread(gwas_file)
  gwas=gwas %>% mutate(SNP=paste0(gsub("chr","",chr),":",pos))


  merged_df <- gwas %>%
    dplyr::rename(gwas_REF = A1, gwas_EFFECT = A2, gwas_beta = b, gwas_se = se, gwas_pval = pval) %>%
    inner_join(pqtl_summary_cis %>%
                dplyr::rename(pqtl_REF = ALLELE0, pqtl_EFFECT = ALLELE1, pqtl_beta = BETA, pqtl_se = SE, pqtl_pval = pvalue,pqtl_log10p=LOG10P,pqtl_FDR=FDR),
              by = "SNP") %>% dplyr::select(rsid,SNP, pqtl_REF, pqtl_EFFECT, pqtl_beta, 
              pqtl_se, pqtl_pval,pqtl_log10p,pqtl_FDR, 
              gwas_REF, gwas_EFFECT, gwas_beta, gwas_se, gwas_pval)


  merged_df<- merged_df %>%
    mutate(
      flip = ifelse(pqtl_EFFECT != gwas_EFFECT & pqtl_EFFECT == gwas_REF, TRUE, FALSE),
      pqtl_EFFECT = ifelse(flip, pqtl_REF, pqtl_EFFECT),
      pqtl_REF = ifelse(flip, pqtl_EFFECT, pqtl_REF),
      pqtl_beta = ifelse(flip, -pqtl_beta, pqtl_beta)
    )

  ## now clump 
  message("Clumping.. ")

  clump_df=data.frame(rsid=merged_df$rsid,pval=merged_df$pqtl_FDR)

  plink_bin=genetics.binaRies::get_plink_binary()
  path_to_binaries="/rds/general/user/ah3918/projects/single_cell_eqtl/live/MRC/REFERENCE_DATASETS/clumping_ref/EUR/EUR"
  r2_cutoff=0.01
  pvalue_cutoff=0.05

  #check if 

  suppressMessages(capture.output(snps_tokeep<-ieugwasr::ld_clump_local(clump_df,clump_r2=r2_cutoff,clump_p=pvalue_cutoff,plink_bin=plink_bin,
  bfile=path_to_binaries,clump_kb=1e6),file=nullfile()))

  mr_input=merged_df %>% filter(rsid %in% snps_tokeep$rsid)

  mr_input=MendelianRandomization::mr_input(snps=mr_input$rsid,bx=mr_input$pqtl_beta,
  bxse=mr_input$pqtl_se,
  by=mr_input$gwas_beta,
  byse=mr_input$gwas_se)

  MR_list[[i]]=mr_input 
  
  
  }
  }
}


names(MR_list)=mrdf$gene_trait

MR_results=lapply(MR_list,function(x){
    if(class(x)!="MRInput"){
        return(x)
    }else{
        mr_res=MendelianRandomization::mr_allmethods(x,method="ivw")
        pvalues=mr_res@Values$`P-value`
        names(pvalues)=mr_res@Values$Method
        ivw_beta <- mr_res@Values %>%
        filter(Method == "IVW") %>%
        pull(`Estimate`)
        n_ivs=length(x@snps)
        return(c(pvalues,n_ivs=n_ivs,ivw_beta=ivw_beta))
    }
})

final_df=as.data.frame(do.call(rbind,MR_results))
final_df=final_df %>% select(IVW,ivw_beta,n_ivs)
final_df$gene_trait<- names(MR_list)
final_df$analysis=final_df$IVW 
final_df$analysis=sapply(final_df$analysis,function(x){
    if(x=="No pQTL file" | x=="No significant pQTLs"){
        return(x)}
        else{
            return("MR")

        }
        
})


final_df=final_df %>% mutate(IVW=as.numeric(IVW),ivw_beta=as.numeric(ivw_beta),n_ivs=as.numeric(n_ivs))
mrdf=mrdf %>% mutate(gene_trait=paste0(gene,"_",GWAS)) %>% left_join(final_df,by="gene_trait")

#count number of proteins assessed 
tmpdf=mrdf %>% filter(analysis!="No pQTL file")
print("Filtering for proteins assessed..")
print(paste0("Number of genes assessed: ",length(unique(tmpdf$gene)),", representing ",length(unique(tmpdf$gene_trait))," unique gene-trait combinations"))

tmpdf=mrdf %>% filter(analysis!="No significant pQTLs" & analysis!="No pQTL file")

print("Filtering for proteins assessed and FDR < 0.05..")
print(paste0("Number of genes assessed: ",length(unique(tmpdf$gene)),", representing ",length(unique(tmpdf$gene_trait))," unique gene-trait combinations"))

print("Filtering for IVW <0.05..")
tmpdf=mrdf %>% filter(analysis!="No significant pQTLs" & analysis!="No pQTL file") %>% filter(IVW<0.05)
print(paste0("Number of genes assessed: ",length(unique(tmpdf$gene)),", representing ",length(unique(tmpdf$gene_trait))," unique gene-trait combinations")) 

tmpdf=tmpdf %>% select(celltype,gene,GWAS,IVW_pval,IVW_beta,IVW,ivw_beta,n_ivs) %>%
rename(`IVW p-value eQTL`=IVW_pval,`IVW beta eQTL`=IVW_beta,`IVW p-value pQTL`=IVW,`IVW beta pQTL`=ivw_beta,`Number of IVs pQTL`=n_ivs)
saveRDS(tmpdf,paste0(pqtl_analysis_dir,"/MR_results_UKB_PPP.rds"))

