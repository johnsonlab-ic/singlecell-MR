

# -============================================================================================================
# -=-=-=-=-=-=-=- STEP 1. Read inputs -=-=-=-=-=-=-===
# -============================================================================================================

###define inputs / output 


mateqtlouts="mateqtlouts.rds"
GWAS="AD_2022_COLOC_ready.txt"
maf_file="MAF_mat.csv"
gene_loc_file="Oligodendrocyte_gene_locations.csv"
indiv_numbers_file="n_indivs_pseudobulk.txt"

### define parameters
GWASsignif=5e-8
GWAS_window=1e6
GWAS_type="cc"
GWAS_name=gsub("_2022_COLOC_ready.txt","",GWAS)
#minimum number of SNPs for a region to be considered
min_snps=100


##load helper functions from cellCOLOC_source.r
source("cellCOLOC_source.r")


## read in inputs - GWAS first
GWAS<-as.data.frame(suppressWarnings(data.table::fread(GWAS)))
message(paste0(Sys.time(),": Checking column names.."))

#helper function - please check cellCOLOC_source.r
check_colnames(GWAS,c("rsid","pval","b","pos","chr","se","A2","A1","case.prop","N_total"))

if(any(GWAS$pval<GWASsignif)==FALSE){
    stop(paste0("No associations below specified p-value: ",GWASsignif,". Change 'GWASsignif' if you want to proceed."))
}

##Filter out zero values
GWAS=GWAS %>% dplyr::filter(!b %in% 0 & !se %in% 0 & !pval %in% 0)
GWAS=GWAS %>% dplyr::arrange(rsid, pval) %>% dplyr::distinct(rsid, .keep_all = TRUE)

#this is a helper function; see cellCOLOC_source.r
GWAS<-select_regions(GWAS,window=GWAS_window,pval=GWASsignif)

## read in MatrixEQTL outputs
mateqtlouts=readRDS(mateqtlouts)

#this is a helper function; see cellCOLOC_source.r
mateqtlouts=preprocess_mateqtlouts(mateqtlouts,GWAS)

##read in auxiliary files
gene_locs<-read.table(gene_loc_file)

#this is also a helper function to check which genes overlap with the GWAS regions
genestokeep=lapply(processed_gwas,get_genes_per_region,gene_locs=gene_locs)

indiv_numbers=read.table(indiv_numbers_file)
maf_df=as.data.frame(data.table::fread(maf_file))
snp_locs=as.data.frame(data.table::fread(snp_locs_file))

# -============================================================================================================
# -=-=-=-=-=-=-=- STEP 2. Generating cellCOLOC objects and run COLOC -=-=-=-=-=-=-===
# -============================================================================================================


# create cellCOLOC object as defined in cellCOLOC_source.r. 
# this object contains a list of each cell-type, each containing a list of regions. 
# each subsequent region is a list in itself, containing data frames of beta, se, pvalue, FDR, and MAF information for each SNP in the region;
# each row is a SNP, each column is a gene. This allows us to iterate over each gene, for each region, in each cell-type against the 
# corresponding GWAS data (since COLOC is binary)


cellCOLOC_obj=createCellCOLOCobject(mateqtlouts, processed_gwas, 
genestokeep,indiv_numbers = indiv_numbers,
maf_df=maf_df,gwas_type=GWAS_type,snp_locs=snp_locs,gene_locs=gene_locs)

cellCOLOC_obj@gwas <- remove_small_regions(cellCOLOC_obj@gwas,min_snps=min_snps)
cellTypes_names <- names(cellCOLOC_obj@cellTypes)
for (cellType in cellTypes_names) {
cellCOLOC_obj@cellTypes[[cellType]]@regions <- remove_small_regions(cellCOLOC_obj@cellTypes[[cellType]]@regions,min_snps=min_snps)
}
saveRDS(cellCOLOC_obj,paste0(GWAS_name,"_cellCOLOC_obj.rds"))



## run COLOC - see cellCOLOC_source.r for the function definition
coloc_df=run_coloc_all(cellCOLOC_obj)
coloc_df=coloc_df[order(coloc_df$PP.H4,decreasing=T),]
write.table(coloc_df,paste0(GWAS_name,"_COLOC_results.txt"))
message(paste0(Sys.time(),": cellCOLOC run complete."))
