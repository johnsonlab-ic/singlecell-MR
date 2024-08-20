
##load helper functions from cellCOLOC_source.r. Mostly for cellCOLOC_obj s4 structure and cellMR function
source("../4_colocalisation/cellCOLOC_source.r")


#define inputs /parameters

GWAS_name="AD_2022"
cellCOLOC_obj=readRDS(paste0(GWAS_name,"_cellCOLOC_obj.rds"))
coloc_res=read.table(paste0(GWAS_name,"_COLOC_results.txt"))
eqtl_FDR_cutoff=0.05
pph4_cutoff=0.8
r2_cutoff=0.01
use_coloc_lead_snp=TRUE
plink_bin=genetics.binaRies::get_plink_binary()
path_to_binaries="/path_to_binaries/EUR/EUR"


mr_results=run_cellMR(cellCOLOC_obj,use_coloc_lead_snp=use_coloc_lead_snp,
coloc_df,
plink_bin=plink_bin,
path_to_binaries=path_to_binaries,
pph4_cutoff=pph4_cutoff,
eqtl_FDR_cutoff=0.05,
r2_cutoff=0.01)


if(length(mr_results)>1){
    write.table(mr_results,paste0(GWAS_name,"_MR_results.txt"))
}

mr_res_PCA_999 = run_cellMR_IVPCA(cellCOLOC_obj = cellCOLOC_obj, 
coloc_res = coloc_df, pph4_cutoff = pph4_cutoff, 
eqtl_FDR_cutoff = 0.05, percentage_variance = 0.999)


write.table(mr_res_PCA_999, paste0(GWAS_name, "_MR_IVPCA_0.999_results.txt"))