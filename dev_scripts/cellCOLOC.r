



library(tidyverse)
library(argparse)
source("/rds/general/user/ah3918/projects/roche/live/ALEX/SCRIPTS/GITHUB/eQTL_scripts/nextflow_pipelines/singlecell/cellCOLOC_scripts/cellCOLOC_source.r")

parser <- ArgumentParser()


##### input files
###########################################################


parser$add_argument("--GWAS_path", help = "Expression matrix path")
parser$add_argument("--mateqtlouts_path", help = "Gene locations path")
parser$add_argument("--MatrixEQTL_IO_path", help = "Gene locations path")
parser$add_argument("--gwas_type", help = "either 'cc' or 'quant'")
parser$add_argument("--clumped", help = "Ge",action="store_true")
parser$add_argument("--clumped_file",help="clumped file path")
args <- parser$parse_args(commandArgs(TRUE))

GWAS_path=args$GWAS_path
mateqtlouts_path=args$mateqtlouts_path
mateqtl_io_path=args$MatrixEQTL_IO_path
gwas_type=args$gwas_type
clump=args$clumped

if(clump!=TRUE){
    preprocessgwas=TRUE
    run_name=basename(GWAS_path)
    preprocessgwas_file=paste0(run_name,"_clumped.rds")
}else{
    preprocessgwas=FALSE
    preprocessgwas_file=args$clumped_file
    run_name=basename(preprocessgwas_file)
}




cellCOLOC(mateqtlouts=mateqtlouts_path,GWAS=GWAS_path,
  gene_loc_file=list.files(pattern="gene_locations",mateqtl_io_path,full.names=TRUE)[1],
  maf_file=paste0(mateqtl_io_path,"/MAF_mat.csv"),
  snp_locs_file=paste0(mateqtl_io_path,"/snp_chromlocations.csv"),
  indiv_numbers_file=paste0(mateqtl_io_path,"//n_indivs_pseudobulk.txt"),
  GWASsignif=5e-8,
  GWAS_window=1e6,
  preprocess_GWAS=TRUE,
  GWAS_type=gwas_type,
  GWAS_name=run_name,
  use_coloc_lead_snp=TRUE,
  cellMR=TRUE,
  pph4_cutoff=0.8,
)

