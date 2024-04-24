
library(argparse)
library(dplyr)
parser <- ArgumentParser()


parser$add_argument("--vcf_file", help = "genotype vcf file") # change this line
parser$add_argument("--ncores", help = "number of cores for parallel processing of vcf file")
parser$add_argument("--script_dir", help = "Directory containing helper functions")
args <- parser$parse_args(commandArgs(TRUE))

script_dir=args$script_dir
input_vcf_file=args$vcf_file # change this line
ncores=as.numeric(args$ncores)

message(script_dir)
source(file.path(script_dir, "genotype_funcs.r")) 
#this function will;
# 1. Filter out SNPs with MAF < 0.05, and keep autosomal SNPs only
# 2. Convert VCF to a .gds file, which is used by SeqArray
# 3. Create a genotype matrix (0,1,2), SNP X Individuals
# 4. Create auxiliary files (snp_locations,MAF/Allele freqs)
# It will also automatically output hg38 positions for the SNPs (and liftover from hg19 if necessary),
# as well as convert CHR:POS to Rsids. 


get_genotype_matrix(vcf=input_vcf_file,
                    gds_file="genotype.gds",
                    preprocess=TRUE,
                    parallel=ncores,
                    minmaf=0.05,
                    autosomalonly=TRUE)