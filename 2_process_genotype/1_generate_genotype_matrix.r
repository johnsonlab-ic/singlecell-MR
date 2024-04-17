
suppressMessages(library(argparse))
source("genotype_funcs.r")

parser$add_argument("--genotypefile", help = "genotype vcf file")

input_vcf_file=args$genotypefile

#this function will;
# 1. Filter out SNPs with MAF < 0.05, and keep autosomal SNPs only
# 2. Convert VCF to a .gds file, which is used by SeqArray
# 3. Create a genotype matrix (0,1,2), SNP X Individuals
# 4. Create auxiliary files (snp_locations,MAF/Allele freqs)
# It will also automatically output hg38 positions for the SNPs (and liftover from hg19 if necessary),
# as well as convert CHR:POS to Rsids. 


get_genotype_matrix(vcf=input_vcf_file,
                    preprocess=TRUE,
                    outdir="MatrixEQTL_IO/",
                    minmaf=0.05,
                    autosomalonly=TRUE)