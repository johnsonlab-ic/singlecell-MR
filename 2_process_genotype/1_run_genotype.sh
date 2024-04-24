#PBS -l select=1:ncpus=20:mem=240gb
#PBS -l walltime=02:00:00


vcf_dir=/rds/general/user/ah3918/projects/roche/live/ALEX/PROCESSED_DATA/PROCESSED_GENOTYPE/FINAL/
basedir=/rds/general/user/ah3918/projects/roche/ephemeral/SC_ANALYSIS/
outdir=${basedir}/MatrixEQTL_IO/

module load anaconda3/personal
source activate OSIRIS 

cd $outdir
###run genotype script

Rscript ${basedir}/singlecell-MR/2_process_genotype/generate_genotype_matrix.r \
--vcf_file ${vcf_dir}/final_geno_440samples_renamedsamples_sorted.vcf.gz \
--ncores 10 \
--script_dir ${basedir}/singlecell-MR/2_process_genotype/ \
