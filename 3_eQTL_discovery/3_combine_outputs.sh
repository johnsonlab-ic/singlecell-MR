



#PBS -l select=1:ncpus=20:mem=240gb
#PBS -l walltime=04:00:00

#for the full dataset
basedir=/rds/general/user/ah3918/projects/roche/ephemeral/SC_ANALYSIS/
cov_file="/rds/general/user/ah3918/projects/roche/live/ALEX/METADATA/MASTER_COVARIATE_MATRIX.txt"
full_sample_list="/rds/general/user/ah3918/projects/roche/live/ALEX/METADATA/full_samples.txt"
inputdir=${basedir}/MatrixEQTL_IO/
outdir=${basedir}/MatrixEQTL_IO/full_matrixeqtl_outputs 

module load anaconda3/personal
source activate OSIRIS 

cd $outdir
Rscript ${basedir}/singlecell-MR/3_eQTL_discovery/combine_outputs.r \
--outpath $outdir \


##now for the controls
outdir=${basedir}/MatrixEQTL_IO/controls_matrixeqtl_outputs 
cd $outdir
Rscript ${basedir}/singlecell-MR/3_eQTL_discovery/combine_outputs.r \
--outpath $outdir \