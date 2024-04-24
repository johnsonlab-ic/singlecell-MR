
#PBS -l select=1:ncpus=20:mem=240gb
#PBS -l walltime=24:00:00


basedir=/rds/general/user/ah3918/projects/roche/ephemeral/SC_ANALYSIS/
cov_file="/rds/general/user/ah3918/projects/roche/live/ALEX/METADATA/MASTER_COVARIATE_MATRIX.txt"
control_sample_list="/rds/general/project/roche/live/ALEX/METADATA/control_samples.txt"
inputdir=${basedir}/MatrixEQTL_IO/
outdir=${basedir}/MatrixEQTL_IO/controls_matrixeqtl_outputs 


module load anaconda3/personal
source activate OSIRIS 

mkdir -p $outdir
cd $inputdir

input files
for file in *pseudobulk.csv
do
  celltype="${file%_pseudobulk.csv}"
  echo "Processing file for cell type: $celltype"
  cd $outdir  # Change working directory to the new output directory
  Rscript ${basedir}/singlecell-MR/3_eQTL_discovery/run_eqtl_discovery.r \
    --script_dir ${basedir}/singlecell-MR/3_eQTL_discovery/ \
    --exp_mat ${inputdir}/${celltype}_pseudobulk.csv --geno_mat ${inputdir}/genotype_012mat.csv \
    --exp_loc ${inputdir}/${celltype}_gene_locations.csv --geno_loc ${inputdir}/snp_chromlocations.csv \
    --cov_file $cov_file \
    --filter_genotype_matrix --filter_pseudobulk_matrix \
    --exp_mat_thresh_percent 5 \
    --trans_eqtls false \
    --covs_to_include Age Sex PMI Sample_Source \
    --specify_samples_file $control_sample_list
  cd $inputdir  # Change back to the main output directory
done