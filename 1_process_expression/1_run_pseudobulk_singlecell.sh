
#PBS -l select=1:ncpus=40:mem=480gb
#PBS -l walltime=16:00:00

seuratdir="/rds/general/user/ah3918/projects/roche/live/ALEX/ANALYSIS/SEURAT_OBJECTS/"
outdir="/rds/general/user/ah3918/projects/roche/ephemeral/SC_ANALYSIS/MatrixEQTL_IO/"


cd $outdir
###run pseudobulk script

Rscript pseudobulk_singlecell.r \
--input_seurat_files ${seuratdir}/BRYOIS_92_AD_final_annotated_decontXfiltered.rds \
    ${seuratdir}/BRYOIS_92_MS_final_annotated_decontXfiltered.rds \
    ${seuratdir}/2023-08-23_MRC_60_final_annotated_decontXfiltered.rds \
    ${seuratdir}/2023-08-23_ROCHE_MPD_92_final_annotated_decontXfiltered.rds \
    ${seuratdir}/2023-08-23_MATTHEWS_final_annotated_decontXfiltered.rds \
--min_cells 10 \
--indiv_column Individual_ID \
--celltype_column CellType \
--assay decontXcounts \


