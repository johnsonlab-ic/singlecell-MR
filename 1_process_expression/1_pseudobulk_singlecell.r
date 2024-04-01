


suppressMessages(library(argparse))
source("expression_helper_funcs.r")

parser$add_argument("--input_seurat_files", help = "path to input seurat object",nargs="+")
parser$add_argument("--min_cells", help = "minimum cells to pseudobulk")
parser$add_argument("--indiv_column", help = "Column name containing Individual IDs")
parser$add_argument("--celltype_column", help = "Column name containing cell-type labels")
parser$add_argument("--assay", help = "Assay to use to pseudobulk. Normally, this is the 'counts' slot.")
args <- parser$parse_args(commandArgs(TRUE))


celltype_colname=as.character(args$celltype_column)
indiv_colname=as.character(args$indiv_column)
min.cells=as.numeric(args$min_cells)
assay=args$assay
input_seurat_list_files=args$input_seurat_files


###read in objects. Make sure that all objects have the same genes and celltypes. Assumes that all of these have the same assay
message("First, making sure that all Seurat objects have the same genes and celltypes..")

so=readRDS(input_seurat_list_files[[1]])
Seurat::DefaultAssay(so)=assay
commongenes=rownames(so)
celltypes_tokeep=unique(so[[]][,celltype_colname])
for(i in 2:length(input_seurat_list_files)){
    so=readRDS(input_seurat_list_files[[i]])
    Seurat::DefaultAssay(so)=assay
    genes=rownames(so)
    celltypes=unique(so[[]][,celltype_colname])
    celltypes_tokeep=intersect(celltypes,celltypes_tokeep)
    commongenes=intersect(genes,commongenes)

}

n_seurat_objs=length(input_seurat_list_files)
message(paste0(n_seurat_objs," seurat objects were provided. Reading in Seurat obj 1 .."))
seuratobj<-readRDS(input_seurat_list_files[[1]])
seuratobj=seuratobj[commongenes,]

celltypelist<-Seurat::SplitObject(seuratobj,split.by=celltype_colname)
celltypelist=celltypelist[celltypes_tokeep]
agg_count_list_full<-pseudobulk_counts(celltypelist,min.cells=min.cells,indiv_col=indiv_colname,assay=assay)

for(i in 2:length(input_seurat_list_files)){
message(paste0(" Reading in Seurat obj ",i,".."))
    seuratobj<-readRDS(input_seurat_list_files[[i]])
    seuratobj=seuratobj[commongenes,]
    celltypelist<-Seurat::SplitObject(seuratobj,split.by=celltype_colname)
    agg_count_list<-pseudobulk_counts(celltypelist,min.cells=min.cells,indiv_col=indiv_colname,assay=assay)
    agg_count_list<-agg_count_list[names(agg_count_list) %in% names(agg_count_list_full)]
    agg_count_list<-agg_count_list[names(agg_count_list_full)]
    agg_count_list_full<-Map(cbind,agg_count_list_full,agg_count_list)
}
nullvec=unlist(lapply(agg_count_list_full,is.null))
dropped_celltypes=names(nullvec[which(nullvec==TRUE)])
agg_count_list_full=Filter(Negate(is.null),agg_count_list_full)





agg_count_list_full=lapply(agg_count_list_full,function(df){df[, !is.na(names(df))]})

for(i in 1:length(agg_count_list_full)){
        write.table(agg_count_list_full[[i]],paste0(names(agg_count_list[i]),"_aggregated_counts.csv"))
    }

agg_count_list_full=normalize_pseudobulk(agg_count_list_full)


if(length(dropped_celltypes)>0){
    message(paste0(dropped_celltypes," celltype was dropped due to no individuals passing min.cells criteria."))
}


for(i in 1:length(agg_count_list_full)){
    write.table(agg_count_list_full[[i]],paste0(names(agg_count_list[i]),"_pseudobulk.csv"))
}