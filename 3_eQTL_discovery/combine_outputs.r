library(SCGsuite)
library(argparse)

parser <- ArgumentParser()
# parser$add_argument("--mateqtlouts_list", help = "output path for MatrixEQTL")
parser$add_argument("--outpath", help = "output path for MatrixEQTL")


args <- parser$parse_args(commandArgs(TRUE))

# mateqtlouts_list=args$mateqtlouts_list
outpath=args$outpath

print(list.files(outpath,pattern="_cis_MatrixEQTLout.rds"))

summarize_outputs_cis=function(outpath){
    
   mateqtlouts_list=list.files(outpath,pattern="_cis_MatrixEQTLout.rds") 
  cellnames<-gsub("_cis_MatrixEQTLout.rds","",mateqtlouts_list)
  mateqtlouts_list=list.files(outpath,pattern="_cis_MatrixEQTLout.rds",full.names=TRUE) 
  message("Reading in MatrixEQTL outputs..")
  mateqtlouts<-lapply(mateqtlouts_list,readRDS)
  names(mateqtlouts)<-cellnames

  message("Saving MatrixEQTL outputs as list object..")
  saveRDS(mateqtlouts,"mateqtlouts.rds")

  mateqtlouts<-lapply(mateqtlouts,function(x){
    x<-x[x$FDR<0.2,]
    return(x)
  })
  saveRDS(mateqtlouts,"mateqtlouts_0.2FDR.rds")
}


summarize_outputs_cis(outpath)