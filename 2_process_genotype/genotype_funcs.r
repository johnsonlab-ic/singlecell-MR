

chrpos_to_rsid=function(chromlocs,build=c("hg37","hg38")){

  # chromlocs is a character vector, every element is
  # "chr1:10000","chr2:1201313" etc...
  original_names<-chromlocs
  chromlocs<-strsplit(chromlocs,":")
  chromlocs<-as.data.frame(do.call(rbind,chromlocs))
  chromlocs<-chromlocs[,1:2]
  colnames(chromlocs)<-c("Chr","Pos")

  if(length(grep("chr",chromlocs$Chr))!=0){
      chronly<-strsplit(as.character(chromlocs$Chr),"r")
      chronly<-as.data.frame(do.call(rbind,chronly))
        chromlocs$Chr<-chronly$V2
  }
  

  final<-data.frame(chrom=chromlocs[,1],position=chromlocs$Pos)
  final<-final[order(final$chrom),]
  final$paste<-paste0(final$chrom,":",final$position)



  grSNPS<-GenomicRanges::makeGRangesFromDataFrame(final,
  seqnames.field="chrom",
  start.field="position",
  end.field="position")


  if(build=="hg37"){
    

    snp<-SNPlocs.Hsapiens.dbSNP144.GRCh37::SNPlocs.Hsapiens.dbSNP144.GRCh37
  } else if(build=="hg38"){

    check_version=grep("155",system.file(package="SNPlocs.Hsapiens.dbSNP155.GRCh38"))
    if(length(check_version)==1){
      snp<-SNPlocs.Hsapiens.dbSNP155.GRCh38::SNPlocs.Hsapiens.dbSNP155.GRCh38
    }else{
      snp<-SNPlocs.Hsapiens.dbSNP151.GRCh38::SNPlocs.Hsapiens.dbSNP151.GRCh38
    }
  }

  rsids<-BSgenome::snpsByOverlaps(snp,grSNPS)
  rsids<-as.data.frame(rsids)
  rsids$paste<-paste0("chr",rsids$seqnames,":",rsids$pos)
  chromlocs$paste<-paste0("chr",chromlocs$Chr,":",chromlocs$Pos)

  chromlocs$original_names<-original_names
  chromlocs$rsid<-rsids$RefSNP_id[match(chromlocs$paste,rsids$paste)]
  chromlocs<-chromlocs %>% dplyr::mutate(rsid=dplyr::coalesce(rsid,original_names))


  return(chromlocs$rsid)
}



check_snps=function(chrompos_mat){

  nsnps<-nrow(chrompos_mat)
  chrompos_mat$old_snp=chrompos_mat$annot
  chromlocs<-chrompos_mat[,c("chrom","position","annot")]
  colnames(chromlocs)<-c("Chr","Pos","old_snp")
  chromlocs<-as.data.frame(chromlocs)
  chromlocs$Chr<-as.data.frame(do.call(rbind,strsplit(as.character(chromlocs$Chr),"r")))[,2]
  chromlocs$paste<-paste0(chromlocs$Chr,":",chromlocs$Pos)
  
  chromlocs_GR<-GenomicRanges::makeGRangesFromDataFrame(chromlocs,
      seqnames.field="Chr",
      start.field="Pos",
      end.field="Pos")
  GenomeInfoDb::seqlevelsStyle(chromlocs_GR)="NCBI"

  #check hg19 and hg38 and see which overlaps. If low overlap with hg38, liftover.

  hg19snp<-SNPlocs.Hsapiens.dbSNP144.GRCh37::SNPlocs.Hsapiens.dbSNP144.GRCh37    
  hg19rsids<-BSgenome::snpsByOverlaps(hg19snp,chromlocs_GR)
  check_version=grep("155",system.file(package="SNPlocs.Hsapiens.dbSNP155.GRCh38"))
    if(length(check_version)==1){
      hg38snp<-SNPlocs.Hsapiens.dbSNP155.GRCh38::SNPlocs.Hsapiens.dbSNP155.GRCh38
    }else{
      hg38snp<-SNPlocs.Hsapiens.dbSNP151.GRCh38::SNPlocs.Hsapiens.dbSNP151.GRCh38
    }

  hg38rsids<-BSgenome::snpsByOverlaps(hg38snp,chromlocs_GR)

  if(nrow(as.data.frame(hg19rsids))>nrow(as.data.frame(hg38rsids))){
   message("Detected genome build as hg19. Lifting over..")
    
    GenomeInfoDb::seqlevelsStyle(chromlocs_GR)="UCSC"
    chain_fpath=system.file("data", "hg19ToHg38.over.chain", package="SCGsuite")
    chain<-rtracklayer::import.chain(chain_fpath)

    liftover_GR<-rtracklayer::liftOver(chromlocs_GR,chain)
    liftover_GR<-unlist(liftover_GR)
    
    GenomeInfoDb::seqlevelsStyle(chromlocs_GR)="NCBI"
    GenomeInfoDb::seqlevelsStyle(liftover_GR)="NCBI"
    
    hg38rsids<-BSgenome::snpsByOverlaps(hg38snp,liftover_GR)

    hg38rsids<-as.data.frame(hg38rsids)
    hg19rsids<-as.data.frame(hg19rsids)
    hg19rsids$paste<-paste0(hg19rsids$seqnames,":",hg19rsids$pos)
    chromlocs<-chromlocs[match(hg19rsids$paste,chromlocs$paste),]
    chromlocs$rsid<-hg19rsids$RefSNP_id
    chromlocs<-chromlocs[match(hg38rsids$RefSNP_id,chromlocs$rsid),]
    chromlocs$pos<-hg38rsids$pos
    chromlocs<-chromlocs[,c("rsid","Chr","pos","old_snp")]
    colnames(chromlocs)<-c("annot","chrom","position","old_snp")
    chromlocs$chrom<-paste0("chr",chromlocs$chrom)
    chromlocs$paste<-paste0(chromlocs$chrom,":",chromlocs$position)
    
    chromlocs<-chromlocs[!duplicated(chromlocs$paste),]
    chromlocs<-chromlocs[complete.cases(chromlocs),]
    
    return(chromlocs)



   
  }else{

      rsids=chrpos_to_rsid(chromlocs$paste,build="hg38")
      chrompos_mat$annot=rsids
      return(chrompos_mat)

      
  }

}


get_genotype_matrix=function(vcfs,
  gds_file="merge_test_seqArray.gds",
  outdir=".",
  parallel=FALSE,
  preprocess=FALSE,
  autosomalonly=TRUE,
  minmaf=0.05){



  if(preprocess==TRUE){
    SeqArray::seqVCF2GDS(vcfs, gds_file,parallel=parallel)
  }

  genofile<<-SeqArray::seqOpen(gds_file)

  if(autosomalonly==TRUE){
    SeqArray::seqSetFilterChrom(genofile,1:22)
  }

  SeqArray::seqSetFilterCond(genofile,maf=minmaf)

  #get snp info
  annot<-SeqArray::seqGetData(genofile,"annotation/id")

  #get sample names
  sample.id<-SeqArray::seqGetData(genofile,"sample.id")

  #get genotype dosage matrix (0,1,2)
  geno_mat<-t(SeqArray::seqGetData(genofile,"$dosage"))

  #create SNP locations file
  position<-SeqArray::seqGetData(genofile,"position")
  chrom<-SeqArray::seqGetData(genofile,"chromosome")
  chrompos_mat<-data.frame(annot,chrom,position)
  chrompos_mat<-chrompos_mat[!duplicated(chrompos_mat$annot),]

  #clean / reformat
  colnames(geno_mat)<-sample.id
  geno_mat<-as.data.frame(geno_mat)
  geno_mat$snp<-annot
  geno_mat<-geno_mat[!duplicated(geno_mat$snp),]
  rownames(geno_mat)<-geno_mat$snp

  #
  chrompos_mat$chrom<-paste0("chr",chrompos_mat$chrom)
  snpnumber<-length(rownames(geno_mat))

  message("Checking snps for build and converting to rsids..")

  
  chrompos_mat=check_snps(chrompos_mat)
  geno_mat=geno_mat[match(chrompos_mat$old_snp,geno_mat$snp),]
  geno_mat$snp=chrompos_mat$annot
  rownames(geno_mat)=geno_mat$snp
  message(paste0("A total of ",snpnumber," snps were kept."))

  data.table::fwrite(geno_mat,"genotype_012mat.csv")
  

  #add alt. Allele Freq
  allele<-SeqArray::seqGetData(genofile,"allele")
  allele<-strsplit(allele,",")
  allele<-as.data.frame(do.call(rbind,allele))


    
  #check if MAF info is available. If not, calculate it
  maf_avail<-tryCatch({
        out<-SeqArray::seqGetData(genofile,"annotation/info/MAF")
    },
     error=function(e){
        out<-"error"
         return(out)
    }
  )

  if(length(maf_avail)>1){
      maf<-SeqArray::seqGetData(genofile,"annotation/info/MAF")
      }else if(maf_avail=="error"){
      maf<- SeqArray::seqAlleleFreq(genofile, minor=TRUE)
  }
  
  af_df<-data.frame(ref=allele$V1,
    alt=allele$V2,
    maf=maf,
    snp=annot)
    af_df$ref<-as.character(af_df$ref)
    af_df$alt<-as.character(af_df$alt)
    af_df$snp<-as.character(af_df$snp)


    af_df<-af_df[match(chrompos_mat$old_snp,af_df$snp),]
    af_df$snp<-chrompos_mat$annot



      #remove extra columns
  chrompos_mat=chrompos_mat[,c("annot","chrom","position")]
  data.table::fwrite(chrompos_mat,"snp_chromlocations.csv")



  data.table::fwrite(af_df,"MAF_mat.csv")

  SeqArray::seqClose(genofile)
}
