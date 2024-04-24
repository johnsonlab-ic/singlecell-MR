
##this function creates a gene locations file (hg38) from the expression matrix

get_gene_locations=function(exp_mat){

  yourgenes<-as.character(rownames(exp_mat))


  # Create df and granges of all genes. These differ actually
  ucsc_genes<-TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  allgenes<-as.data.frame(org.Hs.eg.db::org.Hs.egSYMBOL)

  ## get ucsc sequences
  allgranges<-suppressMessages(GenomicFeatures::genes(ucsc_genes,single.strand.genes.only=FALSE,columns="gene_id"))
  allgranges<-unlist(allgranges)
  gene_ids<-names(allgranges)
  allgranges<-data.frame(allgranges)
  allgranges$gene_id<-gene_ids


  ##remove random sequences
  toremove<-c("alt","fix","random")
  seqnames_toremove<-c(grep(paste(toremove,collapse = "|"),allgranges$seqnames))
  allgranges<-allgranges[-seqnames_toremove,]

  ##filter symbol df
  commongenes<-intersect(allgenes$gene_id,allgranges$gene_id)
  allgenes<-allgenes[allgenes$gene_id %in% commongenes,]

  ##create final chrompos_mat
  allgenes$start<-allgranges$start[match(allgenes$gene_id,allgranges$gene_id)]
  allgenes$end<-allgranges$end[match(allgenes$gene_id,allgranges$gene_id)]
  allgenes$chr<-allgranges$seqnames[match(allgenes$gene_id,allgranges$gene_id)]

  ##filter for the genes we have
  allgenes<-allgenes[allgenes$symbol %in% yourgenes,]

  ##reorganize and return
  allgenes<-allgenes[,c("symbol","chr","start","end")]
  colnames(allgenes)<-c("geneid","chr","left","right")
  return(allgenes)


}

##this function takes a list of seurat objects (split by cell type) and pseudobulks them.
pseudobulk_counts=function(seuratlist,min.cells=100,indiv_col="Sample_ID",assay="RNA"){
  agg_count_list<-lapply(seuratlist,function(x){
  Seurat::DefaultAssay(x) <- assay
    metadata <- x[[]]

    counts <- Seurat::GetAssayData(x, slot = "counts")
    unique_ids <- unique(metadata[[indiv_col]])
    indiv_table <- metadata %>% dplyr::count(get(indiv_col))
    indiv_table <- indiv_table %>% dplyr::filter(n > min.cells)
    colnames(indiv_table) <- c("unique_ids", "n")
    unique_ids <- indiv_table$unique_ids
    n_indiv_start <- length(unique_ids)

    if (nrow(indiv_table) == 0) {
      # Return an empty DF if no individuals pass min.cells threshold.
      # This rarely happens when using a single seurat object, but when using a list 
      # of seurat objects it may be common (10 individuals per seurat object for example.)
      emptydf <- data.frame()[1:nrow(counts), ]
      rownames(emptydf) <- rownames(counts)
      return(emptydf)
    } else {
      cells <- rownames(metadata[which(metadata[[indiv_col]] == unique_ids[1]), ])
      counts_2 <- Matrix::rowSums(counts[, cells])
      finalmat <- data.frame(counts_2)
      samplenames <- unique_ids[1]

      if (length(unique_ids) > 1) {
        for (i in 2:length(unique_ids)) {
          cells <- rownames(metadata[which(metadata[[indiv_col]] == unique_ids[[i]]), ])
          sample <- unique_ids[[i]]
          counts_2 <- Matrix::rowSums(counts[, cells])
          finalmat <- cbind(finalmat, counts_2)
          samplenames <- c(samplenames, sample)
        }
      }
      colnames(finalmat) <- samplenames
      indiv_remain <- n_indiv_start - length(samplenames)    

        return(finalmat)
          }
    })
}

##this function takes a list of pseudobulk counts and normalizes them
normalize_pseudobulk=function(agg_count_list){
    agg_count_list<-lapply(agg_count_list,function(x){
    if(length(x)==0){
        return(NULL)
    }
    norm_x<-log2(edgeR::cpm(x)+1)
    return(norm_x)
  })

  return(agg_count_list)

}