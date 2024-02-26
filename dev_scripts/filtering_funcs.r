
#' This function filters the genotype matrix based on one of two criteria.
#' @param genotypemat The genotype matrix in 0,1,2 format. Rows are variants, columns are individuals.
#' @param filter_type The type of filter.
#' @export
filter_genotype_matrix=function(genotypemat,filter_type=c("2.2","2.3")){

  tmp<-genotypemat

  message(paste0("Filtering genotype matrix. Starting with ",nrow(genotypemat)," snps."))

  tmp$counter_0<-rowSums(tmp[1:ncol(tmp)]==0)
  tmp$counter_1<-rowSums(tmp[1:ncol(tmp)]==1)
  tmp$counter_2<-rowSums(tmp[1:ncol(tmp)]==2)

  if(filter_type=="2.2"){

    message("Retaining snps with at least 2 individuals in each of the 2 genotypic categories.")

    option1=tmp[tmp$counter_0>=2 & tmp$counter_1>=2,]
    option2=tmp[tmp$counter_0>=2 & tmp$counter_2>=2,]
    option3=tmp[tmp$counter_1>=2 & tmp$counter_2>=2,]

    intersected_rows<-unique(c(rownames(option1),rownames(option2),rownames(option3)))

  } else if (filter_type=="2.3"){

    message("Retaining snps with at least 2 individuals in each of the 3 genotypic categories.")

    tmp<-tmp[tmp$counter_0>=2 & tmp$counter_1>=2 & tmp$counter_2>=2,]

    intersected_rows<-row.names(tmp)


  }

  genotypemat<-genotypemat[intersected_rows,]

  message(paste0("Filtering complete. ",nrow(genotypemat), " snps retained."))

  return(genotypemat)
}

#' This function filters the expression matrix.
#' @param exp_mat The expression matrix. Rows are genes, columns are individuals.
#' @param minimum_indivs The filtering parameter. Minimum number of individuals a gene has to be expressed in to be included.
#' @export
filter_pseudobulk=function(exp_mat,minimum_indivs){


  exp_mat<-exp_mat[rowSums(exp_mat>0)>=minimum_indivs,]
  return(exp_mat)



}
