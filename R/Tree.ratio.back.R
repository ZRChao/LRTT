Tree.ratio.back = function(tree.ratio){
  p <- length(tree.ratio$otu.dif)
  group <- tree.ratio$group
  otutab <-tree.ratio$alltab[, which(as.numeric(colnames(tree.ratio$alltab))<=p)]
  taxa.index <- tree.ratio$taxa.index
  
  label <- unique(group)[1]
  otudif_TF <- tree.ratio$otu.dif
  otu_pvalue <- tree.ratio$otu.pvalue
  
  difset <- which(tree.ratio$otu.dif == T)
  
  if(length(difset) == 0){
    difset.detected <- NULL
  }
  else{
    nondifset <- setdiff(1:p, difset)

    if(length(nondifset) == 0){
      ratio_sum <- log(otutab[, difset] + 1e-20) - log(rowSums(otutab) + 1e-20)
    }

    if(length(nondifset) == 1 ){
      ratio_sum <- log(otutab[, difset] + 1e-20) - log(otutab[, nondifset] + 1e-20)
    }

    else{
      ratio_sum <- log(otutab[, difset] + 1e-20) - log(rowSums(otutab[, nondifset]) + 1e-20)
    }
    difset_sumpv <-  apply(as.matrix(ratio_sum), 2, function(x)
      return(ifelse(length(unique(x)) == 1, 1, t.test(x[group==label], x[group!=label])$p.value)))

    difset_sumpv.adj <- p.adjust(difset_sumpv, method = "fdr", n = length(difset_sumpv))                      
    otu_pvalue[difset] <- difset_sumpv.adj                       
    difset.detected <- names(which(difset_sumpv.adj <= 0.05))
  }
  if(length(difset.detected) !=0) {
     otudif_TF[as.numeric(difset.detected)] <- T ;
      otudif_TF[!as.numeric(difset.detected)]<- F  
   }
   else{
     otudif_TF[1:p] <- F
   }
   
  res <- data.frame(OTU=1:p, adj.pvalue=otu_pvalue, Differential=otudif_TF)                         
  return(res)
}

#############################################################################################

