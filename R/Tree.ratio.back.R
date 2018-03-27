Tree.ratio.back = function(p, tree.ratio, taxa.index, otutab, group=c(1:(2*N))){

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
      return(ifelse(length(unique(x)) == 1, 1, t.test(x[group==1], x[group==2])$p.value)))

    difset_sumpv.adj <- p.adjust(difset_sumpv, method = "fdr", n = length(difset_sumpv))
    difset.detected <- names(which(difset_sumpv.adj <= 0.05))
  }
  return(difset.detected)
}

#############################################################################################

