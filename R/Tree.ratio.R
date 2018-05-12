#############################################################################
###            Tree.Ratio based on the tree structure                   #####
### simulate the tree structure based on the pacakge ape                #####
### which need define the OTU numbers p and tip.labed is 1:p            #####
### all the edge.length equals to 1                                     #####
#############################################################################

#-----------------------------------------------------------------------#####


Tree.ratio = function(p, tree, taxa.index, all.tab, group){
  result <- list()
  taxa_pv <- rep(0, ncol(taxa.index))
  names(taxa_pv) <- names(sort(colSums(taxa.index)))

  otu_dif <- rep(F, p)
  otu_pvalue <- rep(0, p)
  names(otu_dif) <- as.character(1:p)
  names(otu_pvalue)<-as.character(1:p)

  taxa_leafs <- unique(sort(colSums(taxa.index)))
  label <- unique(group)[1]
  taxa.pvname <- c()
  taxa.dif <- c()
  n = 0
  for(t in taxa_leafs){
    taxa_names <- names(which(colSums(taxa.index) == t))
    for(j in 1:length(taxa_names)){
      n <- n + 1
      parent <- taxa_names[j]
      taxa.pvname[n] <- parent
      childs <- tree$edge[which(tree$edge[, 1] == parent), 2]
      child1 <- all.tab[, colnames(all.tab) == as.character(childs[1])]
      child2 <- all.tab[, colnames(all.tab) == as.character(childs[2])]

      ratio = log(child1 + 1e-20) - log(child2 + 1e-20)
      taxa_pv[n] <- pv <- ifelse(length(unique(ratio)) == 1, 1,
                                 t.test(ratio[group==label], ratio[group!=label])$p.value)
      childs <- which(taxa.index[, colnames(taxa.index) == parent] == 1)
      otu_pvalue[childs] <- pv*length(taxa_names)
      if(pv <= 0.05/length(taxa_names)){
        #temp = all_table[, colnames(all_table) == parent]
        all.tab[, colnames(all.tab) == parent] <- 1
        taxa.dif <- append(taxa.dif, parent)

        if(t != p){
          otu_dif[childs] <- T
        }
        else{
          break
        }
      }
      else{
        all.tab[, which(colnames(all.tab) == parent)] <- child1 + child2
      }
    }
  }
  # names(taxa_pv)<-taxa.pvname
  result$taxa.pvalue <- taxa_pv
  result$otu.dif <- otu_dif
  result$otu.pvalue <- otu_pvalue
  result$alltab <- all.tab
  result$taxa.dif <- taxa.dif
  return(result)
}

###########################################################################################

