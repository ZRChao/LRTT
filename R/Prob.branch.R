## function to simulate the probability along the tree branch

Prob.branch = function(tree, seed, dif.taxa){
  p <- min(tree$edge[, 1] - 1)
  prob1 <- rep(0, nrow(tree$edge))
  prob2 <- rep(0, nrow(tree$edge))
  for(j in (p + 1) : (p + tree$Nnode)){
    set.seed(j*seed)
    pi0 <- runif(1, 0.1, 0.9)
    prob1[which(tree$edge[, 1] == j)] <- c(pi0, 1 - pi0)
  }

  prob1 -> prob2   ###case_probability
  for(i in 1:length(dif.taxa)){
    parent <- dif.taxa[i]
    temp_diff <- which(tree$edge[, 1] == parent)
    prob2[temp_diff] <- prob1[rev(temp_diff)]
  }
  return(cbind(prob1, prob2))
}
