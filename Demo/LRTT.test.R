# test code

library(LRTT)
p=10
tree <- Tree.Sim(p)
plot.phylo(tree, type = "phylogram")

taxa <- Taxa.index(p, tree)
taxa

p.b <- Prob.branch(tree, 2, 15)
cbind(tree$edge, p.b)

p.otu1 <- Prob.mult(p, tree, p.b[, 1])
p.otu2 <- Prob.mult(p, tree, p.b[, 2])
sum(p.otu1)
sum(p.otu2)
dif <- which(p.otu1 != p.otu2)

bit <- BIT.Sim(p,2, N=10, tree, p.b[, 1], p.b[, 2])
dim(bit)

dtm <- DTM.Sim(p, 2, N=10, tree, p.b[, 1], p.b[, 2], 0.1)
dim(dtm)

lnm <- LNM.Sim(p, 2, N=10, dif)
dim(lnm)
t.count <- lnm %*% taxa

ancom <- ANCOM.Sim(p, 2, 10, dif)
dim(ancom)

group <- c(rep(1, 10), rep(2, 10))

tree.res <- Tree.ratio(p, tree, taxa, bit, group)
str(tree.res)
which(tree.res$otu.dif == T)

tree.detected <- Tree.ratio.back(p, tree.res, taxa, bit[, colnames(bit)<=p], group)
tree.detected

