### R code from vignette source 'LRTT.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: preliminaries
###################################################
options(prompt = "> ")


###################################################
### code chunk number 2: LRTT-prelim (eval = FALSE)
###################################################
## # install.packages("devtools")
## library(devtools)
## install_github("ZRChao/LRTT")
## 
## library(help = LRTT)


###################################################
### code chunk number 3: LRTT.Rnw:71-78
###################################################
library(LRTT)
p <- 10
tree <- Tree.Sim(p)

#plot.phylo(tree, type = "phylogram")
taxa.index <- Taxa.index(p, tree)
taxa.index


###################################################
### code chunk number 4: Probability simulation
###################################################
set.seed(2)
dif.taxa <- sample(tree$edge[, 1], 1)
prob <- Prob.branch(tree, seed = 1, dif.taxa)

prob.m1 <- Prob.mult(p, tree, prob[, 1])
prob.m2 <- Prob.mult(p, tree, prob[, 2])

dif.otu <- which(prob.m1 != prob.m2)
dif.otu


###################################################
### code chunk number 5: Four simulation
###################################################
## BIT and DTM will return all count data
data.bit <- BIT.Sim(p, seed = 2, N = 20, tree = tree, prob[, 1], prob[, 2])
# data.dtm <- DTM.Sim(p, seed = 2, N = 20, tree, prob[, 1], prob[, 2], theta = 0.1)
dim(data.bit)
# dim(data.dtm)

## LNM and ANCOM will return only OTU count data
# data.lnm <- LNM.Sim(p, seed = 2, N = 20, dif.otu)
data.ancom <- ANCOM.Sim(p, seed = 2, N = 20, dif.otu)
# dim(data.lnm)
dim(data.ancom)

data.allancom <- cbind(data.ancom %*% taxa.index[, -1], data.ancom)
dim(data.allancom)


###################################################
### code chunk number 6: data given
###################################################
data("Sim.data")
summary(Sim.data)


###################################################
### code chunk number 7: Tree.Ratio
###################################################
data.bit <- as.matrix(Sim.data$BIT)
grouplabel <- Sim.data$group
dim(data.bit)

tree <- Sim.data$tree
p <- min(tree$edge[, 1]) -1
p

taxa.index <- Taxa.index(p, tree)
colnames(data.bit) <- as.character(1:p)
all.bit <- cbind(data.bit %*% taxa.index , data.bit)

tree.results <- Tree.ratio(p, tree, taxa.index, all.tab = all.bit ,
                           group = grouplabel)
str(tree.results)


###################################################
### code chunk number 8: Tree.Ratio.back
###################################################
tree.dectected <- Tree.ratio.back(p, tree.ratio = tree.results,
                      taxa.index, otutab = data.bit, group = grouplabel)
tree.dectected


