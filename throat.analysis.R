###-----------------Ratio test incroporate the phylogeney information apply to throat data----------###
##
##  This code main to compare different methods including (t.test(relative), zig(metagenomeSeq),
## ANCOM(ancom.R), tree.ratio(Tree.Ratio) and NonTree.Ratio (Same methods just like ANCOM) on the throat
## data from Chong WU (https://github.com/ChongWu-Biostat/MiSPU/tree/master/data)

# library nacessary packages and load the data
rm(list = ls())
library(MiSPU)
library(ancom.R)
library(metagenomeSeq)
library(LRTT)
source("~/Rcode/Zig.pv.adj.R")
data("throat.otu.tab")
data("throat.meta")
data("throat.tree")

###------------------------------------------- pre-process ----------------------------------------###
## this is a binary tree which is suitable for our methods
## Because the OTU numbers is a little big, so here take 0.1 as the significant level
threshold <- 0.1
otu.number <- ncol(throat.otu.tab)
is.binary.tree(throat.tree)
colnames(throat.otu.tab) <- as.character(1:otu.number)
throat.group <- which(throat.meta[,3] == "Smoker")


###-------------------- t.test of the relative otu.table data for differential depth---------------###
throat.log.rel.data = log((throat.otu.tab + 1e-20)/rowSums(throat.otu.tab + 1e-20))
throat.t.rel.pvalue <- apply(throat.log.rel.data, 2, function(x)
  return(ifelse(length(unique(x)) == 1, 1, t.test(x[throat.group], x[-throat.group])$p.value)))
throat.t.rel.pvalue.adj <- p.adjust(throat.t.rel.pvalue, method = "fdr", n = otu.number)
throat.t.rel.detected <-  names(which(throat.t.rel.pvalue.adj <= threshold))


###-------------------------- ZIG which is use the quantile normalization -------------------------###
rownames(throat.otu.tab) <- as.character(1:nrow(throat.otu.tab))
throat.zig.group <- AnnotatedDataFrame(data.frame(throat.meta$SmokingStatus))
throat.zig.otu.tab <- newMRexperiment(t(throat.otu.tab), phenoData = throat.zig.group)
throat.zig.otu.norm <- cumNorm(throat.zig.otu.tab, p = cumNormStatFast(throat.zig.otu.tab))
throat.zig.model <- model.matrix(~ 1 + throat.meta$SmokingStatus,
                                 data = pData(throat.zig.otu.norm))
throat.zig.fit <- fitFeatureModel(throat.zig.otu.norm, throat.zig.model)
throat.zig.result <- MRcoefs(throat.zig.fit, number = otu.number)
throat.zig.detected <- rownames(throat.zig.result[which(throat.zig.result[, 4] <= threshold ), ])


###----------------- ANCOM which take logarithm ratio transform for each OTU ---------------------###
## Have to infrom that use ancom.R::ANCOM which in this data will consume several minutes, just be
## PATIENT!!! (Good Luck to you for more OTU numbers data use this function), You can also try to
## use ancomP in python which is much more faster.
## Reference :1) https://www.niehs.nih.gov/research/resources/software/biostatistics/ancom/index.cfm
##            2) https://github.com/alk224/akutils-v1.2/blob/master/akutils_resources/R-instructions_ancom.r
start = Sys.time()
throat.ancom.otu.tab <- cbind(throat.otu.tab, throat.meta$SmokingStatus)
throat.ancom.result <- ANCOM(throat.ancom.otu.tab, sig = 0.05, multcorr = 2)
Sys.time() - start
throat.ancom.detected <- throat.ancom.result$detected

# here you can try to use the python code to do this
# Reference: https://github.com/mortonjt/ancomP
write.table(throat.otu.tab, "/Users/zhouchao/Documents/Tree/OTU/throat.otu.tab.txt",
            row.names = F, col.names = F)
write.table(t(throat.meta$SmokingStatus), "/Users/zhouchao/Documents/Tree/OTU/throat.group.txt",
            row.names = F, col.names = F, sep = ",")
throat.ancom.py.result <- read.table("/Users/zhouchao/Documents/Tree/throat/ancom.result.csv",
                                     sep = ",", head = T)
throat.ancom.py.delected <- which(throat.ancom.py.result[, 3] == "True")


###---------------------------- Tree Ratio test with correction ----------------------------------###
throat.taxa.index <- Taxa_index(p = otu.number, phy_tree = throat.tree)
throat.taxa.table <- as.matrix(throat.otu.tab)%*%throat.taxa.index
throat.all.table <- cbind(as.matrix(throat.otu.tab), throat.taxa.table)
throat.tree.ratio <- Tree_ratio(p = otu.number, tree = throat.tree, all_table = throat.all.table,
                               group = throat.group, taxa_index = throat.taxa.index)
throat.taxa.diff <- throat.tree.ratio$taxa.dif
throat.tree.back <- Tree_ratio_back(p = otu.number, tree_ratio = throat.tree.ratio,
                                   otutab = throat.otu.tab, group = throat.group)
throat.tree.back.adj <- p.adjust(throat.tree.back, method = "fdr", n = length(throat.tree.back))
throat.tree.ratio.detected <- names(which(throat.tree.back.adj < 0.1))


###------------finally we can get the differential OTU detected by this methods--------------------###
## we can use the FWER
## (Family-Wise Error Rate) to compare the type I error of this methods. So, we permutate the label of
## the Smoking Status 1000 times to calculate it, which show that ANCOM always have some detected which 
## is the type I error, so we would like to infrom you be careful to the results from ANCOM.
