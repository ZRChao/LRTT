\name{Tree.ratio}
\alias{Tree.ratio}
\title{Tree.ratio}
\usage{
Tree.ratio(p, tree, taxa.index, alltab, group = c(1:N))
}
\arguments{
\item{p}{the number of OTU}

\item{tree}{the tree structure, phylo object}
\item{taxa.index}{matrix return by Taxa.index}
\item{alltab}{the count data matrix}

\item{group}{the class of each sample}
}
\value{
A list data
\item{alltab}{data after test by prune the tree}
\item{taxa.pvalue}{pvalue of the internal nodes}
\item{taxa.dif}{differential internal nodes}
\item{otu.dif}{differential OTU before corrected}

}
\description{
Tree ratio test based on the tree.
}
\seealso{
\code{\link{Taxa.index}}
}
\examples{
library(MiSPU)
data(throat.otu.tab)
data(throat.tree)
data(throat.meta)

p <- ncol(throat.otu.tab)
# throat.taxa.index <- Taxa.index(p, throat.tree)
data(throat.taxa.index )
colnames(throat.otu.tab) <- as.character(1:p)
throat.taxa.tab <- throat.otu.tab \%*\% throat.taxa.index
  
throat.alltab <- cbind(throat.taxa.tab,  throat.otu.tab)
group <- throat.meta$SmokingStatus

result <- Tree.ratio(p, throat.tree, throat.taxa.index, throat.alltab, group)

throat.detected <- Tree.ratio.back(p, throat.tree, throat.taxa.index, results, group)

}
