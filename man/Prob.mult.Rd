\name{Prob.mult}
\alias{Prob.mult}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probability of the leaves
}
\description{
Given the probability on each branch, calculate the probability on the leaves along the tree branch, and then we can calculate the differential leaves.
}
\usage{
Prob.mult(p, tree, prob)
}
\arguments{
  \item{p}{OTU numbers}
  \item{tree}{tree structure, object of phylo type}
  \item{prob}{probability on the branches, a vector with length of branches}
}
\details{
On the data simulation, after you can set the probability along the branch, you need judge the probability on the leaves, then you can judge which one is differntial.
}
\value{
Return one vector of the leaves probability
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Chao ZHOU
}
\note{
Sum of the results must to 1
}


