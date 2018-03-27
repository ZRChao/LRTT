\name{DTM.Sim}
\alias{DTM.Sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dirichlet multinomial tree distribution simulation
}
\description{
Function to do simulation with dirichlet multinomial tree distribution based on the tree structure.
}
\usage{
DTM.Sim(p, seed = 1, N, tree, pi1, pi2, theta = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{OTU numbers}
  \item{seed}{random seeds}
  \item{N}{sample size}
  \item{tree}{tree structure with phylo types}
  \item{pi1, pi2}{the parameters of each branch along the tree which both is one vector with length of the branches}
  \item{theta}{dispersion parameter}
}
\details{
Like BIT.Sim()
}
\value{
Return one count matrix of the OTU and internal nodes.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Chao ZHOU
}
\note{
%%  ~~further notes~~
}


\seealso{
\link{BIT.Sim}
}
