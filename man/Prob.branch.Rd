\name{Prob.branch}
\alias{Prob.branch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probability along the tree branch
}
\description{
Each internal node will have his children's probability sum to 1.
}
\usage{
Prob.branch(tree, seed, dif.taxa)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{tree which is phylo class}
  \item{seed}{random seed}
  \item{dif.taxa}{different internal nodes}
}
\details{
Along the tree structure, simulate the probability on each nodes which satisfied that sum to 1 on each internal node.
}
\value{
Return two vector of two groups
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Prob.mult}}
}
\examples{
}
