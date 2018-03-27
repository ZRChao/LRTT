\name{BIT.Sim}
\alias{BIT.Sim}
\title{Binomial Tree Simulation

}
\description{
Function to do simulation with multinomial distribution based on the tree structure.
}
\usage{
BIT.Sim(p, seed = 2, N, tree, prob1, prob2)
}

\arguments{
  \item{p}{OTU numbers}
  \item{seed}{random seeds}
  \item{N}{sample size}
  \item{tree}{tree structure with phylo types}
  \item{prob1, prob2}{the probability of each branch along the tree which both is one vector with length of the branches }
}
\details{
You can random set the probability of the branches, finally use function Prob.mult() to see which OTU on the leafs is differential one.
}
\value{
Return one count matrix of the OTU and internal nodes.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Chao ZHOU
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{Prob.mult}
}
\examples{}
