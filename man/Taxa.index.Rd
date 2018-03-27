\name{Taxa.index}
\alias{Taxa.index}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Relationship between leaves and internal nodes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Taxa.index(p, tree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{OTU numbers}
  \item{tree}{tree structure, object of phylo type}
}
\details{
This function have the same function with caper::clade.members.list
}
\value{
One matrix with dimension p*s, s is the number of the internal nodes
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
p <- 6
tree <- Tree.Sim(p)
otu.tab <- matrix(1:100, p)
colnames(otu.tab) <- as.character(1:p)
taxa <- taxa.index(p, tree)

##get the count on the internal nodes
taxa.otu <- taxa * otu.tab
}
