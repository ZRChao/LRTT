\name{Sim.data}
\alias{Sim.data}
\docType{data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Four Simulation type data based on one tree
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
When given one tree, we can simulation count data with this structure, here is one sample case for you to do some test for the function on \strong{LRTT}.
}
\usage{
data(Sim.data)
}
\format{
  The format is:
 chr "Sim.data"
}

%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{tree}{one tree sturucture with phylo class}
   \item{BIT}{Binominal Tree Distribution}
   \item{DTM}{Dirichlet Multinomial Tree Distribution}
   \item{LNM}{Logistical Normal Multinomial Distribution}
   \item{ANOCM}{Each OTU follow Poission Distribution}
   \item{dif.otu}{Differntial OTU on the leaves}
   \item{dif.taxa}{Differential internal nodes}
   \item{prob.branch}{probability along the branch which corrseponds to the tree edges}
   \item{group}{the subject group lable}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
data(Sim.data)
tree <- Sim.data$tree
data.bit <- Sim.data$BIT
p <- dim(BIT)[2]

##
taxa <- Taxa.index(p, tree)
taxa.count <- taxa\%*\%data.bit
}
\keyword{datasets}
