\name{ANCOM.Sim}
\alias{ANCOM.Sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ANCOM simulation}
\description{
Simulation data for ANCOM situation
}
\usage{
ANCOM.Sim(p, seed, N, dif)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{OTU numbers}
  \item{seed}{random seeds}
  \item{N}{sample size}
  \item{dif}{the differential OTU you can set}
}
\details{
The simulation here is the same from ANCOM papers which set each OTU distributed from Poission distribution.
}
\value{
Return one 2N*p matrix
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Mandal, Siddhartha, et al. "Analysis of composition of microbiomes: a novel method for studying microbial composition." Microbial ecology in health and disease 26.1 (2015): 27663.
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
data <- ANCOM.Sim(p = 50, N = 50, dif = c(1:10))
dim(data)
}

