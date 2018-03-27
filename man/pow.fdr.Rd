\name{pow.fdr}
\alias{pow.fdr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power and fdr calculation
}
\description{
Give the p.value vector and true positive sets, return the power and fdr
}
\usage{
pow.fdr(p.value, dif, cutf = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p.value}{p.value adjusted or not}
  \item{dif}{true positive}
  \item{cutf}{significant level and default value is 0.05}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Return power and fdr of given p.value
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
set.seed(1)
p.value <- runif(10, 0, 0.1)
true <- c(1,2,3)
pow.fdr(p.value, true, 0.05)
}
