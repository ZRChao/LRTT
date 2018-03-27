\name{LNM.Sim}
\alias{LNM.Sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Logistical Normal Multinomial Distribution Simulation
}
\description{

}
\usage{
LNM.Sim(p, seed, N, dif)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{OTU numbers}
  \item{seed}{random seeds}
  \item{N}{sample size}
  \item{dif}{the differential OTU}
}

\details{
Return one count matrix of leaves and internal nodes which colnames corresponds to the tree.
}
\value{
One count matrix
}
\references{
Xia, Fan, et al. "A logistic normal multinomial regression model for microbiome compositional data analysis." Biometrics 69.4 (2013): 1053-1063.

}
\author{
Chou ZHOU
}

\examples{
p <- 10
dif <- c(2,3)

data <- LNM.Sim(p, seed =1, N =10 , dif)
dim(data)

OTU.dif <- data[, colnames(daya) == dif]
boxplot(OTU.dif[1:N, 1], OTU.dif[-c(1:N), 1])
}
