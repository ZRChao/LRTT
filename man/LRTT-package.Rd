\name{LRTT-package}
\alias{LRTT-package}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Differential Abundance Analysis of Microbiome Data Incorporating Phylogeny
}
\description{
\strong{LRTT} simplified of Log Ratio Tree Test which do differential abundance based on the
tree structure. The package here provides four data simulation function include with tree structure or not and test based on the tree with several other function to do analysis based on the tree.

    The complete list of functions can be displayed with library(help = LRTT).

    More information on LRTT can be found :
    \url{https://github.com/ZRChao/LRTT}
}

\details{
      Package: LRTT

      Type: R package

      Version: 0.1.0

      Date: 2018-03-25

      License: GPL-2
}

\author{
Chao ZHOU & Tao WANG

Maintainer : Chao ZHOU <Supdream8@sjtu.edu.cn>
}
\references{
Chao ZHOU, Tao WANG. (2018) Differential Abundance Analysis of Microbiome Data Incorporating Phylogeny (Under review).

Mandal, Siddhartha, et al. "Analysis of composition of microbiomes: a novel method for studying microbial composition." Microbial ecology in health and disease 26.1 (2015): 27663.


Xia, Fan, et al. "A logistic normal multinomial regression model for microbiome compositional data analysis." Biometrics 69.4 (2013): 1053-1063.


Torben Tvedebrink (2010). Overdispersion in allelic counts and theta-correction in forensic genetics. Theoretical Population Biology, 78(3), 200-210.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
library(MiSPU)
load(throat.otu.tab)
load(throat.tree)
load(throat.meta)

p = ncol(throat.otu.tab)
throat.taxa.index <- Taxa.index(p, throat.tree)
throat.alltab <- cbind(throat.taxa.index*throat.otu.tab,  throat.otu.tab)
group <- throat.meta$SmokingState

result <- Tree.Ratio(p, throat.tree, throat.taxa.index, throat.alltab, group)
throat.detected <- Tree.Ratio.back(p, throat.tree, throat.taxa.index, results, group)
}
