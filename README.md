## Log Ratio Tree Test
Differential Abundance Analysis for Microbiome data Incorporating Phylogeny. We call this methods *LRTT*（ Log Ratio Tree Test）.

## Installation

```R
# Install the development version from GitHub
devtools::install_github("ZRChao/LRTT")
```
## Contents

### Simulation
* `BIT.Sim.R` for Multinomial(Binomial) Tree distribution
* `DTM.Sim.R` for Dirichlet Multinomial Tree distribution
* `LNM.Sim.R` for Logistical Normal Multinomial distribution 
* `ANCOM.Sim.R` for Poission distribution (parameters set follow as ANCOM paper)

### Tree relate function
* `Tree.Sim.R` : Tree simulation
* `Taxa.index.R` : relationship between leafs and internal nodes
* `Prob.mult.R` :calculate the probability of each leafs by multiple each probability along the branch
* `Tree.ratio.R` : Tree ratio test based on the tree structure
* `Tree.ratio.back.R`: Correct steps of the tree ratio test.

### Other function
* `pow.fdr.R` : calculate power and fdr with p.value 
* `Zig.pv.adj.R` : metagenomeSeq::fitFeature
*

### Require 
* ancom.R::ANCOM, metagenomeSeq::fitFeature, mvtnorm, gtools, rdirimult 

### Example
 
#You can see this in real data application of throat.R 

```R
library(MiSPU)
data(throat.otu.tab)
data(throat.tree)
data(throat.meta)
data(throat.taxa.index )

p <- ncol(throat.otu.tab)
# throat.taxa.index <- Taxa.index(p, throat.tree)
colnames(throat.otu.tab) <- as.character(1:p)

throat.taxa.tab <- as.matrix(throat.otu.tab) %*% throat.taxa.index
throat.alltab <- cbind(throat.taxa.tab,  throat.otu.tab)
group <- throat.meta$SmokingStatus


result <- Tree.ratio(p, throat.tree, throat.taxa.index, throat.alltab, group)

throat.detected <- Tree.ratio.back(p, throat.tree, throat.taxa.index, results, group)
```
