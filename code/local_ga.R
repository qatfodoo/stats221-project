library(GA)
source("ga.R")

## Try GA
gen.alg <- ga(type=c("binary"), fitness=Fitness, nBits=n.feat * (max.lev + 1),
              maxiter=150, suggestions=Suggestions(10))