source("dat_format.R")
source("MLE.R")
source("SIR_fit.R")

# Genetic Algorithm to optimize feature selection

n.feat <- 21
max.lev <- 5

# Chromosome to levels
ChromToLevel <- function(chrom) {
  
  ## Format sequence for each feature
  seqs <- sapply(1:n.feat, function(n) {
    chrom[(1 + (n - 1) * (max.lev + 1)):(n * (max.lev + 1))]
  }, simplify=F)
  sums <- unlist(lapply(seqs, function(x) {sum(x)}))
  
  if (!all(sums == 1)) {
    return(-100) # Unvalid format of chromosome
  }
  
  ## Get chromosome encoded levels
  lev <- unlist(lapply(seqs, function(x) {return(match(1, x) - 1)}))
  
  return(lev)
}

# Clustering result of a chromosome
ClustRes <- function(chrom) {
  
  ## Get chromosome encoded levels
  lev <- ChromToLevel(chrom)
  cat.age <- lev[length(lev)]
  age <- cat.age != 0
  lev <- lev[1:(length(lev) - 1)]
  fn <- feature.names[lev != 0]
  lev.feat <- lev[lev != 0]
  
  ## Compute categorical data and fit MLE
  feat <- ComputeCatDf(feature.names, rep(2, length(feature.names)), age=T, cat.age=2)
  obs <- feat[2:length(feat)] # Remove sdr_id col
  # theta0 uniform
  theta0 <- PriorTheta(c(rep(2, length(feature.names)), 4))
  # G0 uniform
  G0 <- matrix(data=0.5, nrow=dim(obs)[1], 2)
  # Compute estimates
  gom.est <- gomMLE(obs, G0, theta0, epsilon=1e-2, verbose=F)
  G.hat <- data.frame(sdr.id=feat$sdr.id, g.hat=gom.est$G.hat)
  colnames(G.hat) <- c("sdr.id", "g.low", "g.high")
  clust <- data.frame(sdr.id=G.hat$sdr.id, high=as.numeric(G.hat$g.high >= 0.5))
  clust <- clust[order(clust$sdr.id), ] # Make sure order is right
  hclust.res <- EvalHardClust(clust) # Metrics from clustering
  fclust.res <- EvalFuzzClust(G.hat)
  
  return(list(hclust.res=hclust.res, fclust.res=fclust.res, g=G.hat))
  
}

# Fitness of a chromosome (using ChromRes within did not work)
Fitness <- function(chrom) {
  
  ## Get chromosome encoded levels
  lev <- ChromToLevel(chrom)
  cat.age <- lev[length(lev)]
  age <- cat.age != 0
  lev <- lev[1:(length(lev) - 1)]
  fn <- feature.names[lev != 0]
  lev.feat <- lev[lev != 0]
  
  ## Compute categorical data and fit MLE
  feat <- ComputeCatDf(feature.names, rep(2, length(feature.names)), age=T, cat.age=2)
  obs <- feat[2:length(feat)] # Remove sdr_id col
  # theta0 uniform
  theta0 <- PriorTheta(c(rep(2, length(feature.names)), 4))
  # G0 uniform
  G0 <- matrix(data=0.5, nrow=dim(obs)[1], 2)
  # Compute estimates
  gom.est <- gomMLE(obs, G0, theta0, epsilon=1e-2, verbose=F)
  G.hat <- data.frame(sdr.id=feat$sdr.id, g.hat=gom.est$G.hat)
  colnames(G.hat) <- c("sdr.id", "g.low", "g.high")
  clust <- data.frame(sdr.id=G.hat$sdr.id, high=as.numeric(G.hat$g.high >= 0.5))
  clust <- clust[order(clust$sdr.id), ] # Make sure order is right
  fclust.res <- EvalFuzzClust(G.hat)
  
  return(as.numeric(fclust.res$cl.met))
  
}

# Compute initial chrom of certain length with max level
InitialChrom <- function(n.feat=19, max.lev=5) {
  if (max.lev < 2) {
    stop("Need at least max level of 2.")
  }
  seq0 <- rep(0, max.lev + 1)
  seq0[3] <- 1 # 2 levels
  res <- rep(seq0, n.feat)
}

# Compute random suggestions
Suggestions <- function(num=100, n.feat=21, max.lev=5) {
  if (max.lev < 2) {
    stop("Need at least max level of 2.")
  }
  sugg <- matrix(0, ncol=n.feat * (max.lev + 1), nrow=num)
  for (i in 1:num) {
    row <- c()
    for (j in 1:n.feat) {
      seq0 <- rep(0, (max.lev + 1))
      lev <- sample(1:(max.lev + 1), 1)
      seq0[lev] <- 1
      row <- c(row, seq0)
    }
    sugg[i, ] <- row
  }
  
  return(sugg)
}