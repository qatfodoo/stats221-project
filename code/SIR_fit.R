library(SDMTools)

# Assess classification results with actual SIR model fits

sir.fit <- read.csv(file="../dat/SIR_fit.csv")
r.sq <- data.frame(sdr.id=sir.fit$sdr_id, rsq=sir.fit$r_squared)
r.sq <- r.sq[order(r.sq$sdr.id), ]

# Evaluate a hard clustering results given SIR fits
EvalHardClust <- function(clust) {
  
  # Obtain relevant classification in right order
  cl <- clust[which(clust$sdr.id %in% r.sq$sdr.id), ]
  if (!all(r.sq$sdr.id == cl$sdr.id)) { # Check id considered are the same (ordering)
    stop("SDR ids should match with SIR fit results.")
  }
  
  rsq.low <- r.sq$rsq[cl$high == 0]
  rsq.high <- r.sq$rsq[cl$high == 1]
  # Remove NA fit values
  rsq.low <- rsq.low[!is.na(rsq.low)]
  rsq.high <- rsq.high[!is.na(rsq.high)]
  # Cluster centroids
  low.centr <- mean(rsq.low)
  high.centr <- mean(rsq.high)
  
  # Inter-centroidal separation
  inter <- abs(high.centr - low.centr)
  # Intra-cluster variance
  var.low <- var(rsq.low - low.centr)
  var.high <- var(rsq.high - high.centr)
  intra <- (var.low + var.high) / 2 # Mean intra variance
  
  res <- list(cl.met=(inter / intra), low.centr=low.centr, high.centr=high.centr,
              var.low=var.low, var.high=var.high)
    
}

# Evaluate a fuzzy clustering results given SIR fits
EvalFuzzClust <- function(g.hat) {
  
  # Remove NA sdrid
  rsq.clean <- r.sq[which(!is.na(r.sq$rsq)), ]
  # Obtain relevant classification in right order
  g <- g.hat[which(g.hat$sdr.id %in% rsq.clean$sdr.id), ]
  if (!all(rsq.clean$sdr.id == g$sdr.id)) { # Check id considered are the same (ordering)
    stop("SDR ids should match with SIR fit results.")
  }
  
  # Cluster centroids
  low.centr <- wt.mean(rsq.clean$rsq, g$g.low)
  high.centr <- wt.mean(rsq.clean$rsq, g$g.high)
  
  # Inter-centroidal separation
  inter <- abs(high.centr - low.centr)
  # Intra-cluster variance
  var.low <- wt.var(rsq.clean$rsq, g$g.low)
  var.high <- wt.var(rsq.clean$rsq, g$g.high)
  intra <- (var.low + var.high) / 2 # Mean intra variance
  
  res <- list(cl.met=(inter / intra), low.centr=low.centr, high.centr=high.centr,
              var.low=var.low, var.high=var.high)
  
}

# Plot GoM versus R-sq, by country