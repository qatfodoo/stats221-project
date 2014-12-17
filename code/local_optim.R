source("dat_format.R")
source("MLE.R")
source("SIR_fit.R")

## Local try of MLE

# 2 level per cat
feat <- ComputeCatDf(feature.names, rep(2, length(feature.names)), age=T, cat.age=2)
obs <- feat[2:length(feat)] # Remove sdr_id col

# theta0 uniform
theta0 <- PriorTheta(c(rep(2, length(feature.names)), 4))

# G0 uniform
G0 <- matrix(data=0.5, nrow=dim(obs)[1], 2)

# Compute estimates
t1.sim <- as.numeric(Sys.time())
gom.est <- gomMLE(obs, G0, theta0, epsilon=1e-2)
t2.sim <- as.numeric(Sys.time())
dt.sim <- (t2.sim - t1.sim) / 60 # dt in min
print(paste("MLE elapsed time (min)", dt.sim, sep=": "))

print(gom.est) # Print param estimates
save(file='./out/bicat_gom_est.RData', gom.est)

# Analyze results
G.hat <- data.frame(sdr.id=feat$sdr.id, g.hat=gom.est$G.hat)
G.hat <- G.hat[order(G.hat$sdr.id), ] # Make sure order is right
colnames(G.hat) <- c("sdr.id", "g.low", "g.high")
clust <- data.frame(sdr.id=G.hat$sdr.id, high=as.numeric(G.hat$g.high >= 0.5))
clust <- clust[order(clust$sdr.id), ] # Make sure order is right

hclust.res <- EvalHardClust(clust) # Metrics from clustering
print(paste("Hard clust metric: ", hclust.res$cl.met, sep=""))
fclust.res <- EvalFuzzClust(G.hat)
print(paste("Fuzzy clust metric: ", fclust.res$cl.met, sep=""))
