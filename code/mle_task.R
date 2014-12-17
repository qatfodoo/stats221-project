source("MLE.R")
source("dat_format.R")
source("ga.R")

## Slurm job to compute Ebola risk region MLE

# 10 spread profile schemas per machines
N.sch <- 10
if (Sys.getenv("SLURM_JOB_ID") != "") { # Divide computation per tasks
  
  job.id <- as.numeric(Sys.getenv("SLURM_JOB_ID"))
  print(paste("Job id", job.id, sep=": "))
  task.id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  print(paste("Task id", task.id, sep=": "))
  
  gc() # garbage collection
  load('../dat/schem_sugg.RData')
  
  clust.met <- c()
  
  for (i in (1 + (task.id - 1) * N.sch):(task.id * N.sch)) {
    chrom <- sugg[i, ]
    clust.out <- ClustRes(chrom)
    
    save(clust.out, file=paste("./out/clustres_", task.id, ".dat", sep=""))
    
    clust.met <- c(clust.met, clust.out$clust.res$cl.met) # Record clustering metrics
    
  }

  t1.sim <- as.numeric(Sys.time())
  
  t2.sim <- as.numeric(Sys.time())
  dt.sim <- (t2.sim - t1.sim) / 60 # dt in min
  print(paste(paste("MLE elapsed time (min), task", task.id, sep=" "), dt.sim, sep=": "))
  
}
