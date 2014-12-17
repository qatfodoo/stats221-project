
## Computing categorical data from covariates

data.cov <- read.table(file="../dat/merged_covariate_df.csv", header=TRUE, sep='\t')

# Colnames to be considered
c.name <- colnames(data.cov)
region.names <- c.name[1:7]
feature.names <- c.name[match('iwi', c.name):(match('age09', c.name) - 1)]

# Transform standard feature into categorical data based on quantiles
CatTransf <- function(array, ncat=2) {
  
  quant.val <- seq(0, 1, length.out=(ncat + 1))
  quant.val <- quant.val[2:(length(quant.val) - 1)]
  
  quant <- sapply(quant.val, function(x) { quantile(array, x) })
  
  res <- c()
  for (i in 1:length(array)) {
    quant.add <- c(quant, array[i])
    quant.add <- sort(quant.add)
    
    level <- match(array[i], quant.add) - 1
    
    res <- c(res, level)
  }
  
  return(res)
}

# Compute age groups based on "young" and "old" groups
AgeTransf <- function(young.ar, old.ar, ncat=3) {
  
  quant.val <- seq(0, 1, length.out=(ncat + 1))
  quant.val <- quant.val[2:(length(quant.val) - 1)]
  y.quant <- sapply(quant.val, function(x) { quantile(young.ar, x) })
  o.quant <- sapply(quant.val, function(x) { quantile(old.ar, x) })
  
  res <- c()
  for (i in 1:length(young.ar)) {
    yquant.add <- c(y.quant, young.ar[i])
    yquant.add <- sort(yquant.add, decreasing=T)
    oquant.add <- c(o.quant, old.ar[i])
    oquant.add <- sort(oquant.add)
    
    ylevel <- match(young.ar[i], yquant.add) - 1
    olevel <- match(old.ar[i], oquant.add) - 1
    level <- ncat * olevel + ylevel # old pop principal in lexigocraphic order
    
    res <- c(res, level)
  }
  
  return(res)
  
}

# Given array of categories and level, compute categorized observations
ComputeCatDf <- function(cat.arr, lev.arr, age=T, cat.age=2) {
  
  df <- data.frame(sdr.id=data.cov$sdr_id)
  if (length(cat.arr) != length(lev.arr)){
    stop("Length of both arrays should be the same")
  }
  
  for (i in 1:length(cat.arr)) {
    array <- data.cov[, cat.arr[i]]
    cat.lev <- CatTransf(array, ncat=lev.arr[i])
    
    df[, cat.arr[i]] <- cat.lev
  }
  
  if (age) {
    y.arr <- data.cov$age09 + data.cov$age1019
    o.arr <- data.cov$age5059 + data.cov$age6069 + data.cov$age7079 +
      data.cov$age8089 + data.cov$age90hi
    
    cat.age <- AgeTransf(y.arr, o.arr, ncat=cat.age)
    df[, 'pop.age'] <- cat.age
    
  }
  
  return(df)
}

# Given level number per feature, compte uniform prior param
PriorTheta <- function(lev.arr) {
  
  theta0 <- sapply(lev.arr, function(n) {
    th.df <- data.frame(level=0:(n - 1), high=rep(1 / n, n), low=rep(1 / n, n))
    return(as.vector(th.df))
  }, simplify=F)
  
  return(theta0)
}

## Test on data
# Check ages add to 100
age <- data.cov[, c('age09', 'age1019', 'age2029', 'age3039', 'age4049', 'age5059', 
                    'age6069', 'age7079', 'age8089', 'age90hi')]
sum(age[1, ]) # adds to 100

cat.floor <- CatTransf(data.cov$bad_floor, ncat=2)
print(cat.floor)

y.arr <- data.cov$age09 + data.cov$age1019
o.arr <- data.cov$age5059 + data.cov$age6069 + data.cov$age7079 +
  data.cov$age8089 + data.cov$age90hi
cat.age <- AgeTransf(y.arr, o.arr, ncat=2)
print(cat.age)

