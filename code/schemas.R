source("ga.R")
library(ggplot2)

## Some data visualization, and compute schemas to be tried on Odyssey

data.cov <- read.table(file="../dat/merged_covariate_df.csv", header=TRUE, sep='\t')

i <- 14
feat <- feature.names[i]
print(feat)
cov <- data.cov[feat]
ggplot(cov, aes(x=small_house)) + # Change name after x = 
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=sqrt(var(cov)) / 10,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot

# Compute 100 schemas to be used on Odyssey
n.feat <- 21
max.lev <- 5
num <- 100
sugg <- matrix(0, ncol=n.feat * (max.lev + 1), nrow=num)
# Bicat schema
seq0 <- rep(0, max.lev + 1)
seq0[3] <- 1 # 2 levels
bicat <- rep(seq0, n.feat)
sugg[1, ] <- bicat
# Empirical schema
emp.lev <- c(3, 3, 3, 3, 4, 4, 4, 3, 4, 3, 3, 4, 5, 5, 4, 4, 3, 4, 4, 5, 3)
chrom.emp <- c()
for (i in 1:length(emp.lev)) {
  seq.emp <- rep(0, max.lev + 1)
  seq.emp[emp.lev[i] + 1] <- 1
  chrom.emp <- c(chrom.emp, seq.emp)
}
sugg[2, ] <- chrom.emp

for (i in 3:num) {
  row <- c()
  for (j in 1:n.feat) {
    seq0 <- rep(0, (max.lev + 1))
    lev <- sample(1:(max.lev + 1), 1)
    seq0[lev] <- 1
    row <- c(row, seq0)
  }
  sugg[i, ] <- row
}

save(sugg, file='../dat/schem_sugg.RData')
