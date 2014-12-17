source("ga.R")
library(ggplot2)

## Some data visualization, and compute schemas to be tried on Odyssey

data.cov <- read.table(file="../dat/merged_covariate_df.csv", header=TRUE, sep='\t')

for (i in 1:length(feature.names)) {
  feat <- feature.names[i]
  cov <- data.cov[feat]
  jpeg(file = paste("./out/dist_", feat, sep=""))
  ggplot(cov, aes(x=feat)) + 
    geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   binwidth=.5,
                   colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
  dev.off()
}