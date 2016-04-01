# install and library all necessary package

library(qtl)
library(ggplot2)
library(reshape2)


## Read in cross data
setwd("~/Desktop/Brassica_project/QTL_mapping/map_construction/one_map/qtl/")
cross <- read.cross(format = c("csv"), file="order_from_blast.csv",
                  genotypes=c("N","H","D","X","Y"), estimate.map = FALSE)

# check order 
cross.rf <- est.rf(cross)
checkAlleles(cross.rf)  

plot.rf(cross, col.scheme = "redblue", alternate.chrid = T)



