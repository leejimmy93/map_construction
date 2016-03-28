library(qtl)
library(ggplot2)
library(reshape2)

###################### Read in cross data & data initial check ########################################
setwd("~/Desktop/Brassica_project/QTL_mapping/map_construction/")
cross <- read.cross(file="151212_rqrl_19_chrom_original.csv",
                   genotypes=c("N","H","D","X","Y"))
         
## Take a look at of the cross data
summary(cross) # some markers are on identical positions    
names(cross$geno) 
pull.map(cross, chr = c("A02","A06","C01","C09"))

# Check FnP map, these four chr have more than one linkage group, C2_1 already was renamed.
# go to original .csv file, change the second groups' names to x_1 

## import the modified cross data
cross_m <- read.cross(file = "~/Desktop/Brassica_project/QTL_mapping/map_construction/151212_rqrl_19_chrom.csv", 
                      genotypes=c("N","H","D","X","Y"))

summary(cross_m) 

################### Genotype data checking ##################################################################

# 1) segregation distortion, Null hypothesis: no significant difference 
test <- geno.table(cross_m)
test[test$P.value < 0.01,] # which critical value should be used if not 0.05 or 0.01??? 

#### may remove later depends on the figure result
#          chr missing AA AB BB not.BB not.AA      P.value
# FPBPN0006 A03      11 48 52 27      0      0 3.867612e-03
# FPBPN0205 A03       6 37 44 51      0      0 1.480172e-04
# FPBRN0434 A03       1 43 75 19      0      0 8.057372e-03
# FPBPN0364 A05       6 18 68 46      0      0 2.478752e-03
# FPBPN0981 A10       9 20 64 45      0      0 7.837486e-03
# FPBPN1513 C02       6 43 70 19      0      0 9.991100e-03
# FPBPN0103 C09       4 54 39 41      0      0 2.346263e-06

# 2) compare individual's genotype similarity. Assumption: no two individuals with very similary genotypes
cg <- comparegeno(cross_m) 
dim(cg) 
hist(cg, breaks = 200,
     xlab="Proportion of identical genotypes")
rug(cg)
which(cg >0.9, arr.ind = TRUE) # no individuals with high genotype similarities

# 3) Check marker order
# 3.1) Check pairwise recombination fraction 
cross_m <- est.rf(cross_m) 
checkAlleles(cross_m) 

# No apparent problem, make a plot to see the data
plot.rf(cross_m, col.scheme = "redblue", alternate.chrid = T)
# looks like they are on the right chromosomes but wrong order within each chromosome

# estimate the genetic map
nm <- est.map(cross_m)
plot(nm, alternate.chrid = T) 
pull.map(cross_m, chr = c("A03","A05")) # longer length was plot than what was indicated in the data? 
pull.map(cross_m)

################ Remove the bad genotypes from a matrix in a rqtl cross object: assign "NA" to double recombinants
################ Written by Julin ###################################################################

cross.drop.marker <- cross_m

for (chr in names(cross.drop.marker$geno)) { # for each chromosome in cross genotype data
  my.chr <- get(chr,cross.drop.marker$geno) # return the genotype data, including data & map
  print(paste(chr,"NA before",sum(is.na(my.chr$data)))) 
  if(ncol(my.chr$data) > 3) { 
    my.chr$data[,2:(ncol(my.chr$data)-1)] <- sapply(2:(ncol(my.chr$data)-1),function(i) {
      apply(my.chr$data[,(i-1):(i+1)],1,function(gt) {
        if (any(is.na(gt))) return(gt[2]) #technically should be looking at the next genotyped marker.
        if ( (length(unique(gt)) == 2) & (gt[1] == gt[3])) return(NA)
        if ( length(unique(gt))  == 3) return(NA)
        return(gt[2])
      })
    })
  }
  cross.drop.marker$geno <- within(cross.drop.marker$geno,assign(chr,my.chr))
  print(paste(chr,"NA after",sum(is.na(get(chr,cross.drop.marker$geno)$data))))
}

write.csv(pull.geno(cross.drop.marker),file="F2_415_Dropped_geno.csv")
## Try to improve the map

# calculate new map and replace map
newmap.drop.marker <- est.map(cross.drop.marker,verbose=T,error.prob=.01) # why going through 
# the last section of code doesn't assign new map 
cross.drop.marker <- replace.map(cross.drop.marker, newmap.drop.marker)

############################ draw graph for new cross data ###############################

plot.map(cross_m, alternate.chrid = T) # the old genetic map
plot.map(cross_m,newmap.drop.marker, alternate.chrid = T) # genetic map comparison

geno.image(cross.drop.marker, alternate.chrid = T) # grid with color pixels for different genotypes
par(mfrow=c(1,2)) # create three rows and two colomns of plots with mfrow(1,2)
geno.image(cross_m, chr = "A03") 
geno.image(cross.drop.marker, chr = "A03") # compare grid w/ & w/o bad genotype info

# compare rf graph compare w/ & w/o bad genotype data
plot.rf(cross_m, col.scheme = "redblue", alternate.chrid = T) 
plot.rf(cross.drop.marker, col.scheme = "redblue", alternate.chrid = T) # why no difference??? 
par(mfrow=c(1,1)) # reset to 1:1

##################### drop markers with double crossover ##################
##### dosen't seem work 

geno <- pull.geno(cross_m)
dim(geno)
all.marker <- colnames(geno)

# calculate the number of double crossover for each marker
double.crossover.count <- sapply(all.marker, function(marker){
  sum(is.na(geno_new[,marker]))-sum(is.na(geno[,marker]))
}
)

double.crossover.count <- as.data.frame(double.crossover.count)
double.crossover.count
plot(density(double.crossover.count$double.crossover.count))
marker.drop <- rownames(subset(double.crossover.count, double.crossover.count>10))
marker.drop

# drop marker with double crossover more than 
cross.drop.marker.drop <- drop.markers(cross.drop.marker, marker.drop)
summary(cross.drop.marker.drop)
set.seed(5678)
cross.drop.marker.drop <- orderMarkers(cross.drop.marker.drop,  
                          window=4, use.ripple = T, maxit=4000, 
                          error.prob=0.0001)

map.new <- est.map(cross.drop.marker.drop, error.prob = 0.001)
plot.map(map.new)
plot.rf(cross.drop.marker.drop, col.scheme = "redblue", alternate.chrid = T)

# pairwise comparison of rf and LOD score before & after dropping marker
par(mfrow=c(2,1))
plot.rf(cross.drop.marker, chr = "A03", col.scheme = "redblue", alternate.chrid = T)
plot.rf(cross.drop.marker.drop, chr = "A03", col.scheme = "redblue", alternate.chrid = T)

# try ripple again
set.seed(1234)
cross.drop.marker.drop <- orderMarkers(cross.drop.marker.drop,  
                                       window=4, use.ripple = T, maxit=4000, 
                                       error.prob=0.0001)
map.new2 <- est.map(cross.drop.marker.drop, error.prob = 0.001)
summary.map(cross.drop.marker.drop)
plot.map(map.new2)
plot.rf(cross.drop.marker.drop, col.scheme = "redblue", alternate.chrid = T)

plot.rf(cross.drop.marker, chr = "A03", col.scheme = "redblue", alternate.chrid = T)
plot.rf(cross.drop.marker.drop, chr = "A03", col.scheme = "redblue", alternate.chrid = T)
pull.map(cross.drop.marker.drop, chr = "A03")

############################### go through each chomosome separately #######################
rf.full <- pull.rf(cross.drop.marker)

# go through each chromosome separately 
####### A01
plot.rf(cross.drop.marker, chr = "A01", col.scheme = "redblue") # A01 is OK??? 
geno.image(cross.drop.marker, chr = "A01")

A01 <- markernames(cross.drop.marker, chr = "A01")
plot(rf.full, A01[2], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)

####### A02
plot.map(cross.drop.marker, chr = "A02")
plot.rf(cross.drop.marker, chr = "A02", col.scheme = "redblue") 
set.seed(100)
rip1 <- ripple(cross.drop.marker,chr = "A02",window=8,n.cluster=2)
cross.drop.marker2 <- switch.order(cross.drop.marker, "A02", rip1[2,])
plot.rf(cross.drop.marker2, chr = "A02", col.scheme = "redblue") 
set.seed(201)
rip1 <- ripple(cross.drop.marker2,chr = "A02",window=8,n.cluster=2)
summary(rip1)
cross.drop.marker3 <- switch.order(cross.drop.marker2, "A02", rip1[2,])
plot.rf(cross.drop.marker3, chr = "A02", col.scheme = "redblue") 
geno.image(cross.drop.marker3, chr = "A02")

A02 <- markernames(cross.drop.marker, chr = "A02")
plot(rf.full, A02[1], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)

######## A03 
plot.map(cross.drop.marker3, chr = "A03")
plot.rf(cross.drop.marker3, chr = "A03", col.scheme = "redblue")

## looks like chr assignment is OK, order within chromosome doesn't look right
rf.A03 <- pull.rf(cross.drop.marker3, chr = "A03")
A03 <- markernames(cross.drop.marker3, chr='A03')
A03

par(mfrow=c(1,1))
plot(rf.A03, A03[28], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)
find.markerpos(cross.drop.marker3,"FPBPN0164")

# A04, no improvement, OK????? 
plot.map(cross.drop.marker3, chr = "A04")
plot.rf(cross.drop.marker3, chr="A04", col.scheme = "redblue")
rip1 <- ripple(cross.drop.marker3, chr = "A04", window = 7)

# A05 
plot.map(cross.drop.marker3, chr = "A05")
plot.rf(cross.drop.marker3, chr = "A05", col.scheme = "redblue")


cross.drop.marker2 <- est.rf(cross.drop.marker2)
plot.rf(cross.drop.marker2, col.scheme = "redblue", alternate.chrid = T)

#not working well
plot.map(newmap.drop.marker,cross.drop.marker2)

#does not help
rip1 <- ripple(cross.drop.marker2,chr = "A03",window=7,n.cluster=2)
summary(rip1)

#try to work through this with ripple

cross.drop.marker <- replace.map(cross.drop.marker,newmap.drop.marker)

A01.rip1 <- ripple(cross = cross.drop.marker, chr="A01",  window = 7, n.cluster = 2, verbose = TRUE)
summary(A01.rip1) #no change
plot.rf(cross.drop.marker,chr="A01")

A02.rip <- ripple(cross = cross.drop.marker, chr="A02",  window = 7, n.cluster = 2, verbose = TRUE)
summary(A02.rip)
plot.rf(cross.drop.marker,chr="A02")
cross.drop.marker.rip <- switch.order(cross.drop.marker,"A02",A02.rip[2,])
plot.rf(cross.drop.marker.rip,"A02")
A02.rip <- ripple(cross = cross.drop.marker.rip, 
                   chr="A02",  window = 7, n.cluster = 2, verbose = TRUE)
summary(A02.rip) # no further change
plot.map(cross.drop.marker,cross.drop.marker.rip,chr="A02")

A03.rip <- ripple(cross = cross.drop.marker, chr="A03",  window = 7, n.cluster = 2, verbose = TRUE)
summary(A03.rip)
plot.rf(cross.drop.marker,chr="A03")
cross.drop.marker.rip <- switch.order(cross.drop.marker.rip,"A03",A03.rip[2,])
plot.rf(cross.drop.marker.rip,"A03")
A03.rip <- ripple(cross = cross.drop.marker.rip, 
                  chr="A03",  window = 7, n.cluster = 2, verbose = TRUE)
summary(A03.rip) # no further change
plot.map(cross.drop.marker,cross.drop.marker.rip,chr="A03")

A04.rip <- ripple(cross = cross.drop.marker, chr="A04",  window = 7, n.cluster = 2, verbose = TRUE)
summary(A04.rip)
plot.rf(cross.drop.marker,chr="A04")
# no change

A05.rip <- ripple(cross = cross.drop.marker, chr="A05",  window = 7, n.cluster = 2, verbose = TRUE)
summary(A05.rip)
plot.rf(cross.drop.marker,chr="A05")
plot.map(cross.drop.marker,chr="A05")
#no change

A06.rip <- ripple(cross = cross.drop.marker, chr="A06",  window = 7, n.cluster = 2, verbose = TRUE)
summary(A06.rip)
plot.rf(cross.drop.marker,chr="A06")
cross.drop.marker.test <- switch.order(cross.drop.marker.rip,"A06",c(2:4,6,8,1,5,7,9,10,16,18,17,15,19:24,14:11,25:30))
plot.rf(cross.drop.marker.test,"A06")
A06.rip <- ripple(cross = cross.drop.marker.test, 
                  chr="A06",  window = 7, n.cluster = 2, verbose = TRUE)
summary(A06.rip) 
cross.drop.marker.test <- switch.order(cross.drop.marker.test,"A06",A06.rip[2,])
plot.rf(cross.drop.marker.test,"A06")
cross.drop.marker.test <- switch.order(cross.drop.marker.test,"A06",c(1:8,10:23,9,24:30))
plot.rf(cross.drop.marker.test,"A06")
A06.rip <- ripple(cross = cross.drop.marker.test, 
                  chr="A06",  window = 7, n.cluster = 2, verbose = TRUE)
summary(A06.rip) 
cross.drop.marker.test <- switch.order(cross.drop.marker.test,"A06",A06.rip[2,])
plot.rf(cross.drop.marker.test,"A06")
A06.rip <- ripple(cross = cross.drop.marker.test, 
                  chr="A06",  window = 7, n.cluster = 2, verbose = TRUE)
cross.drop.marker.test <- switch.order(cross.drop.marker.test,"A06",A06.rip[2,])
plot.rf(cross.drop.marker.test,"A06")
A06.rip <- ripple(cross = cross.drop.marker.test, 
                  chr="A06",  window = 7, n.cluster = 2, verbose = TRUE)
summary(A06.rip) 
plot.map(cross.drop.marker,cross.drop.marker.test,chr="A06")
cross.drop.marker.rip <- cross.drop.marker.test

A07.rip <- ripple(cross = cross.drop.marker, chr="A07",  window = 7, n.cluster = 2, verbose = TRUE)
summary(A07.rip)
plot.rf(cross.drop.marker,chr="A07")
cross.drop.marker.rip <- switch.order(cross.drop.marker.rip,"A07",A07.rip[2,])
plot.rf(cross.drop.marker.rip,"A07")
A07.rip <- ripple(cross = cross.drop.marker.rip, 
                  chr="A07",  window = 7, n.cluster = 2, verbose = TRUE)
summary(A07.rip) # no further change
plot.map(cross.drop.marker,cross.drop.marker.rip,chr="A07")



