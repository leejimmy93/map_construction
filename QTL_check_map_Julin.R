library(qtl)
library(ggplot2)
library(reshape2)

###################### Read in cross data & data initial check ########################################
setwd("~/Desktop/Brassica_project/QTL_mapping/map_costruction/")
cross <- read.cross(file="~/Desktop/Brassica_project/QTL_mapping/map_costruction/151212_rqrl_19_chrom_original.csv",
                   genotypes=c("N","H","D","X","Y"))
         
## Take a look at of the cross data
summary(cross) # some markers are on identical positions    
names(cross$geno) 
pull.map(cross, chr = c("A02","A06","C01","C09"))

# Check FnP map, these four chr have more than one linkage group, C2_1 already was renamed.
# go to original .csv file, change the second groups' names to x_1 

## import the modified cross data
cross_m <- read.cross(file = "151212_rqrl_19_chrom.csv", 
                      genotypes=c("N","H","D","X","Y"))

summary(cross_m) 

################### Genotype data checking ##################################################################

# 1) segregation distortion, Null hypothesis: no significant difference 
test <- geno.table(cross_m)
test[test$P.value < 0.01,] # which critical value should be used if not 0.05 or 0.01??? 

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
nm <- est.map(cross_m, error.prob = 0.001)
plot(nm, alternate.chrid = T) 
pull.map(cross_m, chr = c("A03","A05")) # longer length was plot than what was indicated in the data? 

###########################################################################
## Focus on A03 for a bit to visuzlize problems. 

plot.rf(cross_m,chr="A03", col.scheme = "redblue")

A03 <- cross_m$geno$A03$data # assign chr A03 geno data into A03

dimnames(A03) <- list(F2=1:nrow(A03),marker=markernames(cross_m,"A03")) # create a list 
# with "F2" and "marker" info inside, and assign these as the dimension names of object A03

# distance calculation 
A03.dist <- dist(A03) # compute the distance between the rows "individuals here"
A03.clust <- hclust(A03.dist) # hierachical cluster analysis based on distance between individuals

A03.melt <- melt(A03) # make A03 to long form 138 * 39
A03.melt$F2 <- factor(A03.melt$F2,levels=A03.clust$order) # encode distance order between individuals as F2 orders

# using ggplot to make plot for A03 markers
pl <- ggplot(A03.melt,aes(x=F2,y=marker,fill=as.factor(value)))
pl <- pl + geom_tile()
pl <- pl + ggtitle("F2 genotypes")
pl + scale_fill_manual(values=c("magenta","grey50","green"),labels=c("N","H","D"))

## flag markers that are flanked by recombination on both sides.
## Perhaps distinguish between being flanked by similar or different genotypes.

# calculate recombination
A03.recomb <- sapply(2:(ncol(A03)-1),function(i) { # apply the function to all markers except for the 1st and last one
  tmp <- A03[,(i-1):(i+1)] # assign three marker data (i-1), i, (i+1) to tmp
  apply(tmp,1,function(F2.gt) { # apply function F2.gt to tmp by row
    F2.gt <- as.factor(F2.gt) # define F2.gt as factor then assign to F2.gt 
    if (any(is.na(F2.gt))) return(NA) # if "NA", return "NA"
    if (nlevels(F2.gt) == 1) return("no_recombination") 
    if ( (nlevels(F2.gt) == 2) & (F2.gt[1] == F2.gt[3])) return("double_recomb_same") # double recombination same
    if ( (nlevels(F2.gt) == 2) & (F2.gt[1] != F2.gt[3])) return("single_recomb") # single recombination
    if (nlevels(F2.gt) == 3) return("double_recomb_different") # double recombination same
  })
})
#### Check tomorrow, how to define double recomination??????? counting crossover, calc.errorlod #########################
A03.recomb
class(A03.recomb)
dim(A03.recomb)
colnames(A03.recomb)
rownames(A03.recomb)

### try calc.errorlod, problem
map <- est.map(cross_m, error.prob = 0.01)
cross_mm <- replace.map(cross_m, map)
corss_mm <- calc.errorlod(cross_mm) # problem
top <- top.errorlod(cross_mm, cutoff = 5) 
top
plot.geno(cross_mm, A03, top$id[top$chr==A03], cutoff = 5)

############### I added the above to check A03.recomb ################################################################
# transform to a dataframe with marker and recomination frequency 
colnames(A03.recomb) <- colnames(A03)[2:(ncol(A03)-1)] # assign colnames to A03.recomb (except the 1st and last one)

A03.recomb.count <- as.data.frame(matrix(nrow=ncol(A03.recomb),ncol=length(na.omit((unique(as.character(A03.recomb))))), 
                                         dimnames=list(marker=colnames(A03.recomb),result=na.omit(unique(as.character(A03.recomb))))))
# the above make a dataframe with marker # as row #, recombination type as col #, also assign dimnames: marker names as marker
# recomination type as result

A03.recomb.count$marker <- rownames((A03.recomb.count)) # add another colomn "marker"

for(marker in colnames(A03.recomb)) { # for every marker in A03.recomb
  tmp <- as.data.frame(table(A03.recomb[,marker])) # pull out recomb type, make a frequency table and define as a data.frame
  A03.recomb.count[marker,as.character(tmp$Var1)] <- tmp$Freq # assign frequency for each reocom to A03.recom.count
}

A03.recomb.count$percent_double_same <- apply(A03.recomb.count[,1:4],1, # apply to each row of A03.recomb.count
                                              function(x) x["double_recomb_same"]/sum(x,na.rm=TRUE)*100) 
# x calculate the percentage of double recomb same out of all 

A03.recomb.count$percent_double_all <- apply(A03.recomb.count[,1:4],1,
                                             function(x) sum(x[c("double_recomb_same","double_recomb_different")],na.rm=TRUE) / 
                                               sum(x,na.rm=TRUE)*100) 
# calcuate the percentage of double recomb out of all, regardless of same of diff

head(A03.recomb.count)

with(A03.recomb.count,plot(percent_double_all,percent_double_same)) 

write.csv(A03.recomb.count[order(A03.recomb.count$percent_double_same,decreasing = TRUE),
                           c("marker","percent_double_same")],file="B.napus_A03marker.csv")


## Now lets do this for the entire genome

#Remove the bad genotypes from a matrix in a rqtl cross object: assign "NA" to double recombnants

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

newmap.drop.marker <- est.map(cross.drop.marker,verbose=T,error.prob=.01)
plot.map(cross_m,cross.drop.marker)
################# drop marker ##########################################
# explore cross with some markers dropped
plot.rf(cross.drop.marker)
summary(cross.drop.marker)
summary(cross) #??????????? Check this tomorrow!!! 
########################################################################
plot.map(cross.drop.marker,newmap.drop.marker)

cross.drop.marker2 <- orderMarkers(cross.drop.marker,window=7,use.ripple = TRUE,verbose=T)
cross.drop.marker2 <- est.rf(cross.drop.marker2)
plot.rf(cross.drop.marker2)

#not working well
plot.map(newmap.drop.marker,cross.drop.marker2)

plot.rf(newmap.drop.marker)

#does not help
rip1 <- ripple(cross.drop.marker2,"A03",window=7,n.cluster=2)
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



