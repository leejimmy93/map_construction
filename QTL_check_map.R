# install and library all necessary package
install.packages("qtl")

library(qtl)
library(ggplot2)
library(reshape2)


## Read in cross data
setwd("~/Desktop/Brassica_project/QTL_mapping/SNP/")
cross <- read.cross(format = c("csv"), file="temp.csv",
                  genotypes=c("N","H","D","X","Y"), estimate.map = FALSE)

# below is just for practice... 
# plot.map(cross)
# class(cross)
# names(cross)

cross <- jittermap(cross)

cross <- est.rf(cross)

plot.rf(cross)

## Focus on A03 for a bit to visuzlize problems.

plot.rf(cross,chr="A03")

A03 <- cross$geno$A03$data

dimnames(A03) <- list(F2=1:nrow(A03),marker=markernames(cross,"A03"))

A03.dist <- dist(A03)
A03.clust <- hclust(A03.dist)

A03.melt <- melt(A03)
head(A03.melt)
A03.melt$F2 <- factor(A03.melt$F2,levels=A03.clust$order)

pl <- ggplot(A03.melt,aes(x=F2,y=marker,fill=as.factor(value)))
pl <- pl + geom_tile()
pl <- pl + ggtitle("F2 genotypes")
pl + scale_fill_manual(values=c("magenta","grey50","green"),labels=c("N","H","D"))

## flag markers that are flanked by recombination on both sides.
## Perhaps distinguish between being flanked by similar or different genotypes.

A03.recomb <- sapply(2:(ncol(A03)-1),function(i) {
  tmp <- A03[,(i-1):(i+1)]
  apply(tmp,1,function(F2.gt) {
    F2.gt <- as.factor(F2.gt)
    if (any(is.na(F2.gt))) return(NA)
    if (nlevels(F2.gt) == 1) return("no_recombination")
    if ( (nlevels(F2.gt) == 2) & (F2.gt[1] == F2.gt[3])) return("double_recomb_same")
    if ( (nlevels(F2.gt) == 2) & (F2.gt[1] != F2.gt[3])) return("single_recomb")
    if (nlevels(F2.gt) == 3) return("double_recomb_different")
  })
})

colnames(A03.recomb) <- colnames(A03)[2:(ncol(A03)-1)]

A03.recomb.count <- as.data.frame(matrix(nrow=ncol(A03.recomb),ncol=length(na.omit((unique(as.character(A03.recomb))))),
                                         dimnames=list(marker=colnames(A03.recomb),result=na.omit(unique(as.character(A03.recomb))))))

A03.recomb.count$marker <- rownames((A03.recomb.count))

for(marker in colnames(A03.recomb)) {
  tmp <- as.data.frame(table(A03.recomb[,marker]))
  A03.recomb.count[marker,as.character(tmp$Var1)] <- tmp$Freq
}

A03.recomb.count$percent_double_same <- apply(A03.recomb.count[,1:4],1,
                                              function(x) x["double_recomb_same"]/sum(x,na.rm=TRUE)*100)

A03.recomb.count$percent_double_all <- apply(A03.recomb.count[,1:4],1,
                                             function(x) sum(x[c("double_recomb_same","double_recomb_different")],na.rm=TRUE) / 
                                               sum(x,na.rm=TRUE)*100)

head(A03.recomb.count)

with(A03.recomb.count,plot(percent_double_all,percent_double_same))

write.csv(A03.recomb.count[order(A03.recomb.count$percent_double_same,decreasing = TRUE),
                           c("marker","percent_double_same")],file="~/desktop/B.napus_A03marker.csv")



## Now lets do this for the entire genome

#Remove the bad genotypes from a matrix in a rqtl cross object

cross.drop.marker <- cross

for (chr in names(cross.drop.marker$geno)) {
  my.chr <- get(chr,cross.drop.marker$geno)
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

write.csv(pull.geno(cross.drop.marker),file="~/Desktop/F2_415_Dropped_geno.csv")

## Try to improve the map

newmap.drop.marker <- est.map(cross.drop.marker,verbose=T,error.prob=.01)

plot.map(cross,cross.drop.marker)

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
summary(A06.rip) 
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



