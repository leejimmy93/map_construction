# install everything
install.packages("onemap")
install.packages("tkrplot", type="source")
library(tcltk)
library(tkrplot)
library(onemap)

############################## formating data for my analysis ############
setwd("~/Desktop/Brassica_project/QTL_mapping/one_map/")
cross_geno_data <- t(read.csv("F2_415_dropp_geno_Ruijuan.csv", header = T))
cross_geno_data

# formating data for onemap 
cross_geno_data[is.na(cross_geno_data)] <- "-"
cross_geno_data <- gsub("1", "A", cross_geno_data)
cross_geno_data <- gsub("2", "H", cross_geno_data)
cross_geno_data <- gsub("3", "B", cross_geno_data)
cross_geno_data

cross_geno_data <- as.data.frame(cross_geno_data)
ncol(cross_geno_data)
ncol(cross_geno_data)
nrow(cross_geno_data)
write.table(cross_geno_data, file = "cross_geno_data.txt")
# change file format in linux 
# cat cross_geno_data.txt | sed 's/"//g' | awk '{first = $1; $1 = ""; print $0}' | sed 's/ //g' > tmp

cross_pheno_data <- t(read.csv("F2_415_pheno.csv", header = T))
cross_pheno_data
cross_pheno_data <- as.data.frame(cross_pheno_data)
cross_pheno_data
# change file format in linux 
# cat cross_pheno_data.txt | sed 's/\"\.\"/\"\-\"/g'  | sed 's/"//g' | awk '{print "*" $0}' > pheno.txt 

########### real analysis start from here ################################
data <- read.mapmaker(file = "input_for_onemap.txt")
data

# estimate two-point rf
twopts.f2 <- rf.2pts(data, LOD = 6, max.rf = 0.25)
?rf.2pts

# assign marker to linkage groups
mark.all.f2 <- make.seq(twopts.f2, "all")
?make.seq
LGs.f2 <- group(mark.all.f2, LOD = 6, max.rf = 0.25)
LGs.f2

#################### LG1 ##############################
LG1.f2 <- make.seq(LGs.f2, 1)
LG1.f2

# use touchdown
LG1.f2.ord <- order.seq(input.seq = LG1.f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 6,
                        draw.try = T, wait = 1,
                        touchdown = T)

LG1.f2.ord

# safe mode order
LG1.f2.safe <- make.seq(LG1.f2.ord, "safe")
LG1.f2.safe
# all order
LG1.f2.final <- make.seq(LG1.f2.ord, "force")
LG1.f2.final

# ripple
ripple.seq(LG1.f2.final, ws = 5, LOD = 6) 

################# LG2 ##################################
LG2.f2 <- make.seq(LGs.f2, 2)
LG2.f2

# use touchdown
LG2.f2.ord <- order.seq(input.seq = LG2.f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 6,
                        draw.try = T, wait = 1,
                        touchdown = T)

LG2.f2.ord

# safe mode order
LG2.f2.safe <- make.seq(LG2.f2.ord, "safe")
LG2.f2.safe
# all order
LG2.f2.final <- make.seq(LG2.f2.ord, "force")
LG2.f2.final

# ripple
ripple.seq(LG2.f2.final, ws = 5, LOD = 6) 

######################## LG3 ########################

LG3.f2 <- make.seq(LGs.f2, 3)
LG3.f2

# use touchdown
LG3.f2.ord <- order.seq(input.seq = LG3.f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 6,
                        draw.try = T, wait = 1,
                        touchdown = T)

LG3.f2.ord

# safe mode order
LG3.f2.safe <- make.seq(LG3.f2.ord, "safe")
LG3.f2.safe
# all order
LG3.f2.final <- make.seq(LG3.f2.ord, "force")
LG3.f2.final

# ripple
ripple.seq(LG3.f2.final, ws = 3, LOD = 6) 

######################## LG4 ######################

LG4.f2 <- make.seq(LGs.f2, 4)
LG4.f2

# use touchdown
LG4.f2.ord <- order.seq(input.seq = LG4.f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 6,
                        draw.try = T, wait = 1,
                        touchdown = T)

LG4.f2.ord

# safe mode order
LG4.f2.safe <- make.seq(LG4.f2.ord, "safe")
LG4.f2.safe
# all order
LG4.f2.final <- make.seq(LG4.f2.ord, "force")
LG4.f2.final

# ripple
ripple.seq(LG4.f2.final, ws = 5, LOD = 6)

###################### LG5 ########################

LG5.f2 <- make.seq(LGs.f2, 5)
LG5.f2

# use touchdown
LG5.f2.ord <- order.seq(input.seq = LG5.f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 6,
                        draw.try = T, wait = 1,
                        touchdown = T)

LG5.f2.ord

# safe mode order
LG5.f2.safe <- make.seq(LG5.f2.ord, "safe")
LG5.f2.safe
# all order
LG5.f2.final <- make.seq(LG5.f2.ord, "force")
LG5.f2.final

################# LG6 ##############################
LG6.f2 <- make.seq(LGs.f2, 6)
LG6.f2

# use touchdown
LG6.f2.ord <- order.seq(input.seq = LG6.f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 6,
                        draw.try = T, wait = 1,
                        touchdown = T)

LG6.f2.ord

# safe mode order
LG6.f2.safe <- make.seq(LG6.f2.ord, "safe")
LG6.f2.safe
# all order
LG6.f2.final <- make.seq(LG6.f2.ord, "force")
LG6.f2.final

################# LG7 ##############################
LG7.f2 <- make.seq(LGs.f2, 7)
LG7.f2

# use touchdown
LG7.f2.ord <- order.seq(input.seq = LG7.f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 6,
                        draw.try = T, wait = 1,
                        touchdown = T)

LG7.f2.ord

# safe mode order
LG7.f2.safe <- make.seq(LG7.f2.ord, "safe")
LG7.f2.safe
# all order
LG7.f2.final <- make.seq(LG7.f2.ord, "force")
LG7.f2.final

################# LG8 ##############################
LG8.f2 <- make.seq(LGs.f2, 8)
LG8.f2

# use touchdown
LG8.f2.ord <- order.seq(input.seq = LG8.f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 6,
                        draw.try = T, wait = 1,
                        touchdown = T)

LG8.f2.ord

# safe mode order
LG8.f2.safe <- make.seq(LG8.f2.ord, "safe")
LG8.f2.safe
# all order
LG8.f2.final <- make.seq(LG8.f2.ord, "force")
LG8.f2.final

################# LG9 ##############################
LG9.f2 <- make.seq(LGs.f2, 9)
LG9.f2

# use touchdown
LG9.f2.ord <- order.seq(input.seq = LG9.f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 6,
                        draw.try = T, wait = 1,
                        touchdown = T)

LG9.f2.ord

# safe mode order
LG9.f2.safe <- make.seq(LG9.f2.ord, "safe")
LG9.f2.safe
# all order
LG9.f2.final <- make.seq(LG9.f2.ord, "force")
LG9.f2.final

################# LG10 ##############################
LG10.f2 <- make.seq(LGs.f2, 10)
LG10.f2

# use touchdown
LG10.f2.ord <- order.seq(input.seq = LG10.f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 6,
                        draw.try = T, wait = 1,
                        touchdown = T)

LG10.f2.ord

# safe mode order
LG10.f2.safe <- make.seq(LG10.f2.ord, "safe")
LG10.f2.safe
# all order
LG10.f2.final <- make.seq(LG10.f2.ord, "force")
LG10.f2.final

################# LG11/A09 ##############################
LG11.f2 <- make.seq(LGs.f2, 11)
LG11.f2

# use touchdown
LG11.f2.ord <- order.seq(input.seq = LG11.f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 6,
                        draw.try = T, wait = 1,
                        touchdown = T)

LG11.f2.ord

# safe mode order
LG11.f2.safe <- make.seq(LG11.f2.ord, "safe")
LG11.f2.safe
# all order
LG11.f2.final <- make.seq(LG11.f2.ord, "force")
LG11.f2.final

################# LG12/A10 ##############################
LG12.f2 <- make.seq(LGs.f2, 12)
LG12.f2

# use touchdown
LG12.f2.ord <- order.seq(input.seq = LG12.f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 6,
                        draw.try = T, wait = 1,
                        touchdown = T)

LG12.f2.ord

# safe mode order
LG12.f2.safe <- make.seq(LG12.f2.ord, "safe")
LG12.f2.safe
# all order
LG12.f2.final <- make.seq(LG12.f2.ord, "force")
LG12.f2.final

################# LG13/C01 ##############################
LG13.f2 <- make.seq(LGs.f2, 13)
LG13.f2

# use touchdown
LG13.f2.ord <- order.seq(input.seq = LG13.f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 6,
                        draw.try = T, wait = 1,
                        touchdown = T)

LG13.f2.ord

# safe mode order
LG13.f2.safe <- make.seq(LG13.f2.ord, "safe")
LG13.f2.safe
# all order
LG13.f2.final <- make.seq(LG13.f2.ord, "force")
LG13.f2.final

################# LG14/C01-1 ##############################
LG14.f2 <- make.seq(LGs.f2, 14)
LG14.f2

# use touchdown
LG14.f2.ord <- order.seq(input.seq = LG14.f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 6,
                        draw.try = T, wait = 1,
                        touchdown = T)

LG14.f2.ord

# safe mode order
LG14.f2.safe <- make.seq(LG14.f2.ord, "safe")
LG14.f2.safe
# all order
LG14.f2.final <- make.seq(LG14.f2.ord, "force")
LG14.f2.final

################# LG15/C02 ##############################
LG15.f2 <- make.seq(LGs.f2, 15)
LG15.f2

# use touchdown
LG15.f2.ord <- order.seq(input.seq = LG15.f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 6,
                        draw.try = T, wait = 1,
                        touchdown = T)

LG15.f2.ord

# safe mode order
LG15.f2.safe <- make.seq(LG15.f2.ord, "safe")
LG15.f2.safe
# all order
LG15.f2.final <- make.seq(LG15.f2.ord, "force")
LG15.f2.final

################# LG16/C02-2 ##############################
LG16.f2 <- make.seq(LGs.f2, 16)
LG16.f2

# use touchdown
LG16.f2.ord <- order.seq(input.seq = LG16.f2, n.init = 5,
                         subset.search = "twopt",
                         twopt.alg = "rcd", THRES = 6,
                         draw.try = T, wait = 1,
                         touchdown = T)

LG16.f2.ord

# safe mode order
LG16.f2.safe <- make.seq(LG16.f2.ord, "safe")
LG16.f2.safe
# all order
LG16.f2.final <- make.seq(LG16.f2.ord, "force")
LG16.f2.final

################# LG17/C02-3 ##############################
LG17.f2 <- make.seq(LGs.f2, 17)
LG17.f2

# use touchdown
LG17.f2.ord <- order.seq(input.seq = LG17.f2, n.init = 2,
                         subset.search = "twopt",
                         twopt.alg = "rcd", THRES = 6,
                         draw.try = T, wait = 1,
                         touchdown = T)

LG17.f2.ord

# safe mode order
LG17.f2.safe <- make.seq(LG17.f2.ord, "safe")
LG17.f2.safe
# all order
LG17.f2.final <- make.seq(LG17.f2.ord, "force")
LG17.f2.final

################# LG18/C003 ##############################
LG18.f2 <- make.seq(LGs.f2, 18)
LG18.f2

# use touchdown
LG18.f2.ord <- order.seq(input.seq = LG18.f2, n.init = 5,
                         subset.search = "twopt",
                         twopt.alg = "rcd", THRES = 6,
                         draw.try = T, wait = 1,
                         touchdown = T)

LG18.f2.ord

# safe mode order
LG18.f2.safe <- make.seq(LG18.f2.ord, "safe")
LG18.f2.safe
# all order
LG18.f2.final <- make.seq(LG18.f2.ord, "force")
LG18.f2.final

################# LG19/C04 ##############################
LG19.f2 <- make.seq(LGs.f2, 19)
LG19.f2

# use touchdown
LG19.f2.ord <- order.seq(input.seq = LG19.f2, n.init = 5,
                         subset.search = "twopt",
                         twopt.alg = "rcd", THRES = 6,
                         draw.try = T, wait = 1,
                         touchdown = T)

LG19.f2.ord

# safe mode order
LG19.f2.safe <- make.seq(LG19.f2.ord, "safe")
LG19.f2.safe
# all order
LG19.f2.final <- make.seq(LG19.f2.ord, "force")
LG19.f2.final

################# LG20/C05 ##############################
LG20.f2 <- make.seq(LGs.f2, 20)
LG20.f2

# use touchdown
LG20.f2.ord <- order.seq(input.seq = LG20.f2, n.init = 5,
                         subset.search = "twopt",
                         twopt.alg = "rcd", THRES = 6,
                         draw.try = T, wait = 1,
                         touchdown = T)

LG20.f2.ord

# safe mode order
LG20.f2.safe <- make.seq(LG20.f2.ord, "safe")
LG20.f2.safe
# all order
LG20.f2.final <- make.seq(LG20.f2.ord, "force")
LG20.f2.final

################# LG21/C06 ##############################
LG21.f2 <- make.seq(LGs.f2, 21)
LG21.f2

# use touchdown
LG21.f2.ord <- order.seq(input.seq = LG21.f2, n.init = 5,
                         subset.search = "twopt",
                         twopt.alg = "rcd", THRES = 6,
                         draw.try = T, wait = 1,
                         touchdown = T)

LG21.f2.ord

# safe mode order
LG21.f2.safe <- make.seq(LG21.f2.ord, "safe")
LG21.f2.safe
# all order
LG21.f2.final <- make.seq(LG21.f2.ord, "force")
LG21.f2.final

################# LG22/C07 ##############################
LG22.f2 <- make.seq(LGs.f2, 22)
LG22.f2

# use touchdown
LG22.f2.ord <- order.seq(input.seq = LG22.f2, n.init = 5,
                         subset.search = "twopt",
                         twopt.alg = "rcd", THRES = 6,
                         draw.try = T, wait = 1,
                         touchdown = T)

LG22.f2.ord

# safe mode order
LG22.f2.safe <- make.seq(LG22.f2.ord, "safe")
LG22.f2.safe
# all order
LG22.f2.final <- make.seq(LG22.f2.ord, "force")
LG22.f2.final

################# LG23/C08 ##############################
LG23.f2 <- make.seq(LGs.f2, 23)
LG23.f2

# use touchdown
LG23.f2.ord <- order.seq(input.seq = LG23.f2, n.init = 5,
                         subset.search = "twopt",
                         twopt.alg = "rcd", THRES = 6,
                         draw.try = T, wait = 1,
                         touchdown = T)

LG23.f2.ord

# safe mode order
LG23.f2.safe <- make.seq(LG23.f2.ord, "safe")
LG23.f2.safe
# all order
LG23.f2.final <- make.seq(LG23.f2.ord, "force")
LG23.f2.final

################# LG24/C09 ##############################
LG24.f2 <- make.seq(LGs.f2, 24)
LG24.f2

# use touchdown
LG24.f2.ord <- order.seq(input.seq = LG24.f2, n.init = 5,
                         subset.search = "twopt",
                         twopt.alg = "rcd", THRES = 6,
                         draw.try = T, wait = 1,
                         touchdown = T)

LG24.f2.ord

# safe mode order
LG24.f2.safe <- make.seq(LG24.f2.ord, "safe")
LG24.f2.safe
# all order
LG24.f2.final <- make.seq(LG24.f2.ord, "force")
LG24.f2.final

################# LG25/C09-1 ##############################
LG25.f2 <- make.seq(LGs.f2, 25)
LG25.f2

# use touchdown
LG25.f2.ord <- order.seq(input.seq = LG25.f2, n.init = 5,
                         subset.search = "twopt",
                         twopt.alg = "rcd", THRES = 6,
                         draw.try = T, wait = 1,
                         touchdown = T)

LG25.f2.ord

# safe mode order
LG25.f2.safe <- make.seq(LG25.f2.ord, "safe")
LG25.f2.safe
# all order
LG25.f2.final <- make.seq(LG25.f2.ord, "force")
LG25.f2.final

###################################################
########### drop marker with double crossover ###########
data <- read.mapmaker(file = "F2_415_dropp_geno_Ruijuan.txt")
data

# estimate two-point rf
twopts.f2 <- rf.2pts(data, LOD = 3, max.rf = 0.5)
?rf.2pts

# assign marker to linkage groups
mark.all.f2 <- make.seq(twopts.f2, "all")
?make.seq
LGs.f2 <- group(mark.all.f2, LOD = 3, max.rf = 0.5)
LGs.f2

#################### LG1 ##############################
LG1.f2 <- make.seq(LGs.f2, 1)
LG1.f2

# use touchdown
LG1.f2.ord <- order.seq(input.seq = LG1.f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 6,
                        draw.try = T, wait = 1,
                        touchdown = T)

LG1.f2.ord

# safe mode order
LG1.f2.safe <- make.seq(LG1.f2.ord, "safe")
LG1.f2.safe
# all order
LG1.f2.final <- make.seq(LG1.f2.ord, "force")
LG1.f2.final

# ripple
ripple.seq(LG1.f2.final, ws = 5, LOD = 3) 



