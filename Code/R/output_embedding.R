
setwd("/Users/Runze/Documents/GitHub/LLG/Code/R")
# setwd("E:/GitHub/LLG/Code/R")
# setwd("/cis/home/rtang/LLG/Code/R")

isSVD = 0
m = 5
dataName = "desikan"
set.seed(12345)


source("function_collection.R")
tmpList = read_data(dataName, DA=F, newGraph=F)
A_all = tmpList[[1]]
n = tmpList[[2]]
M = tmpList[[3]]
rm(tmpList)

# dVec = 1:n
# nD = length(dVec)

A_sum = add(A_all)

source("getElbows.R")
source("USVT.R")

sampleVec = sample.int(M, m)
A_bar = add(A_all[sampleVec])/m

A_bar_diag_aug = diag_aug(A_bar)

# ZG
nElbow = 3
evalVec = ase(A_bar, ceiling(n*3/5), isSVD)[[1]]
dZG = getElbows(evalVec, n=nElbow, plot=F)[[nElbow]]

A.ase = ase(A_bar_diag_aug, dZG, isSVD)
d = dZG

xHat <- A.ase[[3]] %*% diag(sqrt(A.ase[[1]]))

for (j in 1:5) {
  s = ""
  for (i in 1:n) {
    s = paste0(s,",",xHat[i,j])
  }
  s = substr(s,2,nchar(s))
  write(s,file=paste0("../../Result/eigenvector_dim", j, ".csv"))
}



