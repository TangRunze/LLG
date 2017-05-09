
setwd("/Users/Runze/Documents/GitHub/LLG/Code/R")
# setwd("E:/GitHub/LLG/Code/R")
# setwd("/cis/home/rtang/LLG/Code/R")

isSVD = 0
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

A_bar = add(A_all)/M

A_bar_diag_aug = diag_aug(A_bar)

# ZG
nElbow = 3
evalVec = ase(A_bar, ceiling(n*3/5), isSVD)[[1]]
dZG = getElbows(evalVec, n=nElbow, plot=F)[[nElbow]]

A.ase = ase(A_bar_diag_aug, dZG, isSVD)
d = dZG

print(dZG)

xHat <- A.ase[[3]] %*% diag(sqrt(A.ase[[1]]))

for (j in 1:8) {
  s = ""
  for (i in 1:n) {
    s = paste0(s,",",xHat[i,j])
  }
  s = substr(s,2,nchar(s))
  write(s,file=paste0("../../Result/eigenvector_dim", j, ".csv"))
  write(s,file=paste0("../Python/data/eigenvector_dim", j, ".csv"))
}

nameVec <- read.table("../../Data/desikan.txt")
fileName <- "../../Result/eigenvector_vertex_name.csv"
write("Vertices selected based on embeddings", file=fileName)
for (j in 1:8) {
  write(paste0("Dimension ", j, ":"), file=fileName, append = TRUE)
  tmp <- order(xHat[,j])
  s <- paste0("Smallest 5 vertices")
  for (i in 1:5) {
    s = paste0(s, ", ", nameVec[tmp[i], 1], " ", xHat[tmp[i], j])
  }
  write(s, file=fileName, append = TRUE)
  s <- paste0("Largest 5 vertices")
  for (i in 1:5) {
    s = paste0(s, ", ", nameVec[tmp[n - i + 1], 1], " ", xHat[tmp[n - i + 1], j])
  }
  write(s,file=fileName, append = TRUE)
}
