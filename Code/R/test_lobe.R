
setwd("/Users/Runze/Documents/GitHub/LLG/Code/R")
# setwd("E:/GitHub/LLG/Code/R")
# setwd("/cis/home/rtang/LLG/Code/R")

isSVD = 0
dataName = "desikan"

# colorVec <- c("#af8dc3", "#f7f7f7", "#7fbf7b")
colorVec <- c("#d8b365", "#f5f5f5", "#5ab4ac")

source("function_collection.R")
tmpList = read_data(dataName, DA=F, newGraph=F)
A_all = tmpList[[1]]
n = tmpList[[2]]
M = tmpList[[3]]
rm(tmpList)

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

xHat0 <- A.ase[[3]] %*% diag(sqrt(A.ase[[1]]))


data <- read.csv("../../Data/desikan_lobe.csv", header = F)
clName <- data[(n + 1):dim(data)[1],][, 1]
vertexName <- data[1:n, 1]
cl <- data[1:n, 2]
cl[36:70] <- cl[36:70] + max(cl)

nv <- !(cl == 5 | cl == 10)
cl <- cl[nv]
cl[cl > 5] <- cl[cl > 5] - 1
xHat <- xHat0[nv, ]

require(mclust)
ari <- rep(0, 8)
ari_std <- rep(0, 8)
for (i in 1:8) {
  model <- Mclust(xHat[, 1:5], G = i)
  ari[i] <- adjustedRandIndex(cl, model$classification)
  
  xHat1 <- xHat/sqrt(rowSums(xHat^2))
  model1 <- Mclust(xHat1[, 1:5], G = i)
  ari_std[i] <- adjustedRandIndex(cl, model1$classification)
}

df <- data.frame(ari=c(ari, ari_std), d=rep(1:8, time=2), which=rep(1:2, each=8))
gg <- ggplot(df, aes(d, ari, linetype=factor(which))) + 
  geom_point(size=1.5) +
  geom_line() +
  scale_linetype_manual(name="",values=c(1,2),
                        labels=c("regular", "projected to unit sphere"))+
  xlab("dimension of latent positions") + ylab("ari")
