
setwd("/Users/Runze/Documents/GitHub/LLG/Code/R")
# setwd("E:/GitHub/LLG/Code/R")
# setwd("/cis/home/rtang/LLG/Code/R")

isSVD = 0
dataName = "desikan"

colorVec <- c("#af8dc3", "#f7f7f7", "#7fbf7b")
# colorVec <- c("#d8b365", "#f5f5f5", "#5ab4ac")

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

xHat <- A.ase[[3]] %*% diag(sqrt(A.ase[[1]]))


data <- read.csv("../../Data/desikan_lobe.csv", header = F)
clName <- data[(n + 1):dim(data)[1],][, 1]
vertexName <- data[1:n, 1]
cl <- data[1:n, 2]
cl0 <- cl
# cl[36:70] <- cl[36:70] + 10
cl[1:35] <- cl[1:35]*2 - 1
cl[36:70] <- cl[36:70]*2
nv <- order(cl)

require(ggplot2)

xHat <- sign(xHat)
xHatLobe <- matrix(rep(0, length(unique(cl))*d), ncol = d)
for (i in 1:10) {
  nv <- (cl == i)
  xHatLobe[i, ] <- colMeans(xHat[nv, ])
}

# df <- data.frame(value=c(xHatLobe),
#                  v=rep(sapply(unique(cl), function(i) {
#                    if (i < 10) {
#                      return(paste0("lh, ", clName[i]))
#                    } else {
#                      return(paste0("rh, ", clName[i - 10]))
#                    }
#                  }), times=d),
#                  d=rep(sapply(1:d, function(i) {
#                    if (i < 10) {
#                      return(paste0("0", i))
#                    } else {
#                      return(paste0(i))
#                    }
#                  }), each=length(unique(cl))))
df <- data.frame(value=c(xHatLobe),
                 v=rep(sapply(1:10, function(i) {
                   if (i %% 2 == 1) {
                     return(paste0(clName[(i + 1)/2], ", lh"))
                   } else {
                     return(paste0(clName[i/2], ", rh"))
                   }
                 }), times=d),
                 d=rep(sapply(1:d, function(i) {
                   if (i < 10) {
                     return(paste0("0", i))
                   } else {
                     return(paste0(i))
                   }
                 }), each=length(unique(cl))))

gg <- ggplot(df, aes(d, v)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2(low = colorVec[1], mid = colorVec[2], high = colorVec[3]) + 
  xlab("dimension") + ylab("vertex")

ggsave("../../Draft/eigenvector_lobe_threshold.pdf",
       plot=gg+theme(text=element_text(size=10,family="Times")),
       width=6,height=6)
