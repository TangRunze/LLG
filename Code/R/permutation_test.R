test1 <- function(xHat, cl) {
  diffWithin <- c()
  diffCross <- c()
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (cl[i] == cl[j]) {
        # diffWithin <- c(diffWithin, norm(xHat[i, ] - xHat[j, ], "2")^2)
        diffWithin <- c(diffWithin, norm(xHat[i, ] - xHat[j, ], "2"))
      } else {
        # diffCross <- c(diffCross, norm(xHat[i, ] - xHat[j, ], "2")^2)
        diffCross <- c(diffCross, norm(xHat[i, ] - xHat[j, ], "2"))
      }
    }
  }
  return(mean(diffWithin) - mean(diffCross))
}



setwd("/Users/Runze/Documents/GitHub/LLG/Code/R")
# setwd("E:/GitHub/LLG/Code/R")
# setwd("/cis/home/rtang/LLG/Code/R")

source("function_collection.R")
source("getElbows.R")
require(ggplot2)

set.seed(12345)


indDim <- 1:8
nIter <- 100
switchVec <- 1:10

isSVD <- 0
dataName <- "desikan"

###### Get xHat ######
tmpList <- read_data(dataName, DA=F, newGraph=F)
A_all <- tmpList[[1]]
n <- tmpList[[2]]
M <- tmpList[[3]]
rm(tmpList)
A_sum <- add(A_all)
A_bar <- add(A_all)/M
A_bar_diag_aug <- diag_aug(A_bar)
nElbow <- 3
evalVec <- ase(A_bar, ceiling(n*3/5), isSVD)[[1]]
dZG <- getElbows(evalVec, n=nElbow, plot=F)[[nElbow]]
A.ase <- ase(A_bar_diag_aug, dZG, isSVD)
xHat0 <- A.ase[[3]] %*% diag(sqrt(A.ase[[1]]))
xHat0 <- xHat0[, indDim]

###### Get lobe ID ######
data <- read.csv("../../Data/desikan_lobe.csv", header = F)
cl0 <- data[1:n, 2]

###### Get adjacent matrix in Desikan atlas ######
A0 <- as.matrix(read.csv("../../Data/adjacent_list.csv", header = F))

clTmp <- cl0
clTmp[(n/2 + 1):n] = clTmp[(n/2 + 1):n] + 10
t0 <- test1(xHat0, clTmp)

tVec <- matrix(rep(0, length(switchVec)*nIter), ncol = nIter)
for (switchID in 1:length(switchVec)) {
  nSwitch <- switchVec[switchID]
  for (iIter in 1:nIter) {
    
    print(c(nSwitch, iIter))
    
    ###### Switch the vertices ######
    cl <- cl0
    A <- A0
    xHat <- xHat0
    for (iSwitch in 1:nSwitch) {
      # Select pairs to switch
      flag <- F
      while (flag == F) {
        nv1 <- sample(1:n, 1)
        nv2 <- 1:n
        tmp <- (A[nv1, ] == 1) & (cl[nv2] != cl[nv1])
        if (sum(tmp) > 0) {
          nv2 <- nv2[tmp]
          nv2 <- nv2[sample(1:length(nv2), 1)]
          flag = T
        }
      }
      
      # Switch
      P <- diag(n)
      P[nv1, nv1] <- 0
      P[nv2, nv2] <- 0
      P[nv1, nv2] <- 1
      P[nv2, nv1] <- 1
      A <- P %*% A %*% P
      
      tmp <- cl[nv1]
      cl[nv1] <- cl[nv2]
      cl[nv2] <- tmp
      
      tmp <- xHat[nv1, ]
      xHat[nv1, ] <- xHat[nv2, ]
      xHat[nv2, ] <- tmp
    }
    
    clTmp <- cl0
    clTmp[(n/2 + 1):n] = clTmp[(n/2 + 1):n] + 10
    tVec[switchID, iIter] <- test1(xHat, clTmp)
  }
}


# df <- data.frame(value=c(t0, t(tVec)), flip=c(0, rep(switchVec, each=nIter)))
df <- data.frame(value=c(rep(t0, 3), t(tVec)), flip=c(rep(0, 3), rep(switchVec, each=nIter)))
# gg <- ggplot(data = df, aes(x=factor(flip), y=value))+
#   geom_boxplot(aes(fill=factor(flip)), notch = T)+
#   labs(title = paste0("dimension ", min(indDim), " to dimension ", max(indDim)),
#        x = "number of flips", y = "within lobes - cross lobes, 2-norm square", fill = "")
# ggsave(paste0("../../Draft/boxplot_flip_2norm^2_", min(indDim), "_", max(indDim), ".pdf"),
#        plot=gg+theme(text=element_text(size=10,family="Times")),
#        width=6, height=4)
gg <- ggplot(data = df, aes(x=factor(flip), y=value))+
  geom_violin(aes(fill=factor(flip)))+
  labs(title = paste0("dimension ", min(indDim), " to dimension ", max(indDim)),
       x = "number of flips", y = "within lobes - cross lobes, 2-norm", fill = "")
ggsave(paste0("../../Draft/boxplot_flip_2norm_", min(indDim), "_", max(indDim), ".pdf"),
       plot=gg+theme(text=element_text(size=10,family="Times")),
       width=6, height=4)



