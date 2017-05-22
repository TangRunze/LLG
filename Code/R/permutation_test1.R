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

uniform_select <- function(tmpA) {
  total <- sum(tmpA)
  r <- sample.int(total, 1)
  ind <- 1:((dim(tmpA)[1])^2)
  ind <- ind[c(tmpA) == 1]
  nVec <- ind[r] %% (dim(tmpA)[1])
  nVec <- c((ind[r] - nVec)/(dim(tmpA)[1]), nVec)
  if (nVec[2] > 0) {
    nVec[1] = nVec[1] + 1
  } else {
    nVec[2] = dim(tmpA)[1]
  }
  return(list(nVec[2], nVec[1]))
}


setwd("/Users/Runze/Documents/GitHub/LLG/Code/R")
# setwd("E:/GitHub/LLG/Code/R")
# setwd("/cis/home/rtang/LLG/Code/R")

source("function_collection.R")
source("getElbows.R")
require(ggplot2)

set.seed(12345)


indDim <- 1:8
nIter <- 1000
switchVec <- 1:10
# nIter <- 10
# switchVec <- 1:3

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
    iSwitch <- 1
    while (iSwitch <= nSwitch) {
      tmpA <- matrix(1, n, n)
      for (c in unique(cl)) {
        nv <- (cl == c)
        tmpA[nv, nv] <- 0
      }
      tmpA <- tmpA*A0
      nList <- uniform_select(tmpA)
      n1 <- nList[[1]]
      n2 <- nList[[2]]
      
      c1 <- cl[n1]
      c2 <- cl[n2]
      
      A1 <- matrix(0, n, n)
      nv1 <- (cl == c1)
      nv2 <- (cl == c2)
      A1[nv1, nv2] <- A0[nv1, nv2]
      A1[n1, ] <- 0
      A1[, n2] <- 0
      tmpA <- matrix(1, n, n)
      for (c in unique(cl)) {
        nv <- (cl == c)
        tmpA[nv, nv] <- 0
      }
      tmpA <- tmpA*A1
      if (sum(tmpA) > 0) {
        nList <- uniform_select(tmpA)
        n3 <- nList[[1]]
        n4 <- nList[[2]]
        
        if ((n1 != n3) | (n2 != n4)) {
          cl[n2] <- c1
          cl[n3] <- c2
        }
        iSwitch <- iSwitch + 1
      }
    }
    cl[(n/2 + 1):n] = cl[(n/2 + 1):n] + 10
    tVec[switchID, iIter] <- test1(xHat0, cl)
  }
}


# df <- data.frame(value=c(t0, t(tVec)), flip=c(0, rep(switchVec, each=nIter)))
# df <- data.frame(value=c(t0, t0 + 1e-2, t0 - 1e-2, t(tVec)), flip=c(0, 0, 0, rep(switchVec, each=nIter)))
# gg <- ggplot(data = df, aes(x=factor(flip), y=value))+
#   geom_boxplot(aes(fill=factor(flip)), notch = T)+
#   labs(title = paste0("dimension ", min(indDim), " to dimension ", max(indDim)),
#        x = "number of flips", y = "within lobes - cross lobes, 2-norm square", fill = "")
# ggsave(paste0("../../Draft/boxplot_flip_2norm^2_", min(indDim), "_", max(indDim), ".pdf"),
#        plot=gg+theme(text=element_text(size=10,family="Times")),
#        width=6, height=4)
# gg <- ggplot(data = df, aes(x=factor(flip), y=value, fill=factor(flip)))+
#   geom_violin(draw_quantiles = T)+
#   # geom_boxplot(aes(fill=factor(flip)), notch = T, width = 0.2)+
#   stat_summary(fun.y=mean, geom="point", size=2, show.legend = F)+
#   labs(title = paste0("dimension ", min(indDim), " to dimension ", max(indDim)),
#        x = "number of flips", y = "within lobes - cross lobes, 2-norm", fill = "")

df <- data.frame(value=c(t(tVec)), flip=rep(switchVec, each=nIter))
df0 <- data.frame(yi = t0)
gg <- ggplot(data = df, aes(x=factor(flip), y=value, fill=factor(flip)))+
  geom_violin(draw_quantiles = T)+
  geom_hline(data = df0, aes(yintercept = yi, linetype = factor(yi)), show.legend = TRUE) +
  scale_linetype_manual(name = "true lobe assigment", values = "dashed", labels = "0") +
  guides(fill=guide_legend(title="number of flips"))+
  theme(legend.position="bottom")+
  # geom_boxplot(aes(fill=factor(flip)), notch = T, width = 0.2)+
  # stat_summary(fun.y=mean, geom="point", size=2, show.legend = F)+
  labs(title = paste0("dimension ", min(indDim), " to dimension ", max(indDim)),
       x = "number of flips", y = "within lobes - cross lobes, 2-norm", fill = "")
ggsave(paste0("../../Draft/violinplot_new_flip_2norm_", min(indDim), "_", max(indDim), ".pdf"),
       plot=gg+theme(text=element_text(size=10,family="Times")),
       width=6, height=6)


pVec <- rep(0, length(switchVec))
for (i in 1:length(switchVec)) {
  pfun <- ecdf(tVec[i, ])
  pVec[i] <- pfun(t0)
}
df <- data.frame(value=pVec, flip=switchVec)
gg <- ggplot(data = df, aes(x=factor(flip), y=value))+
  geom_point(size=2)+
  # geom_boxplot(aes(fill=factor(flip)), notch = T, width = 0.2)+
  labs(title = paste0(""),
       x = "number of flips", y = "p-value", fill = "")
ggsave(paste0("../../Draft/pvalue_new_flip_2norm_", min(indDim), "_", max(indDim), ".pdf"),
       plot=gg+theme(text=element_text(size=10,family="Times")),
       width=6, height=4)
