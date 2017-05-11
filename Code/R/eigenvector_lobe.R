
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

xHat <- A.ase[[3]] %*% diag(sqrt(A.ase[[1]]))


data <- read.csv("../../Data/desikan_lobe.csv", header = F)
clName <- data[(n + 1):dim(data)[1],][, 1]
vertexName <- data[1:n, 1]
cl <- data[1:n, 2]
cl0 <- cl
cl[36:70] <- cl[36:70] + 10
nv <- order(cl)

require(ggplot2)

# df <- data.frame(value=c(xHat),
#                  v=rep(sapply(1:n, function(i) {if (i < 10) {return(paste0("0", i))} else {return(paste0(i))}}), times=d),
#                  d=rep(sapply(1:d, function(i) {if (i < 10) {return(paste0("0", i))} else {return(paste0(i))}}), each=n))
df <- data.frame(value=c(xHat[nv, ]),
                 v=rep(sapply(1:n, function(i) {
                   if (cl[nv[i]] < 10) {
                     return(paste0("lh, ", clName[cl0[nv[i]]], ", ", vertexName[nv[i]]))
                   } else {
                     return(paste0("rh, ", clName[cl0[nv[i]]], ", ", vertexName[nv[i]]))
                   }
                 }), times=d),
                 d=rep(sapply(1:d, function(i) {
                   if (i < 10) {
                     return(paste0("0", i))
                   } else {
                     return(paste0(i))
                   }
                 }), each=n))

# gg <- ggplot(df, aes(d, v)) + 
#   geom_tile(aes(fill = value), colour = "white") +
#   # scale_fill_gradient(low = "white", high = "steelblue")
#   scale_fill_gradient(low = "white", high = "grey10") + 
#   xlab("dimension") + ylab("vertex")

gg <- ggplot(df, aes(d, v)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2(low = colorVec[1], mid = colorVec[2], high = colorVec[3]) + 
  xlab("dimension") + ylab("vertex")

ggsave("../../Draft/eigenvector.pdf",
       plot=gg+theme(text=element_text(size=10,family="Times")),
       width=6,height=10)






xHatLobe <- matrix(rep(0, length(unique(cl))*d), ncol = d)
i <- 0
for (c in unique(cl)) {
  i <- i + 1
  nv <- (cl == c)
  xHatLobe[i, ] <- colMeans(xHat[nv, ])
}

df <- data.frame(value=c(xHatLobe),
                 v=rep(sapply(unique(cl), function(i) {
                   if (i < 10) {
                     return(paste0("lh, ", clName[i]))
                   } else {
                     return(paste0("rh, ", clName[i - 10]))
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
  # scale_fill_gradient(low = "white", high = "steelblue") +
  scale_fill_gradient(low = "white", high = "grey10") +
  xlab("dimension") + ylab("vertex")

gg <- ggplot(df, aes(d, v)) + 
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient2(low = colorVec[1], mid = colorVec[2], high = colorVec[3]) + 
  xlab("dimension") + ylab("vertex")

ggsave("../../Draft/eigenvector_lobe.pdf",
       plot=gg+theme(text=element_text(size=10,family="Times")),
       width=6,height=6)
