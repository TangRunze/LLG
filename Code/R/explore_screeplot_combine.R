rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/LLG/Code/R")

dataNameVec = c("JHU", "desikan", "CPAC200")
dataNameDisplayVec = c("JHU", "Desikan", "CPAC200")

isSVD <- 0

source("function_collection.R")
require(ggplot2)
df <- list()
df_LR <- list()
pp <- list()
for (iData in 1:length(dataNameVec)) {
  dataName = dataNameVec[iData]
  tmpList = read_data(dataName, DA=F, newGraph=F)
  A_all = tmpList[[1]]
  n = tmpList[[2]]
  M = tmpList[[3]]
  rm(tmpList)
  Abar = add(A_all)/M
  AbarDiagAug = diag_aug(Abar)
  eigenResult = eigen(AbarDiagAug)$values^2
  eigenResult = 1 - cumsum(sort(abs(eigenResult), decreasing = T))/sum(abs(eigenResult))
  
  df[[iData]] <- data.frame(re = c(eigenResult), d = 1:n,
                            which = "raw data",
                            dataname = dataNameDisplayVec[[iData]])
  
  source("getElbows.R")
  nElb = 3
  dMax = ceiling(n*3/5)
  evalVec = ase(Abar, dMax, isSVD)[[1]]
  dHat = getElbows(evalVec, n=nElb, plot=F)[[nElb]]
  
  A.ase = ase(diag_aug(Abar), dHat, isSVD)
  if (dHat == 1) {
    Ahat = A.ase[[1]] * A.ase[[3]] %*% t(A.ase[[2]])
  } else {
    Ahat <- A.ase[[3]][,1:dHat] %*% diag(A.ase[[1]][1:dHat]) %*% t(A.ase[[2]][,1:dHat])
  }
  eigenResult = eigen(Ahat)$values^2
  eigenResult = 1 - cumsum(sort(abs(eigenResult), decreasing = T))/sum(abs(eigenResult))
  
  df_LR[[iData]] <- data.frame(re = c(eigenResult), d = 1:n,
                               which = "low-rank data",
                               dataname = dataNameDisplayVec[[iData]])
  
  
  error_by_dim_df <- rbind(df[[iData]], df_LR[[iData]])
  
  label_y <- with(error_by_dim_df, .75*max(re)+.25*min(re))
  
  pp[[iData]] <- ggplot(error_by_dim_df,aes(x=d,y=re,linetype=factor(which)))+
    # facet_wrap(~dataname, nrow=1)+
    scale_linetype_manual(name="",values=c(1, 2))+
    geom_line()+
    xlab("dimension d")+ylab("relative error")+
    theme(panel.grid.major = element_line(colour="grey95"),
          panel.grid.minor = element_blank())+
    theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
    # theme(legend.position="bottom")
    theme(legend.position="none")+
    ggtitle(dataNameDisplayVec[[iData]])
  
}

library(gridExtra)
library(grid)

pdf(paste0("../../Draft/screeplot_ratio_all.pdf"), onefile=FALSE, family="Times", width=6, height=2.5)
pp_all <- grid_arrange_shared_legend2(list(pp[[1]], pp[[2]], pp[[3]]), 1, 3)
dev.off()
