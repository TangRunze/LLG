rm(list = ls())
# setwd("E:/GitHub/LLG/Code/R")
setwd("/Users/Runze/Documents/GitHub/LLG/Code/R")
# setwd("/cis/home/rtang/LLG/Code/R")

pp <- list()

isSVD <- 0

lSize = .8
legendSize = 1.5

# mVec = c(5, 10, 20, 30, 40, 50)
mVec = c(1, 5, 10)
# mVec = 1
dZG3Mean = rep(0, length(mVec))
dZG3L = rep(0, length(mVec))
dZG3U = rep(0, length(mVec))
dUSVT07Mean = rep(0, length(mVec))
dUSVT07L = rep(0, length(mVec))
dUSVT07U = rep(0, length(mVec))

# dataName = "CPAC200"
dataName = "Desikan"
# dataName = "JHU"


library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)

for (iM in 1:length(mVec)) {
  m = mVec[iM]
  
  # Eigen-decomposition
  if (isSVD == 0) {
    fileName = paste("../../Result/result_", dataName, "_reorder_brute_", "m_", m, "_eig.RData", sep="")
  } else {
    fileName = paste("../../Result/result_", dataName, "_reorder_brute_", "m_", m, "_svd.RData", sep="")
  }
  load(fileName)
  
  dZG3Mean[iM] = mean(dim_ZG)
  dZG3L[iM] = mean(dim_ZG) - 1.96*sd(dim_ZG)/sqrt(length(dim_ZG))
  dZG3U[iM] = mean(dim_ZG) + 1.96*sd(dim_ZG)/sqrt(length(dim_ZG))
  dUSVT07Mean[iM] = mean(dim_USVT)
  dUSVT07L[iM] = mean(dim_USVT) - 1.96*sd(dim_USVT)/sqrt(length(dim_USVT))
  dUSVT07U[iM] = mean(dim_USVT) + 1.96*sd(dim_USVT)/sqrt(length(dim_USVT))
  
  if (iM == 1) {
    errorAbarMean = array(0, dim=c(length(mVec), 1))
    errorAbarLower = array(0, dim=c(length(mVec), 1))
    errorAbarUpper = array(0, dim=c(length(mVec), 1))
    errorPhatEIGMean = array(0, dim=c(length(mVec), n))
    errorPhatEIGLower = array(0, dim=c(length(mVec), n))
    errorPhatEIGUpper = array(0, dim=c(length(mVec), n))
  }
  
  error_A_bar = error_A_bar^2*n*(n-1)
  error_P_hat = error_P_hat^2*n*(n-1)
  
  errorAbarMean[iM] = rep(mean(error_A_bar))
  errorAbarLower[iM] = errorAbarMean[iM] - sqrt(var(error_A_bar))/sqrt(length(error_A_bar))*1.96
  errorAbarUpper[iM] = errorAbarMean[iM] + sqrt(var(error_A_bar))/sqrt(length(error_A_bar))*1.96
  errorPhatEIGMean[iM,] = rowMeans(error_P_hat)
  errorPhatEIGLower[iM,] = errorPhatEIGMean[iM,] - 
    sqrt(apply(error_P_hat, 1, var))/sqrt(dim(error_P_hat)[2])*1.96
  errorPhatEIGUpper[iM,] = errorPhatEIGMean[iM,] + 
    sqrt(apply(error_P_hat, 1, var))/sqrt(dim(error_P_hat)[2])*1.96
  
}

errorPhatUSVT07 = rep(0, length(mVec))
errorPhatZG3 = rep(0, length(mVec))
for (iM in 1:length(mVec)) {
  x = dUSVT07Mean[iM]
  x1 = floor(x)
  y1 = errorPhatEIGMean[iM, x1]
  x2 = ceiling(x)
  y2 = errorPhatEIGMean[iM, x2]
  errorPhatUSVT07[iM] = (y2-y1)/(x2-x1)*(x-x1)+y1
  
  x = dZG3Mean[iM]
  x1 = floor(x)
  y1 = errorPhatEIGMean[iM, x1]
  x2 = ceiling(x)
  y2 = errorPhatEIGMean[iM, x2]
  errorPhatZG3[iM] = (y2-y1)/(x2-x1)*(x-x1)+y1
}

error_by_dim_df <- rbind(
  data.frame(mse=errorAbarMean,lci=errorAbarLower,uci=errorAbarUpper,
             which="Abar",m=mVec,d=n),
  data.frame(mse=errorAbarMean,lci=errorAbarLower,uci=errorAbarUpper,
             which="Abar",m=mVec,d=1),
  data.frame(mse=c(errorPhatEIGMean),lci=c(errorPhatEIGLower),uci=c(errorPhatEIGUpper),
             which="Phat",m=rep(mVec,n),d=rep(1:n,each=length(mVec)))) %>%
  mutate(m=factor(paste0("M=",m),sapply(mVec, function(m) {paste0("M=", m)})))

dim_selection_df <- rbind(
  data.frame(mse=errorPhatZG3,lci=errorPhatZG3,uci=errorPhatZG3,
             which="ZG 3rd",m=mVec,d=dZG3Mean),
  data.frame(mse=errorPhatUSVT07,lci=errorPhatUSVT07,uci=errorPhatUSVT07,
             which="USVT c=0.7",m=mVec,d=dUSVT07Mean)) %>%
  mutate(m=factor(paste0("M=",m),sapply(mVec, function(m) {paste0("M=", m)})))


label_y <- with(error_by_dim_df, .75*max(mse)+.25*min(mse))

gg <- ggplot(error_by_dim_df,aes(x=d,y=mse,linetype=factor(which),shape=factor(which)))+
  facet_wrap(~m, nrow=1)+
  #   geom_point(data=subset(dim_selection_df,which=="ZG 3rd"),size=2,colour="red")+
  #   geom_point(data=subset(dim_selection_df,which=="USVT c=0.7"),size=2,colour="blue")+
  geom_point(data=dim_selection_df,size=1.5)+
  scale_linetype_manual(name="",values=c(1,2,0,0))+
  scale_shape_manual(name="",values=c(-1,-1,15,17))+
  #   geom_point(dim_selection_df,aes(shape=which))+
  geom_line()+
  geom_linerange(aes(ymin=lci,ymax=uci),linetype=1,alpha=.5,size=.5)+
  #   geom_vline(data=dim_selection_df,
  #              aes(xintercept=value,color=which,linetype=variable))+
  #   scale_linetype_manual(name="",values=c(1,2,3,4))+
  #   geom_text(data=dim_selection_df %>% filter(variable=="mean"),
  #             aes(x=value+n/30,y=label_y,linetype=variable,label=which,color=which),angle=90)+
  #   scale_color_discrete(guide=FALSE)+
  scale_y_continuous(limits = c(0, 0.50))+
  xlab("")+ylab("MSE")+
  # theme(strip.text.x = element_text(size=20,face="bold"))+
  # theme(axis.text=element_text(size=15),
  #       axis.title=element_text(size=20,face="bold"))+
  theme(panel.grid.major = element_line(colour="grey95"),
        panel.grid.minor = element_blank())+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
  theme(legend.position="none")+
  ggtitle(paste0(dataName, ", N=", n))

ggsave("../../Draft/corr_data_MSE_desikan_reorder.pdf",
       plot=gg+theme(text=element_text(size=10,family="Times")),
       width=3.5,height=2)
