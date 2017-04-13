rm(list = ls())

pp <- list()

lSize = .8
legendSize = 1.5

################## JHU #########################
dataName = "JHU"

mVec = c(1, 5, 10)

# 2nd and 3rd elbow
if (dataName == "JHU") {
  dZG3Mean = c(15.81, 13.59, 12.95)
  dZG3L = c(15.41929, 13.25443, 12.6837)
  dZG3U = c(16.20071, 13.92557, 13.2163)
} else if (dataName == "Desikan") {
  dZG3Mean = c(15.95, 11.09, 9.87)
  dZG3L = c(15.34967, 10.70389, 9.629296)
  dZG3U = c(16.55033, 11.47611, 10.110704)
} else if (dataName == "CPAC200") {
  dZG3Mean = c(60.91, 47.64, 42.54)
  dZG3L = c(59.17358, 46.85278, 41.98464)
  dZG3U = c(62.64642, 48.42722, 43.09536)
}


# USVT
if (dataName == "JHU") {
  dUSVT07Mean = c(4.61, 11.61, 18.18)
  dUSVT07L = c(4.465526, 11.39987, 17.95997)
  dUSVT07U = c(4.754474, 11.82013, 18.40003)
} else if (dataName == "Desikan") {
  dUSVT07Mean = c(6.36, 12.68, 19.55)
  dUSVT07L = c(6.198711, 12.47551, 19.35032)
  dUSVT07U = c(6.521289, 12.88449, 19.74968)
} else if (dataName == "CPAC200") {
  dUSVT07Mean = c(10.08, 18.86, 27.77)
  dUSVT07L = c(9.732408, 18.51539, 27.44055)
  dUSVT07U = c(10.427592, 19.20461, 28.09945)
}


library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)

for (iM in 1:length(mVec)) {
  m = mVec[iM]
  
  # Eigen-decomposition
  fileName = paste("../../Result/result_", dataName, "_brute_", "m_", m, "_eig.RData", sep="")
  load(fileName)
  
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
             which="Phat",m=rep(mVec,n),d=rep(1:n,each=3))) %>%
  mutate(m=factor(paste0("M=",m),c("M=1","M=5","M=10")))

dim_selection_df <- rbind(
  data.frame(mse=errorPhatZG3,lci=errorPhatZG3,uci=errorPhatZG3,
             which="ZG 3rd",m=mVec,d=dZG3Mean),
  data.frame(mse=errorPhatUSVT07,lci=errorPhatUSVT07,uci=errorPhatUSVT07,
             which="USVT c=0.7",m=mVec,d=dUSVT07Mean)) %>%
  mutate(m=factor(paste0("M=",m),c("M=1","M=5","M=10")))


label_y <- with(error_by_dim_df, .75*max(mse)+.25*min(mse))

gg <- ggplot(error_by_dim_df,aes(x=d,y=mse,linetype=factor(which),shape=factor(which)))+
  facet_wrap(~m)+
  geom_point(data=dim_selection_df,size=1.5)+
  scale_linetype_manual(name="",values=c(1,2,0,0))+
  scale_shape_manual(name="",values=c(-1,-1,15,17))+
  geom_line()+
  geom_linerange(aes(ymin=lci,ymax=uci),linetype=1,alpha=.5,size=.5)+
  xlab("")+ylab("MSE")+
  theme(panel.grid.major = element_line(colour="grey95"),
        panel.grid.minor = element_blank())+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
  theme(legend.position="bottom")+
  ggtitle(paste0(dataName, ", N=", n))

pp[[1]]=gg

ggsave("../../Draft/corr_data_MSE_jhu_poster.pdf",
  plot=gg+theme(text=element_text(size=10,family="CM Roman")),
    width=5.5,height=2.5)


################## Desikan #########################
# dataName = "CPAC200"
dataName = "Desikan"
# dataName = "JHU"

mVec = c(1, 5, 10)

# 2nd and 3rd elbow
if (dataName == "JHU") {
  dZG3Mean = c(15.81, 13.59, 12.95)
  dZG3L = c(15.41929, 13.25443, 12.6837)
  dZG3U = c(16.20071, 13.92557, 13.2163)
} else if (dataName == "Desikan") {
  dZG3Mean = c(15.95, 11.09, 9.87)
  dZG3L = c(15.34967, 10.70389, 9.629296)
  dZG3U = c(16.55033, 11.47611, 10.110704)
} else if (dataName == "CPAC200") {
  dZG3Mean = c(60.91, 47.64, 42.54)
  dZG3L = c(59.17358, 46.85278, 41.98464)
  dZG3U = c(62.64642, 48.42722, 43.09536)
}


# USVT
if (dataName == "JHU") {
  dUSVT07Mean = c(4.61, 11.61, 18.18)
  dUSVT07L = c(4.465526, 11.39987, 17.95997)
  dUSVT07U = c(4.754474, 11.82013, 18.40003)
} else if (dataName == "Desikan") {
  dUSVT07Mean = c(6.36, 12.68, 19.55)
  dUSVT07L = c(6.198711, 12.47551, 19.35032)
  dUSVT07U = c(6.521289, 12.88449, 19.74968)
} else if (dataName == "CPAC200") {
  dUSVT07Mean = c(10.08, 18.86, 27.77)
  dUSVT07L = c(9.732408, 18.51539, 27.44055)
  dUSVT07U = c(10.427592, 19.20461, 28.09945)
}


library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)

for (iM in 1:length(mVec)) {
  m = mVec[iM]
  
  # Eigen-decomposition
  fileName = paste("../../Result/result_", dataName, "_brute_", "m_", m, "_eig.RData", sep="")
  load(fileName)
  
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
             which="Phat",m=rep(mVec,n),d=rep(1:n,each=3))) %>%
  mutate(m=factor(paste0("M=",m),c("M=1","M=5","M=10")))

dim_selection_df <- rbind(
  data.frame(mse=errorPhatZG3,lci=errorPhatZG3,uci=errorPhatZG3,
             which="ZG 3rd",m=mVec,d=dZG3Mean),
  data.frame(mse=errorPhatUSVT07,lci=errorPhatUSVT07,uci=errorPhatUSVT07,
             which="USVT c=0.7",m=mVec,d=dUSVT07Mean)) %>%
  mutate(m=factor(paste0("M=",m),c("M=1","M=5","M=10")))


label_y <- with(error_by_dim_df, .75*max(mse)+.25*min(mse))

gg <- ggplot(error_by_dim_df,aes(x=d,y=mse,linetype=factor(which),shape=factor(which)))+
  facet_wrap(~m)+
  geom_point(data=dim_selection_df,size=1.5)+
  scale_linetype_manual(name="",values=c(1,2,0,0))+
  scale_shape_manual(name="",values=c(-1,-1,15,17))+
  geom_line()+
  geom_linerange(aes(ymin=lci,ymax=uci),linetype=1,alpha=.5,size=.5)+
  xlab("")+ylab("MSE")+
  theme(panel.grid.major = element_line(colour="grey95"),
        panel.grid.minor = element_blank())+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
  theme(legend.position="bottom")+
  ggtitle(paste0(dataName, ", N=", n))

pp[[2]]=gg

ggsave("../../Draft/corr_data_MSE_desikan_poster.pdf",
  plot=gg+theme(text=element_text(size=10,family="CM Roman")),
    width=5.5,height=2.5)

################## CPAC200 #########################
dataName = "CPAC200"
# dataName = "Desikan"
# dataName = "JHU"

mVec = c(1, 5, 10)

# 2nd and 3rd elbow
if (dataName == "JHU") {
  dZG3Mean = c(15.81, 13.59, 12.95)
  dZG3L = c(15.41929, 13.25443, 12.6837)
  dZG3U = c(16.20071, 13.92557, 13.2163)
} else if (dataName == "Desikan") {
  dZG3Mean = c(15.95, 11.09, 9.87)
  dZG3L = c(15.34967, 10.70389, 9.629296)
  dZG3U = c(16.55033, 11.47611, 10.110704)
} else if (dataName == "CPAC200") {
  dZG3Mean = c(60.91, 47.64, 42.54)
  dZG3L = c(59.17358, 46.85278, 41.98464)
  dZG3U = c(62.64642, 48.42722, 43.09536)
}


# USVT
if (dataName == "JHU") {
  dUSVT07Mean = c(4.61, 11.61, 18.18)
  dUSVT07L = c(4.465526, 11.39987, 17.95997)
  dUSVT07U = c(4.754474, 11.82013, 18.40003)
} else if (dataName == "Desikan") {
  dUSVT07Mean = c(6.36, 12.68, 19.55)
  dUSVT07L = c(6.198711, 12.47551, 19.35032)
  dUSVT07U = c(6.521289, 12.88449, 19.74968)
} else if (dataName == "CPAC200") {
  dUSVT07Mean = c(10.08, 18.86, 27.77)
  dUSVT07L = c(9.732408, 18.51539, 27.44055)
  dUSVT07U = c(10.427592, 19.20461, 28.09945)
}


library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)

for (iM in 1:length(mVec)) {
  m = mVec[iM]
  
  # Eigen-decomposition
  fileName = paste("../../Result/result_", dataName, "_brute_", "m_", m, "_eig.RData", sep="")
  load(fileName)
  
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
             which="Phat",m=rep(mVec,n),d=rep(1:n,each=3))) %>%
  mutate(m=factor(paste0("M=",m),c("M=1","M=5","M=10")))

dim_selection_df <- rbind(
  data.frame(mse=errorPhatZG3,lci=errorPhatZG3,uci=errorPhatZG3,
             which="ZG 3rd",m=mVec,d=dZG3Mean),
  data.frame(mse=errorPhatUSVT07,lci=errorPhatUSVT07,uci=errorPhatUSVT07,
             which="USVT c=0.7",m=mVec,d=dUSVT07Mean)) %>%
  mutate(m=factor(paste0("M=",m),c("M=1","M=5","M=10")))


label_y <- with(error_by_dim_df, .75*max(mse)+.25*min(mse))

gg <- ggplot(error_by_dim_df,aes(x=d,y=mse,linetype=factor(which),shape=factor(which)))+
  facet_wrap(~m)+
  geom_point(data=dim_selection_df,size=1.5)+
  scale_linetype_manual(name="",values=c(1,2,0,0))+
  scale_shape_manual(name="",values=c(-1,-1,15,17))+
  geom_line()+
  geom_linerange(aes(ymin=lci,ymax=uci),linetype=1,alpha=.5,size=.5)+
  xlab("dimension")+ylab("MSE")+
  theme(panel.grid.major = element_line(colour="grey95"),
        panel.grid.minor = element_blank())+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
  theme(legend.position="bottom")+
  ggtitle(paste0(dataName, ", N=", n))

pp[[3]]=gg

ggsave("../../Draft/corr_data_MSE_CPAC200_poster.pdf",
  plot=gg+theme(text=element_text(size=10,family="CM Roman")),
    width=5.5,height=2.5)

source("function_collection.R")

library(gridExtra)
library(grid)

grid_arrange_shared_legend2(list(pp[[1]], pp[[2]], pp[[3]]), 3, 1)

