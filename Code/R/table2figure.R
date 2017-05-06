rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/LLG/Code/R")

RE_ZG_Mean = list()
RE_ZG_UB = list()
RE_ZG_LB = list()
RE_USVT_Mean = list()
RE_USVT_UB = list()
RE_USVT_LB = list()

dataNameVec = c("JHU", "Desikan", "CPAC200")
mVec = c(1, 5, 10)

for (dataName in dataNameVec) {
  for (iM in 1:length(mVec)) {
    m = mVec[iM]
    
    # Eigen-decomposition
    if (dataName == "Desikan") {
      fileName = paste("../../Result/result_desikan_brute_", "m_", m, "_eig.RData", sep="")
    } else {
      fileName = paste("../../Result/result_", dataName, "_brute_", "m_", m, "_eig.RData", sep="")
    }
    load(fileName)
    
    error_A_bar = error_A_bar^2*n*(n-1)
    error_P_hat = error_P_hat^2*n*(n-1)
    
    error_P_hat_ZG = error_P_hat_ZG^2*n*(n-1)
    error_P_hat_USVT = error_P_hat_USVT^2*n*(n-1)
    
    RE_ZG_Mean[[dataName]][iM] = mean(error_P_hat_ZG/error_A_bar)
    RE_ZG_LB[[dataName]][iM] = RE_ZG_Mean[[dataName]][iM] -
      sqrt(var(error_P_hat_ZG/error_A_bar))/sqrt(length(error_A_bar))*1.96
    RE_ZG_UB[[dataName]][iM] = RE_ZG_Mean[[dataName]][iM] +
      sqrt(var(error_P_hat_ZG/error_A_bar))/sqrt(length(error_A_bar))*1.96
    
    RE_USVT_Mean[[dataName]][iM] = mean(error_P_hat_USVT/error_A_bar)
    RE_USVT_LB[[dataName]][iM] = RE_USVT_Mean[[dataName]][iM] -
      sqrt(var(error_P_hat_USVT/error_A_bar))/sqrt(length(error_A_bar))*1.96
    RE_USVT_UB[[dataName]][iM] = RE_USVT_Mean[[dataName]][iM] +
      sqrt(var(error_P_hat_USVT/error_A_bar))/sqrt(length(error_A_bar))*1.96
  }
}

library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)

df_ZG <- list()
df_USVT <- list()
df_Abar <- list()
for (dataName in dataNameVec) {
  df_ZG[[dataName]] <- data.frame(re = c(RE_ZG_Mean[[dataName]]),
                                  lci = c(RE_ZG_LB[[dataName]]),
                                  uci = c(RE_ZG_UB[[dataName]]),
                                  which = "ZG 3rd", m = mVec, dataname = dataName)
  df_USVT[[dataName]] <- data.frame(re = c(RE_USVT_Mean[[dataName]]),
                                    lci = c(RE_USVT_LB[[dataName]]),
                                    uci = c(RE_USVT_UB[[dataName]]),
                                    which = "USVT c=0.7", m = mVec, dataname = dataName)
  df_Abar[[dataName]] <- rbind(
    data.frame(re=1,lci=0,uci=0,which="Abar",m=mVec[1],dataname = dataName),
    data.frame(re=1,lci=0,uci=0,which="Abar",m=mVec[length(mVec)],dataname = dataName))
}

error_by_dim_df <- rbind(
  df_ZG[[dataNameVec[1]]], df_ZG[[dataNameVec[2]]], df_ZG[[dataNameVec[3]]],
  df_USVT[[dataNameVec[1]]], df_USVT[[dataNameVec[2]]], df_USVT[[dataNameVec[3]]],
  df_Abar[[dataNameVec[1]]], df_Abar[[dataNameVec[2]]], df_Abar[[dataNameVec[3]]]) %>%
  mutate(dataname = factor(dataname, dataNameVec))

label_y <- with(error_by_dim_df, .75*max(re)+.25*min(re))

gg <- ggplot(error_by_dim_df,aes(x=m,y=re,linetype=factor(which),shape=factor(which)))+
  facet_wrap(~dataname, nrow=1)+
  geom_point(data=error_by_dim_df,size=1.5)+
  scale_linetype_manual(name="",values=c(0,0,1),
                        labels=c("ZG 3rd","USVT c=0.7",expression(bar(A))))+
  scale_shape_manual(name="",values=c(17,15,-1),
                     labels=c("ZG 3rd","USVT c=0.7",expression(bar(A))))+
  scale_x_continuous(breaks=mVec)+
  geom_line()+
  # geom_hline(yintercept = 0, linetype="dashed")+
  # geom_hline(yintercept = 0)+
  xlab("M")+ylab("RE")+
  theme(panel.grid.major = element_line(colour="grey95"),
        panel.grid.minor = element_blank())+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
  theme(legend.position="bottom")

ggsave("../../Draft/corr_data_REdiff.pdf",
       plot=gg+theme(text=element_text(size=10,family="Times")),
       width=5.5,height=3)

