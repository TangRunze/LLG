rm(list = ls())
# setwd("E:/GitHub/LLG/Code/R")
setwd("/Users/Runze/Documents/GitHub/LLG/Code/R")
# setwd("/cis/home/rtang/LLG/Code/R")

RE_ZG_Mean = list()
RE_ZG_UB = list()
RE_ZG_LB = list()
RE_USVT_Mean = list()
RE_USVT_UB = list()
RE_USVT_LB = list()

# dataName = "JHU"
dataName = "desikan"
# dataName = "CPAC200"

# mVec = c(1, 5, 10, 50)
mVec = c(1, 5, 10)

for (iM in 1:length(mVec)) {
  m = mVec[iM]
  
  # Eigen-decomposition
  # fileName = paste("../../Result/result_", dataName, "_reorder_brute_", "m_", m, "_eig.RData", sep="")
  fileName = paste("../../Result/result_", dataName, "_brute_", "m_", m, "_eig.RData", sep="")
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
  
  print(sqrt(var(error_P_hat_ZG/error_A_bar))/sqrt(length(error_A_bar))*1.96)
  print(sqrt(var(error_P_hat_USVT/error_A_bar))/sqrt(length(error_A_bar))*1.96)
}

