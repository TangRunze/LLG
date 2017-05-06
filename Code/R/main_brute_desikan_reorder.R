
# args <- commandArgs(trailingOnly = TRUE)

setwd("/Users/Runze/Documents/GitHub/LLG/Code/R")
# setwd("E:/GitHub/LLG/Code/R")
# setwd("/cis/home/rtang/LLG/Code/R")

# m = as.numeric(args[1])
# m = 1
# isSVD = 1

mVec = c(1, 5, 10)
# mVec = 20

for (m in mVec) {
  for (isSVD in 0) {
    
    print(c(m, isSVD))
    
    nIter = 200
    nCores = 2
    
        # dataName = "CPAC200"
    dataName = "desikan"
    #         dataName = "JHU"
    #     dataName = "slab907"
    #     dataName = "slab1068"
    #     dataName = "Talairach"
    
    
    source("function_collection.R")
    tmpList = read_data(dataName, DA=F, newGraph=F)
    A_all = tmpList[[1]]
    n = tmpList[[2]]
    M = tmpList[[3]]
    rm(tmpList)
    
    #     nZeroVec <- sapply(1:M, function(ind) {sum(A_all[[ind]] == 0)})
    #     hist(nZeroVec)
    
    dVec = 1:n
    nD = length(dVec)
    
    A_sum = add(A_all)
    
    error_P_hat = matrix(0, nD, nIter)
    error_A_bar = matrix(0, nD, nIter)
    
    require(parallel)
    
    # ptm <- proc.time()
    # proc.time() - ptm
    
    # out <- mclapply(1:nIter, function(x) sapply(dVec, function(d) dim_brute(M, m, d, A_all, A_sum)),
    #                 mc.cores=nCores)
    # out = array(unlist(out), dim = c(2, nD, nIter))
    # error_A_bar = out[1,,]
    # error_P_hat = out[2,,]
    
    out <- mclapply(1:nIter, function(x) dim_brute_reorder2(M, m, dVec, A_all, A_sum, isSVD), 
                    mc.cores=nCores)
    out = array(unlist(out), dim = c(2*(nD+3), nIter))
    
    error_A_bar = out[1,]
    error_P_hat = out[2:(nD+1),]
    dim_ZG = out[nD+2,]
    dim_USVT = out[nD+3,]
    error_P_hat_ZG = rep(0, length(dim_ZG))
    error_P_hat_USVT = rep(0, length(dim_USVT))
    for (i in 1:length(dim_ZG)) {
      error_P_hat_ZG[i] = error_P_hat[dim_ZG[i], i]
      error_P_hat_USVT[i] = error_P_hat[dim_USVT[i], i]
    }
    
    if (isSVD) {
      fileName = paste("../../Result/result_", dataName, "_brute_", "m_", m, "_svd.RData", sep="")
    } else {
      fileName = paste("../../Result/result_", dataName, "_brute_", "m_", m, "_eig.RData", sep="")
    }
    
    save(error_A_bar, error_P_hat, error_P_hat_ZG, error_P_hat_USVT, 
         dim_ZG, dim_USVT, n, M, m, dVec, nIter, file=fileName)
    
    error_A_bar = out[nD+3 + 1,]
    error_P_hat = out[(nD+3 + 2):(nD+3 + nD+1),]
    dim_ZG = out[(nD+3 + nD+2),]
    dim_USVT = out[(nD+3 + nD+3),]
    error_P_hat_ZG = rep(0, length(dim_ZG))
    error_P_hat_USVT = rep(0, length(dim_USVT))
    for (i in 1:length(dim_ZG)) {
      error_P_hat_ZG[i] = error_P_hat[dim_ZG[i], i]
      error_P_hat_USVT[i] = error_P_hat[dim_USVT[i], i]
    }
    
    if (isSVD) {
      fileName = paste("../../Result/result_", dataName, "_reorder_brute_", "m_", m, "_svd.RData", sep="")
    } else {
      fileName = paste("../../Result/result_", dataName, "_reorder_brute_", "m_", m, "_eig.RData", sep="")
    }
    
    save(error_A_bar, error_P_hat, error_P_hat_ZG, error_P_hat_USVT, 
         dim_ZG, dim_USVT, n, M, m, dVec, nIter, file=fileName)
  }
}