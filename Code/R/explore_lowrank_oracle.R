load_abar <- function(dataName){
  tmpList <- read_data(dataName, DA=F, newGraph=F)
  A_all <- tmpList[[1]]
  n <- tmpList[[2]]
  M <- tmpList[[3]]
  rm(tmpList)
  add(A_all)/M
}

dataNameVec = c("JHU", "desikan", "CPAC200")
dataNameDisplayVec = c("JHU", "Desikan", "CPAC200")

isSVD <- 0

source("function_collection.R")
require(ggplot2)
for (iData in 1:length(dataNameVec)) {
  dataName = dataNameVec[iData]
  Abar = load_abar(dataName)
  n <- nrow(Abar)
  eval_sq = eigen(diag_aug(Abar))$values^2
  rankDError = sum(eval_sq) - cumsum(sort(eval_sq, decreasing = T))
  abarError <- sum(Abar*(1-Abar))

  plot(rankDError*10/abarError)
  
  yMax = max(eval_sq)
  yMin = min(eval_sq)
  df <- data.frame(eval=eval_sq, k=1:n)
  label_y <- with(df, .75*yMax+.25*yMin)
  
  pp_scree[[iData]] <- ggplot(df,aes(x=k,y=eval))+
    geom_line()+
    scale_linetype_manual(name="",values=c("longdash","dotted","dotdash"))+
    xlab("dimension d") + ylab("relative error")+
    theme(panel.grid.major = element_line(colour="grey95"),
          panel.grid.minor = element_blank())+
    theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
    theme(legend.position="none")+
    ggtitle(dataNameDisplayVec[[iData]])
  
  ggsave(paste0("../../Draft/screeplot_ratio_", dataName, ".pdf"),
         pp_scree[[iData]]+theme(text=element_text(size=10,family="Times")),
         width=2, height=2)
  
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
  eval_sq = eigen(Ahat)$values^2
  eval_sq = 1 - cumsum(sort(abs(eval_sq), decreasing = T))/sum(abs(eval_sq))
  
  yMax = max(eval_sq)
  yMin = min(eval_sq)
  df <- data.frame(eval=eval_sq, k=1:n)
  label_y <- with(df, .75*yMax+.25*yMin)
  
  pp_scree[[iData]] <- ggplot(df,aes(x=k,y=eval))+
    geom_line()+
    scale_linetype_manual(name="",values=c("longdash","dotted","dotdash"))+
    xlab("dimension d") + ylab("relative error")+
    theme(panel.grid.major = element_line(colour="grey95"),
          panel.grid.minor = element_blank())+
    theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
    theme(legend.position="none")+
    ggtitle(dataNameDisplayVec[[iData]])
  
  ggsave(paste0("../../Draft/screeplot_ratio_", dataName, "_lowrank.pdf"),
         pp_scree[[iData]]+theme(text=element_text(size=10,family="Times")),
         width=2, height=2)
  
}

