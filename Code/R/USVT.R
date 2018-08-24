require("irlba")
# A is the adjacency matrix, Abar, or igraph object
# minS is the smallest number of singular values to use
# m is the number of graphs used to compute A
# usvt <- function(A,minS=1, m=1){
#   if(class(A)=="igraph"){
#     A <- get.adjacency(A)
#   }
#   
#   
#   n <- nrow(A)
#   
#   # the threshold
#   tau <- sqrt(n/m)
#   s <- minS-1
#   usv <- NULL
#   # Collect singular values and vectors until you reach a singular value less than tau
#   repeat{
#     s <- s+1
#     usv.new <- irlba.r(A,s,usv=usv)
#     d <- usv.new$d[s]
#     if(d < tau){break}
#     else{ usv <- usv.new }
#   }
#   
#   return(usv)
# }
# 
# 
# # this code allows you to repeatedly collect more and more dimensions of without recomputing the old dimensions
# irlba.r <- function(A,s,usv=NULL){
#   if(is.null(usv)){
#     return(irlba(A,s))
#   }
#   s.old <- length(usv$d)
#   if(s.old < s){
#     n <- nrow(A)
#     U <- cbind(usv$u, scale(matrix(rnorm(n*(s-s.old)),nrow=n),center=F,scale=T)/sqrt(n-1))
#     D <- c(usv$d,rep(1,s-s.old))
#     return(irlba(A,s,s,v=U,tol = 1e-5))
#   }
#   else{
#     return(usv)
#   }
# }



usvt <- function(A, minS=1, m=1, sigma  = 0.7){
  if(class(A)=="igraph"){
    A <- get.adjacency(A)
  }
  
  
  n <- nrow(A)
  
  # the threshold
  tau <- sigma * sqrt(n/m)
  usv <- svd(A)
  s = sum(usv$d>=tau)
  usv$d = usv$d[1:s]
  return(usv)
}