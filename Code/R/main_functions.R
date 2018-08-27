
library(irlba)


source("getElbows.R")
source("USVT.R")

# ASE using SVD or eigen-decomposition.
decompose <- function(A, dim, is_svd = TRUE){
  n <- nrow(A)
  if (is_svd) {
    if (dim < n / 2){
      A.svd <- irlba(A, nu = dim, nv = dim)
    } else {
      A.svd <- svd(A, dim, dim)
      A.svd$d <- A.svd$d[1:dim]
    }
    values <- A.svd$d
    lvectors <- A.svd$u
    rvectors <- A.svd$v
  }
  else {
    if (dim < n / 2){
      shift <- n * max(A)

      A.eig <- irlba(A, nu = dim, nv = dim, shift = shift)

      values <- A.eig$d - shift
      lvectors <- A.eig$u
      rvectors <- A.eig$v
    }
    else{
      A.eig <- eigen(A)
      A.eig$values <- A.eig$values[1:dim]
      A.eig$vectors <- A.eig$vectors[,1:dim]

      values <- A.eig$values
      lvectors <- A.eig$vectors
      rvectors <- lvectors
    }

  }
  return(list(values = values, rvectors = rvectors, lvectors = lvectors))
}

add_dims <- function(mat, decomposition, dims){
  mat + with(decompose,
    lvectors[, dims] %*%
      Matrix::diag(values[dims]) %*%
      t(rvectors[, dims]))
}

low_rank_approx_dec <- function(dec){
  ndim <- length(dec$values)
  if ( ndim > 1 ){
    with(dec, lvectors %*% diag(values) %*% t(rvectors))
  } else{
    with(dec, outer(as.vector(lvectors), as.vector(rvectors)) * values)
  }
}


low_rank_approx <- function(A, ndim, is_svd = TRUE){
  if (is_svd){
    usv <- irlba(A, ndim)
    if ( ndim > 1 ){
      with(usv, u %*% diag(d) %*% t(v))
    } else{
      with(usv, outer(as.vector(u), as.vector(v)) * d)
    }
  } else {
    require(rARPACK)
    A.eig <- eigs_sym(matrix(A, ncol=dim(A)[1]), ndim, which = "LA")
    if ( ndim > 1 ){
      with(A.eig, vectors %*% diag(values) %*% t(vectors))
    } else{
      with(A.eig, outer(as.vector(vectors), as.vector(vectors)) * values)
    }
  }
}

get_dim_zg <- function(eval, is_svd = TRUE, elbow = 3){
  getElbows(eval, n = elbow, plot = F)[[elbow]]
}

get_dim_usvt <- function(eval, n, m, t = 0.7){
  tau <- t*sqrt(n/m)
  sum(eval > tau)
}

spectral_diag_aug <- function(dec, ...){
  with(dec, rowSums((lvectors %*% Matrix::Diagonal(x = values)) * rvectors))
}

compute_phat <- function(alist, abar = NULL,
  dim = "ZG",
  diag_aug = FALSE, diag_aug_sec = FALSE,
  threshold = TRUE, is_svd = FALSE, ...){

  m <- length(alist)
  if (is.null(abar)){
    abar <- Reduce("+", alist) / m
  }

  n <- nrow(abar)

  # First diagonal augmentation
  if (diag_aug){
    diag(abar) <- rowSums(abar) / (n - 1)
  }

  # select dimension
  if ( !is.numeric(dim) ){
    # Automatic so first compute big decompositions
    abar_dec <- decompose(abar, n, is_svd)

    # Zhu and Ghodsi
    if (dim == "ZG"){
      dim <- get_dim_zg(abar_dec$values, is_svd, ...)
    }
    # Universal Singular Value Thresholding
    if ( dim == "USVT" ){
      dim <- get_dim_usvt(abar_dec$values, n, m, ...)
    }

    # Truncate the decomposition
    abar_dec <- with(abar_dec, 
      list(values = values[1:dim], lvectors = lvectors[, 1:dim],
        rvectors = rvectors[, 1:dim]))
  } else{
    # Compute trunctated decomposition
    abar_dec <- decompose(abar, dim, is_svd)
  }

  # Second diagonal augmentation
  if( diag_aug_sec ){
    diag(abar) <- spectral_diag_aug(abar_dec)
    abar_dec <- decompose(abar, dim, is_svd)
  }


  phat <- low_rank_approx_dec(abar_dec)

  if (threshold){
    diag(phat) <- 0
    phat[phat > 1] <- 1
    phat[phat < 0] <- 0
  }
  list(phat = phat, dim = dim)

}


low_rank_from_dec <- function(dec, dim){
  d <- seq(dim)
  if (dim != 1){
    with(dec,
      lvectors[, d] %*% Matrix::diag(values[d]) %*% t(rvectors[, d]))
  } else{
     with(dec,
      outer(lvectors[, 1], rvectors[, d]) %*%  values[1])
   }
}

