#' @title Determing the number of factors in approximate factor models.
#' @description The function is adapted by Bai & Ng method to minimize information criterion with two penalties in static factor models.
#' @param Y p by n matrix of raw data, where p is the dimensionality, n is the sample size.
#' @return estimated number of factors based on two infomation criterion using Bai & Ng method.
#' @examples 
#' \dontrun{
#' p=100
#' n=100
#' Y<-array(rnorm(p*n),dim=c(p,n))
#' K<-FM(Y)
#' }
#' @export
FM<-function (Y) 
{
  p = nrow(Y)
  n = ncol(Y)
  Y <- Y - t(t(apply(t(Y), 2, mean))) %*% matrix(1, 1, n)
  rmax = 10
  IC = array(0, c(2, rmax))
  frob = array(0, c(rmax))
  for (k in 1:rmax) {
    V <- eigen(t(Y) %*% Y)$vectors
    V = as.matrix(V)
    Dd <- eigen(t(Y) %*% Y)$values
    Dd = as.vector(Dd)
    #eigenvectors corresponding to the K largest eigenvalues of Y'Y
    F <- V[, 1:k]
    LamPCA = Y %*% F/n
    uhat = Y - LamPCA %*% t(F)
    #two penalties
    frob[k] = sum(diag(uhat %*% t(uhat)))/(p * n)
    first = log((p * n)/(p + n)) * (p + n)/(p * n)
    second = log(min(p, n)) * (p + n)/(p * n)
    IC[1, k] = log(frob[k]) + k * first
    IC[2, k] = log(frob[k]) + k * second
  }
  K1 = which.min(IC[1, ])
  K2 = which.min(IC[2, ])
  result <- list(K1 = K1, K2 = K2, IC = IC)
  return(result)
}

#' @title Determining the Number of Factors in the general dynamic factor model.
#' @description The function is used to minimize information criterion with two penalties in the general dynamic factor model. Comparing to the static factor model, parameters of information 
#' criterion in this model is not fixed, changing the fixed p,n and c to vectors.
#' @param Y p by n matrix of raw data, where p is the dimensionality, n is the sample size.
#' @return estimated number of factors based on two information criterion in the General Dynamic Factor Model.
#' @importFrom stats sd
#' @examples 
#' \dontrun{
#' p=100
#' n=100
#' Y<-array(rnorm(p*n),dim=c(p,n))
#' K<-DFM(Y)
#' }
#' @export
DFM<-function (Y) 
{
  p = nrow(Y)
  n = ncol(Y)
  Y <- Y - t(t(apply(t(Y), 2, mean))) %*% matrix(1, 1, n)
  c = seq(0.05, 5, length = 100)
  re = 20
  rmax = 10
  IC = array(0, c(2, re, rmax, 100))
  first = array(1, c(re))
  second = array(1, c(re))
  pi = array(1, c(re))
  ni = array(1, c(re))
  for (i in 1:re) {
    pi[i] = min(i * floor(p/re) + min(p, 5), p)
    ni[i] = min(i * floor(n/re) + min(n, 5), n)
    if (i == re) {
      pi[i] = p
      ni[i] = n
    }
    Yi = Y[1:pi[i], 1:ni[i]]
    frob = array(0, c(rmax))
    for (k in 1:min(pi[i], ni[i], rmax)) {
      V <- eigen(t(Yi) %*% Yi)$vectors
      V = as.matrix(V)
      Dd <- eigen(t(Yi) %*% Yi)$values
      Dd = as.vector(Dd)
      F <- V[, 1:k]
      LamPCA = Yi %*% F/ni[i]
      uhat = Yi - LamPCA %*% t(F)
      frob[k] = sum(diag(uhat %*% t(uhat)))/(pi[i] * ni[i])
      first[i] = log((pi[i] * ni[i])/(pi[i] + ni[i])) * 
        (pi[i] + ni[i])/(pi[i] * ni[i])
      second[i] = log(min(pi[i], ni[i])) * (pi[i] + ni[i])/(pi[i] * 
                                                             ni[i])
      for (l in 1:100) {
        IC[1, i, k, l] = log(frob[k]) + c[l] * k * first[i]
        IC[2, i, k, l] = log(frob[k]) + c[l] * k * second[i]
      }
    }
  }
  rhat = array(0, c(2, re, 100))
  for (i in 1:re) {
    for (l in 1:100) {
      m = min(pi[i], ni[i], rmax)
      temp1 = which.min(IC[1, i, 1:m, l])
      rhat[1, i, l] = temp1
      temp2 = which.min(IC[2, i, 1:m, l])
      rhat[2, i, l] = temp2
    }
  }
  Sc1 = array(0, c(100))
  Sc2 = array(0, c(100))
  for (l in 1:100) {
    # i using sd, with l fixed
    Sc1[l] = sd(rhat[1, , l])
    Sc2[l] = sd(rhat[2, , l])
  }
  c1vec = which(Sc1 == 0)
  ctemp1 = c1vec[1]
  K1 = rhat[1, 1, ctemp1]
  
  c2vec = which(Sc2 == 0)
  ctemp2 = c2vec[1]
  K2 = rhat[2, 1, ctemp2]
  
  result <- list(K1 = K1, K2 = K2)
  return(result)
}
