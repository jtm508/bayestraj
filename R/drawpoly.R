#' drawpoly
#'
#' Draw the polynomial degree for time
#'
#' @param y: Vector, outcomes
#' @param X: Matrix, design matrix
#' @param time_index: Integer, column of X corresponding to time
#' @param z: Vector, current variable selection matrix
#'
#' @export

drawpoly = function(y,X,time_index,z) {
  D = dim(X)[2] - time_index + 1
  marg.lik = rep(0,D+1)
  z[time_index:length(z)] = 0
  marg.lik[1] = marg_lik(y, X[,z %in% c(-1,1),drop=FALSE])
  for (d in 1:D) {
    z[time_index+d-1] = 1
    marg.lik[d+1] = marg_lik(y, X[,z %in% c(-1,1),drop=FALSE])
  }
  p = logp(marg.lik)
  d_new = sample(0:D,1,prob=p)
  z_new = rep(0,D)
  z_new[1:d_new] = 1
  return(z_new)
}

#funcion to calcuulate p_i/sum(p) i=1..k when we are given lop(p)
logp = function(logx){
  a = max(logx)
  logx.a  = logx - a
  exp.a = exp(logx.a)
  denom = sum(exp.a)
  p = exp.a / denom
  return(p)
}