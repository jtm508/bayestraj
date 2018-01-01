#' drawz
#'
#' Draw a component of z, the vectod determining which variables are included in the model.
#'
#' @param z: Matrix, current variable selection matrix
#' @param j: Integer, class flag
#' @param b: Integer, elemtn of z[j,] to draw
#' @param y: Vector, outcomes
#' @param X: Matrix, design matrix
#' @param index: Vector: List of class assignments
#' @param marg.lik.c: Float: Current marginal likelihood
#'
#' @export

drawz = function(z,j,b,y,X,index,marg.lik.c) {
  zp = z[j,]
  zp[b] = 1 - zp[b]
  marg.lik.p = marg_lik(y[index==j], X[index==j,zp==1,drop=FALSE])
  r = (marg.lik.p - marg.lik.c) * (-1)^(zp[b]==0)
  z.b = rbinom(1, 1, 1 / (1 + exp(-r)))
  if (z.b == zp[b]) {marg.lik.c = marg.lik.p}
  return(list(z.b,marg.lik.c))
}