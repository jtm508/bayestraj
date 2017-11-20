#' logsumexp
#'
#' Compute log(sum(exp(X))) while avoiding numerical underflow
#'
#' @param X: Vector

logsumexp = function(X){
  a = max(X)
  y = a + log(sum(exp(X-a)))
  return(y)
}