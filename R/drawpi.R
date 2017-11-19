#' drawpi
#'
#' Draw pi from posterior distribution
#'
#' @param g: Vector, group memberships
#' @param alpha: Vector, parameters of Dirichlet prior
#' @param K: Integer, number of groups
#'
#' @importFrom MCMCpack rdirichlet
#'
#' @export

drawpi = function(g, alpha, K) {
  pi = rdirichlet(1, alpha + table(factor(g, levels = 1:K)))
  return(pi)
}