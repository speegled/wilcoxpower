#' Calculates the effect size given two prob distributions on finite set
#'
#' @param p1 A vector of probs.
#' @param p2 A second vector of probs.
#' @return The Vargha and Delaney's A effect size of a population with probability distributions given by p1 and p2.
#' @examples
#' compute_effsize(c(.5, .5), c(.4, .6))

compute_effsize <- function(p1, p2) {
  p1 <- p1/sum(p1)
  p2 <- p2/sum(p2)
  m <- outer(p1, p2)
  s1 <- sum(diag(m)) * 1/2
  m[upper.tri(m, diag = TRUE)] <- 0
  return(sum(m) + s1)
}
