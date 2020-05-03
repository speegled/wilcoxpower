#' Estimates power given a method of sampling from two distributions
#'
#' @param sample_size Number of samples in each population
#' @param dsamp1 Method of sampling from population 1
#' @param dsamp2 Method of sampling from population 2
#' @param alpha Significance level of the hypothesis test
#' @param any_ties Set to FALSE if there are no possible ties in the populations
#' @param num_replicate Number of replications in simulated power. Default 10000
#' @return Estimate of power from wilcox.test
#' @examples
#' estimate_power(sample_size = 40)
#'

estimate_power <- function(sample_size,
                           dsamp1,
                           dsamp2,
                           alpha = .05,
                           any_ties = TRUE,
                           num_replicate = 10000) {
  if(any_ties) {
    s <- replicate(num_replicate, {
      d1 <- dsamp1(sample_size)
      d2 <- dsamp2(sample_size)
      r <- rank(c(d1, d2))
      n.x <- as.double(length(d1))
      n.y <- as.double(length(d2))
      statistic <- sum(r[seq_along(d1)]) - n.x * (n.x + 1)/2
      nties <- table(r)
      z <- statistic - n.x * n.y/2
      sigma <- sqrt((n.x * n.y/12) * ((n.x + n.y + 1) - sum(nties^3 - nties)/((n.x + n.y) * (n.x + n.y - 1))))
      z <- (z - sign(z) /2)/sigma
      return(2 * min(pnorm(z), pnorm(z, lower.tail = FALSE)))
    })
    return(mean(s < alpha))
  } else {
    s <- replicate(num_replicate, {
      d1 <- dsamp1(sample_size)
      d2 <- dsamp2(sample_size)
      r <- rank(c(d1, d2))
      n.x <- as.double(length(d1))
      n.y <- as.double(length(d2))
      statistic <- sum(r[seq_along(d1)]) - n.x * (n.x + 1)/2
      z <- statistic - n.x * n.y/2
      sigma <- sqrt((n.x * n.y/12) * ((n.x + n.y + 1)))
      z <- (z - sign(z) /2)/sigma
      return(2 * min(pnorm(z), pnorm(z, lower.tail = FALSE)))
      # suppressWarnings(wilcox.test(d1, d2))$p.value
    })
    return(mean(s < alpha))
  }
}
