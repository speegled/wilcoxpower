#' Computes power for wilcox.test for ordinal responses with a small number of levels.
#'
#' This function will estimate the power for wilcox.test for several different scenarios. It
#' will estimate power for the probability distributions of the populations provided
#' by the user, for a balanced and an unbalanced probability distribution of the
#' population, and it will also compute the Noether power estimate.
#'
#' If p1 and p2 are not provided, then the function will choose semi-random values
#' for p1 and p2 that have the desired effect size. If p1 and p2 are provided without
#' an effect size, then the effect size is taken to be that implied by p1 and p2. If p1, p2
#' and an effect size are all given, but are not consistent, then p1 and p2 are modified to
#' have the given effect size.
#'
#' Other arguments can be included that will be passed on include num_replicate, which has
#' default of 10000 and alpha, which has a default of .05.
#'
#' @param sample_size Number of samples in each population
#' @param effect_size Vargha and Delaney's A: P(X > Y) + 1/2 P(X == Y)
#' @param num_levels The number of levels in the response
#' @param p1 Optional probability distribution in population 1. Should be length num_levels. Defaults to random.
#' @param p2 Optional probability distribution in population 2. Should be length num_levels. Defaults to random.
#' @param power_only Return only the powers
#' @return Estimate of power from wilcox.test in four scenarios: the p1 and p2 given as arguments, unbalanced p1 and p2, balanced p1 and p2, and a noether estimate.
#'
#' If `power_only` is FALSE, then it returns p1, p2 the effect size, and any messages
#' from the constrained optimization problem of finding p1 and p2 for each of the first
#' three power estimates.
#'
#'
#' @examples
#' power_wilcox_test(sample_size = 40)
#' power_wilcox_test(sample_size = 30, p1 = c(1, .2), p2 = c(1, .01))
#' @export

power_wilcox_test <- function(sample_size,
                          effect_size,
                          num_levels,
                          p1,
                          p2,
                          power_only = TRUE,
                          ...) {
  if(missing(effect_size) && !missing(p1) && !missing(p2)) {
    effect_size <- compute_effsize(p1, p2)
  }
  if(missing(num_levels) && !missing(p1) && missing(p2)) {
    num_levels <- length(p1)
  }
  if(missing(num_levels) && !missing(p1) && !missing(p2)) {
    if(length(p1) == length(p2)) {
      num_levels <- length(p1)
    } else{
      stop("length of p1 and p2 must be equal")
    }
  }
  probs <- find_probs(effect_size, num_levels, p1, p2)
    levs <- length(probs$p1)
    dsamp1 <- function(N) {
      sample(1:levs, N, replace = TRUE, prob = probs$p1)
    }
    dsamp2 <- function(N) {
      sample(1:levs, N, replace = TRUE, prob = probs$p2)
    }
    est <- estimate_power(sample_size = sample_size,
                   dsamp1 = dsamp1,
                   dsamp2 = dsamp2,
                   any_ties = TRUE,
                   ...)
    ret <- probs
    probs2 <- find_probs(effect_size, num_levels, p1, p2, max_ties = TRUE)
    dsamp1 <- function(N) {
      sample(1:num_levels, N, replace = TRUE, prob = probs2$p1)
    }
    dsamp2 <- function(N) {
      sample(1:num_levels, N, replace = TRUE, prob = probs2$p2)
    }
    est2 <- estimate_power(sample_size = sample_size,
                           dsamp1 = dsamp1,
                           dsamp2 = dsamp2,
                           any_ties = TRUE,
                           ...)
    probs3 <- find_probs(effect_size, num_levels, p1, p2, unbalanced = TRUE)
    dsamp1 <- function(N) {
      sample(1:num_levels, N, replace = TRUE, prob = probs3$p1)
    }
    dsamp2 <- function(N) {
      sample(1:num_levels, N, replace = TRUE, prob = probs3$p2)
    }
    est3 <- estimate_power(sample_size = sample_size,
                           dsamp1 = dsamp1,
                           dsamp2 = dsamp2,
                           any_ties = TRUE,
                           ...)
    est4 <- pnorm(-1 * sqrt(sample_size * 2 * 3 * (effect_size - .5)^2) - qnorm(.05),
                  lower.tail = FALSE)
    if(power_only) {
      return(list(estimated_power = est,
                  balanced_power = est2,
                  unbalanced_power = est3,
                  noether_power = est4))
    }
    return(list(estimated_power = est,
           estimated_power_probs = probs,
           balanced_power = est2,
           balanced_power_probs = probs2,
           unbalanced_power = est3,
           unbalanced_power_probs = probs3,
           noether_power = est4))
}

