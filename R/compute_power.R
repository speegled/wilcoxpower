compute_power <- function(sample_size,
                          effect_size,
                          num_levels,
                          p1,
                          p2,
                          power_only = TRUE,
                          conservative = FALSE, ...) {
  if(conservative) {
    browser()
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
      sample(1:num_levels, N, replace = TRUE, prob = probs2$p1)
    }
    dsamp2 <- function(N) {
      sample(1:num_levels, N, replace = TRUE, prob = probs2$p2)
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
            probs = probs,
           balanced_power = est2,
           probs2 = probs2,
           unbalanced_power = est3,
           probs3 = probs3,
           noether_power = est4))
}

