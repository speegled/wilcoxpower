#' Create distributions with desired effect size
#'
#' @param effect_size A numeric value between 0 and 1, Vargha and Delaney's A
#' @param num_levels The number of levels in the responses.
#' @param p1 Optional nonnegative values for the distribution of population 1. These will be normalized to sum to 1..
#' @param p2 Optional nonnegative values for the distribution of population 2. Will be normalized to sum to 1.
#' @return Returns population probabilities associated with the levels to achieve desired effect size. These are not unique.
#' @examples
#' find_probs(.6, 3)
#' find_probs(.6, 3, rep(1, 3))

find_probs <- function(effect_size, num_levels, p1, p2) {
  if(missing(p2) && !missing(p1)) {
    if(num_levels != length(p1)) {
      num_levels <- length(p1)
      warning("number of levels reset to length of p1")
    }
    p2 <- runif(1:num_levels, 0, 1)
    dat <- list(num_levels = num_levels,
                      goal = effect_size,
                      p1 = p1)
    to_min <- function(par, dat) {
      (compute_effsize(dat$p1, par[1:dat$num_levels]) - dat$goal)^2
    }
    probs <- constrOptim(p2,
                         f = to_min,
                         ui = diag(num_levels),
                         ci = 0,
                         method = "Nelder-Mead",
                         dat = dat)
    if(probs$value < .0001) {
      p1 <- p1/sum(p1)
      p2 <- probs$par/sum(probs$par)
      return(list(p1 = p1,
                  p2 = p2,
                  effect_size = compute_effsize(p1, probs$par),
                  message = probs$message))
    } else {
      warning("Effect size not possible with provided p1")
    }
   }
  if(missing(num_levels) && missing(p1) && missing(p1)) {
    num_levels <- 5
    warning("Number of levels set to 5")
  }
  if(missing(p1)) {
    p1 <- runif(1:num_levels, 0, 1)
  }
  if(missing(p2)) {
    p2 <- runif(1:num_levels, 0, 1)
  }

  to_min <- function(par, dat) {
    (compute_effsize(par[1:dat$num_levels], par[(dat$num_levels+1):(2 * dat$num_levels)]) - dat$goal)^2
  }

  dat <- data.frame(num_levels = num_levels,
                    goal = effect_size)
  probs <- constrOptim(c(p1, p2),
                       f = to_min,
                       ui = diag(2 * num_levels),
                       ci = 0,
                       method = "Nelder-Mead",
                       dat = dat)
  p1 <- probs$par[1:num_levels]
  p2 <-  probs$par[(num_levels + 1):(2 * num_levels)]
  return(list(p1 = p1/sum(p1),
              p2 = p2/sum(p2),
              effect_size = compute_effsize(p1, p2),
              message = probs$message))
}
