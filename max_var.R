ss <- replicate(10000, {
  d1 <- sample(1:2, 50, replace = T, prob = c(.25, .75))
  d2 <- sample(1:2, 50, T, c(25, 75))
\
})
mean(ss < .05)

mnk <- function(n, k, x) {
  1/4 * (n - k)^2 * x^2 + (n - k) * (n - k/2) * x * (1 - x) + (n - k/2)^2 * (1 - x)^2
}

sum(mnk(5, 0:5, .3) * dbinom(0:5, 5,prob = 2 - 1.4 - .3)) - 5^2 * .7^2
x <- 0.3
A <- 0.7
n <- 5
sqrt(n^3 * x * (1 - x)/4 + n^2 * (2 - 2*A - x) * (2*A - 1 + x)/4 + n *(n-1)*  (sum(mnk(n, 0:n, x) * dbinom(0:n, n, prob = (2 - 2*A - x))) - n^2 * A^2) )

fn <- function(x, n, A) {
  -1 * sqrt(n^3 * x * (1 - x)/4 + n^2 * (2 - 2*A - x) * (2*A - 1 + x)/4 +
         n *(n-1)*  (sum(mnk(n, 0:n, x) * dbinom(0:n, n, prob = (2 - 2*A - x))) - n^2 * A^2) )
}

fn(0.4, 5, .7)

optim(par = c(.5), fn = fn, lower = 0, upper = 1, n = 7, A = 0.8, method = "Brent")












replicate(10000, {
  d1 <- sample(1:2, 5, T, c(3, 7))
  d2 <- sample(1:2, 5, T, c(7, 3))
  wilcox.test(d1, d2)$statistic
})
ss <- .Last.value
mean(ss)
hist(ss)
sd(ss)

library(tidyverse)
n <- 10
out <- data.frame(sd = 0, effect_size = 0)
num_levels <- 2
for(i in 1:100) {
  # aa <- find_probs(.4, num_levels = 3)
  # p1 <- aa$p1
  # p2 <- aa$p2
  p1 <- runif(num_levels)
  p2 <- runif(num_levels)
  ss <- replicate(4000, {
    d1 <- sample(1:num_levels, n, replace = T, prob = p1)
    d2 <- sample(1:num_levels, n, T, p2)
    r <- rank(c(d1, d2))
    n.x <- as.double(length(d1))
    n.y <- as.double(length(d2))
    statistic <- sum(r[seq_along(d1)]) - n.x * (n.x + 1)/2
    return(statistic)
    # nties <- table(r)
    # z <- statistic - n.x * n.y/2
    # sigma <- sqrt((n.x * n.y/12) * ((n.x + n.y + 1) - sum(nties^3 - nties)/((n.x + n.y) * (n.x + n.y - 1))))
    # z <- (z - sign(z) /2)/sigma
    # return(2 * min(pnorm(z), pnorm(z, lower.tail = FALSE)))
    # suppressWarnings(wilcox.test(d1, d2))$statistic
  })
  # if(sd(ss) < 10.6) {
  #   print(p1)
  #   print(p2)
  #   print(sd(ss))
  # }
  out <- bind_rows(out, list(
    sd = sd(ss),
    effect_size = wilcoxpower:::compute_effsize(p1, p2)))
}

power_wilcox_test(sample_size = 20, p1 = c(.75, .1, .15), p2 = c(.58, .06, .36), num_levels = 3)

p1 <- c(.42, .32, .26)
p2 <- c(.26, .32, .42)
wilcoxpower:::compute_effsize(p1, p2)
sd(ss)
hist(ss)
e1071::skewness(ss)
mean(ss)
qwilcox(.05, 10, 10, lower.tail = F)
sum(ss <= 24)

p1 <- c(.75, .1, .15)
p2 <- c(.58, .06, .36)

out <- rbind(out, c(0, 1))
plot(out)
lines(-pp, seq(0.05, 0.95, length.out = 20))
out
warnings()

sapply(seq(0.05, 0.95, length.out = 20), function(A) {
  optim(par = c(.5), fn = fn, lower = 0, upper = 1, n = n, A = A, method = "Brent")$value
})
pp <- .Last.value
pp[18:20] <- pp[3:1]


optim(par = c(.5), fn = fn, lower = 0, upper = 1, n = n, A = .4, method = "Brent")$value




ui <- diag(6)
ui
ui <- rbind(ui, c(1,1,1,0,0,0))
ui <- rbind(ui, c(-1,-1,-1,0,0,0))
ui <- rbind(ui, c(0,0,0,1,1,1))
ui <- rbind(ui, c(0,0,0,-1,-1,-1))

ui <- matrix(c(1,0,1,-1,0,1,1,-1), ncol = 2)
ci <- c(0,0,0,-1)
ui
fn <- function(par, num_levels) {
  p1vec <- par[1:num_levels]
  p2vec <- par[(num_levels + 1): 2 * num_levels]
  pp <- outer(p1vec, p2vec)
  p2 <- sum(diag(pp))
  p1 <- sum(pp[lower.tri(pp)])
  p3 <- 1 - p1 - p2
  if(p3 > 1 || p3 < 0) {
    stop("error in p3")
  }
  -1 * (p1 + 1/4 * p2 - (p1 + 1/2 * p2)^2)
}
theta <- c(.3, .3)
constrOptim(theta = c(.3,.3), f = fn, ui = ui, ci = ci, method = "Nelder-Mead")
ui %*% theta - ci


hin <- function(x) {
  h <- rep(NA, 6)
  h[1] <- x[1]
  h[2] <- x[2]
  h[3] <- x[3]
  h[4] <- x[4]
  h[5] <- x[5]
  h[6] <- x[6]
  h
}

heq <- function(x) {
  h <- rep(NA, 3)
  h[1] <- x[1] + x[2] + x[3] - 1
  h[2] <- x[4] + x[5] + x[6] - 1
  h[3] <- 1/2 * (x[3] * x[6] + x[2] * x[5] + x[1] * x[4]) + x[3]*(x[5] + x[4]) + x[2] * x[4] - 0.7
  h
}
heq(c(.3, 0, .7, 0, 1, 0))
fn <- function(par, num_levels = 3) {
  p1vec <- par[1:num_levels]
  p2vec <- par[(num_levels + 1): (2 * num_levels)]
  pp <- outer(p1vec, p2vec)
  p2 <- sum(diag(pp))
  p1 <- sum(pp[lower.tri(pp)])
  p3 <- 1 - p1 - p2
  # if(p3 > 1 || p3 < 0) {
  #   stop("error in p3")
  # }
  nties <- p2
  -1 * (p1 + 1/4 * p2 - (p1 + 1/2 * p2)^2)
}
#/sqrt(1 * 2 * 5/24 - sum(nties^3 - nties)/48
sqrt(n * (n + 1) * (2 * n + 1)/24 - sum(NTIES^3 -
                      NTIES)/48)

p1 <- c(.3, 0, .7)
p2 <- c(0, 1, 0)

ss_big_var <- replicate(10000, {
  d1 <- sample(1:3, 30, T, p1)
  d2 <- sample(1:3, 30, T, p2)
  suppressWarnings(wilcox.test(d1, d2))$statistic
})

p1 <- c(.2, .3, .5)

p2 <- c(.5, .3, .2)
ss_smaller_var <- replicate(10000, {
  d1 <- sample(1:3, 30, T, p1)
  d2 <- sample(1:3, 30, T, p2)
  suppressWarnings(wilcox.test(d1, d2))$statistic
})

sort(unique(ss_big_var))
mean(ss_big_var > qwilcox(.025, 30, 30, lower.tail = F))
mean(ss_big_var > 73)
mean(ss_smaller_var > qwilcox(.025, 30, 30, lower.tail = F))
mean(ss_smaller_var > 73)
sd(ss_smaller_var)
sd(ss_big_var)
mean(ss_smaller_var)
mean(ss_big_var)
hist(ss_big_var)
hist(ss_smaller_var)



qwilcox(.025, 10, 10, lower.tail = F)
unique(ss) %>% sort()
hist(ss)
length(unique(ss))
qwilcox(.975, 30, 30)

mean(ss >= 569)
2 * pwilcox(569, 30, 30, lower.tail = F)
var(ss)


mean(ss > 72)
p1 <- c(.3, 0, .7)
p1 <- c(.2, .3, .5)
p2 <- c(0, 1, 0)
p2 <- c(.5, .3, .2)
wilcoxpower:::compute_effsize(c(.3, 0, .7), c(0, 1, 0))
aa <- find_probs(0.7, 3)
aa
library(alabama)
fn(c(aa$p1, aa$p2))
heq(c(aa$p1, aa$p2))
theta <- c(aa$p1, aa$p2)
constrOptim.nl(par = theta, fn = fn, hin = hin, heq = heq)

power_wilcox_test(sample_size = 10, p1 = c(.3, 0.01, .7), p2 = c(0.01, 1, 0.01))
fn(c(.2, .3, .5, .5, .3, .2))

#'
#' OK, .3 0 .7, 0 1 0 has the highest variance of the test statistic, but does not
#' have the smallest power. Must be because of the ties correction. Look into that tomorrow.
#'
#' Nope. It's because smaller variance only leads to higher power when the expected value of
#' the test statistic is greater than the critical value. There is also some weirdness with
#' the granularity of the big variance set, but I don't think that is causing anything.
#'
p1 <- c(.2, .3, .5)
p2 <- c(.5, .3, .2)
p1 <- c(.3, 0, 7)
p2 <- c(0, 1, 0)
dd <- replicate(1000, {
  d1 <- sample(1:3, 2, T, p1)
  d2 <- sample(1:3, 15, T, p2)
  c(as.numeric(wilcox.test(d1[1], d2)$statistic), as.numeric(wilcox.test(d1[2], d2)$statistic))
})
cov(dd[1,], dd[2,])
e1071::skewness(ss_big_var)
library(fGarch)
snormFit(ss_smaller_var)
hist(ss_big_var, probability = T)
curve(dsnorm(x, snormFit(ss_big_var)$par), add = T)



hist(ss_smaller_var, probability = T)
curve(dsnorm(x, mean = 624.8, sd = 57.5, xi = 0.9), add = T)
psnorm(550, 624.8, 57.5, 0.9)
mean(ss_smaller_var < 550)



mean(ss_smaller_var)
#install.packages("fGarch")

d1 <- sample(1:3, 10000, T, prob = p1)
d2 <- sample(1:3, 10000, T, p2)
cov(d1, d2)

p1 <- c(.3, 0, .7)
p2 <- c(.5, 1, 0)
d1 <- sample(1:3, 10000, T, prob = p1)
d2 <- sample(1:3, 10000, T, p2)
cov(d1, d2)



aa <- find_probs(.7, 3)
p1 <- aa$p1
p2 <- aa$p2

ss <- replicate(10000, {
  d1 <- sample(1:3, 30, T, p1)
  d2 <- sample(1:3, 30, T, p2)
  suppressWarnings(wilcox.test(d1, d2))$statistic
})
snormFit(ss)$par

curve(dsnorm(x, xi = 1), from = -2, to = 2)








num_levels <- 3

mcov <- function(p1vec, p2vec, p1, p2) {
  # computes covariance of W statistic when p1vec is the same, but p2vec is different
  # p2 is the probability of a tie
  # p1 is the probability that p1vec is better than p2vec
  if(missing(p1) || missing(p2)) {
    pp <- outer(p1vec, p2vec)
    p2 <- sum(diag(pp))
    p1 <- sum(pp[lower.tri(pp)])
  }
  p2p2 <- outer(p2vec, p2vec)
  probxy_1 <- sum(sapply(nrow(p2p2):1, function(x) {
    p1vec[x] * sum(p2p2[-(x:nrow(p2p2)),-(x:nrow(p2p2))])
  }))
  probxy_12 <- sum(sapply(nrow(p2p2):1, function(x) {
    p1vec[x] * (sum(p2p2[x, -(x:nrow(p2p2))]) + sum(p2p2[-(x:nrow(p2p2)), x]))
  }))
  probxy_14 <- sum(p1vec * diag(p2p2))

  cov <- 1 * probxy_1 + 1/2 * probxy_12 + 1/4 * probxy_14 - (p1 + 1/2 * p2)^2
}


hin <- function(x) {
  h <- rep(NA, 2 * num_levels)
  for(i in 1:(2 * num_levels)) {
    h[i] <- x[i]
  }
  h
}

heq <- function(x) {
  h <- rep(NA, 3)
  p1 <- x[1:num_levels]
  p2 <- x[(num_levels + 1):(2 * num_levels)]
  h[1] <- sum(p1) - 1
  h[2] <- sum(p2) - 1
  m <- outer(p1, p2)
  s1 <- sum(diag(m)) * 1/2
  m[upper.tri(m, diag = TRUE)] <- 0
  h[3] <- sum(m) + s1 - 0.6
  h
}

fn <- function(par, num_levels = 3) {
  # browser()
  p1vec <- par[1:num_levels]
  p2vec <- par[(num_levels + 1): (2 * num_levels)]
  pp <- outer(p1vec, p2vec)
  p2 <- sum(diag(pp))
  p1 <- sum(pp[lower.tri(pp)])
  p3 <- 1 - p1 - p2
  # if(p3 > 1 || p3 < 0) {
  #   stop("error in p3")
  # }
  nties <- p2
  p2p2 <- outer(p2vec, p2vec)
  probxy_1 <- sum(sapply(nrow(p2p2):1, function(x) {
    p1vec[x] * sum(p2p2[-(x:nrow(p2p2)),-(x:nrow(p2p2))])
  }))
  probxy_12 <- sum(sapply(nrow(p2p2):1, function(x) {
    p1vec[x] * (sum(p2p2[x, -(x:nrow(p2p2))]) + sum(p2p2[-(x:nrow(p2p2)), x]))
  }))
  probxy_14 <- sum(p1vec * diag(p2p2))

  cov <- 1 * probxy_1 + 1/2 * probxy_12 + 1/4 * probxy_14 - (p1 + 1/2 * p2)^2

  -1 * (-1 * n^2 * (p1 + 1/4 * p2 - (p1 + 1/2 * p2)^2) - choose(n, 2) * mcov(p1vec, p2vec) - choose(n, 2) * mcov(p2vec, p1vec))
    #- choose(n, 2) * mcov(p2vec, p1vec)
}

fn <- function(par, num_levels = 3) {
  p1vec <- par[1:num_levels]
  p2vec <- par[(num_levels + 1): (2 * num_levels)]
  Fx <- stepfun(1:length(p1vec), c(0, cumsum(p1vec)))
  Fy <- stepfun(1:length(p2vec), c(0, cumsum(p2vec)))
  A <- wilcoxpower:::compute_effsize(p1vec, p2vec)
  1 * (n^2 * (A - A^2 - 1/4 * sum(p1vec * p2vec)) +
    n^2 * (n - 1) * (sum(p1vec * (Fy(0:(num_levels-1)) + p2vec/2)^2) +
                       sum(p2vec * (1 - Fx(1:num_levels)  + p1vec/2) ^2) -  2 * A^2))
}




fn(par = c(pp$par[1:3], pp$par[4:6])) %>% abs() %>% sqrt()
alabama::constrOptim.nl(par = theta, fn = fn, hin = hin, heq = heq)
pp <- .Last.value
pp$par
n <- 60
aa <- wilcoxpower::find_probs(effect_size = 0.6, num_levels = 3)
theta <- c(aa$p1, aa$p2)
fn(c(aa$p1, aa$p2))
aa$p1
wilcoxpower:::compute_effsize(c(.155, .708, .137), c(.74, 0, .26))
p1 = c(1 - A, 0.0001, A)
p2 = c(0.0001, 1, 0.0001)
p1 = c(2 * (1 - A), 1 - 2 *(1 - A), 0.0001)
p2 = c(1, .0001, .0001)
A <- .6
wilcoxpower::power_wilcox_test(sample_size = n,  num_levels = 3, p1 = pp$par[1:3], p2 = pp$par[4:6], one_value = F)
ss <- replicate(10000, {
  d1 <- sample(1:3, 60, T, c(.3, 0, .7))
  d2 <- sample(1:3, 60, T,c(0, 1, 0) )
  wilcox.test(d1, d2)$statistic
})
var(ss)
fn(c(.3, 0, .7, 0, 1, 0))

p1 <- c(.2, .3, .5)
p2 <- c(.5, .4, .1)
Fx <- stepfun(1:length(p), c(0, cumsum(p)))
Fy <- stepfun(1:length(p), c(0, cumsum(q)))
wilcoxpower:::compute_effsize(p1, p2)

p2p2 <- outer(p2, p2)
probxy_1 <- sum(sapply(nrow(p2p2):1, function(x) {
  p1[x] * sum(p2p2[-(x:nrow(p2p2)),-(x:nrow(p2p2))])
}))
probxy_12 <- sum(sapply(nrow(p2p2):1, function(x) {
  p1[x] * (sum(p2p2[x, -(x:nrow(p2p2))]) + sum(p2p2[-(x:nrow(p2p2)), x]))
}))
probxy_14 <- sum(p1 * diag(p2p2))

cov <- 1 * probxy_1 + 1/2 * probxy_12 + 1/4 * probxy_14





p <- c(.3, 0, .7)
q <- c(.5, .3, .2)
A <- wilcoxpower:::compute_effsize(p, q)
Fx <- stepfun(1:length(p), c(0, cumsum(p)))
Fy <- stepfun(1:length(p), c(0, cumsum(q)))

fn(par = c(p, q))
n <- 20
num_levels <- length(p)
n^2 * (sum(p* (Fy(0:(num_levels - 1)) + q/4))  - A^2) + n^2 * (n - 1) * (
  sum(p * (Fy(0:(num_levels-1)) + q/2)^2) + sum(q * (1 - Fx(1:num_levels)  + p/2) ^2) -  2 * A^2
)
sqrt(.Last.value)
ss <- replicate(10000, {
  d1 <- sample(1:3, n, T, p)
  d2 <- sample(1:3, n, T, q)
  wilcox.test(d1, d2)$statistic
})
var(ss)

table(ss)/10000
sum(p * (Fy(0:(num_levels-1)) + q/2)^2) - A^2


sum(q * (Fx(0:(num_levels - 1)) + p/2) ^2) - A^2

ss <- replicate(10000, {
  d1 <- sample(1:3, n, T, p)
  d2 <- sample(1:3, n, T, q)
  d3 <- sample(1:3, n, T, q)
  wilcox.test(d1, d2)$statistic
})








