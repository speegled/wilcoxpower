% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_probs.R
\name{find_probs}
\alias{find_probs}
\title{Create distributions with desired effect size}
\usage{
find_probs(
  effect_size,
  num_levels,
  p1,
  p2,
  max_ties = FALSE,
  unbalanced = FALSE
)
}
\arguments{
\item{effect_size}{A numeric value between 0 and 1, Vargha and Delaney's A}

\item{num_levels}{The number of levels in the responses.}

\item{p1}{Optional nonnegative values for the distribution of population 1. These will be normalized to sum to 1..}

\item{p2}{Optional nonnegative values for the distribution of population 2. Will be normalized to sum to 1.}
}
\value{
Returns population probabilities associated with the levels to achieve desired effect size. These are not unique.
}
\description{
Create distributions with desired effect size
}
\examples{
find_probs(.6, 3)
find_probs(.6, 3, rep(1, 3))
}
