# wilcoxpower
Power computation for `wilcox.test` in R

This R package will help with power computations for one and two sample tests using `wilcox.test` in R for **finite** Likert-like populations.
Currently only the two-sample helper functions are implemented. 
The functions provide a **relatively conservative** estimate of power given the sample size, number of levels in the responses, and the anticipated effect size.
Exact power computations are not possible unless the population distribution is known, which is typically not the case.
