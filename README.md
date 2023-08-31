
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rpsftmPDT

<!-- badges: start -->
<!-- badges: end -->

The goal of rpsftmPDT is to evaluate the treatment effect on overall
survival (OS) in the presence of treatment crossover and other
post-discontinuation therapies (PDTs). It provides an estimated
treatment effect and a variance estimate if needed based on
nonparametric bootstrap.

## Installation

You can install the development version of rpsftmPDT from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("siyi48/rpsftmPDT")
```

## Usage

The two main functions `rpsft.cox` and `rpsft.ipcw` provide two
different approach to adjust for treatment crossover and PDTs and assess
the treatment effect on the OS.

To use the two main functions, users need to identify three periods of
the patients may experience, including the progression-free survival
period, the treatment crossover period, and the PDT period if the
patients progress. A basic example which shows a common usage of the
function `rpsft.cox` is as follows.

``` r
library(rpsftmPDT)
## basic example code for the function `rpsft.cox`
x <- cbind(dat$x1, dat$x2)
a <- dat$a
t.pfs <- dat$t.pfs
t.co <- dat$t.co
t.os <- dat$t.os
cen.time <- dat$c
delta.pfs <- dat$delta.pfs
delta.dp <- dat$delta.dp
delta.pdt <- dat$delta.pdt
delta.co <- dat$delta.co
delta.os <- dat$delta.os
mat.pdt <- x
# if no variance estimate is needed
res.rpsftcox <- rpsft.cox(t.pfs, t.co, t.os, delta.os,
                          cen.time, a, mat.pdt, delta.co, delta.pdt,
                          include.pdt = TRUE, const = 1,
                          tau.lower = -3, tau.upper = 1,
                          grid.length = 100)
res.rpsftcox
#> $tau.est
#> [1] -0.6969697
# if the bootstrap variance estimate is included
res.rpsftcox <- rpsft.cox(t.pfs, t.co, t.os, delta.os,
                          cen.time, a, mat.pdt, delta.co, delta.pdt,
                          include.pdt = TRUE, const = 1,
                          tau.lower = -3, tau.upper = 1,
                          grid.length = 100, var.est = TRUE, B = 100)
res.rpsftcox
#> $tau.est
#> [1] -0.6969697
#> 
#> $var.est
#> [1] 0.0195938
```

Also, we provide a basic example to apply the function `rpsft.ipcw` is
as follows.

``` r
## basic example code for the function `rpsft.ipcw`
x <- cbind(dat$x1, dat$x2)
a <- dat$a
t.pfs <- dat$t.pfs
t.co <- dat$t.co
t.os <- dat$t.os
cen.time <- dat$c
delta.pfs <- dat$delta.pfs
delta.dp <- dat$delta.dp
delta.pdt <- dat$delta.pdt
delta.co <- dat$delta.co
delta.os <- dat$delta.os
mat.pdt <- cbind(x, t.co) # better to include the survival time before the initiation of the PDTs to formulate the PDT-indicator model
# if no variance estimate is needed
res.rpsftipcw <- rpsft.ipcw(t.pfs, t.co, t.os, delta.os, 
                            cen.time, a, mat.pdt, delta.co, delta.pdt,
                            include.pdt = TRUE, const = 1,
                            tau.lower = -3, tau.upper = 1,
                            grid.length = 100)
res.rpsftipcw
#> $tau.est
#> [1] -0.8585859
# if the bootstrap variance estimate is included
res.rpsftipcw <- rpsft.ipcw(t.pfs, t.co, t.os, delta.os,
                            cen.time, a, mat.pdt, delta.co, delta.pdt,
                            include.pdt = TRUE, const = 1,
                            tau.lower = -3, tau.upper = 1,
                            grid.length = 100, var.est = TRUE, B = 100)
res.rpsftipcw
#> $tau.est
#> [1] -0.8585859
#> 
#> $var.est
#> [1] 0.05056157
```
