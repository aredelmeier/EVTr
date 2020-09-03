
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EVTr

<!-- badges: start -->

[![R build
status](https://github.com/aredelmeier/EVTr/workflows/R/badge.svg)](https://github.com/aredelmeier/EVTr/actions?query=workflow%3AR-CMD-check)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

Injured athletes or performers can cost teams, sports associations, and
coaches a lot of money. Therefore, it can be useful to predict the
number of seriously injured athletes and the length of these injuries at
the start of a season or year. It turns out that individuals (especially
professional athletes) are injured all the time and often get better
within the day or the week. However, we are interested in *long term
injuries* that last a season or many years since these are the most
costly. The problem is that these types of injuries happen the least.
Similar to massive earthquakes and stock market crashes, how do we model
these types of scarce and infrequent injuries?

We use extreme value theory - and specifically the peaks-over-threshold
method - to model this type of process. To begin, we do various
simulation studies to study how various distributional assumptions,
censoring, and data manipulations can affect estimation. We assume that
individuals live in two states - they can either be injured or healthy.
Individuals start healthy (or in a healthy state) and then remain
healthy for a while until they switch to being injured (or in an injured
state). Individuals continue alternating between these two states until
our study is over (or *censoring* occurs). Then we use extreme value
theory to study these injury periods. Specifically, we ask

  - Does censoring matter? Can we assume the censored injury was not
    censored and use MLE normally?
  - Will the censored injury be the longest injury and therefore should
    censoring be taken into account?
  - What is the proportion of injuries being censored as appose to
    healthy periods? Will the distribution of the injuries/healthy
    periods play a roll in this proportion?

This package accompanies a [masters
thesis](https://escholarship.mcgill.ca/concern/theses/fn107132j) called
“Using extreme value theory to model severe injuries of circus
artists”.

## Installation

The development version from [GitHub](https://github.com/) can be
installed with:

``` r
# install.packages("devtools")
devtools::install_github("aredelmeier/EVTr")
```

To install the package and *also* download the vignette use:

``` r
# install.packages("devtools")
devtools::install_github("aredelmeier/EVTr", build_vignettes = TRUE)
```

Then, you can access the vignette using:

``` r
vignette("understanding_EVTr",package = "EVTr")
```

Once the package is loaded, find help on the most important functions in
the package using:

``` r
library(EVTr) # load package
help(EVTr)
```

To find specific help on one of the functions use (with correct name of
function):

``` r
help(mle)
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(EVTr)

df <- gen_data(censor = 10,
         xi = 1,
         n = 5,
         num_inj = 5,
         rate_exp = 1,
         ne = NULL,
         specific = c("delete_censored_obs"),
         seed = 1)
```

Then you can easily calculate the non-censored genearlized Pareto
distribution MLE with:

``` r
# non-censored MLE
mle(data = df, method = "MLE")
#> [1] 0.21929
```

Or the censored MLE with:

``` r
# censored MLE
mle(data = df, method = "CensMLE")
#> [1] 0.21929
```
