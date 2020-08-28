## ----setup--------------------------------------------------------------------
library(EVTr)

## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 3
)

## -----------------------------------------------------------------------------
gen_data(censor = 10,
         xi = 1,
         n = 5,
         num_inj = 5,
         rate_exp = 1,
         ne = NULL,
         specific = c("delete_censored_obs"),
         seed = 1)

