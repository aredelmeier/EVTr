## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 3
)

## ----setup--------------------------------------------------------------------
library(EVTr)

## -----------------------------------------------------------------------------
df <- gen_data(censor = 10,
         xi = 1,
         n = 5,
         num_inj = 5,
         rate_exp = 1,
         ne = NULL,
         specific = c("delete_censored_obs"),
         seed = 1)

head(df, 6)

## -----------------------------------------------------------------------------
# non-censored MLE
mle(data = df, method = "MLE")


## -----------------------------------------------------------------------------
# censored MLE
mle(data = df, method = "CensMLE")


## -----------------------------------------------------------------------------
# non-censored MLE
df_censored <- gen_data(censor = 10,
         xi = 1,
         n = 5,
         num_inj = 5,
         rate_exp = 1,
         ne = NULL,
         specific = c("keep_censored_obs"),
         seed = 1)

mle(data = df_censored, method = "MLE")


## -----------------------------------------------------------------------------
# censored MLE
mle(data = df_censored, method = "CensMLE")


