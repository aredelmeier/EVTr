library(EVTr)
library(data.table)

# helpful to plot multiple ggplots on one page
# comes from here: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Comes from package "QRM" but has been updated
xiplot <- function (data, models = 30, start = 15, end = 500, reverse = TRUE,
                    ci = 0.95, auto.scale = TRUE, labels = TRUE, table = FALSE,
                    ...) {
  if (is.timeSeries(data))
    data <- series(data)
  data <- as.numeric(data)
  qq <- 0
  if (ci)
    qq <- qnorm(1 - (1 - ci)/2)
  x <- trunc(seq(from = min(end, length(data)), to = start,
                 length = models))
  gpd.dummy <- function(nex, data) {
    out <- fit.GPD(data = data, nextremes = nex, information = "expected")
    c(out$threshold, out$par.ests[1], out$par.ses[1])
  }
  mat <- apply(as.matrix(x), 1, gpd.dummy, data = data)
  mat <- rbind(mat, x)
  dimnames(mat) <- list(c("threshold", "shape", "se", "exceedances"),
                        NULL)
  thresh <- mat[1, ]
  y <- mat[2, ]
  yrange <- range(y)
  if (ci) {
    u <- y + mat[3, ] * qq
    l <- y - mat[3, ] * qq
    yrange <- range(y, u, l)
  }
  index <- x
  if (reverse)
    index <- -x
  if (auto.scale) {
    plot(x=index, y=y, ylim = yrange, type = "l", xlab = "",
         ylab = "", axes = FALSE, ...)
  }
  else {
    plot(index, y, type = "l", xlab = "", ylab = "", axes = FALSE,
         ...)
  }
  #axis(1, at = index, labels = paste(x), tick = FALSE)
  axis(2)
  axis(1, at = index, labels = paste(format(signif(thresh, 3))), tick = FALSE)
  box()
  if (ci) {
    lines(index, u, lty = 2, col = 2)
    lines(index, l, lty = 2, col = 2)
  }
  if (labels) {
    labely <- "Shape (xi)"
    if (ci)
      labely <- paste(labely, " (CI, p = ", ci, ")", sep = "")
    title(xlab = "Threshold", ylab = labely)
    #mtext("Threshold", side = 3, line = 3)
  }
  if (table)
    print(mat)
  return(invisible(list(x = index, y = y, upper = u, lower = l)))
}



###############################################################
######### Study 1: Delete Censored Obs  #######################
###############################################################

n.sim = 10
n = 10
num_inj = 30
xi = 0.5
censor = c(12, 36)
rate_exp = 1
method = c("MLE", "CensMLE")

results_MLE <- matrix(NA, nrow = n.sim, ncol = length(censor))

results <- list()
for(met in method) {
  for(i in 1:n.sim) {
    row_MLE = NULL
    row_CensMLE = NULL

    for (cens in censor) {
      row_MLE <- c(row_MLE, simulation(censor = cens,
                                       xi = xi,
                                       n = n,
                                       num_inj = num_inj,
                                       rate_exp = rate_exp,
                                       method = met,
                                       specific = "keep_censored_obs",
                                       seed = i))
    }
    results_MLE[i, ] <- row_MLE

  }
  results[[met]] <- data.table(results_MLE)
}

results_MLE_dt <- rbindlist(results, idcol = TRUE)

colnames(results_MLE_dt) <- c("method", as.character(censor))
results_MLE_dt[, n := 1:.N, by = method]

df <- data.table::melt(results_MLE_dt, id.vars = c("method", "n"), variable.name = "censor")
setkeyv(df, "method")


boxplot(I(value - 0.5) ~ censor + method, data = df, main = "Simulation Study Mar 21")
abline(h=0,col="red")

