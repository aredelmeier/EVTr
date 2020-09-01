
#' Calculate excess above the threshold. Threshold calculated using ne = count of excesses above the threshold.
#'
#' @param data Dataframe. Contains at least the columns "ID", "Injury_Length", and "Censored".
#'
#' @param ne Integer. Count of excesses above the threshold.
#'
#' @param theshold Integer. If \code{ne} is not supplied, \code{threshold} can be used to specify the threshold
#' to use.
#'
#' @return Dataframe. Of observations above the threshold (decided by \code{ne}) and the excess above the threshold.
#'
#' @export
#'
#' @import QRM
#'
#' @author Annabelle Redelmeier
#'
#' @examples
#' data <- data.frame(ID = c(rep(1, 3), rep(2, 2), rep(3, 4), 5), Injury_Length = rexp(10))
#'
#' excess(data, ne = 1)


excess <- function(data, ne = NULL, threshold = NULL) {

  Injury_Length <- Injury_Length_before <- threshold <- NULL

  if (!is.data.frame(data)) {
    stop("data must be a data.frame.")
  }
  data <- data.table(data)

  if (!("Injury_Length" %in% names(data))) {
    stop("data must include a column called 'Injury_Length' with the data.")
  }

  if (!is.null(ne)) {
    thresh <- QRM::findthreshold(data = data$Injury_Length, ne = ne)

    if (length(thresh) > 1) {
      stop("Please include a larger threshold so that not all data is below threshold.")
    }

  } else if (!is.null(threshold)) {
    thresh <- threshold
  }
  else {
    stop("Please include either 'ne' for the number above the threshold or 'threshold' for the
         actual threshold to be used.")
  }

  exc <- data[Injury_Length >= thresh]
  exc[, excess := Injury_Length - thresh]
  exc[, Injury_Length_before := Injury_Length]
  exc[, Injury_Length := excess]
  exc[, threshold := thresh]
  exc[, excess := NULL]

  return(exc)
}
