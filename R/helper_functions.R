
#' Calculate excess above the threshold. Threshold calculated using ne = count of excesses above the threshold.
#'
#' @param data Dataframe. Contains at least the columns "ID", "Injury_Length", and "Censored".
#'
#' @param ne Integer. Count of excesses above the threshold.
#'
#' @return Dataframe. Of observations above the threshold (decided by \code{ne}) and the excess above the threshold.
#'
#' @export
#'
#' @author Annabelle Redelmeier
#'
#' @examples
#' data <- data.frame(ID = c(rep(1, 3), rep(2, 2), rep(3, 4), 5), Injury_Length = rexp(10))
#'
#' excess(data, ne = 1)


excess <- function(data, ne){
  data <- data.table(data)

  thresh <- QRM::findthreshold(data = data$Injury_Length, ne = ne)

  exc <- data[Injury_Length >= thresh]
  exc[, excess := Injury_Length - thresh]
  exc[, Injury_Length_before := Injury_Length]
  exc[, Injury_Length := excess]
  exc[, threshold := thresh]
  exc[, excess := NULL]

  return(exc)
}
