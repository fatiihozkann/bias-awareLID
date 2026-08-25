#' CFA residualization with Bartlett scores
#'
#' Fits a confirmatory factor model (MLR, FIML), extracts Bartlett factor
#' scores, and returns per-item regression residuals on the full score matrix.
#' If \code{model} is \code{NULL}, a one-factor model over all items is used.
#'
#' @param data data frame or matrix of item responses (numeric).
#' @param model lavaan model syntax, or \code{NULL} for one factor.
#' @return list with \code{residuals} (matrix) and \code{fit} (lavaan object).
#' @export
cfa_residualize <- function(data, model = NULL) {
  X <- as.matrix(data); items <- colnames(X)
  if (is.null(model)) model <- paste0("F1 =~ ", paste(items, collapse = " + "))
  fit <- lavaan::cfa(model, data = as.data.frame(X), ordered = FALSE,
                     estimator = "MLR", missing = "fiml")
  fs <- lavaan::lavPredict(fit, method = "Bartlett")
  fsc <- stats::complete.cases(fs)
  Xres <- X
  for (j in seq_along(items)) {
    u <- !is.na(X[, j]) & fsc
    if (sum(u) > 10) {
      Xres[u, j] <- stats::resid(stats::lm(X[u, j] ~ fs[u, ]))
      Xres[!u, j] <- NA
    }
  }
  list(residuals = Xres, fit = fit)
}
