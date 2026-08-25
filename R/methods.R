#' @export
print.biasaware_lid <- function(x, ...) {
  s <- x$settings
  cat("Bias-aware LID analysis\n")
  cat(sprintf("  Groups: %s | B = %d, q = %.2f, tau = %.2f, floor = %.2f, seed = %d\n",
      paste(names(s$groups), as.integer(s$groups), sep = "=", collapse = ", "),
      s$B, s$q, s$tau, s$floor, s$seed))
  cat(sprintf("  BH-significant: %d | qualified E0: %d | beta* = %.2f | final: %d\n",
      x$n_bh, x$n_E0, x$beta_star, x$n_final))
  if (nrow(x$removed))
    cat("  Penalty-removed:", paste(x$removed$item_i, x$removed$item_j,
        sep = "-", collapse = ", "), "\n")
  cat(sprintf("  Sign reversals: E0 %d, final %d\n",
      x$reversals_E0, x$reversals_final))
  invisible(x)
}

#' @export
summary.biasaware_lid <- function(object, ...) {
  e <- object$edges[object$edges$final, ]
  e <- e[order(-e$wTO), ]
  cat("Final bias-aware edges (", nrow(e), "):\n", sep = "")
  print(utils::head(e[, c("item_i", "item_j", "wTO", "group_disparity",
                          "sign_reversal")], 25), row.names = FALSE)
  invisible(e)
}
