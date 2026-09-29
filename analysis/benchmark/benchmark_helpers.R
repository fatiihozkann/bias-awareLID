canonical_pair <- function(item_i, item_j) {
  paste(sort(c(item_i, item_j)), collapse = "--")
}

pair_matrix <- function(pairs, p) {
  result <- matrix(FALSE, p, p)

  for (pair in pairs) {
    result[pair[1], pair[2]] <- TRUE
    result[pair[2], pair[1]] <- TRUE
  }

  diag(result) <- FALSE
  result
}

signed_q3_star <- function(residuals) {
  correlation <- make_psd(
    safe_cor(residuals, method = "pearson"),
    corr = TRUE
  )

  off_diagonal_mean <- mean(
    correlation[upper.tri(correlation, diag = FALSE)]
  )

  signed <- correlation - off_diagonal_mean
  diag(signed) <- 0
  signed
}

unsigned_adjacency <- function(signed_q3) {
  adjacency <- abs(signed_q3)
  diag(adjacency) <- 0
  adjacency
}

precision_recall_f1 <- function(detected, truth) {
  detected <- unique(detected)
  truth <- unique(truth)

  true_positive <- length(intersect(detected, truth))
  false_positive <- length(setdiff(detected, truth))
  false_negative <- length(setdiff(truth, detected))

  precision <- if (true_positive + false_positive > 0) {
    true_positive / (true_positive + false_positive)
  } else {
    NA_real_
  }

  recall <- if (true_positive + false_negative > 0) {
    true_positive / (true_positive + false_negative)
  } else {
    NA_real_
  }

  f1 <- if (
    is.finite(precision) &&
    is.finite(recall) &&
    precision + recall > 0
  ) {
    2 * precision * recall / (precision + recall)
  } else {
    NA_real_
  }

  false_discovery_rate <- if (
    true_positive + false_positive > 0
  ) {
    false_positive / (true_positive + false_positive)
  } else {
    0
  }

  data.frame(
    TP = true_positive,
    FP = false_positive,
    FN = false_negative,
    precision = precision,
    recall = recall,
    F1 = f1,
    FDR = false_discovery_rate
  )
}

matrix_detected_pairs <- function(selection, items) {
  indices <- which(
    selection & upper.tri(selection, diag = FALSE),
    arr.ind = TRUE
  )

  if (nrow(indices) == 0L) {
    return(character(0))
  }

  apply(
    indices,
    1,
    function(index) {
      canonical_pair(
        items[index[1]],
        items[index[2]]
      )
    }
  )
}

pair_indices_to_names <- function(pairs, items) {
  vapply(
    pairs,
    function(pair) {
      canonical_pair(
        items[pair[1]],
        items[pair[2]]
      )
    },
    character(1)
  )
}

sign_reversal_matrix <- function(signed_group_matrices) {
  p <- nrow(signed_group_matrices[[1]])
  result <- matrix(FALSE, p, p)

  for (i in seq_len(p - 1L)) {
    for (j in seq.int(i + 1L, p)) {
      values <- vapply(
        signed_group_matrices,
        function(group_matrix) group_matrix[i, j],
        numeric(1)
      )

      result[i, j] <- min(values) < 0 && max(values) > 0
      result[j, i] <- result[i, j]
    }
  }

  diag(result) <- FALSE
  result
}

classify_bias_aware_pairs <- function(
  E0,
  final_selection,
  items
) {
  removed <- E0 & !final_selection

  list(
    qualified = matrix_detected_pairs(E0, items),
    retained = matrix_detected_pairs(final_selection, items),
    group_specific = matrix_detected_pairs(removed, items)
  )
}
