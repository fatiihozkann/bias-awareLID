# Targeted-benchmark data generator
# Design:
#   N = 900
#   p = 10
#   3 groups of 300
#   X1-X2 and X3-X4: consistent rho = .30
#   X5-X6 and X7-X8: group-specific rho = .40
#   X9 and X10: no planted residual dependence

generate_targeted_benchmark_data <- function(
  dtype = c("continuous", "polytomous", "dichotomous"),
  seed = NULL,
  affected_group = NULL
) {
  dtype <- match.arg(dtype)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  N <- 900L
  p <- 10L
  n_groups <- 3L
  group_sizes <- rep(300L, n_groups)
  items <- paste0("X", seq_len(p))

  if (is.null(affected_group)) {
    affected_group <- sample.int(n_groups, 1L)
  }

  if (!affected_group %in% seq_len(n_groups)) {
    stop("affected_group must be 1, 2, or 3")
  }

  loadings <- rep(0.70, p)

  consistent_pairs <- list(
    c(1L, 2L),
    c(3L, 4L)
  )

  group_specific_pairs <- list(
    c(5L, 6L),
    c(7L, 8L)
  )

  planted_pairs <- c(
    consistent_pairs,
    group_specific_pairs
  )

  group_data <- vector("list", n_groups)

  for (group_index in seq_len(n_groups)) {
    residual_covariance <- diag(1, p)

    for (pair in consistent_pairs) {
      residual_covariance[pair[1], pair[2]] <- 0.30
      residual_covariance[pair[2], pair[1]] <- 0.30
    }

    if (group_index == affected_group) {
      for (pair in group_specific_pairs) {
        residual_covariance[pair[1], pair[2]] <- 0.40
        residual_covariance[pair[2], pair[1]] <- 0.40
      }
    }

    residual_covariance <- make_psd(
      residual_covariance,
      corr = TRUE
    )

    latent_trait <- rnorm(group_sizes[group_index])

    residuals <- MASS::mvrnorm(
      n = group_sizes[group_index],
      mu = rep(0, p),
      Sigma = residual_covariance
    )

    responses <- outer(
      latent_trait,
      loadings
    ) + residuals

    responses <- as.data.frame(responses)
    names(responses) <- items
    responses$group <- factor(group_index)

    group_data[[group_index]] <- responses
  }

  data <- do.call(rbind, group_data)
  rownames(data) <- NULL

  if (dtype == "polytomous") {
    for (item in items) {
      data[[item]] <- to_polytomous5(data[[item]])
    }
  }

  if (dtype == "dichotomous") {
    for (item in items) {
      data[[item]] <- to_dichotomous(data[[item]])
    }
  }

  truth_long <- do.call(
    rbind,
    lapply(seq_len(n_groups), function(group_index) {
      do.call(
        rbind,
        lapply(seq_along(planted_pairs), function(pair_index) {
          pair <- planted_pairs[[pair_index]]
          truth_type <- if (pair_index <= 2L) {
            "consistent"
          } else {
            "group_specific"
          }

          rho_planted <- if (truth_type == "consistent") {
            0.30
          } else if (group_index == affected_group) {
            0.40
          } else {
            0.00
          }

          data.frame(
            item_i = items[pair[1]],
            item_j = items[pair[2]],
            truth_type = truth_type,
            group = as.character(group_index),
            rho_planted = rho_planted,
            affected_group = affected_group,
            stringsAsFactors = FALSE
          )
        })
      )
    })
  )

  list(
    data = data,
    items = items,
    affected_group = affected_group,
    consistent_pairs = consistent_pairs,
    group_specific_pairs = group_specific_pairs,
    planted_pairs = planted_pairs,
    truth_long = truth_long,
    design = data.frame(
      N = N,
      p = p,
      n_groups = n_groups,
      dtype = dtype,
      affected_group = affected_group,
      stringsAsFactors = FALSE
    )
  )
}
