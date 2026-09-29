suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
})

prepare_pisa_bootstrap_data <- function(
  data_file = "~/Downloads/CY08MSP_STU_QQQ.SAV"
) {
  items <- paste0(
    "ST267Q0",
    seq_len(8),
    "JA"
  )

  data <- read_sav(data_file)

  plausible_values <- grep(
    "^PV[0-9]+MATH$",
    names(data),
    value = TRUE
  )

  if (length(plausible_values) == 0L) {
    stop("No mathematics plausible-value columns found")
  }

  plausible_value_data <- as.data.frame(
    lapply(
      data[, plausible_values],
      as.numeric
    )
  )

  data$MATH_MEAN <- rowMeans(
    plausible_value_data,
    na.rm = TRUE
  )

  gifted_cutoff <- quantile(
    data$MATH_MEAN,
    0.90,
    na.rm = TRUE
  )

  gifted_data <- data[
    data$MATH_MEAN >= gifted_cutoff,
    ,
    drop = FALSE
  ]

  analysis_data <- gifted_data |>
    select(
      all_of(items),
      ST004D01T,
      ESCS
    ) |>
    rename(
      Gender = ST004D01T
    )

  for (item in items) {
    analysis_data[[item]] <- as.numeric(
      analysis_data[[item]]
    )

    analysis_data[[item]][
      analysis_data[[item]] < 0 |
      analysis_data[[item]] > 90
    ] <- NA
  }

  analysis_data$Gender <- as.numeric(
    analysis_data$Gender
  )

  analysis_data$ESCS <- as.numeric(
    analysis_data$ESCS
  )

  analysis_data$Gender[
    analysis_data$Gender < 0
  ] <- NA

  analysis_data$ESCS[
    analysis_data$ESCS < -90
  ] <- NA

  analysis_data <- analysis_data |>
    filter(
      !is.na(Gender),
      !is.na(ESCS)
    ) |>
    filter(
      if_all(
        all_of(items),
        ~ !is.na(.)
      )
    ) |>
    mutate(
      SES = ntile(
        ESCS,
        3
      ),
      SES_label = case_when(
        SES == 1 ~ "Low",
        SES == 2 ~ "Middle",
        TRUE ~ "High"
      ),
      Gender_label = ifelse(
        Gender == 1,
        "Female",
        "Male"
      ),
      Group = paste(
        Gender_label,
        SES_label,
        sep = "_"
      )
    )

  analysis_data$Group <- factor(
    analysis_data$Group
  )

  list(
    data = analysis_data,
    items = items,
    group_variable = "Group",
    groups = analysis_data$Group,
    gifted_cutoff = gifted_cutoff,
    plausible_values = plausible_values
  )
}
