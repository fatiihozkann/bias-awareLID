suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
})

prepare_timss_bootstrap_data <- function(
  data_file = "~/Downloads/bsgusam7.sav",
  sampling_seed = 2025L
) {
  items <- paste0(
    "BSBM19",
    LETTERS[1:9]
  )

  data <- read_sav(data_file) |>
    select(
      all_of(items),
      ITSEX,
      BSBG04
    ) |>
    rename(
      Gender = ITSEX,
      Books = BSBG04
    )

  for (variable in c(items, "Gender", "Books")) {
    data[[variable]][data[[variable]] < 0] <- NA
  }

  data <- data |>
    filter(
      !is.na(Gender),
      !is.na(Books)
    ) |>
    mutate(
      SES = case_when(
        Books <= 2 ~ "Low",
        Books <= 3 ~ "Middle",
        TRUE ~ "High"
      ),
      Gender_label = ifelse(
        Gender == 1,
        "Female",
        "Male"
      ),
      Group = paste(
        Gender_label,
        SES,
        sep = "_"
      )
    ) |>
    filter(
      if_all(
        all_of(items),
        ~ !is.na(.)
      )
    )

  set.seed(sampling_seed)

  analysis_data <- data |>
    group_by(Group) |>
    slice_sample(
      n = 800,
      replace = FALSE
    ) |>
    ungroup()

  item_data <- analysis_data |>
    select(all_of(items)) |>
    mutate(
      across(
        everything(),
        as.numeric
      )
    )

  analysis_data[, items] <- item_data

  list(
    data = analysis_data,
    items = items,
    group_variable = "Group",
    groups = factor(analysis_data$Group),
    sampling_seed = sampling_seed
  )
}
