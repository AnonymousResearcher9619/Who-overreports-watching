library(dplyr)
library(survey)
library(haven)

# =========================================================================
# Weighted frequency table of the 4-category misreporting typology
# (misreport_4k), pooled across all debates/waves (full long-format
# data_pooled, N = 6944 rows), respondent-clustered via RADIOMETER_ID2.

# IMPORTANT: Before running this script, the variable `misreport_4k` must
# already exist in the environment/data. It is created in the script
# "Rscript03 - Table 3 - Transformations and overreporting desriptives".
# Therefore, that script must be run first, or the prepared data containing
# `misreport_4k` must otherwise be loaded into the current R session.
# =========================================================================

get_weighted_misreport_pooled <- function(data, var_name) {
  
  d <- data |>
    dplyr::filter(!is.na(weight)) |>
    dplyr::mutate(
      var = haven::as_factor(.data[[var_name]])
    ) |>
    dplyr::filter(!is.na(var))
  
  n_unweighted <- nrow(d)
  
  des <- survey::svydesign(
    ids = ~RADIOMETER_ID2,
    weights = ~weight,
    data = d
  )
  
  tab <- survey::svymean(~var, des, na.rm = TRUE) * 100
  pct <- setNames(as.numeric(tab), names(tab))
  
  get_pct <- function(label) {
    idx <- grepl(label, names(pct), fixed = TRUE)
    if (!any(idx)) return(NA_real_)
    as.numeric(pct[idx][1])
  }
  
  list(
    n = n_unweighted,
    overreporting         = get_pct("overreporting"),
    validated_watchers    = get_pct("validated watchers"),
    validated_nonwatchers = get_pct("validated nonwatchers"),
    underreporting        = get_pct("underreporting")
  )
}

res_pooled <- get_weighted_misreport_pooled(data_pooled, "misreport_4k")

# ---- Derived rows (same logic/labels as the per-debate table) ----
# row 5 = overreporters / (overreporters + validated nonwatchers)
# row 6 = overreporters / (validated watchers + overreporters)
row5_overrep_among_nonwatchers <- res_pooled$overreporting /
  (res_pooled$overreporting + res_pooled$validated_nonwatchers) * 100

row6_overrep_among_selfreport <- res_pooled$overreporting /
  (res_pooled$validated_watchers + res_pooled$overreporting) * 100

fmt <- function(x, d = 1) formatC(x, digits = d, format = "f")

latex_body <- c(
  "\\begin{table}[H]",
  "\\centering",
  "\\caption{Descriptives -- dependent variable}",
  "\\label{tab:descriptives_dependent}",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{lllc}",
  "\\toprule",
  "\\multicolumn{4}{l}{\\textbf{Dependent variable}} \\\\",
  "Type of respondent & Exposure self-report & Exposure passive measure & Pooled data \\\\",
  "\\midrule",
  paste0(
    "1. validated watchers & Yes & Yes & ",
    fmt(res_pooled$validated_watchers),
    " \\\\"
  ),
  paste0(
    "2. overreporters & Yes & No & ",
    fmt(res_pooled$overreporting),
    " \\\\"
  ),
  paste0(
    "3. underreporters & No & Yes & ",
    fmt(res_pooled$underreporting),
    " \\\\"
  ),
  paste0(
    "4. validated nonwatchers & No & No & ",
    fmt(res_pooled$validated_nonwatchers),
    " \\\\"
  ),
  "\\midrule",
  "Total \\% & & & 100.0 \\\\",
  paste0(
    "Total N & & & ",
    res_pooled$n,
    " \\\\"
  ),
  "\\midrule",
  paste0(
    "5 overreporters among validated non-watchers & ",
    "5=2/(2+4); in \\% & & ",
    fmt(row5_overrep_among_nonwatchers),
    " \\\\"
  ),
  paste0(
    "6 overreporters among self-report watchers & ",
    "6=2/(1+2); in \\% & & ",
    fmt(row6_overrep_among_selfreport),
    " \\\\"
  ),
  "\\bottomrule",
  "\\end{tabular}%",
  "}",
  "\\vspace{0.5em}",
  "\\begin{minipage}{\\textwidth}",
  "\\footnotesize",
  "\\raggedright",
  "Source: pooled CPES 2021--2023 (N = 6944) -- weighted by gender, age, education, social status, media consumption and internet use.",
  "\\end{minipage}",
  "\\end{table}"
)

cat(latex_body, sep = "\n")