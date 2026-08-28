# =========================================================================
# MISREPORTING BY DEBATE: ALL RESPONDENTS VS. EXCL. SPEEDERS & STRAIGHTLINERS
#
# For each of the six debates, computes the weighted percentage of each
# misreport_4k category (overreporting, validated watchers, validated
# nonwatchers, underreporting), plus the two derived "overreporters
# among X" rates, both:
#   (a) "All"  - full sample
#   (b) "EXCL" - excluding speeders (speeding == "x < 2.decile") and
#                straightliners (straightlining == "1")
#
# DATA STRUCTURE: data_pooled is long format, one row per respondent
# (RADIOMETER_ID2) per debate they were asked about. Since each debate's
# estimate is computed on the subset with debate == d, and a respondent
# cannot be asked about the same debate twice, that subset already has
# exactly one row per respondent - so svydesign(ids = ~1, ...) is used
# (no clustering needed), consistent with the mean-plot, weighted-
# Spearman, and per-debate misreport-frequency scripts.
#
# 'misreport_4k' must already exist in data_pooled (created by the
# misreport-typology script Rscript03). Weighted by 'weight' (survey::svymean, no
# additional CI/SE computation - the target table only needs point
# estimates).
# =========================================================================

library(dplyr)
library(survey)
library(haven)
library(purrr)
library(tibble)

# =========================================================================
# 1. DEBATE METADATA
# =========================================================================

debate_key <- c("2021-CT", "2021-NOVA", "2023-NOVA1",
                "2023-CT", "2023-PRIMA", "2023-NOVA2")

debate_labels_display <- c(
  "2021-CT" = "2021-CT", "2021-NOVA" = "2021-NOVA", "2023-NOVA1" = "2023-NOVA",
  "2023-CT" = "2023-CT",  "2023-PRIMA" = "2023-PRIMA", "2023-NOVA2" = "2023-NOVA2"
)

# =========================================================================
# 2. HELPER: WEIGHTED MISREPORT TYPOLOGY FOR ONE DEBATE (ALL OR EXCL.)
# =========================================================================

compute_misreport_stats <- function(data, debate_id, exclude_speeders_straightliners = FALSE) {
  
  d <- data |>
    filter(debate == debate_id, !is.na(weight), !is.na(misreport_4k))
  
  if (exclude_speeders_straightliners) {
    d <- d |> filter(!(trimws(as.character(speeding)) == "x < 2.decile" |
                         trimws(as.character(straightlining)) == "1"))
  }
  
  d <- d |> mutate(var = haven::as_factor(misreport_4k))
  
  n_unweighted <- nrow(d)
  
  des <- survey::svydesign(
    ids     = ~1,
    weights = ~weight,
    data    = d
  )
  
  tab <- survey::svymean(~var, des, na.rm = TRUE) * 100
  pct <- setNames(as.numeric(tab), names(tab))
  
  get_pct <- function(label) {
    idx <- grepl(label, names(pct), fixed = TRUE)
    if (!any(idx)) return(NA_real_)
    as.numeric(pct[idx][1])
  }
  
  tibble::tibble(
    n                      = n_unweighted,
    overreporting          = get_pct("overreporting"),
    validated_watchers     = get_pct("validated watchers"),
    validated_nonwatchers  = get_pct("validated nonwatchers"),
    underreporting         = get_pct("underreporting")
  )
}

cat("=== Diagnostic: speeding / straightlining raw values (check for stray whitespace) ===\n")
cat("Unique speeding values (data_pooled):\n")
print(unique(data_pooled$speeding))
cat("Unique straightlining values (data_pooled):\n")
print(unique(data_pooled$straightlining))

# =========================================================================
# 3. COMPUTE FOR ALL SIX DEBATES, BOTH SAMPLES
# =========================================================================

results_all  <- purrr::map_dfr(debate_key, function(d) {
  compute_misreport_stats(data_pooled, d, exclude_speeders_straightliners = FALSE) |>
    mutate(debate = d, sample = "All")
})

results_excl <- purrr::map_dfr(debate_key, function(d) {
  compute_misreport_stats(data_pooled, d, exclude_speeders_straightliners = TRUE) |>
    mutate(debate = d, sample = "EXCL")
})

results_all_excl <- bind_rows(results_all, results_excl) |>
  mutate(
    overrep_among_selfreport   = overreporting / (validated_watchers + overreporting) * 100,
    overrep_among_nonwatchers  = overreporting / (overreporting + validated_nonwatchers) * 100
  )

print(results_all_excl)

# =========================================================================
# 4. LATEX TABLE EXPORT
# =========================================================================

fmt <- function(x, d = 1) {
  if (is.na(x)) return("")
  formatC(x, digits = d, format = "f")
}

get_val <- function(debate_id, sample_type, field) {
  row <- results_all_excl |> filter(debate == debate_id, sample == sample_type)
  if (nrow(row) == 0) return(NA_real_)
  row[[field]][1]
}

# Build one complete data row as a SINGLE-LINE string: "Label & v1 & v2 & ... & v12 \\"
# (no embedded newlines inside a row - avoids any risk of cells rendering
# on separate lines / rows being visually stacked in the compiled PDF).
build_row <- function(label, field, digits = 1, is_n = FALSE) {
  cells <- purrr::map_chr(debate_key, function(d) {
    all_v  <- get_val(d, "All",  field)
    excl_v <- get_val(d, "EXCL", field)
    all_txt  <- if (is_n) (if (is.na(all_v))  "" else as.character(all_v))  else fmt(all_v, digits)
    excl_txt <- if (is_n) (if (is.na(excl_v)) "" else as.character(excl_v)) else fmt(excl_v, digits)
    paste(all_txt, excl_txt, sep = " & ")
  })
  paste0(label, " & ", paste(cells, collapse = " & "), " \\\\")
}

latex_table <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\scriptsize",
  "\\caption{Overreporting and underreporting of exposure to televised election debates -- without and with speeders and straightliners correction}",
  "\\setlength{\\tabcolsep}{2.5pt}",
  "\\begin{tabular}{l*{12}{c}}",
  "\\toprule",
  paste(
    "&",
    paste(
      sapply(c("2021 CT", "2021 NOVA", "2023 NOVA", "2023 CT", "2023 PRIMA", "2023 NOVA2"),
             function(h) paste0("\\multicolumn{2}{c}{", h, "}")),
      collapse = " & "
    ),
    "\\\\"
  ),
  paste(
    "Type of respondent",
    paste(rep("All & EXCL", 6), collapse = " & "),
    sep = " & "
  ) |> paste0(" \\\\"),
  "\\midrule",
  build_row("Validated watchers",     "validated_watchers"),
  build_row("Overreporters",          "overreporting"),
  build_row("Underreporters",         "underreporting"),
  build_row("Validated nonwatchers",  "validated_nonwatchers"),
  "\\midrule",
  build_row("N", "n", is_n = TRUE),
  "\\midrule",
  build_row("Overreporters among passive non-watchers (\\%)", "overrep_among_nonwatchers"),
  build_row("Overreporters among self-report watchers (\\%)",   "overrep_among_selfreport"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{minipage}{0.95\\textwidth}",
  "\\footnotesize",
  "Note: EXCL excludes speeders and straightliners.",
  "Source: CEPS 2021 and CPEPS 2023. Weighted.",
  "\\end{minipage}",
  "\\end{table}"
)

cat(latex_table, sep = "\n")
writeLines(latex_table, "misreport_by_debate_all_vs_excl_pooled.tex")

# =========================================================================
# END OF SCRIPT
# Outputs produced:
#   - misreport_by_debate_all_vs_excl_pooled.tex   (LaTeX table, ready for \input{})
#   - results_all_excl (in-memory tibble with underlying numbers, printed above)
# =========================================================================