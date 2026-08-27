# =========================================================================

# DEBATE REACH: NIELSEN VS. PASSIVE MEASURE (3-MINUTE THRESHOLD) VS. SELF-REPORT 
#
# For each debate, computes the weighted percentage (with 95% CI) of:
# (a) "Passive measure (at least 3 minutes)": deb_passive_secs >= 180
# (b) "Self-report (saw at least part of)": selfreport %in%
# c("Yes, I saw a part of it", "Yes, I saw the whole debate")
# for comparison with Nielsen reach 3 minutes.
# DATA STRUCTURE: data_pooled is long format, one row per respondent
# (RADIOMETER_ID2) per debate they were asked about. Since each debate's
# estimate is computed on the subset with debate == d, and a respondent
# cannot be asked about the same debate twice, each row represents one
# respondent within a single debate's estimate. No clustering or repeated-
# measures structure is therefore used within a single debate's estimate.
#
# CONFIDENCE INTERVALS: computed via survey::svyciprop(..., method =
# "mean"), which gives a Wald-type CI based on Taylor-linearized
# variance (i.e. proportion +/- 1.96 * design-based SE), matching the
# "Taylor linearization (Wald method)" noted in the table source line.
# =========================================================================

library(dplyr)
library(survey)
library(tibble)
library(purrr)

# =========================================================================
# 1. DEBATE METADATA
# =========================================================================

debate_key <- c("2021-CT", "2021-NOVA", "2023-NOVA1",
                "2023-CT", "2023-PRIMA", "2023-NOVA2")
debate_labels_display <- setNames(
  c("2021-CT", "2021-NOVA", "2023-NOVA", "2023-CT", "2023-PRIMA", "2023-NOVA2"),
  debate_key
)
nielsen_reach <- setNames(
  c(20.8, 21.4, 35.8, 32.0, 29.5),
  c("2021-CT", "2021-NOVA", "2023-NOVA1", "2023-CT", "2023-NOVA2")
)
debate_ids_display <- c("2021-CT", "2021-NOVA", "2023-NOVA1", "2023-CT", "2023-NOVA2")

# =========================================================================
# 2. HELPER: WEIGHTED % + WALD (TAYLOR-LINEARIZED) 95% CI FOR ONE DEBATE
# =========================================================================

get_debate_reach <- function(data, debate_id) {
  d <- data |>
    filter(debate == debate_id, !is.na(weight)) |>
    mutate(
      passive_3min = as.numeric(deb_passive_secs >= 180),
      selfreport_any = as.numeric(
        selfreport %in% c("Yes, I saw a part of it", "Yes, I saw the whole debate")
      )
    ) |>
    filter(!is.na(passive_3min), !is.na(selfreport_any))
  design <- survey::svydesign(
    ids     = ~1,
    weights = ~weight,
    data    = d
  )
  
  # Wald (Taylor-linearized) CI: method = "mean" uses the design-based
  
  # SE from svymean() and a normal-approximation (Wald) interval, rather
  
  # than the default logit-transformed CI of svyciprop().
  
  ci_passive <- survey::svyciprop(~passive_3min,   design, method = "mean", level = 0.95)
  ci_self    <- survey::svyciprop(~selfreport_any, design, method = "mean", level = 0.95)
  tibble::tibble(
    debate       = debate_id,
    n            = nrow(d),
    passive_pct  = as.numeric(ci_passive) * 100,
    passive_lo   = confint(ci_passive)[1] * 100,
    passive_hi   = confint(ci_passive)[2] * 100,
    self_pct     = as.numeric(ci_self) * 100,
    self_lo      = confint(ci_self)[1] * 100,
    self_hi      = confint(ci_self)[2] * 100
  )
}

# =========================================================================
# 3. COMPUTE FOR ALL REQUESTED DEBATES
# =========================================================================

reach_results <- purrr::map_dfr(debate_ids_display, ~ get_debate_reach(data_pooled, .x))
print(reach_results)

# =========================================================================
# 4. LATEX TABLE EXPORT
# =========================================================================

fmt_pct   <- function(x, d = 1) formatC(x, digits = d, format = "f")
fmt_ci    <- function(lo, hi, d = 1) paste0(fmt_pct(lo, d), "; ", fmt_pct(hi, d))

rows_txt <- purrr::map_chr(debate_ids_display, function(id) {
  r <- reach_results |> filter(debate == id)
  paste(
    debate_labels_display[id], "&",
    fmt_pct(nielsen_reach[id]), "&",
    paste0(fmt_pct(r$passive_pct), " [", fmt_ci(r$passive_lo, r$passive_hi), "]"), "&",
    paste0(fmt_pct(r$self_pct),    " [", fmt_ci(r$self_lo,    r$self_hi),    "]"),
    "\\\\"
  )
})

latex_table <- c(
  "\\begin{table}[H]",
  "\\centering",
  "\\caption{Debate reach: Nielsen ratings vs. passive measure vs. self-report (3-minute threshold)}",
  "\\label{tab:reach_comparison}",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{lccc}",
  "\\toprule",
  "Debate name & Nielsen reach & Passive measures & Self-report \\\\",
  " & (at least 3 minutes) & (at least 3 minutes) [95\\% CI] & (saw at least part of) [95\\% CI] \\\\",
  "\\midrule",
  rows_txt,
  "\\bottomrule",
  "\\end{tabular}%",
  "}",
  "",
  "\\vspace{0.5em}",
  "",
  "\\begin{minipage}{\\textwidth}",
  "\\footnotesize",
  "\\raggedright",
  "Source: ATO: personal communication (for Nielsen reach); CPES 2021 (N = 804), CPES 2023 (N = 1326, 1358); weighted by gender, age, education, social status, media consumption and internet use (for passive measures and self-reports).",
  "Confidence intervals were calculated using Taylor linearization (Wald method)",
  "\\end{minipage}",
  "",
  "\\end{table}"
)

cat(latex_table, sep = "\n")
writeLines(latex_table, "reach_comparison_table_pooled.tex")

# =========================================================================
# END OF SCRIPT
# Outputs produced:
#   - reach_comparison_table_pooled.tex   (LaTeX table, ready for \input{})
#   - reach_results (in-memory tibble with underlying numbers, printed above)
# =========================================================================
