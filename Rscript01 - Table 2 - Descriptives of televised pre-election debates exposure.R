# =========================================================================
# DESCRIPTIVE STATISTICS OF DEBATE EXPOSURE: SELF-REPORT VS. PASSIVE
# MEASUREMENT, BY DEBATE (WEIGHTED)
#
# For each of the six debates, this script computes the weighted
# distribution of self-reported exposure (selfreport) and the weighted
# mean/SD and threshold-category distribution of passively measured
# exposure (deb_passive_secs = passive measure in seconds /
# deb_passive_p = passive measure proportion of debate seen).
#
# This is the first descriptive step of the misreporting analysis,
# characterising exposure under each measure separately before the two
# are compared to classify respondents as accurate, over-, or
# under-reporters.
#
# data_pooled is in long format: one row per respondent per debate they
# were asked about (see main text for pooled sample composition in
# Chapter 4.2). Because every statistic below is computed separately
# for one debate at a time, and a respondent can only be asked about a
# given debate once, each debate's subset is effectively one row per
# respondent. Hence, svydesign(ids = ~1, ...) is appropriate here
# (no clustering is needed).
#
# This differs from the modelling scripts, which pool observations
# across debates and therefore cluster on respondent ID.
# =========================================================================

library(dplyr)
library(survey)

# ---- 0. Debate metadata --------------------------------------------------
# debate_key: the six levels of data_pooled$debate, used both to filter
# and to fix the column order. durations: real debate length (mins),
# used as a reference row for table. debate_labels_display: cosmetic short labels
# for the table header only.

debate_key <- c("2021-CT", "2021-NOVA", "2023-NOVA1",
                "2023-CT", "2023-PRIMA", "2023-NOVA2")

durations <- setNames(c(109, 89, 89, 116, 84, 89), debate_key)
debate_labels_display <- setNames(
  c("21-CT", "21-NOVA", "23-NOVA1", "23-CT", "23-PRIMA", "23-NOVA2"),
  debate_key
)
debate_order <- debate_key

# ---- 1. Derived variables ------------------------------------------------
# deb_passive_min: passive measure in minutes (for a readable M/SD row).
# deb_passive_cat: four-category version of passive exposure, combining
# absolute seconds watched and share of the debate's duration watched
# (deb_passive_p, already in data_pooled): no exposure / under 3 minutes
# / at least 3 minutes but under 10% of duration / 10%+ of duration.
# Both are local to this script (not stored in data_pooled.rds).

data_pooled <- data_pooled %>%
  mutate(
    deb_passive_min = deb_passive_secs / 60,
    deb_passive_cat = case_when(
      deb_passive_secs == 0 ~ "0 secs",
      deb_passive_secs > 0 & deb_passive_secs < 180 ~ "0 < x < 3 min",
      deb_passive_secs >= 180 & deb_passive_p < 0.10 ~ "3 <= x < 10%",
      deb_passive_p >= 0.10 ~ "10% <= x",
      TRUE ~ NA_character_
    ),
    deb_passive_cat = factor(
      deb_passive_cat,
      levels = c("0 secs", "0 < x < 3 min", "3 <= x < 10%", "10% <= x")
    )
  )

# ---- 2. Weighted stats for one debate's subset ---------------------------
# Returns: (a) weighted % distribution of selfreport, (b) weighted M/SD
# of the passive measure in secs and mins, (c) weighted % distribution
# of deb_passive_cat.

get_stats <- function(sub_df) {
  sub_df <- sub_df %>% filter(!is.na(weight))  # svydesign needs non-missing weights
  des <- svydesign(ids = ~1, weights = ~weight, data = sub_df)
  
  self_tab <- svymean(~selfreport, des, na.rm = TRUE) * 100
  self_pct <- setNames(as.numeric(self_tab), names(self_tab))
  
  m_sec  <- svymean(~deb_passive_secs, des, na.rm = TRUE)
  sd_sec <- sqrt(svyvar(~deb_passive_secs, des, na.rm = TRUE))
  m_min  <- svymean(~deb_passive_min, des, na.rm = TRUE)
  sd_min <- sqrt(svyvar(~deb_passive_min, des, na.rm = TRUE))
  
  cat_tab <- svymean(~deb_passive_cat, des, na.rm = TRUE) * 100
  cat_pct <- setNames(as.numeric(cat_tab), names(cat_tab))
  
  list(
    self_pct = self_pct,
    m_sec = as.numeric(m_sec), sd_sec = as.numeric(sd_sec),
    m_min = as.numeric(m_min), sd_min = as.numeric(sd_min),
    cat_pct = cat_pct
  )
}

# ---- 3. Run for each debate -----------------------------------------------
# Filtering to debate == d is what "resolves" the long-format data into
# a one-row-per-respondent subset for the statistics above.

results <- lapply(debate_order, function(d) {
  get_stats(data_pooled %>% filter(debate == d))
})
names(results) <- debate_order

# ---- 4. Formatting helpers -------------------------------------------------

fmt <- function(x, d = 1) formatC(x, digits = d, format = "f")

# svymean() prefixes output names with the variable name (e.g.
# "selfreportNo"), so match a category by substring rather than exact name.
get_named_pct <- function(pct_vec, label) {
  idx <- grepl(label, names(pct_vec), fixed = TRUE)
  if (!any(idx)) return("--")
  fmt(pct_vec[idx][1])
}

# Plain-text labels for MATCHING against svymean() output:
self_levels <- c("No", "Yes, I saw a part of it", "Yes, I saw the whole debate")
cat_levels  <- c("0 secs", "0 < x < 3 min", "3 <= x < 10%", "10% <= x")

# Same four categories, escaped for LaTeX DISPLAY (unescaped "%" would
# make everything after it on the line a comment in Overleaf); aligned
# by position with cat_levels.
cat_levels_display <- c(
  "0 secs",
  "$0 < x < 3$ min",
  "$3 \\leq x < 10\\%$",
  "$10\\% \\leq x $"
)

# ---- 5. Build the LaTeX table ----------------------------------------------
# One column per debate; rows: duration, self-report categories, passive
# measure M/SD (secs, mins), passive-measure threshold categories.

row_duration <- paste("\\textbf{Duration (mins)}",
                      paste(durations[debate_order], collapse = " & "), sep = " & ")

self_rows <- sapply(self_levels, function(lv) {
  vals <- sapply(results, function(res) get_named_pct(res$self_pct, lv))
  paste(lv, paste(vals, collapse = " & "), sep = " & ")
})

secs_M  <- paste("M", paste(fmt(sapply(results, `[[`, "m_sec")), collapse = " & "), sep = " & ")
secs_SD <- paste("SD", paste(fmt(sapply(results, `[[`, "sd_sec")), collapse = " & "), sep = " & ")

min_M   <- paste("M", paste(fmt(sapply(results, `[[`, "m_min")), collapse = " & "), sep = " & ")
min_SD  <- paste("SD", paste(fmt(sapply(results, `[[`, "sd_min")), collapse = " & "), sep = " & ")

cat_rows <- sapply(seq_along(cat_levels), function(i) {
  lv <- cat_levels[i]
  lv_disp <- cat_levels_display[i]
  vals <- sapply(results, function(res) get_named_pct(res$cat_pct, lv))
  paste(lv_disp, paste(vals, collapse = " & "), sep = " & ")
})

latex_table <- c(
  "\\begin{tabular}{lcccccc}",
  "\\hline",
  paste(" &", paste(debate_labels_display[debate_order], collapse = " & "), "\\\\"),
  "\\hline",
  paste0(row_duration, " \\\\"),
  "{\\textbf{Self-report (\\%)}} & & & & & & \\\\",
  paste0(self_rows, " \\\\"),
  "{\\textbf{Passive measure (secs)}} & & & & & & \\\\",
  paste0(secs_M, " \\\\"),
  paste0(secs_SD, " \\\\"),
  "{\\textbf{Passive measure (mins)}} & & & & & & \\\\",
  paste0(min_M, " \\\\"),
  paste0(min_SD, " \\\\"),
  "{\\textbf{Passive measure (thresholds) (\\%)}} & & & & & & \\\\",
  paste0(cat_rows, " \\\\"),
  "\\hline",
  "\\end{tabular}"
)

cat(latex_table, sep = "\n")