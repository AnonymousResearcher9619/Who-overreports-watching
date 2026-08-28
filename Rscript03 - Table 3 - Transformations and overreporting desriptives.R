library(dplyr)
library(survey)
library(haven)
library(purrr)

# =========================================================================
# MISREPORTING TYPOLOGIES
#
# Classifies each respondent-debate observation by comparing self-
# reported exposure (selfreport) against passively measured exposure,
# under three alternative thresholds for what counts as "having watched"
# on the passive measure: at least 3 minutes (180 seconds; the main
# specification), any exposure at all (0 seconds), and at least 10% of
# the debate's duration. For each threshold, a 4-category typology is
# built misreporting variable misreport_4k (overreporting / validated
# watchers / validated nonwatchers / underreporting).
# 
# The second half of the script produces a
# descriptive, weighted per-debate breakdown of the main (3-minute)
# 4-category typology, formatted as a LaTeX table.
# =========================================================================


# ---- Misreporting by the 180-second (3-minute) threshold: main spec -----
# 1 = overreporting   (said "yes", passive < 180s)
# 2 = validated watchers    (said "yes", passive >= 180s)
# 3 = validated nonwatchers (said "no",  passive < 180s)
# 4 = underreporting  (said "no",  passive >= 180s)

data_pooled <- data_pooled |>
  dplyr::mutate(
    misreport_4k = dplyr::case_when(
      deb_passive_secs < 180 &
        selfreport %in% c(
          "Yes, I saw the whole debate",
          "Yes, I saw a part of it"
        ) ~ 1,
      deb_passive_secs >= 180 &
        selfreport %in% c(
          "Yes, I saw the whole debate",
          "Yes, I saw a part of it"
        ) ~ 2,
      deb_passive_secs < 180 &
        selfreport == "No" ~ 3,
      deb_passive_secs >= 180 &
        selfreport == "No" ~ 4,
      TRUE ~ NA_real_
    ),
    misreport_4k = haven::labelled(
      misreport_4k,
      labels = c(
        overreporting = 1,
        `validated watchers` = 2,
        `validated nonwatchers` = 3,
        underreporting = 4
      ),
      label = "Misreporting: 4 categories (3-minute threshold)"
    )
  )

table(data_pooled$misreport_4k, useNA = "ifany")



# ---- Misreporting by the 0-second threshold (robustness check) -----------
# Same logic as above, but "watched" on the passive measure means ANY
# exposure (> 0 seconds) rather than 3+ minutes.

data_pooled <- data_pooled |>
  dplyr::mutate(
    misreport_4k_0sec = dplyr::case_when(
      deb_passive_secs == 0 &
        selfreport %in% c(
          "Yes, I saw the whole debate",
          "Yes, I saw a part of it"
        ) ~ 1,
      deb_passive_secs > 0 &
        selfreport %in% c(
          "Yes, I saw the whole debate",
          "Yes, I saw a part of it"
        ) ~ 2,
      deb_passive_secs == 0 &
        selfreport == "No" ~ 3,
      deb_passive_secs > 0 &
        selfreport == "No" ~ 4,
      TRUE ~ NA_real_
    ),
    misreport_4k_0sec = haven::labelled(
      misreport_4k_0sec,
      labels = c(
        overreporting = 1,
        `seen at least part of debate` = 2,
        `did not see any of the debate` = 3,
        underreporting = 4
      ),
      label = "Misreporting: 4 categories (0-second threshold)"
    )
  )

table(data_pooled$misreport_4k_0sec, useNA = "ifany")


# ---- Misreporting by the 10% threshold (robustness check) ----------------
# Same logic again, but "watched" means at least 10% of that specific
# debate's duration (deb_passive_p, already in data_pooled).

data_pooled <- data_pooled |>
  dplyr::mutate(
    misreport_4k_10pct = dplyr::case_when(
      deb_passive_p < 0.1 &
        selfreport %in% c(
          "Yes, I saw the whole debate",
          "Yes, I saw a part of it"
        ) ~ 1,
      deb_passive_p >= 0.1 &
        selfreport %in% c(
          "Yes, I saw the whole debate",
          "Yes, I saw a part of it"
        ) ~ 2,
      deb_passive_p < 0.1 &
        selfreport == "No" ~ 3,
      deb_passive_p >= 0.1 &
        selfreport == "No" ~ 4,
      TRUE ~ NA_real_
    ),
    misreport_4k_10pct = haven::labelled(
      misreport_4k_10pct,
      labels = c(
        overreporting = 1,
        `seen at least part of debate` = 2,
        `did not see any of the debate` = 3,
        underreporting = 4
      ),
      label = "Misreporting: 4 categories (10% threshold)"
    )
  )

table(data_pooled$misreport_4k_10pct, useNA = "ifany")
table(data_pooled$misreport_debate_10pct, useNA = "ifany")


# =========================================================================
# WEIGHTED FREQUENCY TABLE OF THE MAIN (3-MINUTE) 4-CATEGORY TYPOLOGY,
# PER DEBATE
#
# Descriptive companion to misreport_4k above: how common is each of the
# four categories within each debate, weighted. Computed one debate at a
# time, within which a respondent contributes exactly one row, so
# svydesign(ids = ~1, ...) is used here (no clustering needed) -
# consistent with the mean-plot and weighted-Spearman scripts, which use
# the same reasoning.
# =========================================================================

debate_key <- c("2021-CT", "2021-NOVA", "2023-NOVA1",
                "2023-CT", "2023-PRIMA", "2023-NOVA2")
debate_labels_display <- setNames(debate_key, debate_key)
debate_order <- debate_key

# Computes the weighted % of each misreport_4k category within one
# debate's subset. var_name is passed as a string so the same function
# also works for the _0sec/_10pct variants if needed.
get_weighted_misreport <- function(data, var_name, debate_id) {
  
  d <- data |>
    dplyr::filter(debate == debate_id, !is.na(weight)) |>
    dplyr::mutate(
      var = haven::as_factor(.data[[var_name]])   # turn labelled numeric into factor
    ) |>
    dplyr::filter(!is.na(var))
  
  n_unweighted <- nrow(d)
  
  des <- survey::svydesign(
    ids = ~1,
    weights = ~weight,
    data = d
  )
  
  tab <- survey::svymean(~var, des, na.rm = TRUE) * 100
  pct <- setNames(as.numeric(tab), names(tab))
  
  # svymean() prefixes output names with the variable name, so match a
  # category by substring rather than requiring an exact name match.
  get_pct <- function(label) {
    idx <- grepl(label, names(pct), fixed = TRUE)
    if (!any(idx)) return(NA_real_)
    as.numeric(pct[idx][1])
  }
  
  data.frame(
    debate = debate_id,
    n = n_unweighted,
    overreporting          = get_pct("overreporting"),
    validated_watchers     = get_pct("validated watchers"),
    validated_nonwatchers  = get_pct("validated nonwatchers"),
    underreporting         = get_pct("underreporting")
  )
}

misreport_results <- purrr::map_dfr(
  debate_key,
  ~ get_weighted_misreport(data_pooled, "misreport_4k", .x)
)

misreport_results


# =========================================================================
# LATEX TABLE EXPORT
# =========================================================================
# Table layout: the four misreport_4k categories as rows 1-4, total N
# per debate, then two derived diagnostic rates: overreporting among
# validated non-watchers (row 5) and among self-reported watchers
# (row 6) - i.e. what share of people who could have over-reported
# actually did, from each starting point.
# =========================================================================

fmt <- function(x, d = 1) formatC(x, digits = d, format = "f")

row_watchers     <- misreport_results$validated_watchers[match(debate_order, misreport_results$debate)]
row_overrep      <- misreport_results$overreporting[match(debate_order, misreport_results$debate)]
row_underrep     <- misreport_results$underreporting[match(debate_order, misreport_results$debate)]
row_nonwatchers  <- misreport_results$validated_nonwatchers[match(debate_order, misreport_results$debate)]
row_n            <- misreport_results$n[match(debate_order, misreport_results$debate)]

# row 5 = overreporters / (overreporters + validated nonwatchers)
# row 6 = overreporters / (validated watchers + overreporters)
row5_overrep_among_nonwatchers <- row_overrep / (row_overrep + row_nonwatchers) * 100
row6_overrep_among_selfreport  <- row_overrep / (row_watchers + row_overrep) * 100

latex_rows <- c(
  paste("1. validated watchers & Yes & Yes &",
        paste(fmt(row_watchers), collapse = " & "), "\\\\"),
  paste("2. overreporters & Yes & No &",
        paste(fmt(row_overrep), collapse = " & "), "\\\\"),
  paste("3. underreporters & No & Yes &",
        paste(fmt(row_underrep), collapse = " & "), "\\\\"),
  paste("4. validated nonwatchers & No & No &",
        paste(fmt(row_nonwatchers), collapse = " & "), "\\\\"),
  "\\midrule",
  paste("Total \\% & & &",
        paste(rep("100.0", 6), collapse = " & "), "\\\\"),
  paste("Total N & & &",
        paste(row_n, collapse = " & "), "\\\\"),
  "\\midrule",
  paste("5 overreporters among passive non-watchers & 5=2/(2+4); in \\% & &",
        paste(fmt(row5_overrep_among_nonwatchers), collapse = " & "), "\\\\"),
  paste("6 overreporters among self-report watchers & 6=2/(1+2); in \\% & &",
        paste(fmt(row6_overrep_among_selfreport), collapse = " & "), "\\\\")
)

latex_table <- c(
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{lllcccccc}",
  "\\toprule",
  "Type of respondent & Exposure self-report & Exposure passive measure & 2021-CT & 2021-NOVA & 2023-NOVA1 & 2023-CT & 2023-PRIMA & 2023-NOVA2 \\\\",
  "\\midrule",
  latex_rows,
  "\\bottomrule",
  "\\end{tabular}%",
  "}"
)

cat(latex_table, sep = "\n")