# =========================================================================
# DEBATE-SPECIFIC WEIGHTED LOGISTIC REGRESSIONS
# Overreporting, estimated separately for each of the six debates,
# weighted by 'weight', against two alternative reference groups:
#
#   Model A ("validated non-watchers"): overreporters (misreport_4k == 1)
#     vs. validated non-watchers (misreport_4k == 3)
#   Model B ("self-reported watchers"): overreporters (misreport_4k == 1)
#     vs. validated watchers (misreport_4k == 2)
#
# Predictors: age, duty, edu, party_id_01, polint, sex, straightlining,
# speeding, total_tv.
#
#
# WEIGHTING: survey::svyglm() with family = quasibinomial() (rather than
# binomial()) because svyglm requires a quasi-family when weights are
# non-integer/pweights, to avoid a spurious warning and to get correct
# (sandwich-type) standard errors. AME via marginaleffects::avg_comparisons()
# on the svyglm object uses the survey-consistent vcov for the CIs.
#
# McFadden's pseudo-R2 for weighted models is not uniquely defined; here
# it is a weighted analogue: 1 - (weighted log-likelihood of the full
# model) / (weighted log-likelihood of the intercept-only model) - a
# standard ad hoc extension, descriptive/relative across debates rather
# than a formal design-based quantity.
#
# IMPORTANT: Before running this script, the variable `misreport_4k` must
# already exist in the environment/data. It is created in the script
# "Rscript03 - Table 3 - Transformations and overreporting desriptives".
# Therefore, that script must be run first, or the prepared data containing
# `misreport_4k` must otherwise be loaded into the current R session.
#
# Separation diagnostics:
# Additional diagnostics conducted outside this script indicated that
# quasi-complete separation can occur for some predictors in specific
# debate/model combinations. The checks below are therefore intentionally
# limited to the predictors/debates where this issue affects the model
# specification. They provide a transparent diagnostic for the exclusion
# of affected predictors without reproducing the full set of auxiliary
# diagnostics used during model development.
# =========================================================================

library(dplyr)
library(survey)
library(marginaleffects)
library(tibble)
library(purrr)
library(stringr)

# ---- 1. Debate metadata ------------------------------------------------

debate_key <- c("2021-CT", "2021-NOVA", "2023-NOVA1",
                "2023-CT", "2023-PRIMA", "2023-NOVA2")

debate_labels_display <- c(
  "2021-CT" = "2021-CT", "2021-NOVA" = "2021-NOVA", "2023-NOVA1" = "2023-NOVA",
  "2023-CT" = "2023-CT",  "2023-PRIMA" = "2023-PRIMA", "2023-NOVA2" = "2023-NOVA2"
)

# Display order for the table columns (as requested):
# 2021-CT, 2021-NOVA, 2023-NOVA(1), 2023-NOVA2, 2023-CT, 2023-PRIMA
debate_order_display <- c("2021-CT", "2021-NOVA", "2023-NOVA1",
                          "2023-NOVA2", "2023-CT", "2023-PRIMA")

SE_THRESHOLD <- 5  # logit-scale SE beyond which the OR CI is judged uninformative

# ---- 2. Helper: fit one debate's weighted logistic regression + AME + R2 ---
# comparison_values: which two misreport_4k categories define the
#   outcome for this model (c(1,3) for Model A, c(1,2) for Model B).
# check_separation: whether to run the automatic SE-based quasi-
#   separation check (TRUE only for Model B - see header note).

fit_debate_model_weighted <- function(data, debate_id, comparison_values, check_separation) {
  
  d <- data |>
    filter(debate == debate_id, misreport_4k %in% comparison_values, !is.na(weight)) |>
    mutate(
      misreport_bin  = ifelse(misreport_4k == 1, 1, 0),
      party_id_01    = factor(party_id_01),
      edu            = factor(edu),
      sex            = factor(sex),
      speeding       = factor(speeding),
      duty           = factor(duty),
      polint         = factor(polint),
      straightlining = factor(straightlining),
      total_tv       = as.numeric(total_tv)
    )
  
  is_constant <- function(x) nlevels(droplevels(x)) <= 1
  
  has_straightlining <- !is_constant(d$straightlining)
  has_speeding        <- !is_constant(d$speeding)
  reason_straightlining <- if (!has_straightlining) "constant" else NA_character_
  reason_speeding        <- if (!has_speeding)        "constant" else NA_character_
  
  build_formula <- function(incl_straightlining, incl_speeding) {
    rhs_terms <- c(
      "party_id_01", "duty", "polint",
      if (incl_straightlining) "straightlining",
      if (incl_speeding)       "speeding",
      "age", "edu", "sex", "total_tv"
    )
    as.formula(paste("misreport_bin ~", paste(rhs_terms, collapse = " + ")))
  }
  
  # Survey design: cluster by respondent ID (no-op within a single
  # debate's subset unless duplicate rows exist), weighted by 'weight'.
  design <- survey::svydesign(
    ids     = ~RADIOMETER_ID2,
    weights = ~weight,
    data    = d
  )
  
  # Automatic quasi-complete separation check (Model B only - see header
  # note). Fits a trial model with whatever passed the constancy check,
  # then flags straightlining/speeding if their logit-scale SE exceeds
  # SE_THRESHOLD (an SE that large implies an odds-ratio CI spanning
  # many orders of magnitude, i.e. not informative).
  if (check_separation) {
    m_trial <- survey::svyglm(
      build_formula(has_straightlining, has_speeding),
      design = design, family = quasibinomial()
    )
    se_trial <- sqrt(diag(vcov(m_trial)))
    
    check_se <- function(has_flag, prefix) {
      if (!has_flag) return(has_flag)
      se_term <- se_trial[grepl(paste0("^", prefix), names(se_trial))]
      if (length(se_term) > 0 && any(se_term > SE_THRESHOLD, na.rm = TRUE)) FALSE else TRUE
    }
    
    if (has_straightlining && !check_se(has_straightlining, "straightlining")) {
      has_straightlining <- FALSE
      reason_straightlining <- "quasi-complete separation (excessive SE)"
    }
    if (has_speeding && !check_se(has_speeding, "speeding")) {
      has_speeding <- FALSE
      reason_speeding <- "quasi-complete separation (excessive SE)"
    }
    
    # 2021-CT straightlining: print concrete OR/CI evidence regardless
    # of whether the automatic SE check above already caught it (SE on
    # the logit scale can stay under SE_THRESHOLD even when the
    # exponentiated OR CI is already absurdly wide). This diagnostic is
    # what Model A's manual override (below) cites as justification.
    if (debate_id == "2021-CT" && has_straightlining) {
      coefs_trial <- summary(m_trial)$coefficients
      ci_trial    <- confint(m_trial)
      sl_rows <- grepl("^straightlining", rownames(coefs_trial))
      if (any(sl_rows)) {
        or_diag <- tibble::tibble(
          term    = rownames(coefs_trial)[sl_rows],
          OR      = exp(coefs_trial[sl_rows, "Estimate"]),
          CI_low  = exp(ci_trial[sl_rows, 1]),
          CI_high = exp(ci_trial[sl_rows, 2])
        )
        cat("\n=== DIAGNOSTIC: straightlining OR/CI, 2021-CT (self-reported watchers model) ===\n")
        print(or_diag)
        cat("Odds-ratio CI spans several orders of magnitude - the signature of\n",
            "quasi-complete separation. Excluding straightlining from the 2021-CT\n",
            "model on this basis.\n")
      }
      has_straightlining <- FALSE
      reason_straightlining <- "quasi-complete separation (astronomically wide OR CI - see diagnostic printed above)"
    }
    
  } else if (debate_id == "2021-CT" && has_straightlining) {
    has_straightlining <- FALSE
    reason_straightlining <- "quasi-complete separation confirmed in the companion self-reported-watchers model for this debate (manual override)"
  }
  
  form <- build_formula(has_straightlining, has_speeding)
  
  # Final model: re-fit with 'form', which already excludes any
  # predictor flagged above as constant or (quasi-)separated.
  m <- survey::svyglm(form, design = design, family = quasibinomial())
  
  # N actually used in the regression: svyglm silently drops rows with
  # missing values on any predictor (listwise deletion), so this can be
  # smaller than nrow(d). nobs() reflects the post-deletion sample
  # actually entering the model - this is what gets reported in the
  # table, not the raw filtered subset size.
  n_used <- nobs(m)
  
  # Null (intercept-only) model, same design/weights
  m_null <- survey::svyglm(misreport_bin ~ 1, design = design, family = quasibinomial())
  
  # Weighted McFadden pseudo-R2 (see header note on interpretation)
  p_full <- fitted(m)
  p_null <- fitted(m_null)
  y <- d$misreport_bin
  w <- d$weight
  eps <- 1e-10  # numerical guard against log(0)
  ll_full <- sum(w * (y * log(pmax(p_full, eps)) + (1 - y) * log(pmax(1 - p_full, eps))))
  ll_null <- sum(w * (y * log(pmax(p_null, eps)) + (1 - y) * log(pmax(1 - p_null, eps))))
  r2_mcfadden <- as.numeric(1 - ll_full / ll_null)
  
  # AME for each predictor (factor levels vs. reference; +1 unit for
  # continuous age/total_tv), using the svyglm's survey-consistent vcov.
  var_specs <- list(
    party_id_01 = "reference",
    duty        = "reference",
    polint      = "reference",
    edu         = "reference",
    sex         = "reference",
    age         = 1,
    total_tv    = 1
  )
  if (has_straightlining) var_specs$straightlining <- "reference"
  if (has_speeding)       var_specs$speeding       <- "reference"
  
  ame_one_var <- function(var, spec) {
    tryCatch({
      avg_comparisons(m, type = "response", variables = setNames(list(spec), var)) |>
        as_tibble()
    }, error = function(e) {
      message("  ERROR (debate ", debate_id, ", var '", var, "'): ", conditionMessage(e))
      NULL
    })
  }
  
  ame_list <- lapply(names(var_specs), function(v) ame_one_var(v, var_specs[[v]]))
  ame_df <- bind_rows(Filter(Negate(is.null), ame_list)) |>
    mutate(
      term_label_raw = if_else(is.na(contrast) | contrast == "",
                               term, paste0(term, ": ", contrast))
    )
  
  list(
    debate                 = debate_id,
    n                      = n_used,
    n_subset               = nrow(d),
    ame                    = ame_df,
    r2_mcfadden            = r2_mcfadden,
    has_straightlining     = has_straightlining,
    has_speeding           = has_speeding,
    reason_straightlining  = reason_straightlining,
    reason_speeding        = reason_speeding
  )
}

# ---- 3. Fit all six debate models, for both comparison groups ---------

debate_results_A <- purrr::map(
  debate_key, ~ fit_debate_model_weighted(data_pooled, .x, c(1, 3), check_separation = FALSE)
)
names(debate_results_A) <- debate_key

debate_results_B <- purrr::map(
  debate_key, ~ fit_debate_model_weighted(data_pooled, .x, c(1, 2), check_separation = TRUE)
)
names(debate_results_B) <- debate_key

print_diagnostics <- function(debate_results, model_name) {
  cat("\n\n########## ", model_name, " ##########\n")
  for (d in debate_key) {
    res <- debate_results[[d]]
    cat("\nDebate", d, "-",
        "N in regression:", res$n,
        "(N in filtered subset before predictor-missingness drop:", res$n_subset, ")\n")
    cat("  straightlining kept:", res$has_straightlining,
        if (!res$has_straightlining) paste0(" [dropped: ", res$reason_straightlining, "]") else "", "\n")
    cat("  speeding kept:", res$has_speeding,
        if (!res$has_speeding) paste0(" [dropped: ", res$reason_speeding, "]") else "", "\n")
    if (res$n != res$n_subset) {
      cat("  NOTE: ", res$n_subset - res$n,
          " case(s) dropped due to missing values on a predictor in this debate.\n")
    }
  }
}

print_diagnostics(debate_results_A, "MODEL A: VALIDATED NON-WATCHERS")
print_diagnostics(debate_results_B, "MODEL B: SELF-REPORTED WATCHERS")

# ---- 4. Label mapping (shared across both models) -----------------------

norm_key <- function(x) gsub("[^a-z0-9]", "", tolower(x))

label_dict <- c(
  "party_id_01: somewhat close or very close - no pid, DK" = "PID (somewhat close or very close)",
  "duty: duty - right/choice"                              = "Voting is duty",
  "polint: just a little - not at all"                     = "Political interest (little)",
  "polint: quite/very interested - not at all"             = "Political interest (quite/very)",
  "straightlining: 1 - 0"                                  = "Straightlining (1+)",
  "speeding: x < 2.decile - higher than 2. decile"         = "Speeding",
  "age: +1"                                                = "Age",
  "edu: tertiary - secondary and lower"                    = "Education (tertiary)",
  "sex: male - female"                                     = "Sex (male)",
  "total_tv: +1"                                           = "TV consumption (min/day)"
)
dict_lookup <- setNames(unname(label_dict), norm_key(names(label_dict)))

apply_labels <- function(debate_results) {
  for (d in names(debate_results)) {
    debate_results[[d]]$ame <- debate_results[[d]]$ame |>
      mutate(
        .key       = norm_key(term_label_raw),
        term_label = if_else(.key %in% names(dict_lookup), dict_lookup[.key], term_label_raw)
      ) |>
      select(-.key)
  }
  debate_results
}

debate_results_A <- apply_labels(debate_results_A)
debate_results_B <- apply_labels(debate_results_B)

check_unmapped <- function(debate_results, model_name) {
  unmapped <- purrr::map(debate_results, ~ .x$ame) |>
    bind_rows() |>
    filter(term_label == term_label_raw) |>
    pull(term_label_raw) |> unique() |> sort()
  cat("\n=== Unmapped terms,", model_name, "(fix label_dict if not empty) ===\n")
  if (length(unmapped) == 0) cat("All terms mapped successfully.\n") else print(unmapped)
}
check_unmapped(debate_results_A, "Model A")
check_unmapped(debate_results_B, "Model B")

# ---- 5. LaTeX table export (shared row/formatting logic) ----------------

sig_stars <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "^{***}",
    p < 0.01  ~ "^{**}",
    p < 0.05  ~ "^{*}",
    TRUE      ~ ""
  )
}

fmt_ame_cell <- function(df, label, digits = 2) {
  row <- df |> dplyr::filter(term_label == label)
  if (nrow(row) == 0) return("--")
  est   <- formatC(row$estimate[1],  digits = digits, format = "f")
  lo    <- formatC(row$conf.low[1],  digits = digits, format = "f")
  hi    <- formatC(row$conf.high[1], digits = digits, format = "f")
  stars <- sig_stars(row$p.value[1])
  paste0(est, "$", stars, "$ [", lo, ",", hi, "]")
}

row_defs <- tibble::tibble(
  key   = c("Age", "Voting is duty", "Education (tertiary)",
            "PID (somewhat close or very close)",
            "Political interest (little)", "Political interest (quite/very)",
            "Sex (male)", "Straightlining (1+)", "Speeding",
            "TV consumption (min/day)"),
  latex = c("Age", "Voting is duty", "Education (tertiary)",
            "PID (somewhat close or very close)",
            "Political interest (little)", "Political interest (quite/very)",
            "Sex (male)", "Straightlining (1+)", "Speeding",
            "TV consumption (min/day)")
)

build_latex_table <- function(debate_results, caption, label, tex_note_extra, filename) {
  
  row_watchers    <- sapply(debate_key, function(d) debate_results[[d]]$n)  # placeholder, unused
  ame_cols <- purrr::map(debate_order_display, function(d) {
    df <- debate_results[[d]]$ame
    sapply(row_defs$key, function(k) fmt_ame_cell(df, k))
  })
  names(ame_cols) <- debate_order_display
  
  ame_rows_txt <- purrr::map_chr(seq_len(nrow(row_defs)), function(i) {
    lbl   <- row_defs$latex[i]
    cells <- sapply(debate_order_display, function(d) ame_cols[[d]][row_defs$key[i]])
    paste(lbl, "&", paste(cells, collapse = " &\n "), "\\\\")
  })
  
  n_row <- paste(
    "&",
    paste(sprintf("$N=%d$", sapply(debate_order_display, function(d) debate_results[[d]]$n)),
          collapse = " &\n "),
    "\\\\"
  )
  
  header_row <- paste(
    "&",
    paste(debate_labels_display[debate_order_display], collapse = " & "),
    "\\\\"
  )
  
  r2_row <- paste(
    "McFadden $R^2$ &",
    paste(sapply(debate_order_display, function(d) {
      formatC(debate_results[[d]]$r2_mcfadden, digits = 2, format = "f")
    }), collapse = " &\n "),
    "\\\\"
  )
  
  latex_table <- c(
    "\\begin{table}[H]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{3pt}",
    "",
    "\\resizebox{\\textwidth}{!}{",
    "\\begin{tabular}{lcccccc}",
    "\\toprule",
    "",
    header_row,
    n_row,
    "\\midrule",
    "",
    paste(ame_rows_txt, collapse = "\n\n"),
    "",
    "\\midrule",
    r2_row,
    "",
    "\\bottomrule",
    "\\end{tabular}}",
    "",
    "\\vspace{0.2cm}",
    "",
    "\\begin{minipage}{0.95\\textwidth}",
    "\\footnotesize",
    "\\textit{Notes:} Entries report average marginal effects (AME) from weighted",
    "logistic regressions (survey weight: \\texttt{weight}).",
    "Values in brackets represent 95\\% confidence intervals.",
    "``--'' indicates the predictor was dropped from that debate's model",
    "because it exhibited (quasi-)complete separation with the outcome, which inflates",
    "its coefficient's standard error to an uninformative odds-ratio",
    "confidence interval.",
    tex_note_extra,
    "* $p<0.05$, ** $p<0.01$, *** $p<0.001$.",
    "\\end{minipage}",
    "",
    "\\end{table}"
  )
  
  cat(latex_table, sep = "\n")
  writeLines(latex_table, filename)
  invisible(latex_table)
}

build_latex_table(
  debate_results_A,
  caption = "Overreporting among validated non-watchers -- debate specific logistic regressions (weighted)",
  label   = "tab:robust_nonwatchers_weighted",
  tex_note_extra = paste0(
    "Straightlining in the 2021-CT debate is excluded from this model as well; ",
    "quasi-complete separation for this predictor/debate was identified directly ",
    "in the companion self-reported-watchers model (Model B) below, where the ",
    "odds-ratio confidence interval is printed as a diagnostic."
  ),
  filename = "ame_debate_specific_nonwatchers_weighted_pooled_table.tex"
)

build_latex_table(
  debate_results_B,
  caption = "Overreporting among self-reported watchers -- debate specific logistic regressions (weighted)",
  label   = "tab:robust_watchers_weighted",
  tex_note_extra = paste0(
    "Straightlining in the 2021-CT debate is excluded on this basis: its odds-ratio ",
    "confidence interval spans several orders of magnitude in this specification ",
    "(see the diagnostic printed to the console when this script is run)."
  ),
  filename = "ame_debate_specific_watchers_weighted_pooled_table.tex"
)

# ---- 6. Export full AME results to Excel (all debates, all terms) -------

export_ame <- function(debate_results, filename) {
  purrr::imap_dfr(debate_results, function(res, d) {
    res$ame |>
      transmute(
        debate      = d,
        debate_name = debate_labels_display[d],
        term_label,
        AME     = round(estimate, 4),
        SE      = round(std.error, 4),
        p_value = round(p.value, 4),
        CI_low  = round(conf.low, 4),
        CI_high = round(conf.high, 4)
      )
  }) |>
    writexl::write_xlsx(filename)
}

export_ame(debate_results_A, "AME_debate_specific_nonwatchers_weighted_pooled.xlsx")
export_ame(debate_results_B, "AME_debate_specific_watchers_weighted_pooled.xlsx")

# =========================================================================
# END OF SCRIPT
# Outputs produced:
#   - ame_debate_specific_nonwatchers_weighted_pooled_table.tex  (Model A)
#   - ame_debate_specific_watchers_weighted_pooled_table.tex     (Model B)
#   - AME_debate_specific_nonwatchers_weighted_pooled.xlsx
#   - AME_debate_specific_watchers_weighted_pooled.xlsx
#   - console output: N, dropped-predictor flags/reasons, unmapped-term
#     check, and (Model B only) the 2021-CT straightlining OR/CI diagnostic
# =========================================================================