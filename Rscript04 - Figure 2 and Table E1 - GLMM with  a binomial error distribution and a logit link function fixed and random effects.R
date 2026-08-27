# =========================================================================
# MAIN SPECIFICATION: MIXED-EFFECTS LOGISTIC MODELS OF OVERREPORTING
#
# Two models test what predicts overreporting of debate exposure
# (claiming to have watched at least part of the debate when the passive
# measure shows under 3 minutes of exposure), each defined against a
# different reference group:
#   Model A ("Validated non-watchers"):
#       overreporters (misreport_4k == 1) vs. validated non-watchers
#       (misreport_4k == 3)
#   Model B ("Self-reported watchers"):
#       overreporters (misreport_4k == 1) vs. validated watchers
#       (misreport_4k == 2)
#
# Both are respondent-clustered random-intercept logistic regressions
# (glmmTMB) with debate as a fixed effect: the random intercept accounts
# for a respondent contributing up to six rows (one per debate) in
# data_pooled's long format, and the debate fixed effect absorbs
# debate-specific differences in overall misreporting rates so the
# predictor effects aren't confounded with which debate is being asked
# about.
#
# IMPORTANT: Before running this script, the variable `misreport_4k` must
# already exist in the environment/data. It is created in the script
# "Rscript03 - Table 3 - Transformations and overreporting desriptives".
# Therefore, that script must be run first, or the prepared data containing
# `misreport_4k` must otherwise be loaded into the current R session.

# =========================================================================

library(glmmTMB)
library(lme4)
library(performance)
library(broom.mixed)
library(texreg)
library(marginaleffects)
library(ggplot2)
library(dplyr)
library(writexl)
library(forcats)
library(pROC)
library(stringr)
library(tibble)
library(purrr)

# ---- 1. Data preparation ---------------------------------------------------
# panel_bin_A: overreporters (1) vs. validated non-watchers (3)
# panel_bin_B: overreporters (1) vs. validated watchers (2)
# misreport_bin: 1 = overreporter, 0 = reference category. Categorical
# predictors are (re-)cast to factor for safety (most are already real
# factors in data_pooled); total_tv and age kept numeric.

panel_bin_A <- data_pooled |>
  filter(misreport_4k %in% c(1, 3)) |>
  mutate(
    misreport_bin  = ifelse(misreport_4k == 1, 1, 0),
    party_id_01    = factor(party_id_01),
    edu            = factor(edu),
    sex            = factor(sex),
    speeding       = factor(speeding),
    debate         = factor(debate),
    duty           = factor(duty),
    polint         = factor(polint),
    straightlining = factor(straightlining),
    total_tv       = as.numeric(total_tv),
    age = as.numeric(age)
  )

panel_bin_B <- data_pooled |>
  filter(misreport_4k %in% c(1, 2)) |>
  mutate(
    misreport_bin  = ifelse(misreport_4k == 1, 1, 0),
    party_id_01    = factor(party_id_01),
    edu            = factor(edu),
    sex            = factor(sex),
    speeding       = factor(speeding),
    debate         = factor(debate),
    duty           = factor(duty),
    polint         = factor(polint),
    straightlining = factor(straightlining),
    total_tv       = as.numeric(total_tv)
  )

# ---- 2. Models --------------------------------------------------------------

m1_mixed <- glmmTMB(
  misreport_bin ~ party_id_01 + duty + polint + straightlining +
    speeding + age + edu + sex + total_tv + debate +
    (1 | RADIOMETER_ID2),
  data   = panel_bin_A,
  family = binomial(link = "logit")
)

m2_mixed <- glmmTMB(
  misreport_bin ~ party_id_01 + duty + polint + straightlining +
    speeding + age + edu + sex + total_tv + debate +
    (1 | RADIOMETER_ID2),
  data   = panel_bin_B,
  family = binomial(link = "logit")
)

summary(m1_mixed)
summary(m2_mixed)

# ---- 3. Convergence checks ---------------------------------------------------

cat("\n=== CONVERGENCE ===\n")
cat("Model A - pdHess:", m1_mixed$sdr$pdHess,
    ifelse(m1_mixed$sdr$pdHess, "[OK]", "[PROBLEM]"), "\n")
cat("Model B - pdHess:", m2_mixed$sdr$pdHess,
    ifelse(m2_mixed$sdr$pdHess, "[OK]", "[PROBLEM]"), "\n")

check_convergence(m1_mixed)
check_convergence(m2_mixed)

# ---- 4. Intraclass correlation (ICC) -----------------------------------------
# Share of total variance (logit scale) attributable to between-
# respondent differences: how much of the variation in overreporting is
# a stable respondent trait vs. debate-to-debate noise.

cat("\n=== ICC ===\n")
print(icc(m1_mixed))
print(icc(m2_mixed))

var_u_A <- as.numeric(VarCorr(m1_mixed)$cond$RADIOMETER_ID2[1])
var_u_B <- as.numeric(VarCorr(m2_mixed)$cond$RADIOMETER_ID2[1])
icc_A <- round(var_u_A / (var_u_A + pi^2 / 3), 3)
icc_B <- round(var_u_B / (var_u_B + pi^2 / 3), 3)
cat("ICC Model A:", icc_A, "\n")
cat("ICC Model B:", icc_B, "\n")

# ---- 5. Random-effects diagnostics -------------------------------------------
# Histogram + Q-Q plot of the estimated random intercepts (Model A), and
# a likelihood-ratio test of whether the random intercept is needed at
# all (boundary test: mixture chi-sq with 0.5 weight on chi-sq_0).

re_A <- ranef(m1_mixed)$cond$RADIOMETER_ID2[, 1]
hist(re_A, breaks = 40,
     main = "Distribution of random intercepts - Model A",
     xlab = "Random intercept (logit scale)",
     col  = "steelblue", border = "white")
qqnorm(re_A, main = "Q-Q plot - Model A")
qqline(re_A, col = "red", lwd = 2)

m_no_re_A <- glm(
  misreport_bin ~ party_id_01 + duty + polint + straightlining +
    speeding + age + edu + sex + total_tv + debate,
  data = panel_bin_A, family = binomial()
)
lrt_stat_A <- as.numeric(-2 * (logLik(m_no_re_A) - logLik(m1_mixed)))
p_lrt_A    <- 0.5 * pchisq(lrt_stat_A, df = 1, lower.tail = FALSE)
cat("\nLRT random intercept, Model A:\n")
cat("  chi2 =", round(lrt_stat_A, 2), "\n")
cat("  p    =", format(p_lrt_A, scientific = TRUE), "\n")
print(AIC(m_no_re_A, m1_mixed))
print(BIC(m_no_re_A, m1_mixed))

m_no_re_B <- glm(
  misreport_bin ~ party_id_01 + duty + polint + straightlining +
    speeding + age + edu + sex + total_tv + debate,
  data = panel_bin_B, family = binomial()
)
lrt_stat_B <- as.numeric(-2 * (logLik(m_no_re_B) - logLik(m2_mixed)))
p_lrt_B    <- 0.5 * pchisq(lrt_stat_B, df = 1, lower.tail = FALSE)
cat("\nLRT random intercept, Model B:\n")
cat("  chi2 =", round(lrt_stat_B, 2), "\n")
cat("  p    =", format(p_lrt_B, scientific = TRUE), "\n")
print(AIC(m_no_re_B, m2_mixed))
print(BIC(m_no_re_B, m2_mixed))

# ---- 6. Predictive performance (AUC) -----------------------------------------
# re.form = NA gives population-average (fixed-effects-only) predictions,
# consistent with the AME below.

pred_A <- predict(m1_mixed, type = "response", re.form = NA)
pred_B <- predict(m2_mixed, type = "response", re.form = NA)

auc_A <- round(as.numeric(auc(roc(panel_bin_A$misreport_bin, pred_A, quiet = TRUE))), 3)
auc_B <- round(as.numeric(auc(roc(panel_bin_B$misreport_bin, pred_B, quiet = TRUE))), 3)
cat("AUC Model A:", auc_A, "\n")
cat("AUC Model B:", auc_B, "\n")


# ---- 7. Average marginal effects (AME) ---------------------------------------
# avg_comparisons() gives AME on the response (probability) scale for
# both factor and continuous predictors, using the model's conditional
# variance-covariance matrix for confidence intervals.

ame_one_var <- function(model, var, spec, model_label) {
  tryCatch({
    avg_comparisons(
      model,
      type      = "response",
      vcov      = vcov(model)$cond,
      variables = setNames(list(spec), var)
    ) |>
      as_tibble() |>
      mutate(model = model_label)
  }, error = function(e) {
    message("  ERROR for variable '", var, "': ", conditionMessage(e))
    NULL
  })
}

# "reference" = AME of each factor level vs. its reference level;
# 1 = AME of a one-unit increase, for continuous predictors (age, total_tv)
var_specs <- list(
  party_id_01    = "reference",
  duty           = "reference",
  polint         = "reference",
  straightlining = "reference",
  speeding       = "reference",
  edu            = "reference",
  sex            = "reference",
  debate         = "reference",
  age            = 1,
  total_tv       = 1
)

cat("\n=== Computing AME - Model A ===\n")
ame_list_A <- lapply(names(var_specs), function(v) {
  cat("  ", v, "...")
  res <- ame_one_var(
    m1_mixed, v, var_specs[[v]],
    "Model A - Overreporting among validated nonwatchers"
  )
  cat(if (!is.null(res)) " OK\n" else " FAILED\n")
  res
})
ame_m1_mixed <- bind_rows(Filter(Negate(is.null), ame_list_A))

cat("\n=== Computing AME - Model B ===\n")
ame_list_B <- lapply(names(var_specs), function(v) {
  cat("  ", v, "...")
  res <- ame_one_var(
    m2_mixed, v, var_specs[[v]],
    "Model B - Overreporting among self-reported watchers"
  )
  cat(if (!is.null(res)) " OK\n" else " FAILED\n")
  res
})
ame_m2_mixed <- bind_rows(Filter(Negate(is.null), ame_list_B))

# Build a raw term label ("term: contrast") for each row, then map it to
# a human-readable label below.
ame_mixed_all <- bind_rows(ame_m1_mixed, ame_m2_mixed) |>
  mutate(
    term_label_raw = if_else(
      is.na(contrast) | contrast == "",
      term,
      paste0(term, ": ", contrast)
    )
  )

cat("\n=== Available term_label_raw values (for label_dict) ===\n")
print(sort(unique(ame_mixed_all$term_label_raw)))

# Human-readable label mapping - match these against the printed output
# above and adjust if factor level text differs.
label_dict <- c(
  "party_id_01: somewhat close or very close - no pid, DK" = "PID (somewhat close or very close)",
  "duty: duty - right/choice"                              = "Voting is duty",
  "polint: just a little - not at all"                     = "Political interest (just a little)",
  "polint: quite/very interested - not at all"             = "Political interest (quite / very interested)",
  "straightlining: 1 - 0"                                  = "Straightlining (1+)",
  "speeding: x < 2.decile - higher than 2. decile"         = "Speeding (x<2.decile)",
  "age: +1"                                                = "Age",
  "edu: tertiary - secondary and lower"                    = "Education (tertiary)",
  "sex: male - female"                                     = "Sex (Male)",
  "total_tv: +1"                                           = "Total TV consumption (mins/day)",
  "debate: 2021-NOVA - 2021-CT"  = "2021-NOVA",
  "debate: 2023-NOVA1 - 2021-CT" = "2023-NOVA1",
  "debate: 2023-CT - 2021-CT"    = "2023-CT",
  "debate: 2023-PRIMA - 2021-CT" = "2023-PRIMA",
  "debate: 2023-NOVA2 - 2021-CT" = "2023-NOVA2"
)

norm_key <- function(x) gsub("[^a-z0-9]", "", tolower(x))

# Mapping via an explicit left_join on the normalized key (rather than
# named-vector indexing): any unmatched row keeps its raw text and shows
# up in the "Unmapped terms" check below, so mismatches are always
# visible rather than silently mis-assigned.
label_dict_df <- tibble::tibble(
  raw_key    = names(label_dict),
  term_label = unname(label_dict)
) |>
  mutate(.key = norm_key(raw_key))

ame_mixed_all <- ame_mixed_all |>
  mutate(.key = norm_key(term_label_raw)) |>
  dplyr::left_join(label_dict_df |> select(.key, term_label), by = ".key") |>
  mutate(term_label = if_else(is.na(term_label), term_label_raw, term_label)) |>
  select(-.key)

cat("\n=== Debate-term mapping check ===\n")
print(
  ame_mixed_all |>
    filter(grepl("^debate", term_label_raw, ignore.case = TRUE)) |>
    distinct(term_label_raw, term_label)
)

cat("\n=== Unmapped terms (fix label_dict if not empty) ===\n")
unmapped <- ame_mixed_all |>
  filter(term_label == term_label_raw) |>
  pull(term_label_raw) |>
  unique() |>
  sort()
if (length(unmapped) == 0) {
  cat("All terms mapped successfully.\n")
} else {
  print(unmapped)
}

ame_m1_mixed <- ame_mixed_all |> filter(model == "Model A - Overreporting among validated nonwatchers")
ame_m2_mixed <- ame_mixed_all |> filter(model == "Model B - Overreporting among self-reported watchers")

# ---- 8. Coefficient plot -----------------------------------------------------

term_order <- c(
  "PID (somewhat close or very close)",
  "Political interest (just a little)",
  "Political interest (quite / very interested)",
  "Voting is duty",
  "Speeding (x<2.decile)",
  "Straightlining (1+)",
  "Age",
  "Education (tertiary)",
  "Sex (Male)",
  "Total TV consumption (mins/day)",
  "2021-NOVA",
  "2023-NOVA1",
  "2023-CT",
  "2023-PRIMA",
  "2023-NOVA2"
)

ame_plot <- ame_mixed_all |>
  filter(term_label %in% term_order) |>
  mutate(term_label = factor(term_label, levels = term_order))

# geom_pointrange() silently drops rows with NA estimate/CI (only a
# console warning) - flagged explicitly here so a missing point is
# never a silent surprise.
na_rows <- ame_plot |> filter(is.na(estimate) | is.na(conf.low) | is.na(conf.high))
if (nrow(na_rows) > 0) {
  cat("\n=== WARNING: terms dropped from the plot due to missing estimate/CI ===\n")
  print(na_rows |> select(model, term_label, estimate, conf.low, conf.high))
}

model_colors <- c(
  "Model A - Overreporting among validated nonwatchers" = "#1b6ca8",
  "Model B - Overreporting among self-reported watchers" = "#c0392b"
)
model_shapes <- c(
  "Model A - Overreporting among validated nonwatchers" = 16,
  "Model B - Overreporting among self-reported watchers" = 17
)

plot_ame <- ggplot(
  ame_plot,
  aes(x = fct_rev(term_label), y = estimate, color = model, shape = model)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_pointrange(
    aes(ymin = conf.low, ymax = conf.high),
    position = position_dodge(width = 0.5),
    size = 0.55,
    linewidth = 0.7
  ) +
  geom_text(
    aes(label = sprintf("%.2f", estimate), y = conf.high + 0.008),
    position = position_dodge(width = 0.5),
    size = 2.3,
    hjust = 0,
    show.legend = FALSE
  ) +
  coord_flip() +
  scale_color_manual(
    values = model_colors,
    guide = guide_legend(title = NULL)
  ) +
  scale_shape_manual(
    values = model_shapes,
    guide = guide_legend(title = NULL)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.28))) +
  labs(
    x = NULL,
    y = "AME",
    color = NULL,
    shape = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.location = "plot",
    legend.margin = margin(0, 0, 2, 0)
  )

plot_ame

ggsave(
  "AME_mixed_ModelA_vs_ModelB_pooled.png",
  plot_ame,
  width = 9,
  height = 7,
  dpi = 300
)

# ---- 9. LaTeX table export --------------------------------------------------
# AME table (estimate + 95% CI + significance stars) in the target
# Overleaf format, plus a model-statistics block (N, respondents, ICC,
# BIC, AUC).

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
  if (is.na(row$estimate[1]) || is.na(row$conf.low[1]) || is.na(row$conf.high[1])) {
    # "--" means the predictor isn't in this model at all; "n/a" means
    # it IS in the model but the estimate/CI came back missing (usually
    # too few observations or quasi-complete separation for that level).
    message("  NOTE: '", label, "' has a missing estimate/CI - shown as 'n/a'. ",
            "Likely cause: too few observations or quasi-complete separation ",
            "for this factor level within this model.")
    return("n/a")
  }
  est   <- formatC(row$estimate[1],  digits = digits, format = "f")
  lo    <- formatC(row$conf.low[1],  digits = digits, format = "f")
  hi    <- formatC(row$conf.high[1], digits = digits, format = "f")
  stars <- sig_stars(row$p.value[1])
  paste0(est, "$", stars, "$ [", lo, ", ", hi, "]")
}

row_defs <- tibble::tibble(
  key   = term_order,
  latex = c(
    "PID (somewhat close or very close)",
    "Political interest (just a little)",
    "Political interest (quite / very interested)",
    "Voting is duty",
    "Speeding ($<$2nd decile)",
    "Straightlining (1+)",
    "Age",
    "Education (tertiary)",
    "Sex (Male)",
    "Total TV consumption (mins/day)",
    "Debate (2021-NOVA)",
    "Debate (2023-NOVA1)",
    "Debate (2023-CT)",
    "Debate (2023-PRIMA)",
    "Debate (2023-NOVA2)"
  )
)

ame_rows <- purrr::map_chr(seq_len(nrow(row_defs)), function(i) {
  key <- row_defs$key[i]
  lbl <- row_defs$latex[i]
  cellA <- fmt_ame_cell(ame_m1_mixed, key)
  cellB <- fmt_ame_cell(ame_m2_mixed, key)
  paste(lbl, "\n&\n", cellA, "\n&\n", cellB, "\n\\\\")
})

n_obs_A  <- nobs(m1_mixed)
n_obs_B  <- nobs(m2_mixed)
n_resp_A <- dplyr::n_distinct(panel_bin_A$RADIOMETER_ID2)
n_resp_B <- dplyr::n_distinct(panel_bin_B$RADIOMETER_ID2)
bic_A    <- round(BIC(m1_mixed), 1)
bic_B    <- round(BIC(m2_mixed), 1)

stat_rows <- c(
  paste("Observations\n&\n", n_obs_A,  "\n&\n", n_obs_B,  "\n\\\\"),
  paste("Respondents\n&\n",  n_resp_A, "\n&\n", n_resp_B, "\n\\\\"),
  paste("ICC\n&\n",          icc_A,    "\n&\n", icc_B,    "\n\\\\"),
  paste("BIC\n&\n",          bic_A,    "\n&\n", bic_B,    "\n\\\\"),
  paste("AUC\n&\n",          auc_A,    "\n&\n", auc_B,    "\n\\\\")
)

latex_table <- c(
  "\\begin{table}[H]",
  "\\centering",
  "\\caption{Main specification (3 mins) - Average marginal effects on the probability of overreporting TV debate watching}",
  "\\label{tab:ame_main}",
  "\\scriptsize",
  "\\setlength{\\tabcolsep}{6pt}",
  "\\resizebox{\\textwidth}{!}{",
  "\\begin{tabular}{lcc}",
  "\\toprule",
  "&",
  "Validated non-watchers &",
  "Self-reported watchers \\\\",
  "\\midrule",
  "",
  paste(ame_rows, collapse = "\n\n"),
  "",
  "\\midrule",
  "",
  paste(stat_rows, collapse = "\n\n"),
  "",
  "\\bottomrule",
  "\\end{tabular}",
  "}",
  "",
  "\\vspace{0.3cm}",
  "\\begin{minipage}{0.97\\textwidth}",
  "\\footnotesize",
  "\\textit{Notes.}",
  "Entries report average marginal effects (AME) from mixed-effects logistic regression models.",
  "Values in brackets are 95\\% confidence intervals.",
  "$^{*}p<0.05$, $^{**}p<0.01$, $^{***}p<0.001$.",
  "\\end{minipage}",
  "\\end{table}"
)

cat(latex_table, sep = "\n")
writeLines(latex_table, "ame_main_table_pooled.tex")

# ---- 10. Additional diagnostics ----------------------------------------------

cat("\n=== VIF - Model A ===\n")
print(performance::check_collinearity(m1_mixed))
cat("\n=== VIF - Model B ===\n")
print(performance::check_collinearity(m2_mixed))

cat("\n=== R2 (marginal / conditional) ===\n")
print(performance::r2_nakagawa(m1_mixed))
print(performance::r2_nakagawa(m2_mixed))

cat("\n=== Overdispersion ===\n")
print(performance::check_overdispersion(m1_mixed))
print(performance::check_overdispersion(m2_mixed))

# Calibration: mean predicted probability vs. observed proportion,
# grouped into deciles of predicted risk.
make_calibration_plot <- function(observed, predicted, title) {
  df <- tibble::tibble(observed = observed, predicted = predicted) |>
    dplyr::mutate(decile = ntile(predicted, 10)) |>
    dplyr::group_by(decile) |>
    dplyr::summarise(
      mean_pred = mean(predicted),
      mean_obs  = mean(observed),
      n = dplyr::n(),
      .groups = "drop"
    )
  
  ggplot2::ggplot(df, ggplot2::aes(x = mean_pred, y = mean_obs)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(size = 2, color = "steelblue") +
    ggplot2::geom_line(color = "steelblue") +
    ggplot2::labs(
      title = title,
      x = "Mean predicted probability (decile)",
      y = "Observed proportion"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}

calib_A <- make_calibration_plot(panel_bin_A$misreport_bin, pred_A, "Calibration - Model A")
calib_B <- make_calibration_plot(panel_bin_B$misreport_bin, pred_B, "Calibration - Model B")

print(calib_A)
print(calib_B)

ggsave("calibration_ModelA_pooled.png", calib_A, width = 6, height = 5, dpi = 300)
ggsave("calibration_ModelB_pooled.png", calib_B, width = 6, height = 5, dpi = 300)

re_summary <- tibble::tibble(
  Model = c("A (validated nonwatchers)", "B (self-report watchers)"),
  `Var(intercept)` = c(var_u_A, var_u_B),
  `SD(intercept)`  = c(sqrt(var_u_A), sqrt(var_u_B)),
  ICC              = c(icc_A, icc_B)
)
print(re_summary)
writexl::write_xlsx(re_summary, "random_effects_summary_pooled.xlsx")

# =========================================================================
# END OF SCRIPT
# Outputs produced:
#   - plot_ame                              (AME coefficient plot, in-memory + PNG)
#   - ame_main_table_pooled.tex             (LaTeX table, ready for \input{})
#   - calibration_ModelA/B_pooled.png       (calibration diagnostic plots)
#   - random_effects_summary_pooled.xlsx    (variance/ICC summary for appendix)
#   - console output: convergence, ICC, LRT, AUC, VIF, R2, overdispersion
# =========================================================================