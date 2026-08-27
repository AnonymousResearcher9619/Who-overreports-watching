# =========================================================================
# ROBUSTNESS CHECK: excluding TV non-watchers (total_tv == 0)
# Same outcome definition and predictors as the main specification
# (misreport_4k, 3-min / 180-second threshold), but estimated on a
# restricted sample that excludes respondents with total_tv == 0, i.e.
# people who report watching no politics in TV at all on a typical day.
# All object names carry an "_exclTV0" suffix so this script can be run 
# in the same session as the main-specification script (and the other robustness
# checks) without overwriting their objects.
#
#   Model A ("Validated non-watchers"):
#       overreporters (misreport_4k == 1) vs. validated non-watchers
#       (misreport_4k == 3), excl. total_tv == 0.
#   Model B ("Self-reported watchers"):
#       overreporters (misreport_4k == 1) vs. validated watchers
#       (misreport_4k == 2), excl. total_tv == 0.
#
# Both models are respondent-clustered random-intercept logistic
# regressions (glmmTMB), with debate as a fixed effect, estimated on the
# long-format panel (one row per respondent per debate). total_tv is
# still included as a continuous predictor.
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

# =========================================================================
# 1. DATA PREPARATION
# =========================================================================
# panel_bin_A_exclTV0: overreporters (1) vs. validated non-watchers (3),
#                       excl. total_tv == 0
# panel_bin_B_exclTV0: overreporters (1) vs. validated watchers (2),
#                       excl. total_tv == 0
# misreport_bin_exclTV0 is coded 1 = overreporter, 0 = reference category.
# =========================================================================

panel_bin_A_exclTV0 <- data_pooled |>
  filter(misreport_4k %in% c(1, 3), total_tv != 0) |>
  mutate(
    misreport_bin_exclTV0 = ifelse(misreport_4k == 1, 1, 0),
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

panel_bin_B_exclTV0 <- data_pooled |>
  filter(misreport_4k %in% c(1, 2), total_tv != 0) |>
  mutate(
    misreport_bin_exclTV0 = ifelse(misreport_4k == 1, 1, 0),
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

cat("N panel_bin_A_exclTV0:", nrow(panel_bin_A_exclTV0),
    "(vs. full-sample equivalent would include total_tv == 0 respondents)\n")
cat("N panel_bin_B_exclTV0:", nrow(panel_bin_B_exclTV0), "\n")

# =========================================================================
# 2. MODELS
# =========================================================================
# Random intercept for RADIOMETER_ID2 accounts for the panel structure:
# respondents contribute one row per debate they were asked about, so
# observations are not independent within a respondent. Fixed effects for
# debate absorb debate-specific differences in overall misreporting rates.
# =========================================================================

m1_mixed_exclTV0 <- glmmTMB(
  misreport_bin_exclTV0 ~ party_id_01 + duty + polint + straightlining +
    speeding + age + edu + sex + total_tv + debate +
    (1 | RADIOMETER_ID2),
  data   = panel_bin_A_exclTV0,
  family = binomial(link = "logit")
)

m2_mixed_exclTV0 <- glmmTMB(
  misreport_bin_exclTV0 ~ party_id_01 + duty + polint + straightlining +
    speeding + age + edu + sex + total_tv + debate +
    (1 | RADIOMETER_ID2),
  data   = panel_bin_B_exclTV0,
  family = binomial(link = "logit")
)

summary(m1_mixed_exclTV0)
summary(m2_mixed_exclTV0)

# =========================================================================
# 3. CONVERGENCE CHECKS
# =========================================================================

cat("\n=== CONVERGENCE (excl. TV non-watchers) ===\n")
cat("Model A - pdHess:", m1_mixed_exclTV0$sdr$pdHess,
    ifelse(m1_mixed_exclTV0$sdr$pdHess, "[OK]", "[PROBLEM]"), "\n")
cat("Model B - pdHess:", m2_mixed_exclTV0$sdr$pdHess,
    ifelse(m2_mixed_exclTV0$sdr$pdHess, "[OK]", "[PROBLEM]"), "\n")

check_convergence(m1_mixed_exclTV0)
check_convergence(m2_mixed_exclTV0)

# =========================================================================
# 4. INTRACLASS CORRELATION (ICC)
# =========================================================================

cat("\n=== ICC (excl. TV non-watchers) ===\n")
print(icc(m1_mixed_exclTV0))
print(icc(m2_mixed_exclTV0))

var_u_A_exclTV0 <- as.numeric(VarCorr(m1_mixed_exclTV0)$cond$RADIOMETER_ID2[1])
var_u_B_exclTV0 <- as.numeric(VarCorr(m2_mixed_exclTV0)$cond$RADIOMETER_ID2[1])
icc_A_exclTV0 <- round(var_u_A_exclTV0 / (var_u_A_exclTV0 + pi^2 / 3), 3)
icc_B_exclTV0 <- round(var_u_B_exclTV0 / (var_u_B_exclTV0 + pi^2 / 3), 3)
cat("ICC Model A:", icc_A_exclTV0, "\n")
cat("ICC Model B:", icc_B_exclTV0, "\n")

# =========================================================================
# 5. RANDOM-EFFECTS DIAGNOSTICS
# =========================================================================

re_A_exclTV0 <- ranef(m1_mixed_exclTV0)$cond$RADIOMETER_ID2[, 1]
hist(re_A_exclTV0, breaks = 40,
     main = "Distribution of random intercepts - Model A (excl. TV non-watchers)",
     xlab = "Random intercept (logit scale)",
     col  = "steelblue", border = "white")
qqnorm(re_A_exclTV0, main = "Q-Q plot - Model A (excl. TV non-watchers)")
qqline(re_A_exclTV0, col = "red", lwd = 2)

# LRT: random intercept vs. no random intercept, Model A
m_no_re_A_exclTV0 <- glm(
  misreport_bin_exclTV0 ~ party_id_01 + duty + polint + straightlining +
    speeding + age + edu + sex + total_tv + debate,
  data = panel_bin_A_exclTV0, family = binomial()
)
lrt_stat_A_exclTV0 <- as.numeric(-2 * (logLik(m_no_re_A_exclTV0) - logLik(m1_mixed_exclTV0)))
p_lrt_A_exclTV0    <- 0.5 * pchisq(lrt_stat_A_exclTV0, df = 1, lower.tail = FALSE)
cat("\nLRT random intercept, Model A (excl. TV non-watchers):\n")
cat("  chi2 =", round(lrt_stat_A_exclTV0, 2), "\n")
cat("  p    =", format(p_lrt_A_exclTV0, scientific = TRUE), "\n")
print(AIC(m_no_re_A_exclTV0, m1_mixed_exclTV0))
print(BIC(m_no_re_A_exclTV0, m1_mixed_exclTV0))

# Same LRT for Model B
m_no_re_B_exclTV0 <- glm(
  misreport_bin_exclTV0 ~ party_id_01 + duty + polint + straightlining +
    speeding + age + edu + sex + total_tv + debate,
  data = panel_bin_B_exclTV0, family = binomial()
)
lrt_stat_B_exclTV0 <- as.numeric(-2 * (logLik(m_no_re_B_exclTV0) - logLik(m2_mixed_exclTV0)))
p_lrt_B_exclTV0    <- 0.5 * pchisq(lrt_stat_B_exclTV0, df = 1, lower.tail = FALSE)
cat("\nLRT random intercept, Model B (excl. TV non-watchers):\n")
cat("  chi2 =", round(lrt_stat_B_exclTV0, 2), "\n")
cat("  p    =", format(p_lrt_B_exclTV0, scientific = TRUE), "\n")
print(AIC(m_no_re_B_exclTV0, m2_mixed_exclTV0))
print(BIC(m_no_re_B_exclTV0, m2_mixed_exclTV0))

# =========================================================================
# 6. PREDICTIVE PERFORMANCE (AUC)
# =========================================================================

pred_A_exclTV0 <- predict(m1_mixed_exclTV0, type = "response", re.form = NA)
pred_B_exclTV0 <- predict(m2_mixed_exclTV0, type = "response", re.form = NA)

auc_A_exclTV0 <- round(as.numeric(auc(roc(panel_bin_A_exclTV0$misreport_bin_exclTV0, pred_A_exclTV0, quiet = TRUE))), 3)
auc_B_exclTV0 <- round(as.numeric(auc(roc(panel_bin_B_exclTV0$misreport_bin_exclTV0, pred_B_exclTV0, quiet = TRUE))), 3)
cat("AUC Model A:", auc_A_exclTV0, "\n")
cat("AUC Model B:", auc_B_exclTV0, "\n")


# =========================================================================
# 7. AVERAGE MARGINAL EFFECTS (AME)
# =========================================================================

# ---- 7a. Helper: AME for a single predictor, with error handling ----
ame_one_var_exclTV0 <- function(model, var, spec, model_label) {
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

# ---- 7b. Predictor specifications ----
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

# ---- 7c. Compute AME: Model A ----
cat("\n=== Computing AME - Model A (excl. TV non-watchers) ===\n")
ame_list_A_exclTV0 <- lapply(names(var_specs), function(v) {
  cat("  ", v, "...")
  res <- ame_one_var_exclTV0(
    m1_mixed_exclTV0, v, var_specs[[v]],
    "Model A - Overreporting among validated nonwatchers (excl. TV non-watchers)"
  )
  cat(if (!is.null(res)) " OK\n" else " FAILED\n")
  res
})
ame_m1_mixed_exclTV0 <- bind_rows(Filter(Negate(is.null), ame_list_A_exclTV0))

# ---- 7d. Compute AME: Model B ----
cat("\n=== Computing AME - Model B (excl. TV non-watchers) ===\n")
ame_list_B_exclTV0 <- lapply(names(var_specs), function(v) {
  cat("  ", v, "...")
  res <- ame_one_var_exclTV0(
    m2_mixed_exclTV0, v, var_specs[[v]],
    "Model B - Overreporting among self-reported watchers (excl. TV non-watchers)"
  )
  cat(if (!is.null(res)) " OK\n" else " FAILED\n")
  res
})
ame_m2_mixed_exclTV0 <- bind_rows(Filter(Negate(is.null), ame_list_B_exclTV0))

# ---- 7e. Build raw term labels and inspect what came back ----
ame_mixed_all_exclTV0 <- bind_rows(ame_m1_mixed_exclTV0, ame_m2_mixed_exclTV0) |>
  mutate(
    term_label_raw = if_else(
      is.na(contrast) | contrast == "",
      term,
      paste0(term, ": ", contrast)
    )
  )

cat("\n=== Available term_label_raw values (excl. TV non-watchers, for label_dict) ===\n")
print(sort(unique(ame_mixed_all_exclTV0$term_label_raw)))

# ---- 7f. Human-readable label mapping ----
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

# Mapping via explicit left_join on the normalized key (robust to
# duplicate/collision issues - see main-specification script notes).
label_dict_df <- tibble::tibble(
  raw_key    = names(label_dict),
  term_label = unname(label_dict)
) |>
  mutate(.key = norm_key(raw_key))

dup_keys <- label_dict_df |> dplyr::filter(duplicated(.key) | duplicated(.key, fromLast = TRUE))
if (nrow(dup_keys) > 0) {
  cat("\n=== WARNING: label_dict has entries that normalize to the SAME key ===\n")
  print(dup_keys)
  stop("Fix label_dict: the entries above collide after normalization and would ",
       "cause rows to be duplicated/mislabeled by the join below.")
}

ame_mixed_all_exclTV0 <- ame_mixed_all_exclTV0 |>
  mutate(.key = norm_key(term_label_raw)) |>
  dplyr::left_join(
    label_dict_df |> select(.key, term_label),
    by = ".key",
    relationship = "many-to-one"
  ) |>
  mutate(term_label = if_else(is.na(term_label), term_label_raw, term_label)) |>
  select(-.key)

cat("\n=== Debate-term mapping check (excl. TV non-watchers) ===\n")
print(
  ame_mixed_all_exclTV0 |>
    filter(grepl("^debate", term_label_raw, ignore.case = TRUE)) |>
    distinct(term_label_raw, term_label)
)

cat("\n=== Unmapped terms (excl. TV non-watchers, fix label_dict if not empty) ===\n")
unmapped_exclTV0 <- ame_mixed_all_exclTV0 |>
  filter(term_label == term_label_raw) |>
  pull(term_label_raw) |>
  unique() |>
  sort()
if (length(unmapped_exclTV0) == 0) {
  cat("All terms mapped successfully.\n")
} else {
  print(unmapped_exclTV0)
}

ame_m1_mixed_exclTV0 <- ame_mixed_all_exclTV0 |> filter(model == "Model A - Overreporting among validated nonwatchers (excl. TV non-watchers)")
ame_m2_mixed_exclTV0 <- ame_mixed_all_exclTV0 |> filter(model == "Model B - Overreporting among self-reported watchers (excl. TV non-watchers)")

# =========================================================================
# 8. COEFFICIENT PLOT
# =========================================================================

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

ame_plot_exclTV0 <- ame_mixed_all_exclTV0 |>
  filter(term_label %in% term_order) |>
  mutate(term_label = factor(term_label, levels = term_order))

na_rows_exclTV0 <- ame_plot_exclTV0 |> filter(is.na(estimate) | is.na(conf.low) | is.na(conf.high))
if (nrow(na_rows_exclTV0) > 0) {
  cat("\n=== WARNING: terms dropped from the plot due to missing estimate/CI (excl. TV non-watchers) ===\n")
  print(na_rows_exclTV0 |> select(model, term_label, estimate, conf.low, conf.high))
}

model_colors_exclTV0 <- c(
  "Model A - Overreporting among validated nonwatchers (excl. TV non-watchers)" = "#1b6ca8",
  "Model B - Overreporting among self-reported watchers (excl. TV non-watchers)" = "#c0392b"
)
model_shapes_exclTV0 <- c(
  "Model A - Overreporting among validated nonwatchers (excl. TV non-watchers)" = 16,
  "Model B - Overreporting among self-reported watchers (excl. TV non-watchers)" = 17
)

plot_ame_exclTV0 <- ggplot(
  ame_plot_exclTV0,
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
    values = model_colors_exclTV0,
    guide = guide_legend(title = NULL)
  ) +
  scale_shape_manual(
    values = model_shapes_exclTV0,
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

plot_ame_exclTV0

ggsave(
  "AME_mixed_ModelA_vs_ModelB_robustness_exclTV0_pooled.png",
  plot_ame_exclTV0,
  width = 9,
  height = 7,
  dpi = 300
)

# =========================================================================
# 9. LATEX TABLE EXPORT
# =========================================================================

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

ame_rows_exclTV0 <- purrr::map_chr(seq_len(nrow(row_defs)), function(i) {
  key <- row_defs$key[i]
  lbl <- row_defs$latex[i]
  cellA <- fmt_ame_cell(ame_m1_mixed_exclTV0, key)
  cellB <- fmt_ame_cell(ame_m2_mixed_exclTV0, key)
  paste(lbl, "\n&\n", cellA, "\n&\n", cellB, "\n\\\\")
})

n_obs_A_exclTV0  <- nobs(m1_mixed_exclTV0)
n_obs_B_exclTV0  <- nobs(m2_mixed_exclTV0)
n_resp_A_exclTV0 <- dplyr::n_distinct(panel_bin_A_exclTV0$RADIOMETER_ID2)
n_resp_B_exclTV0 <- dplyr::n_distinct(panel_bin_B_exclTV0$RADIOMETER_ID2)
bic_A_exclTV0    <- round(BIC(m1_mixed_exclTV0), 1)
bic_B_exclTV0    <- round(BIC(m2_mixed_exclTV0), 1)

stat_rows_exclTV0 <- c(
  paste("Observations\n&\n", n_obs_A_exclTV0,  "\n&\n", n_obs_B_exclTV0,  "\n\\\\"),
  paste("Respondents\n&\n",  n_resp_A_exclTV0, "\n&\n", n_resp_B_exclTV0, "\n\\\\"),
  paste("ICC\n&\n",          icc_A_exclTV0,    "\n&\n", icc_B_exclTV0,    "\n\\\\"),
  paste("BIC\n&\n",          bic_A_exclTV0,    "\n&\n", bic_B_exclTV0,    "\n\\\\"),
  paste("AUC\n&\n",          auc_A_exclTV0,    "\n&\n", auc_B_exclTV0,    "\n\\\\")
)

latex_table_exclTV0 <- c(
  "\\begin{table}[H]",
  "\\centering",
  "\\caption{Robustness check (excl. TV non-watchers) - Average marginal effects on the probability of overreporting TV debate watching}",
  "\\label{tab:ame_robust_exclTV0}",
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
  paste(ame_rows_exclTV0, collapse = "\n\n"),
  "",
  "\\midrule",
  "",
  paste(stat_rows_exclTV0, collapse = "\n\n"),
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
  "Respondents reporting zero total TV consumption (total\\_tv == 0) are excluded from this robustness check.",
  "Values in brackets are 95\\% confidence intervals.",
  "$^{*}p<0.05$, $^{**}p<0.01$, $^{***}p<0.001$.",
  "\\end{minipage}",
  "\\end{table}"
)

cat(latex_table_exclTV0, sep = "\n")
writeLines(latex_table_exclTV0, "ame_robustness_exclTV0_pooled_table.tex")

# =========================================================================
# 10. ADDITIONAL DIAGNOSTICS
# =========================================================================

cat("\n=== VIF - Model A (excl. TV non-watchers) ===\n")
print(performance::check_collinearity(m1_mixed_exclTV0))
cat("\n=== VIF - Model B (excl. TV non-watchers) ===\n")
print(performance::check_collinearity(m2_mixed_exclTV0))

cat("\n=== R2 (marginal / conditional), excl. TV non-watchers ===\n")
print(performance::r2_nakagawa(m1_mixed_exclTV0))
print(performance::r2_nakagawa(m2_mixed_exclTV0))

cat("\n=== Overdispersion, excl. TV non-watchers ===\n")
print(performance::check_overdispersion(m1_mixed_exclTV0))
print(performance::check_overdispersion(m2_mixed_exclTV0))

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

calib_A_exclTV0 <- make_calibration_plot(panel_bin_A_exclTV0$misreport_bin_exclTV0, pred_A_exclTV0, "Calibration - Model A (excl. TV non-watchers)")
calib_B_exclTV0 <- make_calibration_plot(panel_bin_B_exclTV0$misreport_bin_exclTV0, pred_B_exclTV0, "Calibration - Model B (excl. TV non-watchers)")

print(calib_A_exclTV0)
print(calib_B_exclTV0)

ggsave("calibration_ModelA_exclTV0_pooled.png", calib_A_exclTV0, width = 6, height = 5, dpi = 300)
ggsave("calibration_ModelB_exclTV0_pooled.png", calib_B_exclTV0, width = 6, height = 5, dpi = 300)

re_summary_exclTV0 <- tibble::tibble(
  Model = c("A (validated nonwatchers, excl. TV non-watchers)", "B (self-report watchers, excl. TV non-watchers)"),
  `Var(intercept)` = c(var_u_A_exclTV0, var_u_B_exclTV0),
  `SD(intercept)`  = c(sqrt(var_u_A_exclTV0), sqrt(var_u_B_exclTV0)),
  ICC              = c(icc_A_exclTV0, icc_B_exclTV0)
)
print(re_summary_exclTV0)
writexl::write_xlsx(re_summary_exclTV0, "random_effects_summary_exclTV0_pooled.xlsx")

# =========================================================================
# END OF SCRIPT
# Outputs produced:
#   - plot_ame_exclTV0                              (AME coefficient plot, in-memory + PNG)
#   - ame_robustness_exclTV0_pooled_table.tex       (LaTeX table, ready for \input{})
#   - calibration_ModelA/B_exclTV0_pooled.png       (calibration diagnostic plots)
#   - random_effects_summary_exclTV0_pooled.xlsx    (variance/ICC summary for appendix)
#   - console output: convergence, ICC, LRT, AUC, VIF, R2, overdispersion
# =========================================================================