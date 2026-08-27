# =========================================================================
# ROBUSTNESS CHECK (10% threshold), adapted for data_pooled
# Mixed-effects logistic regression models of overreporting TV debate
# watching, average marginal effects (AME), LaTeX table export, and
# model diagnostics.
#
# This mirrors the main-specification script (misreport_4k, 180-second
# threshold) but uses misreport_4k_10pct / misreport_bin_10pct instead,
# i.e. overreporting/underreporting is defined against a threshold of 10%
# of the debate's duration passively watched (deb_passive_p >= 0.10)
# rather than the fixed 3-minute (180-second) threshold used in the main
# specification. All object names carry a "_10pct" suffix so this script
# can be run in the same session as the main-specification and the 0sec
# robustness-check scripts without overwriting their objects.
# To create the misreport_4k variable, run the Rscript03.
#
# Two models are estimated:
#   Model A ("Validated non-watchers"):
#       overreporters (misreport_4k_10pct == 1) vs. validated non-watchers
#       (misreport_4k_10pct == 3)
#       overreporting.
#   Model B ("Self-reported watchers"):
#       overreporters (misreport_4k_10pct == 1) vs. validated watchers
#       (misreport_4k_10pct == 2)
# Both models are respondent-clustered random-intercept logistic
# regressions (glmmTMB), with debate as a fixed effect, estimated on the
# long-format panel (one row per respondent per debate).
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
# panel_bin_A: overreporters (1) vs. validated non-watchers (3)
# panel_bin_B: overreporters (1) vs. validated watchers (2)
# misreport_bin_10pct is coded 1 = overreporter, 0 = reference category.
# All predictors that come in as labelled/character are cast to factor so
# that glmmTMB treats them categorically; age and total_tv is kept numeric.
# =========================================================================

panel_bin_A <- data_pooled |>
  filter(misreport_4k_10pct %in% c(1, 3)) |>
  mutate(
    misreport_bin_10pct = ifelse(misreport_4k_10pct == 1, 1, 0),
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

panel_bin_B <- data_pooled |>
  filter(misreport_4k_10pct %in% c(1, 2)) |>
  mutate(
    misreport_bin_10pct = ifelse(misreport_4k_10pct == 1, 1, 0),
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

# =========================================================================
# 2. MODELS
# =========================================================================
# Random intercept for RADIOMETER_ID2 accounts for the panel structure:
# respondents contribute one row per debate they were asked about, so
# observations are not independent within a respondent. Fixed effects for
# debate absorb debate-specific differences in overall misreporting rates.
# =========================================================================

m1_mixed_10pct <- glmmTMB(
  misreport_bin_10pct ~ party_id_01 + duty + polint + straightlining +
    speeding + age + edu + sex + total_tv + debate +
    (1 | RADIOMETER_ID2),
  data   = panel_bin_A,
  family = binomial(link = "logit")
)

m2_mixed_10pct <- glmmTMB(
  misreport_bin_10pct ~ party_id_01 + duty + polint + straightlining +
    speeding + age + edu + sex + total_tv + debate +
    (1 | RADIOMETER_ID2),
  data   = panel_bin_B,
  family = binomial(link = "logit")
)

summary(m1_mixed_10pct)
summary(m2_mixed_10pct)

# =========================================================================
# 3. CONVERGENCE CHECKS
# =========================================================================

cat("\n=== CONVERGENCE (10pct robustness check) ===\n")
cat("Model A - pdHess:", m1_mixed_10pct$sdr$pdHess,
    ifelse(m1_mixed_10pct$sdr$pdHess, "[OK]", "[PROBLEM]"), "\n")
cat("Model B - pdHess:", m2_mixed_10pct$sdr$pdHess,
    ifelse(m2_mixed_10pct$sdr$pdHess, "[OK]", "[PROBLEM]"), "\n")

check_convergence(m1_mixed_10pct)
check_convergence(m2_mixed_10pct)

# =========================================================================
# 4. INTRACLASS CORRELATION (ICC)
# =========================================================================
# Share of total variance (on the latent logit scale) attributable to
# between-respondent differences, i.e. how much of the variation in
# overreporting is a stable respondent trait vs. debate-to-debate noise.
# =========================================================================

cat("\n=== ICC (10pct robustness check) ===\n")
print(icc(m1_mixed_10pct))
print(icc(m2_mixed_10pct))

var_u_A_10pct <- as.numeric(VarCorr(m1_mixed_10pct)$cond$RADIOMETER_ID2[1])
var_u_B_10pct <- as.numeric(VarCorr(m2_mixed_10pct)$cond$RADIOMETER_ID2[1])
icc_A_10pct <- round(var_u_A_10pct / (var_u_A_10pct + pi^2 / 3), 3)
icc_B_10pct <- round(var_u_B_10pct / (var_u_B_10pct + pi^2 / 3), 3)
cat("ICC Model A:", icc_A_10pct, "\n")
cat("ICC Model B:", icc_B_10pct, "\n")

# =========================================================================
# 5. RANDOM-EFFECTS DIAGNOSTICS
# =========================================================================
# Histogram + Q-Q plot of the estimated random intercepts (Model A), plus
# a likelihood-ratio test of whether the random intercept is needed at
# all (boundary test: mixture chi-sq with 0.5 weight on chi-sq_0).
# =========================================================================

re_A_10pct <- ranef(m1_mixed_10pct)$cond$RADIOMETER_ID2[, 1]
hist(re_A_10pct, breaks = 40,
     main = "Distribution of random intercepts - Model A (10pct check)",
     xlab = "Random intercept (logit scale)",
     col  = "steelblue", border = "white")
qqnorm(re_A_10pct, main = "Q-Q plot - Model A (10pct check)")
qqline(re_A_10pct, col = "red", lwd = 2)

# LRT: random intercept vs. no random intercept, Model A
m_no_re_A_10pct <- glm(
  misreport_bin_10pct ~ party_id_01 + duty + polint + straightlining +
    speeding + age + edu + sex + total_tv + debate,
  data = panel_bin_A, family = binomial()
)
lrt_stat_A_10pct <- as.numeric(-2 * (logLik(m_no_re_A_10pct) - logLik(m1_mixed_10pct)))
p_lrt_A_10pct    <- 0.5 * pchisq(lrt_stat_A_10pct, df = 1, lower.tail = FALSE)
cat("\nLRT random intercept, Model A (10pct check):\n")
cat("  chi2 =", round(lrt_stat_A_10pct, 2), "\n")
cat("  p    =", format(p_lrt_A_10pct, scientific = TRUE), "\n")
print(AIC(m_no_re_A_10pct, m1_mixed_10pct))
print(BIC(m_no_re_A_10pct, m1_mixed_10pct))

# Same LRT for Model B
m_no_re_B_10pct <- glm(
  misreport_bin_10pct ~ party_id_01 + duty + polint + straightlining +
    speeding + age + edu + sex + total_tv + debate,
  data = panel_bin_B, family = binomial()
)
lrt_stat_B_10pct <- as.numeric(-2 * (logLik(m_no_re_B_10pct) - logLik(m2_mixed_10pct)))
p_lrt_B_10pct    <- 0.5 * pchisq(lrt_stat_B_10pct, df = 1, lower.tail = FALSE)
cat("\nLRT random intercept, Model B (10pct check):\n")
cat("  chi2 =", round(lrt_stat_B_10pct, 2), "\n")
cat("  p    =", format(p_lrt_B_10pct, scientific = TRUE), "\n")
print(AIC(m_no_re_B_10pct, m2_mixed_10pct))
print(BIC(m_no_re_B_10pct, m2_mixed_10pct))

# =========================================================================
# 6. PREDICTIVE PERFORMANCE (AUC)
# =========================================================================
# Predictions use re.form = NA (population-average / fixed-effects-only
# predictions), consistent with the AME below being computed on the
# response scale, marginalised over the random effect.
# =========================================================================

pred_A_10pct <- predict(m1_mixed_10pct, type = "response", re.form = NA)
pred_B_10pct <- predict(m2_mixed_10pct, type = "response", re.form = NA)

auc_A_10pct <- round(as.numeric(auc(roc(panel_bin_A$misreport_bin_10pct, pred_A_10pct, quiet = TRUE))), 3)
auc_B_10pct <- round(as.numeric(auc(roc(panel_bin_B$misreport_bin_10pct, pred_B_10pct, quiet = TRUE))), 3)
cat("AUC Model A:", auc_A_10pct, "\n")
cat("AUC Model B:", auc_B_10pct, "\n")

# =========================================================================
# 7. AVERAGE MARGINAL EFFECTS (AME)
# =========================================================================
# avg_comparisons() from marginaleffects gives AME on the response
# (probability) scale for both factor and continuous predictors, using
# the model's own variance-covariance matrix (conditional component of
# glmmTMB) for confidence intervals.
# =========================================================================

# ---- 7a. Helper: AME for a single predictor, with error handling ----
ame_one_var_10pct <- function(model, var, spec, model_label) {
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

# ---- 7c. Compute AME: Model A ----
cat("\n=== Computing AME - Model A (10pct check) ===\n")
ame_list_A_10pct <- lapply(names(var_specs), function(v) {
  cat("  ", v, "...")
  res <- ame_one_var_10pct(
    m1_mixed_10pct, v, var_specs[[v]],
    "Model A - Overreporting among validated nonwatchers (10pct)"
  )
  cat(if (!is.null(res)) " OK\n" else " FAILED\n")
  res
})
ame_m1_mixed_10pct <- bind_rows(Filter(Negate(is.null), ame_list_A_10pct))

# ---- 7d. Compute AME: Model B ----
cat("\n=== Computing AME - Model B (10pct check) ===\n")
ame_list_B_10pct <- lapply(names(var_specs), function(v) {
  cat("  ", v, "...")
  res <- ame_one_var_10pct(
    m2_mixed_10pct, v, var_specs[[v]],
    "Model B - Overreporting among self-reported watchers (10pct)"
  )
  cat(if (!is.null(res)) " OK\n" else " FAILED\n")
  res
})
ame_m2_mixed_10pct <- bind_rows(Filter(Negate(is.null), ame_list_B_10pct))

# ---- 7e. Build raw term labels and inspect what came back ----
ame_mixed_all_10pct <- bind_rows(ame_m1_mixed_10pct, ame_m2_mixed_10pct) |>
  mutate(
    term_label_raw = if_else(
      is.na(contrast) | contrast == "",
      term,
      paste0(term, ": ", contrast)
    )
  )

cat("\n=== Available term_label_raw values (10pct check, for label_dict) ===\n")
print(sort(unique(ame_mixed_all_10pct$term_label_raw)))

# ---- 7f. Human-readable label mapping ----
# Match these strings against the printed output above; adjust if your
# factor level text differs.
# NOTE (vs. the panel_all version): "survey_time" -> "speeding"; edu
# levels are now "secondary and lower" / "tertiary"; debate contrasts
# are now text ("<level> - 2021-CT") instead of numeric ("N - 1").
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

# Mapping is done via an explicit left_join on the normalized key
# (rather than named-vector indexing) - more transparent and robust; see
# the "Debate-term mapping check" diagnostic below.
label_dict_df <- tibble::tibble(
  raw_key    = names(label_dict),
  term_label = unname(label_dict)
) |>
  mutate(.key = norm_key(raw_key))

# Safety check: if two different label_dict entries normalize to the
# SAME .key, a plain left_join would silently match one raw term to
# BOTH entries (a one-to-many join), duplicating that row and assigning
# it two different (wrong) term_labels - exactly the kind of "PRIMA
# shows up as CT, and both colors get mixed up on the plot" symptom.
# Catch that here explicitly, before it can happen silently.
dup_keys <- label_dict_df |> dplyr::filter(duplicated(.key) | duplicated(.key, fromLast = TRUE))
if (nrow(dup_keys) > 0) {
  cat("\n=== WARNING: label_dict has entries that normalize to the SAME key ===\n")
  print(dup_keys)
  stop("Fix label_dict: the entries above collide after normalization and would ",
       "cause rows to be duplicated/mislabeled by the join below.")
}

ame_mixed_all_10pct <- ame_mixed_all_10pct |>
  mutate(.key = norm_key(term_label_raw)) |>
  dplyr::left_join(
    label_dict_df |> select(.key, term_label),
    by = ".key",
    relationship = "many-to-one"  # errors loudly instead of silently
    # duplicating rows if .key isn't
    # unique in label_dict_df
  ) |>
  mutate(term_label = if_else(is.na(term_label), term_label_raw, term_label)) |>
  select(-.key)

cat("\n=== Debate-term mapping check (10pct check) ===\n")
print(
  ame_mixed_all_10pct |>
    filter(grepl("^debate", term_label_raw, ignore.case = TRUE)) |>
    distinct(term_label_raw, term_label)
)

cat("\n=== Unmapped terms (10pct check, fix label_dict if not empty) ===\n")
unmapped_10pct <- ame_mixed_all_10pct |>
  filter(term_label == term_label_raw) |>
  pull(term_label_raw) |>
  unique() |>
  sort()
if (length(unmapped_10pct) == 0) {
  cat("All terms mapped successfully.\n")
} else {
  print(unmapped_10pct)
}

# Split back into per-model tibbles carrying the mapped term_label, for
# use both by the plot and by the LaTeX table builder below.
ame_m1_mixed_10pct <- ame_mixed_all_10pct |> filter(model == "Model A - Overreporting among validated nonwatchers (10pct)")
ame_m2_mixed_10pct <- ame_mixed_all_10pct |> filter(model == "Model B - Overreporting among self-reported watchers (10pct)")

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

ame_plot_10pct <- ame_mixed_all_10pct |>
  filter(term_label %in% term_order) |>
  mutate(term_label = factor(term_label, levels = term_order))

# ggplot's geom_pointrange() silently drops rows with NA in
# estimate/conf.low/conf.high - flag them explicitly here.
na_rows_10pct <- ame_plot_10pct |> filter(is.na(estimate) | is.na(conf.low) | is.na(conf.high))
if (nrow(na_rows_10pct) > 0) {
  cat("\n=== WARNING: terms dropped from the plot due to missing estimate/CI (10pct check) ===\n")
  print(na_rows_10pct |> select(model, term_label, estimate, conf.low, conf.high))
}

model_colors_10pct <- c(
  "Model A - Overreporting among validated nonwatchers (10pct)" = "#1b6ca8",
  "Model B - Overreporting among self-reported watchers (10pct)" = "#c0392b"
)
model_shapes_10pct <- c(
  "Model A - Overreporting among validated nonwatchers (10pct)" = 16,
  "Model B - Overreporting among self-reported watchers (10pct)" = 17
)

plot_ame_10pct <- ggplot(
  ame_plot_10pct,
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
    values = model_colors_10pct,
    guide = guide_legend(title = NULL)
  ) +
  scale_shape_manual(
    values = model_shapes_10pct,
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

plot_ame_10pct

ggsave(
  "AME_mixed_ModelA_vs_ModelB_robustness_10pct_pooled.png",
  plot_ame_10pct,
  width = 9,
  height = 7,
  dpi = 300
)

# =========================================================================
# 9. LATEX TABLE EXPORT
# =========================================================================
# Builds the AME table (est + 95% CI + significance stars) in the target
# Overleaf format, plus a model-statistics block (N, respondents, ICC,
# BIC, AUC).
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
    # Row exists but with a missing estimate/CI - flagged as "n/a"
    # (distinct from "--", which means the predictor isn't in this
    # model at all), and printed to console so it's never silent.
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

# Display labels; order must match term_order 1:1
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

ame_rows_10pct <- purrr::map_chr(seq_len(nrow(row_defs)), function(i) {
  key <- row_defs$key[i]
  lbl <- row_defs$latex[i]
  cellA <- fmt_ame_cell(ame_m1_mixed_10pct, key)
  cellB <- fmt_ame_cell(ame_m2_mixed_10pct, key)
  paste(lbl, "\n&\n", cellA, "\n&\n", cellB, "\n\\\\")
})

n_obs_A_10pct  <- nobs(m1_mixed_10pct)
n_obs_B_10pct  <- nobs(m2_mixed_10pct)
n_resp_A_10pct <- dplyr::n_distinct(panel_bin_A$RADIOMETER_ID2)
n_resp_B_10pct <- dplyr::n_distinct(panel_bin_B$RADIOMETER_ID2)
bic_A_10pct    <- round(BIC(m1_mixed_10pct), 1)
bic_B_10pct    <- round(BIC(m2_mixed_10pct), 1)

stat_rows_10pct <- c(
  paste("Observations\n&\n", n_obs_A_10pct,  "\n&\n", n_obs_B_10pct,  "\n\\\\"),
  paste("Respondents\n&\n",  n_resp_A_10pct, "\n&\n", n_resp_B_10pct, "\n\\\\"),
  paste("ICC\n&\n",          icc_A_10pct,    "\n&\n", icc_B_10pct,    "\n\\\\"),
  paste("BIC\n&\n",          bic_A_10pct,    "\n&\n", bic_B_10pct,    "\n\\\\"),
  paste("AUC\n&\n",          auc_A_10pct,    "\n&\n", auc_B_10pct,    "\n\\\\")
)

latex_table_10pct <- c(
  "\\begin{table}[H]",
  "\\centering",
  "\\caption{Robustness check (10\\%) - Average marginal effects on the probability of overreporting TV debate watching}",
  "\\label{tab:ame_robust_10pct}",
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
  paste(ame_rows_10pct, collapse = "\n\n"),
  "",
  "\\midrule",
  "",
  paste(stat_rows_10pct, collapse = "\n\n"),
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
  "Overreporting/underreporting is defined against a threshold of 10\\% of debate duration passively watched (robustness check to the 3-minute threshold used in the main specification).",
  "Values in brackets are 95\\% confidence intervals.",
  "$^{*}p<0.05$, $^{**}p<0.01$, $^{***}p<0.001$.",
  "\\end{minipage}",
  "\\end{table}"
)

cat(latex_table_10pct, sep = "\n")
writeLines(latex_table_10pct, "ame_robustness_10pct_pooled_table.tex")

# =========================================================================
# 10. ADDITIONAL DIAGNOSTICS
# =========================================================================

# ---- 10a. Multicollinearity (VIF) ----
cat("\n=== VIF - Model A (10pct check) ===\n")
print(performance::check_collinearity(m1_mixed_10pct))
cat("\n=== VIF - Model B (10pct check) ===\n")
print(performance::check_collinearity(m2_mixed_10pct))

# ---- 10b. Marginal / conditional R^2 (Nakagawa) ----
cat("\n=== R2 (marginal / conditional), 10pct check ===\n")
print(performance::r2_nakagawa(m1_mixed_10pct))
print(performance::r2_nakagawa(m2_mixed_10pct))

# ---- 10c. Overdispersion check ----
cat("\n=== Overdispersion, 10pct check ===\n")
print(performance::check_overdispersion(m1_mixed_10pct))
print(performance::check_overdispersion(m2_mixed_10pct))

# ---- 10d. Calibration plots: predicted vs. observed probability (deciles) ----
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

calib_A_10pct <- make_calibration_plot(panel_bin_A$misreport_bin_10pct, pred_A_10pct, "Calibration - Model A (10pct check)")
calib_B_10pct <- make_calibration_plot(panel_bin_B$misreport_bin_10pct, pred_B_10pct, "Calibration - Model B (10pct check)")

print(calib_A_10pct)
print(calib_B_10pct)

ggsave("calibration_ModelA_10pct_pooled.png", calib_A_10pct, width = 6, height = 5, dpi = 300)
ggsave("calibration_ModelB_10pct_pooled.png", calib_B_10pct, width = 6, height = 5, dpi = 300)

# ---- 10e. Random-effects variance summary (for appendix) ----
re_summary_10pct <- tibble::tibble(
  Model = c("A (validated nonwatchers, 10pct)", "B (self-report watchers, 10pct)"),
  `Var(intercept)` = c(var_u_A_10pct, var_u_B_10pct),
  `SD(intercept)`  = c(sqrt(var_u_A_10pct), sqrt(var_u_B_10pct)),
  ICC              = c(icc_A_10pct, icc_B_10pct)
)
print(re_summary_10pct)
writexl::write_xlsx(re_summary_10pct, "random_effects_summary_10pct_pooled.xlsx")

# =========================================================================
# END OF SCRIPT
# Outputs produced:
#   - plot_ame_10pct                               (AME coefficient plot, in-memory + PNG)
#   - ame_robustness_10pct_pooled_table.tex        (LaTeX table, ready for \input{})
#   - calibration_ModelA/B_10pct_pooled.png        (calibration diagnostic plots)
#   - random_effects_summary_10pct_pooled.xlsx     (variance/ICC summary for appendix)
#   - console output: convergence, ICC, LRT, AUC, VIF, R2, overdispersion
# =========================================================================