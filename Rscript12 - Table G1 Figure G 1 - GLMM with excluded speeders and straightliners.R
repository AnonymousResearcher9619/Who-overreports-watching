# =========================================================================
# ROBUSTNESS CHECK: Models C and D
# Same outcome definition as the main specification (misreport_4k, 3-min /
# 180-second threshold), but estimated on a restricted sample that
# excludes speeders (speeding == "x < 2.decile") and straightliners
# (straightlining == 1).
#
#   Model A: Overreporting among passive non-watchers (full sample)
#   Model B: Overreporting among self-reported watchers (full sample)
#   Model C: Overreporting among passive non-watchers (excl. speeders &
#             straightliners)
#   Model D: Overreporting among self-reported watchers (excl. speeders &
#             straightliners)
#
# LABEL HANDLING: the four model labels (label_A..label_D, defined in
# Section 0 below) are the SINGLE SOURCE OF TRUTH for how each model is
# identified throughout this script - when the AME rows are created,
# and again when they are colored/shaped in the plot. Retyping the same
# string in multiple places is exactly how a plotted model can silently
# lose its legend entry/color/shape (an invisible character difference -
# e.g. a typographic dash "-" vs "--", a double space, a stray character
# from copy-paste - makes scale_color_manual()/scale_shape_manual() fail
# to match, and ggplot silently drops that model's points rather than
# erroring). Defining the labels once and reusing the variables
# everywhere makes that entire class of bug impossible.
#
# This script assumes Models A and B (m1_mixed, m2_mixed) and their AME
# results (ame_m1_mixed, ame_m2_mixed) already exist in the environment,
# e.g. from running the main-specification script first (Rscript04). If
# they don't, Section 0 re-estimates them on the full sample so this
# script can also be run standalone - using the same label_A/label_B
# constants, so labels always match regardless of which path is taken.
# =========================================================================

library(glmmTMB)
library(performance)
library(marginaleffects)
library(ggplot2)
library(forcats)
library(dplyr)
library(writexl)
library(pROC)
library(tibble)
library(purrr)

# =========================================================================
# 0. MODEL LABELS (single source of truth)
# =========================================================================

label_A <- "Model A - Overreporting among passive non-watchers"
label_B <- "Model B - Overreporting among self-reported watchers"
label_C <- "Model C - Passive non-watchers (excl. speeders & straightliners)"
label_D <- "Model D - Self-reported watchers (excl. speeders & straightliners)"

need_ab <- !all(c("m1_mixed", "m2_mixed", "ame_m1_mixed", "ame_m2_mixed") %in% ls())

if (need_ab) {
  
  message("Models A/B not found in environment - estimating them now.")
  
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
      total_tv       = as.numeric(total_tv)
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
  
  m1_mixed <- glmmTMB(
    misreport_bin ~ party_id_01 + duty + polint + straightlining +
      speeding + age + edu + sex + total_tv + debate +
      (1 | RADIOMETER_ID2),
    data = panel_bin_A, family = binomial(link = "logit")
  )
  
  m2_mixed <- glmmTMB(
    misreport_bin ~ party_id_01 + duty + polint + straightlining +
      speeding + age + edu + sex + total_tv + debate +
      (1 | RADIOMETER_ID2),
    data = panel_bin_B, family = binomial(link = "logit")
  )
  
  ame_one_var_ab <- function(model, var, spec, model_label) {
    tryCatch({
      avg_comparisons(
        model, type = "response", vcov = vcov(model)$cond,
        variables = setNames(list(spec), var)
      ) |> as_tibble() |> mutate(model = model_label)
    }, error = function(e) {
      message("  ERROR for '", var, "': ", conditionMessage(e))
      NULL
    })
  }
  
  var_specs_ab <- list(
    party_id_01 = "reference", duty = "reference", polint = "reference",
    straightlining = "reference", speeding = "reference",
    edu = "reference", sex = "reference", debate = "reference",
    age = 1, total_tv = 1
  )
  
  ame_m1_mixed <- bind_rows(Filter(Negate(is.null), lapply(
    names(var_specs_ab), function(v) ame_one_var_ab(m1_mixed, v, var_specs_ab[[v]], label_A)
  )))
  ame_m2_mixed <- bind_rows(Filter(Negate(is.null), lapply(
    names(var_specs_ab), function(v) ame_one_var_ab(m2_mixed, v, var_specs_ab[[v]], label_B)
  )))
}

# =========================================================================
# 1. DATA PREPARATION: RESTRICTED SAMPLE (EXCL. SPEEDERS & STRAIGHTLINERS)
# =========================================================================

panel_excl <- data_pooled |>
  filter(!(speeding == 2 | straightlining == 1)) |>
  mutate(
    party_id_01    = factor(party_id_01),
    edu            = factor(edu),
    sex            = factor(sex),
    speeding    = factor(speeding),
    debate         = factor(debate),
    duty           = factor(duty),
    polint         = factor(polint),
    straightlining = factor(straightlining),
    total_tv       = as.numeric(total_tv)
  )

# Sanity check: speeding and straightlining should now be constant
cat("=== Constancy check after exclusion ===\n")
cat("speeding levels in panel_excl:\n")
print(table(panel_excl$speeding, useNA = "ifany"))
cat("straightlining levels in panel_excl:\n")
print(table(panel_excl$straightlining, useNA = "ifany"))
cat("N panel_excl:", nrow(panel_excl), "\n")

panel_c <- panel_excl |>
  filter(misreport_4k %in% c(1, 3)) |>
  mutate(misreport_bin = ifelse(misreport_4k == 1, 1, 0))

panel_d <- panel_excl |>
  filter(misreport_4k %in% c(1, 2)) |>
  mutate(misreport_bin = ifelse(misreport_4k == 1, 1, 0))

cat("N panel_c (1 vs 3, excl.):", nrow(panel_c), "\n")
cat("N panel_d (1 vs 2, excl.):", nrow(panel_d), "\n")

# =========================================================================
# 2. MODELS C AND D
# =========================================================================
# speeding and straightlining are dropped from the formula (constant
# in the restricted sample). Otherwise identical specification to
# Models A/B: respondent random intercept, debate fixed effects.
# =========================================================================

m3_mixed <- glmmTMB(
  misreport_bin ~ party_id_01 + duty + polint +
    age + edu + sex + total_tv + debate +
    (1 | RADIOMETER_ID2),
  data   = panel_c,
  family = binomial(link = "logit")
)

m4_mixed <- glmmTMB(
  misreport_bin ~ party_id_01 + duty + polint +
    age + edu + sex + total_tv + debate +
    (1 | RADIOMETER_ID2),
  data   = panel_d,
  family = binomial(link = "logit")
)

summary(m3_mixed)
summary(m4_mixed)

# =========================================================================
# 3. CONVERGENCE
# =========================================================================

cat("\n=== CONVERGENCE ===\n")
cat("Model C - pdHess:", m3_mixed$sdr$pdHess,
    ifelse(m3_mixed$sdr$pdHess, "[OK]", "[PROBLEM]"), "\n")
cat("Model D - pdHess:", m4_mixed$sdr$pdHess,
    ifelse(m4_mixed$sdr$pdHess, "[OK]", "[PROBLEM]"), "\n")

grad_C <- m3_mixed$sdr$gradient.fixed
grad_D <- m4_mixed$sdr$gradient.fixed
cat("Model C - max|gradient|:", round(max(abs(grad_C)), 6),
    ifelse(max(abs(grad_C)) < 0.001, "[OK]", "[CHECK]"), "\n")
cat("Model D - max|gradient|:", round(max(abs(grad_D)), 6),
    ifelse(max(abs(grad_D)) < 0.001, "[OK]", "[CHECK]"), "\n")

check_convergence(m3_mixed)
check_convergence(m4_mixed)

# =========================================================================
# 4. ICC
# =========================================================================

cat("\n=== ICC ===\n")
print(icc(m3_mixed))
print(icc(m4_mixed))

var_u_C <- as.numeric(VarCorr(m3_mixed)$cond$RADIOMETER_ID2[1])
var_u_D <- as.numeric(VarCorr(m4_mixed)$cond$RADIOMETER_ID2[1])
icc_C <- round(var_u_C / (var_u_C + pi^2 / 3), 3)
icc_D <- round(var_u_D / (var_u_D + pi^2 / 3), 3)
cat("ICC Model C:", icc_C, "\n")
cat("ICC Model D:", icc_D, "\n")

# =========================================================================
# 5. LRT (RANDOM INTERCEPT) + AUC
# =========================================================================

m_no_re_C <- glm(
  misreport_bin ~ party_id_01 + duty + polint +
    age + edu + sex + total_tv + debate,
  data = panel_c, family = binomial()
)
m_no_re_D <- glm(
  misreport_bin ~ party_id_01 + duty + polint +
    age + edu + sex + total_tv + debate,
  data = panel_d, family = binomial()
)

lrt_C <- as.numeric(-2 * (logLik(m_no_re_C) - logLik(m3_mixed)))
lrt_D <- as.numeric(-2 * (logLik(m_no_re_D) - logLik(m4_mixed)))
cat("\nLRT Model C: chi2 =", round(lrt_C, 2),
    "p =", format(0.5 * pchisq(lrt_C, 1, lower.tail = FALSE), scientific = TRUE), "\n")
cat("LRT Model D: chi2 =", round(lrt_D, 2),
    "p =", format(0.5 * pchisq(lrt_D, 1, lower.tail = FALSE), scientific = TRUE), "\n")

cat("\nAIC/BIC:\n")
print(AIC(m_no_re_C, m3_mixed))
print(AIC(m_no_re_D, m4_mixed))
bic_C <- round(BIC(m3_mixed), 1)
bic_D <- round(BIC(m4_mixed), 1)

pred_C <- predict(m3_mixed, type = "response", re.form = NA)
pred_D <- predict(m4_mixed, type = "response", re.form = NA)
auc_C <- round(as.numeric(auc(roc(panel_c$misreport_bin, pred_C, quiet = TRUE))), 3)
auc_D <- round(as.numeric(auc(roc(panel_d$misreport_bin, pred_D, quiet = TRUE))), 3)
cat("\nAUC Model C:", auc_C, "\n")
cat("AUC Model D:", auc_D, "\n")

# =========================================================================
# 6. AME - MODELS C AND D
# =========================================================================
# Same predictors as A/B minus straightlining and speeding (dropped
# because they are constant in the restricted sample).
# =========================================================================

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
    message("  ERROR for '", var, "': ", conditionMessage(e))
    NULL
  })
}

var_specs_excl <- list(
  party_id_01 = "reference",
  duty        = "reference",
  polint      = "reference",
  edu         = "reference",
  sex         = "reference",
  debate      = "reference",
  age         = 1,
  total_tv    = 1
)

cat("\n=== AME Model C ===\n")
ame_m3 <- bind_rows(Filter(Negate(is.null), lapply(
  names(var_specs_excl), function(v) {
    cat("  ", v, "...")
    res <- ame_one_var(m3_mixed, v, var_specs_excl[[v]], label_C)
    cat(if (!is.null(res)) " OK\n" else " FAILED\n"); res
  }
)))

cat("\n=== AME Model D ===\n")
ame_m4 <- bind_rows(Filter(Negate(is.null), lapply(
  names(var_specs_excl), function(v) {
    cat("  ", v, "...")
    res <- ame_one_var(m4_mixed, v, var_specs_excl[[v]], label_D)
    cat(if (!is.null(res)) " OK\n" else " FAILED\n"); res
  }
)))

# =========================================================================
# 7. COMBINE ALL FOUR MODELS AND ADD DISPLAY LABELS
# =========================================================================

norm_key <- function(x) {
  x <- as.character(x)
  x <- tolower(trimws(x))
  gsub("[^a-z0-9]", "", x)
}

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

dict_keys <- norm_key(names(label_dict))

if (anyDuplicated(dict_keys)) {
  dup_keys <- unique(dict_keys[duplicated(dict_keys)])
  stop(
    "Collision in label_dict after norm_key(): ",
    paste(dup_keys, collapse = ", ")
  )
}

dict_lookup <- setNames(unname(label_dict), dict_keys)

add_labels <- function(df) {
  df |>
    mutate(
      term_label_raw = if_else(
        is.na(contrast) | contrast == "",
        term, paste0(term, ": ", contrast)
      ),
      .key       = norm_key(term_label_raw),
      term_label = if_else(.key %in% names(dict_lookup),
                           dict_lookup[.key], term_label_raw)
    ) |>
    select(-.key)
}

ame_all_4 <- bind_rows(
  add_labels(ame_m1_mixed),
  add_labels(ame_m2_mixed),
  add_labels(ame_m3),
  add_labels(ame_m4)
)

cat("\n=== Unmapped terms (fix label_dict if not empty) ===\n")
unmapped <- ame_all_4 |>
  filter(term_label == term_label_raw) |>
  pull(term_label_raw) |> unique() |> sort()
if (length(unmapped) == 0) cat("All terms mapped successfully.\n") else print(unmapped)

# Safety check: every model label actually present in ame_all_4$model
# should be one of label_A..label_D. If not, something upstream (e.g. a
# stale ame_m1_mixed/ame_m2_mixed loaded from an old run) still carries
# an outdated label string, and would silently lose its color/shape in
# the plot below - catch that here instead.
model_values_present <- unique(ame_all_4$model)
expected_labels <- c(label_A, label_B, label_C, label_D)
stray_labels <- setdiff(model_values_present, expected_labels)
if (length(stray_labels) > 0) {
  cat("\n=== WARNING: unexpected model label(s) found in ame_all_4$model ===\n")
  print(stray_labels)
  cat("These will not get a color/shape/legend entry in the plot below.\n",
      "Likely cause: ame_m1_mixed/ame_m2_mixed were created by an earlier\n",
      "run with a different label - re-run Section 0 (or the main-\n",
      "specification script) to regenerate them with label_A/label_B.\n")
}

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

ame_plot_df <- ame_all_4 |>
  filter(term_label %in% term_order) |>
  mutate(
    term_label  = factor(term_label, levels = term_order),
    model_group = case_when(
      model == label_A ~ "A",
      model == label_B ~ "B",
      model == label_C ~ "C",
      model == label_D ~ "D"
    )
  )

# =========================================================================
# 8. AME PLOT (4 MODELS: A, B, C, D)
# =========================================================================
# Designed to remain readable in black-and-white printing: each of the
# four models gets its own POINT SHAPE (circle / triangle / square /
# diamond) in addition to colour, so models are distinguishable by shape
# alone even if colour doesn't render. Filled vs. hollow shapes further
# separate full-sample (A, B: solid) from restricted-sample (C, D:
# hollow) specifications. Point-estimate labels are rounded to 2 decimal
# places.
#
# All keys below are the label_A..label_D variables (not retyped
# strings), so they are guaranteed to match whatever ame_all_4$model
# actually contains - see the Section 0/7 comments on why this matters.
# =========================================================================

model_labels_map <- setNames(
  c(
    "A: Passive non-watchers (full sample)",
    "B: Self-reported watchers (full sample)",
    "C: Passive non-watchers (excl. speeders & straightliners)",
    "D: Self-reported watchers (excl. speeders & straightliners)"
  ),
  c(label_A, label_B, label_C, label_D)
)

model_colors <- setNames(
  c("#1b6ca8", "#c0392b", "#1b6ca8", "#c0392b"),
  c(label_A, label_B, label_C, label_D)
)

# Distinct shape per model (not just per sample), so B&W printouts still
# distinguish all four: filled circle (A), filled triangle (B), hollow
# square (C), hollow diamond (D).
model_shapes <- setNames(
  c(16, 17, 0, 5),  # filled circle, filled triangle, hollow square, hollow diamond
  c(label_A, label_B, label_C, label_D)
)

plot_ame_4 <- ggplot(
  ame_plot_df,
  aes(x = fct_rev(term_label), y = estimate, color = model, shape = model)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_pointrange(
    aes(ymin = conf.low, ymax = conf.high),
    position = position_dodge(width = 0.6), size = 0.5, linewidth = 0.65,
    fill = "white", stroke = 0.9
  ) +
  geom_text(
    aes(label = sprintf("%.2f", estimate), y = conf.high + 0.008),
    position = position_dodge(width = 0.6), size = 2.1, hjust = 0,
    show.legend = FALSE
  ) +
  coord_flip() +
  scale_color_manual(values = model_colors, labels = model_labels_map) +
  scale_shape_manual(values = model_shapes, labels = model_labels_map) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.30))) +
  labs(
    x = NULL, y = "AME", color = "Model", shape = "Model"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "top",
    legend.direction = "vertical"
  )

plot_ame_4

ggsave(
  "AME_4models_full_vs_excl.png",
  plot_ame_4,
  width = 10,
  height = 8,
  dpi = 300
)

# =========================================================================
# 9. LATEX TABLE EXPORT (4-COLUMN FORMAT: A | B | C | D)
# =========================================================================
# Cells for Straightlining/Speeding are "--" for Models C/D since those
# predictors were dropped from those specifications.
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
  est   <- formatC(row$estimate[1],  digits = digits, format = "f")
  lo    <- formatC(row$conf.low[1],  digits = digits, format = "f")
  hi    <- formatC(row$conf.high[1], digits = digits, format = "f")
  stars <- sig_stars(row$p.value[1])
  paste0(est, "$", stars, "$ [", lo, ", ", hi, "]")
}

ame_m1_lbl <- ame_all_4 |> filter(model == label_A)
ame_m2_lbl <- ame_all_4 |> filter(model == label_B)
ame_m3_lbl <- ame_all_4 |> filter(model == label_C)
ame_m4_lbl <- ame_all_4 |> filter(model == label_D)

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

ame_rows_4 <- purrr::map_chr(seq_len(nrow(row_defs)), function(i) {
  key <- row_defs$key[i]
  lbl <- row_defs$latex[i]
  cA <- fmt_ame_cell(ame_m1_lbl, key)
  cB <- fmt_ame_cell(ame_m2_lbl, key)
  cC <- fmt_ame_cell(ame_m3_lbl, key)
  cD <- fmt_ame_cell(ame_m4_lbl, key)
  paste(lbl, "\n&", cA, "\n&", cB, "\n&", cC, "\n&", cD, "\n\\\\")
})

n_obs_A  <- nobs(m1_mixed);  n_obs_B  <- nobs(m2_mixed)
n_obs_C  <- nobs(m3_mixed); n_obs_D  <- nobs(m4_mixed)
n_resp_A <- n_distinct(model.frame(m1_mixed)$RADIOMETER_ID2)
n_resp_B <- n_distinct(model.frame(m2_mixed)$RADIOMETER_ID2)
n_resp_C <- n_distinct(panel_c$RADIOMETER_ID2)
n_resp_D <- n_distinct(panel_d$RADIOMETER_ID2)

var_u_A <- as.numeric(VarCorr(m1_mixed)$cond$RADIOMETER_ID2[1])
var_u_B <- as.numeric(VarCorr(m2_mixed)$cond$RADIOMETER_ID2[1])
icc_A <- round(var_u_A / (var_u_A + pi^2 / 3), 3)
icc_B <- round(var_u_B / (var_u_B + pi^2 / 3), 3)
bic_A <- round(BIC(m1_mixed), 1)
bic_B <- round(BIC(m2_mixed), 1)

pred_A_full <- predict(m1_mixed, type = "response", re.form = NA)
pred_B_full <- predict(m2_mixed, type = "response", re.form = NA)
auc_A <- round(as.numeric(auc(roc(model.frame(m1_mixed)$misreport_bin, pred_A_full, quiet = TRUE))), 3)
auc_B <- round(as.numeric(auc(roc(model.frame(m2_mixed)$misreport_bin, pred_B_full, quiet = TRUE))), 3)

stat_rows_4 <- c(
  paste("Observations\n&", n_obs_A,  "\n&", n_obs_B,  "\n&", n_obs_C,  "\n&", n_obs_D,  "\n\\\\"),
  paste("Respondents\n&",  n_resp_A, "\n&", n_resp_B, "\n&", n_resp_C, "\n&", n_resp_D, "\n\\\\"),
  paste("ICC\n&",          icc_A,    "\n&", icc_B,    "\n&", icc_C,    "\n&", icc_D,    "\n\\\\"),
  paste("BIC\n&",          bic_A,    "\n&", bic_B,    "\n&", bic_C,    "\n&", bic_D,    "\n\\\\"),
  paste("AUC\n&",          auc_A,    "\n&", auc_B,    "\n&", auc_C,    "\n&", auc_D,    "\n\\\\")
)

latex_table_4 <- c(
  "\\begin{table}[H]",
  "\\centering",
  "\\caption{Main specification (3 mins) -- Average marginal effects on the",
  "         probability of overreporting TV debate watching}",
  "\\label{tab:ame_main}",
  "\\scriptsize",
  "\\setlength{\\tabcolsep}{5pt}",
  "",
  "\\resizebox{\\textwidth}{!}{",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  "",
  "& A & B & C & D \\\\",
  "",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  "",
  "&",
  "Passive non-watchers &",
  "Self-reported watchers &",
  "Passive non-watchers &",
  "Self-reported watchers \\\\",
  "",
  "&",
  "(Full sample) &",
  "(Full sample) &",
  "(Excl. speeders \\& straightliners) &",
  "(Excl. speeders \\& straightliners) \\\\",
  "",
  "\\midrule",
  "",
  paste(ame_rows_4, collapse = "\n\n"),
  "",
  "\\midrule",
  "",
  paste(stat_rows_4, collapse = "\n\n"),
  "",
  "\\bottomrule",
  "\\end{tabular}",
  "}",
  "",
  "\\vspace{0.3cm}",
  "",
  "\\begin{minipage}{0.97\\textwidth}",
  "\\footnotesize",
  "\\textit{Notes.}",
  "Entries report average marginal effects (AME) from mixed-effects logistic",
  "regression models with respondent-level random intercepts and debate fixed",
  "effects (six levels).",
  "Values in brackets are 95\\% confidence intervals.",
  "Models~C and~D exclude respondents identified as speeders or straightliners;",
  "survey response quality indicators are therefore omitted from those",
  "specifications (shown as ``--'').",
  "$^{*}p<0.05$, $^{**}p<0.01$, $^{***}p<0.001$.",
  "\\end{minipage}",
  "\\end{table}"
)

cat(latex_table_4, sep = "\n")
writeLines(latex_table_4, "ame_4models_table.tex")

# =========================================================================
# 10. EXPORT AME RESULTS TO EXCEL
# =========================================================================

ame_export_4 <- ame_plot_df |>
  transmute(
    model, term_label,
    AME     = round(estimate, 4),
    SE      = round(std.error, 4),
    p_value = round(p.value, 4),
    CI_low  = round(conf.low, 4),
    CI_high = round(conf.high, 4)
  )

write_xlsx(
  list(
    Model_A  = filter(ame_export_4, model == label_A),
    Model_B  = filter(ame_export_4, model == label_B),
    Model_C  = filter(ame_export_4, model == label_C),
    Model_D  = filter(ame_export_4, model == label_D),
    Combined = ame_export_4
  ),
  "AME_4models_results.xlsx"
)

# =========================================================================
# 11. ADDITIONAL DIAGNOSTICS (MODELS C & D)
# =========================================================================

cat("\n=== VIF - Model C ===\n"); print(performance::check_collinearity(m3_mixed))
cat("\n=== VIF - Model D ===\n"); print(performance::check_collinearity(m4_mixed))

cat("\n=== R2 (marginal / conditional) ===\n")
print(performance::r2_nakagawa(m3_mixed))
print(performance::r2_nakagawa(m4_mixed))

cat("\n=== Overdispersion ===\n")
print(performance::check_overdispersion(m3_mixed))
print(performance::check_overdispersion(m4_mixed))

make_calibration_plot <- function(observed, predicted, title) {
  df <- tibble::tibble(observed = observed, predicted = predicted) |>
    dplyr::mutate(decile = ntile(predicted, 10)) |>
    dplyr::group_by(decile) |>
    dplyr::summarise(mean_pred = mean(predicted), mean_obs = mean(observed),
                     n = dplyr::n(), .groups = "drop")
  
  ggplot2::ggplot(df, ggplot2::aes(x = mean_pred, y = mean_obs)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(size = 2, color = "steelblue") +
    ggplot2::geom_line(color = "steelblue") +
    ggplot2::labs(title = title, x = "Mean predicted probability (decile)", y = "Observed proportion") +
    ggplot2::theme_minimal(base_size = 12)
}

calib_C <- make_calibration_plot(panel_c$misreport_bin, pred_C, "Calibration - Model C")
calib_D <- make_calibration_plot(panel_d$misreport_bin, pred_D, "Calibration - Model D")
print(calib_C); print(calib_D)
ggsave("calibration_ModelC.png", calib_C, width = 6, height = 5, dpi = 300)
ggsave("calibration_ModelD.png", calib_D, width = 6, height = 5, dpi = 300)

re_summary_4 <- tibble::tibble(
  Model = c("A (full)", "B (full)", "C (excl.)", "D (excl.)"),
  `Var(intercept)` = c(var_u_A, var_u_B, var_u_C, var_u_D),
  `SD(intercept)`  = c(sqrt(var_u_A), sqrt(var_u_B), sqrt(var_u_C), sqrt(var_u_D)),
  ICC              = c(icc_A, icc_B, icc_C, icc_D)
)
print(re_summary_4)
writexl::write_xlsx(re_summary_4, "random_effects_summary_4models.xlsx")

# =========================================================================
# END OF SCRIPT
# Outputs produced:
#   - plot_ame_4                        (4-model AME plot, in-memory + PNG)
#   - ame_4models_table.tex             (LaTeX table A|B|C|D, ready for \input{})
#   - AME_4models_results.xlsx          (AME estimates, all 4 models)
#   - calibration_ModelC/D.png          (calibration diagnostic plots)
#   - random_effects_summary_4models.xlsx
#   - console output: convergence, ICC, LRT, AUC, VIF, R2, overdispersion,
#     and a check that every model label present matches label_A..label_D
# =========================================================================