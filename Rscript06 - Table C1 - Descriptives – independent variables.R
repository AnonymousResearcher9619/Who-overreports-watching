library(dplyr)
library(survey)
library(haven)

# =========================================================================
# DESCRIPTIVES OF INDEPENDENT VARIABLES, WEIGHTED, BY WAVE + POOLED
#
# Produces the "Descriptives -- independent variables" table: weighted
# distributions of the study's disposition/control variables, separately
# for each of the three survey waves and for the pooled sample.
#
# data_pooled is long format (one row per respondent per debate); the
# variables described here (party_id_01, polint, duty, sex, edu,
# straightlining, speeding, age, total_tv) are respondent-level, so each
# wave's subset is first reduced to one row per respondent
# (distinct(RADIOMETER_ID2)) before computing statistics - otherwise
# respondents asked about multiple debates would be overweighted.
#
# Waves are identified by the combination of 'year' and 'wave' (the wave
# number is only unique within a given year's survey, not across years).
# =========================================================================

d2021 <- data_pooled |>
  dplyr::filter(year == 2021 & wave == 3) |>
  dplyr::distinct(RADIOMETER_ID2, .keep_all = TRUE)

d2023_w3 <- data_pooled |>
  dplyr::filter(year == 2023 & wave == 3) |>
  dplyr::distinct(RADIOMETER_ID2, .keep_all = TRUE)

d2023_w4 <- data_pooled |>
  dplyr::filter(year == 2023 & wave == 4) |>
  dplyr::distinct(RADIOMETER_ID2, .keep_all = TRUE)

dpooled <- data_pooled  # pooled column: full long-format data, no deduplication

# ---- Diagnostic: category labels and missingness per wave -----------------
# Prints, for each wave, the exact factor levels and missingness of every
# variable described below. Useful for catching label mismatches (e.g.
# differently worded categories across the source projects) or unexpected
# missingness before trusting the table.

diagnose_wave <- function(data, wave_name) {
  cat("\n=====", wave_name, "- N respondents:", nrow(data), "=====\n")
  cat("N missing weight:", sum(is.na(data$weight)), "\n")
  for (v in c("party_id_01", "polint", "duty", "sex", "edu",
              "straightlining", "speeding")) {
    cat("\n--", v, "--\n")
    cat("N missing:", sum(is.na(data[[v]])), " / ", nrow(data), "\n")
    print(table(haven::as_factor(data[[v]]), useNA = "ifany"))
  }
  cat("\nage: N missing =", sum(is.na(data$age)),
      " | total_tv: N missing =", sum(is.na(data$total_tv)), "\n")
}

diagnose_wave(d2021,    "2021 (wave 3 of CEPS 2021)")
diagnose_wave(d2023_w3, "2023 (wave 3 of CPEPS 2023)")
diagnose_wave(d2023_w4, "2023 (wave 4 of CPEPS 2023)")

# ---- Weighted statistics for one wave's subset -----------------------------

get_stats_indep <- function(data) {
  
  data <- data |> dplyr::filter(!is.na(weight))
  
  des <- survey::svydesign(
    ids = ~RADIOMETER_ID2,
    weights = ~weight,
    data = data
  )
  
  pct_var <- function(varname) {
    form <- as.formula(paste0("~", varname))
    tab <- survey::svymean(form, des, na.rm = TRUE) * 100
    setNames(as.numeric(tab), names(tab))
  }
  
  # svymean() prefixes output names with the variable name, so match a
  # category by substring rather than requiring an exact name match.
  get_pct <- function(pct_vec, label) {
    idx <- grepl(label, names(pct_vec), fixed = TRUE)
    if (!any(idx)) return(NA_real_)
    as.numeric(pct_vec[idx][1])
  }
  
  party_pct    <- pct_var("party_id_01")
  polint_pct   <- pct_var("polint")
  duty_pct     <- pct_var("duty")
  sex_pct      <- pct_var("sex")
  edu_pct      <- pct_var("edu")
  straight_pct <- pct_var("straightlining")
  speed_pct    <- pct_var("speeding")
  
  age_m  <- as.numeric(survey::svymean(~age, des, na.rm = TRUE))
  age_sd <- sqrt(as.numeric(survey::svyvar(~age, des, na.rm = TRUE)))
  
  tv_m  <- as.numeric(survey::svymean(~total_tv, des, na.rm = TRUE))
  tv_sd <- sqrt(as.numeric(survey::svyvar(~total_tv, des, na.rm = TRUE)))
  
  list(
    n = nrow(data),
    
    no_pid       = get_pct(party_pct, "no pid, DK"),
    close_pid    = get_pct(party_pct, "somewhat close or very close"),
    
    not_at_all   = get_pct(polint_pct, "not at all"),
    just_little  = get_pct(polint_pct, "just a little"),
    quite_very   = get_pct(polint_pct, "quite/very interested"),
    
    duty_right   = get_pct(duty_pct, "right/choice"),
    duty_duty    = get_pct(duty_pct, "duty"),
    
    age_m = age_m, age_sd = age_sd,
    
    male_pct = get_pct(sex_pct, "male"),
    
    edu_sec = get_pct(edu_pct, "secondary and lower"),
    edu_ter = get_pct(edu_pct, "tertiary"),
    
    straight_0 = get_pct(straight_pct, "0"),
    straight_1 = get_pct(straight_pct, "1"),
    
    speed_high = get_pct(speed_pct, "higher than 2. decile"),
    speed_low  = get_pct(speed_pct, "x < 2.decile"),
    
    tv_m = tv_m, tv_sd = tv_sd
  )
}

res_2021    <- get_stats_indep(d2021)
res_2023_w2 <- get_stats_indep(d2023_w3)
res_2023_w3 <- get_stats_indep(d2023_w4)
res_pooled  <- get_stats_indep(dpooled)

all_res <- list(res_2021, res_2023_w2, res_2023_w3, res_pooled)

# Sanity check: computed N per column should match the known sample sizes.
cat("\n=== Computed N per column (expected: 804 / 1358 / 1326 / 6944) ===\n")
cat("2021 (wave 3):", res_2021$n, "\n")
cat("2023 (wave 3):", res_2023_w2$n, "\n")
cat("2023 (wave 4):", res_2023_w3$n, "\n")
cat("Pooled:", res_pooled$n, "\n")

# ---- LaTeX table export -----------------------------------------------------

fmt <- function(x, d = 1) {
  if (is.na(x)) return("--")
  formatC(x, digits = d, format = "f")
}
row4 <- function(field) paste(sapply(all_res, function(r) fmt(r[[field]])), collapse = " & ")

n_header <- sapply(all_res, function(r) sprintf("$N = %d$", r$n))

latex_table <- c(
  "\\begin{table}[H]",
  "\\centering",
  "\\caption{Descriptives -- independent variables}",
  "\\label{tab:descriptives_independent}",
  "\\begin{adjustbox}{max width=\\textwidth, max height=0.35\\textheight}",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  "\\multicolumn{5}{l}{\\textbf{Independent variables}} \\\\",
  " & \\textbf{2021} & \\textbf{2023} & \\textbf{2023} & \\textbf{Pooled} \\\\",
  " & \\textbf{(wave 3)} & \\textbf{(wave 3)} & \\textbf{(wave 4)} & \\textbf{data} \\\\",
  paste(" &", paste(n_header, collapse = " & "), "\\\\"),
  "\\midrule",
  "\\textbf{Intensity of party identification} (\\%) & & & & \\\\",
  paste0("no pid, DK & ", row4("no_pid"), " \\\\"),
  paste0("somewhat / very close & ", row4("close_pid"), " \\\\"),
  "\\textbf{Political interest} (\\%) & & & & \\\\",
  paste0("not at all & ", row4("not_at_all"), " \\\\"),
  paste0("just a little & ", row4("just_little"), " \\\\"),
  paste0("quite/very interested & ", row4("quite_very"), " \\\\"),
  "\\textbf{Vote is duty (\\%)} & & & & \\\\",
  paste0("right/choice & ", row4("duty_right"), " \\\\"),
  paste0("duty & ", row4("duty_duty"), " \\\\"),
  "\\textbf{Age} & & & & \\\\",
  paste0("M & ", row4("age_m"), " \\\\"),
  paste0("SD & ", row4("age_sd"), " \\\\"),
  paste0("\\textbf{Gender} (\\% male) & ", row4("male_pct"), " \\\\"),
  "\\textbf{Education} (\\%) & & & & \\\\",
  paste0("secondary and lower & ", row4("edu_sec"), " \\\\"),
  paste0("tertiary & ", row4("edu_ter"), " \\\\"),
  "\\textbf{Straightlining} (\\%) & & & & \\\\",
  paste0("0 & ", row4("straight_0"), " \\\\"),
  paste0("1+ & ", row4("straight_1"), " \\\\"),
  "\\textbf{Survey speed} (\\%) & & & & \\\\",
  paste0("$x \\geq 2.\\text{ decile}$ & ", row4("speed_high"), " \\\\"),
  paste0("$x < 2.\\text{ decile}$ & ", row4("speed_low"), " \\\\"),
  "\\textbf{Total TV use (mins/day)} & & & & \\\\",
  paste0("M & ", row4("tv_m"), " \\\\"),
  paste0("SD & ", row4("tv_sd"), " \\\\"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{adjustbox}",
  "\\vspace{0.5em}",
  "\\begin{minipage}{\\textwidth}",
  "\\footnotesize",
  "\\raggedright",
  "Source: CPES 2021 (N = 804), CPES 2023 (N = 1358, 1326), pooled CPES 2021--2023 (N = 6944) -- weighted by gender, age, education, social status, media consumption and internet use.",
  "\\end{minipage}",
  "\\end{table}"
)

cat(latex_table, sep = "\n")
writeLines(latex_table, "descriptives_independent_vars_pooled_table.tex")

# Requires the LaTeX 'adjustbox' package (\usepackage{adjustbox}).