library(dplyr)
library(survey)
library(wCorr)
library(ggplot2)
library(patchwork)
library(purrr)

# =========================================================================
# SELF-REPORT VS. PASSIVE EXPOSURE, BY DEBATE: STRIP/MEAN PLOTS + WEIGHTED
# SPEARMAN CORRELATION
#
# For each debate, this script (a) plots individual and weighted mean
# passive exposure (minutes) by self-report category, with the debate's
# actual duration as a reference line, and (b) computes the weighted
# Spearman correlation between self-report and passive exposure, with
# bootstrap inference. Together these show how closely the two exposure
# measures track each other, debate by debate.
#
# data_pooled is long format (one row per respondent per debate). All
# statistics below are computed on a single debate's subset at a time,
# within which a respondent contributes exactly one row, so ordinary
# (unclustered) survey weighting - svydesign(ids = ~1, ...) - and a
# simple row-level bootstrap are used throughout.
# =========================================================================

# ---- debate keys / labels / reference durations ----
debate_key <- c("2021-CT", "2021-NOVA", "2023-NOVA1",
                "2023-CT", "2023-PRIMA", "2023-NOVA2")

debate_labels_display <- setNames(
  c("2021-CT", "2021-NOVA", "2023-NOVA1", "2023-CT", "2023-PRIMA", "2023-NOVA2"),
  debate_key
)

durations <- setNames(c(109, 89, 89, 116, 84, 89), debate_key)

# ---- Passive exposure in minutes (local to this script) ----
data_pooled <- data_pooled |>
  dplyr::mutate(
    deb_passive_min = deb_passive_secs / 60
  )


# ---- Mean plot for one debate ---------------------------------------------
# Plots individual respondents (jittered) and the weighted mean passive
# exposure (mins) with CIs, grouped by self-report category and ordered
# from lowest to highest mean; the dashed line marks the debate's actual
# duration for scale.

make_debate_plot_mean <- function(
    data,
    group_var,
    outcome_var,
    weight_var,
    title,
    ref_line
) {
  group_var  <- rlang::sym(group_var)
  outcome_var <- rlang::sym(outcome_var)
  weight_var <- rlang::sym(weight_var)
  
  data <- data |>
    dplyr::filter(
      !is.na(!!group_var),
      !is.na(!!outcome_var),
      !is.na(!!weight_var)
    )
  
  theme_debates <- theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 14, face = "plain"),
      plot.title.position = "plot"
    )
  
  # ids = ~1: this data is already one row per respondent (a single
  # debate's subset), so no clustering unit is needed.
  design <- survey::svydesign(
    ids = ~1,
    data = data,
    weights = as.formula(paste0("~", rlang::as_string(weight_var)))
  )
  
  means_df <- survey::svyby(
    ~deb_passive_min,
    as.formula(paste0("~", rlang::as_string(group_var))),
    design,
    survey::svymean,
    vartype = "ci"
  ) |>
    as.data.frame() |>
    dplyr::rename(
      mean = deb_passive_min,
      ci_low = ci_l,
      ci_high = ci_u
    ) |>
    dplyr::arrange(mean)
  
  # Order self-report categories along the x-axis by their weighted mean
  lvl_order <- means_df[[rlang::as_string(group_var)]]
  
  # Wrap long category labels onto multiple lines.
  lvl_order_wrapped <- stringr::str_wrap(lvl_order, width = 12)
  
  data <- data |>
    dplyr::mutate(
      !!group_var := factor(
        stringr::str_wrap(
          .data[[rlang::as_string(group_var)]],
          width = 12
        ),
        levels = lvl_order_wrapped
      )
    )
  
  means_df <- means_df |>
    dplyr::mutate(
      !!group_var := factor(
        stringr::str_wrap(
          .data[[rlang::as_string(group_var)]],
          width = 12
        ),
        levels = lvl_order_wrapped
      )
    )
  
  ggplot2::ggplot() +
    
    # Individual observations
    ggplot2::geom_jitter(
      data = data,
      ggplot2::aes(x = !!group_var, y = !!outcome_var),
      width = 0.15, alpha = 0.25, size = 1, color = "grey40"
    ) +
    
    # Weighted mean points
    ggplot2::geom_point(
      data = means_df,
      ggplot2::aes(x = !!group_var, y = mean),
      color = "black", size = 3
    ) +
    
    # Confidence intervals
    ggplot2::geom_errorbar(
      data = means_df,
      ggplot2::aes(x = !!group_var, ymin = ci_low, ymax = ci_high),
      width = 0.15, color = "black", linewidth = 0.8
    ) +
    
    # Reference line (debate duration)
    ggplot2::geom_hline(
      yintercept = ref_line, linetype = "dashed", color = "black", linewidth = 0.8
    ) +
    
    ggplot2::labs(
      x = "Self-report",
      y = "Passive measurement (mins)",
      title = title
    ) +
    
    ggplot2::coord_cartesian(
      ylim = c(0, 125)
    ) +
    
    theme_debates
}


# ---- One panel per debate --------------------------------------------------

p1 <- make_debate_plot_mean(data_pooled |> dplyr::filter(debate == "2021-CT"),    "selfreport", "deb_passive_min", "weight", title = debate_labels_display["2021-CT"],    ref_line = durations["2021-CT"])
p2 <- make_debate_plot_mean(data_pooled |> dplyr::filter(debate == "2021-NOVA"),  "selfreport", "deb_passive_min", "weight", title = debate_labels_display["2021-NOVA"],  ref_line = durations["2021-NOVA"])
p3 <- make_debate_plot_mean(data_pooled |> dplyr::filter(debate == "2023-NOVA1"), "selfreport", "deb_passive_min", "weight", title = debate_labels_display["2023-NOVA1"], ref_line = durations["2023-NOVA1"])
p4 <- make_debate_plot_mean(data_pooled |> dplyr::filter(debate == "2023-CT"),    "selfreport", "deb_passive_min", "weight", title = debate_labels_display["2023-CT"],    ref_line = durations["2023-CT"])
p5 <- make_debate_plot_mean(data_pooled |> dplyr::filter(debate == "2023-PRIMA"), "selfreport", "deb_passive_min", "weight", title = debate_labels_display["2023-PRIMA"], ref_line = durations["2023-PRIMA"])
p6 <- make_debate_plot_mean(data_pooled |> dplyr::filter(debate == "2023-NOVA2"), "selfreport", "deb_passive_min", "weight", title = debate_labels_display["2023-NOVA2"], ref_line = durations["2023-NOVA2"])


# ---- Collage of all six debates --------------------------------------------

(p1 | p2 | p3) /
  (p4 | p5 | p6)


# ---- Weighted Spearman correlation per debate ------------------------------
# Point estimate via wCorr::weightedCorr() (a genuinely weighted Spearman
# correlation; plain cor()/cor.test() ignores weights). Inference (SE,
# p-value) via a plain weighted row-level bootstrap: since each debate's
# subset already has one row per respondent, resampling rows directly is
# equivalent to resampling respondents - no need to resample by
# RADIOMETER_ID2 and join.

get_weighted_spearman <- function(data, debate_id, n_boot = 2000, seed = 1) {
  
  d <- data |>
    dplyr::filter(
      debate == debate_id,
      !is.na(deb_passive_min),
      !is.na(selfreport),
      !is.na(weight)
    ) |>
    dplyr::mutate(
      # ordered numeric coding of selfreport: No < part < whole
      selfreport_num = as.integer(factor(
        selfreport,
        levels = c("No", "Yes, I saw a part of it", "Yes, I saw the whole debate"),
        ordered = TRUE
      ))
    )
  
  n <- nrow(d)
  
  # Point estimate on the full weighted sample
  rho_hat <- wCorr::weightedCorr(
    x = d$selfreport_num,
    y = d$deb_passive_min,
    weights = d$weight,
    method = "Spearman"
  )
  
  # Row-level weighted bootstrap
  set.seed(seed)
  boot_rho <- replicate(n_boot, {
    boot_d <- d[sample(seq_len(n), n, replace = TRUE), ]
    wCorr::weightedCorr(
      x = boot_d$selfreport_num,
      y = boot_d$deb_passive_min,
      weights = boot_d$weight,
      method = "Spearman"
    )
  })
  
  se_rho <- sd(boot_rho, na.rm = TRUE)
  z_val  <- rho_hat / se_rho
  p_val  <- 2 * stats::pnorm(-abs(z_val))
  
  data.frame(
    debate = debate_id,
    n = n,
    spearman_rho_w = rho_hat,
    se_boot = se_rho,
    p_value = p_val
  )
}

spearman_results_weighted <- purrr::map_dfr(
  debate_key,
  ~ get_weighted_spearman(data_pooled, .x)
) |>
  dplyr::mutate(
    debate_name = debate_labels_display[debate]
  ) |>
  dplyr::select(debate, debate_name, n, spearman_rho_w, se_boot, p_value)

spearman_results_weighted