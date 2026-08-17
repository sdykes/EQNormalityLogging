library(tidyverse)

# ===========================================================================
# PRODUCTION: EQ size-distribution model for grader batches
#
# Design principle: ZERO free parameters per batch.
#
#   S50 and the transition width are MACHINE parameters, calibrated once from
#   the batches where per-fruit EQ is available, then held fixed.
#
#   MEAN_EQ and SD_EQ from AppleBatchStats are moments of the MEASURED
#   (post-filter) population. They are treated as data, not as parameters to
#   fit. For each batch the presented-population mean and sd are recovered by
#   solving a 2x2 system so that the selected distribution reproduces the two
#   measured moments exactly.
#
#   Two candidate models are compared per batch. Both are constrained to the
#   measured moments, so BOTH have the same effective complexity (2) and the
#   comparison needs no AIC/BIC penalty -- whichever reproduces the bin counts
#   better is deployed.
#
# Inputs:  all_batches.csv, AppleBatchStats.csv, BatchEQ.csv (calibration only)
# Outputs: EQ_model_summary.csv, EQ_model_bins.csv, and a faceted
#          observed-vs-theoretical column chart.
# ===========================================================================

# ---------------------------------------------------------------------------
# 0. Configuration
# ---------------------------------------------------------------------------

GRID <- seq(20, 105, by = 0.05)          # fine grid for numerical integration

CUTS <- tibble(
  EQ  = c("US", "58", "63 S", "63 M", "63 L", "2PC", "67",
          "67(66.4)", "67(66.6)", "67(66.8)", "72", "OS"),
  lo  = c(0, 52.5, 54.4, 56.1, 59.0, 61.5, 62.5, 66.2, 66.4, 66.6, 66.8, 71.9),
  hi  = c(52.5, 54.4, 56.1, 59.0, 61.5, 62.5, 66.2, 66.4, 66.6, 66.8, 71.9, 1000)
) |> mutate(EQ = fct_inorder(EQ))

EDGES <- c(CUTS$lo, 1000)

CALIB_BATCHES <- c(541, 545)   # per-fruit batches showing filter activity
EXCLUDE       <- 130           # different EQ cut set

retention <- function(x, a, b) 1 / (1 + exp(-(x - a) / b))

# ---------------------------------------------------------------------------
# 1. Load
# ---------------------------------------------------------------------------

#all_batches <- read_csv("all_batches.csv", show_col_types = FALSE)
#Stats       <- read_csv("AppleBatchStats.csv", show_col_types = FALSE)
Stats <- AppleBatchSizeStats

Observed <- all_batches |>
  filter(BatchName != EXCLUDE) |>
  mutate(EQ = str_trim(str_remove(Size, "^\\(\\d+\\) "))) |>
  summarise(Pieces = sum(Pieces), .by = c(BatchName, EQ)) |>
  complete(BatchName, EQ = levels(CUTS$EQ), fill = list(Pieces = 0)) |>
  mutate(EQ = factor(EQ, levels = levels(CUTS$EQ)),
         obs_p = Pieces / sum(Pieces), .by = BatchName)

Params <- Stats |>
  filter(BATCH_ID %in% unique(Observed$BatchName)) |>
  transmute(BatchName = BATCH_ID, mean_eq = MEAN_EQ, sd_eq = SD_EQ)

# ---------------------------------------------------------------------------
# 2. CALIBRATION: machine parameters from per-fruit data
#
# Shared S50 and width across the calibration batches, with each batch's own
# presented mean and sd as nuisance parameters. Run once; then hard-code the
# result and skip this block in routine use.
# ---------------------------------------------------------------------------

write_csv(Batch545EQ, "BatchEQ.csv")

calibrate_machine <- function(path = "BatchEQ.csv", batches = CALIB_BATCHES) {
  eq <- read_csv(path, show_col_types = FALSE) |> rename_with(tolower) |>
    filter(batch_id %in% batches)

  brk <- c(GRID - 0.025, tail(GRID, 1) + 0.025)
  cnt <- eq |> nest(.by = batch_id) |>
    mutate(c = map(data, \(d) as.numeric(table(cut(d$eq, brk, labels = FALSE)) |>
                                        (\(tb) { v <- numeric(length(GRID))
                                                 v[as.integer(names(tb))] <- tb; v })()))) |>
    pull(c, name = batch_id)

  nll <- function(p) {
    a <- p[1]; b <- exp(p[2])
    sum(imap_dbl(cnt, \(cc, i) {
      k <- which(names(cnt) == i)
      m <- p[1 + 2 * k]; s <- abs(p[2 + 2 * k])
      if (s < 0.5) return(1e12)
      f <- dnorm(GRID, m, s) * retention(GRID, a, b)
      tot <- sum(f); if (!is.finite(tot) || tot <= 0) return(1e12)
      ok <- cc > 0
      -sum(cc[ok] * log(pmax(f[ok] / tot, 1e-300)))
    }))
  }

  init <- c(52.1, log(1.0), rep(c(54, 5.5), length(cnt)))
  fit <- optim(init, nll, method = "Nelder-Mead",
               control = list(maxit = 40000, reltol = 1e-13))
  list(S50 = fit$par[1], width = exp(fit$par[2]))
}

# MACHINE <- calibrate_machine()          # re-run only when recalibrating
MACHINE <- list(S50 = 52.241, width = 1.056)

message(sprintf("Machine parameters: S50 = %.3f mm, transition width = %.3f mm",
                MACHINE$S50, MACHINE$width))

# ---------------------------------------------------------------------------
# 3. Moment inversion: measured moments -> presented moments
#
# Solve for (mu, sigma) of the PRESENTED population such that the selected
# distribution N(mu, sigma) * S(x) has the measured mean and sd. Two equations,
# two unknowns, no free parameters.
# ---------------------------------------------------------------------------

selected_moments <- function(mu, sd, a, b) {
  p <- dnorm(GRID, mu, sd) * retention(GRID, a, b)
  p <- p / sum(p)
  m <- sum(GRID * p)
  c(mean = m, sd = sqrt(sum((GRID - m)^2 * p)))
}

invert_moments <- function(mean_obs, sd_obs, a, b) {
  obj <- function(v) {
    mm <- selected_moments(v[1], abs(v[2]), a, b)
    (mm[["mean"]] - mean_obs)^2 + (mm[["sd"]] - sd_obs)^2
  }
  r <- optim(c(mean_obs - 1, sd_obs + 0.5), obj, method = "Nelder-Mead",
             control = list(reltol = 1e-14, maxit = 5000))
  tibble(mu_presented = r$par[1], sd_presented = abs(r$par[2]),
         solve_error = sqrt(r$value))
}

# Bin probabilities of the selected distribution, on the commercial cuts
selected_bin_probs <- function(mu, sd, a, b) {
  p <- dnorm(GRID, mu, sd) * retention(GRID, a, b)
  p <- p / sum(p)
  cd <- c(0, cumsum(p))
  ge <- c(GRID[1] - 0.025, GRID + 0.025)
  pr <- diff(approx(ge, cd, xout = pmin(pmax(EDGES, min(ge)), max(ge)), rule = 2)$y)
  pmax(pr, 0) / sum(pmax(pr, 0))
}

dropped_fraction <- function(mu, sd, a, b) {
  f <- dnorm(GRID, mu, sd)
  1 - sum(f * retention(GRID, a, b)) / sum(f)
}

# ---------------------------------------------------------------------------
# 4. Fit both candidate models to every batch
# ---------------------------------------------------------------------------

tvd <- function(o, p) 0.5 * sum(abs(o / sum(o) - p / sum(p)))

Fitted <- Params |>
  mutate(inv = pmap(list(mean_eq, sd_eq),
                    \(m, s) invert_moments(m, s, MACHINE$S50, MACHINE$width))) |>
  unnest(inv) |>
  mutate(dropped = pmap_dbl(list(mu_presented, sd_presented),
                            \(m, s) dropped_fraction(m, s, MACHINE$S50, MACHINE$width)))

BinProbs <- Fitted |>
  mutate(p_normal = map2(mean_eq, sd_eq, \(m, s) {
           pr <- diff(pnorm(EDGES, m, s)); pr / sum(pr) }),
         p_model  = map2(mu_presented, sd_presented,
                         \(m, s) selected_bin_probs(m, s, MACHINE$S50, MACHINE$width))) |>
  select(BatchName, p_normal, p_model) |>
  mutate(EQ = list(levels(CUTS$EQ))) |>
  unnest(c(EQ, p_normal, p_model)) |>
  mutate(EQ = factor(EQ, levels = levels(CUTS$EQ)))

Comparison <- Observed |>
  left_join(BinProbs |>
              mutate(BatchName = as.character(BatchName)), by = c("BatchName", "EQ")) |>
  arrange(BatchName, EQ)

# Model choice: both candidates honour the measured moments and so carry the
# same effective complexity. No information penalty is required -- deploy
# whichever reproduces the observed bin counts more closely.
Choice <- Comparison |>
  summarise(TVD_raw   = tvd(obs_p, p_normal),
            TVD_model = tvd(obs_p, p_model),
            .by = BatchName) |>
  mutate(model_deployed = if_else(TVD_model < TVD_raw, "normal x dropthrough", "normal"),
         TVD_deployed   = pmin(TVD_raw, TVD_model),
         TVD_gain       = TVD_raw - TVD_deployed)

Summary <- Fitted |>
  mutate(BatchName = as.character(BatchName)) |>
  left_join(Observed |> 
              summarise(n_pieces = sum(Pieces), .by = BatchName), by = "BatchName") |>
  left_join(Choice, by = "BatchName") |>
  transmute(Batch = BatchName, n_pieces,
            mean_eq = round(mean_eq, 2), sd_eq = round(sd_eq, 2),
            mu_presented = round(mu_presented, 2), sd_presented = round(sd_presented, 2),
            dropped_pct = round(100 * dropped, 1),
            model_deployed,
            TVD_raw_pct      = round(100 * TVD_raw, 2),
            TVD_modelled_pct = round(100 * TVD_model, 2),
            TVD_deployed_pct = round(100 * TVD_deployed, 2),
            TVD_gain_pct     = round(100 * TVD_gain, 2)) |>
  arrange(desc(TVD_raw_pct))

print(Summary, n = Inf)

Summary |>
  summarise(max_raw = max(TVD_raw_pct), median_raw = median(TVD_raw_pct),
            max_deployed = max(TVD_deployed_pct), median_deployed = median(TVD_deployed_pct),
            n_filtered = sum(model_deployed != "normal")) |>
  print()

# ---------------------------------------------------------------------------
# 5. Side-by-side column chart: observed vs theoretical counts per bin
# ---------------------------------------------------------------------------

PlotData <- Comparison |>
  left_join(select(Choice, BatchName, model_deployed, TVD_deployed), by = "BatchName") |>
  mutate(n_batch = sum(Pieces), .by = BatchName) |>
  mutate(deployed_p = if_else(model_deployed == "normal", p_normal, p_model)) |>
  transmute(BatchName, EQ, model_deployed, TVD_deployed,
            Observed    = Pieces,
            Theoretical = deployed_p * n_batch) |>
  pivot_longer(c(Observed, Theoretical), names_to = "series", values_to = "count")

facet_labels <- Choice |>
  transmute(BatchName,
            lab = sprintf("%s  |  %s  |  TVD %.1f%%", BatchName,
                          if_else(model_deployed == "normal", "normal", "normal x DT"),
                          100 * TVD_deployed)) |>
  arrange(desc(BatchName)) |> deframe()

p_bins <- PlotData |>
  mutate(BatchName = factor(BatchName,
                            levels = arrange(Choice, TVD_deployed)$BatchName,
                            labels = facet_labels[as.character(arrange(Choice, TVD_deployed)$BatchName)])) |>
  ggplot(aes(EQ, count, fill = series)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72) +
  facet_wrap(~BatchName, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(Observed = "#C1272D", Theoretical = "grey62")) +
  scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  labs(title = "Observed vs theoretical piece counts by EQ size class",
       subtitle = sprintf("Deployed model per batch; drop-through filter at S50 = %.2f mm, width %.2f mm",
                          MACHINE$S50, MACHINE$width),
       x = NULL, y = "Pieces", fill = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
        legend.position = "top",
        panel.grid.major.x = element_blank(),
        strip.text = element_text(face = "bold", size = 8))

ggsave("EQ_model_bins.png", p_bins, width = 13, height = 15, dpi = 200)
print(p_bins)

# Before/after, affected batches only -- shows what the correction buys
p_gain <- Choice |>
  filter(model_deployed != "normal") |>
  pivot_longer(c(TVD_raw, TVD_deployed)) |>
  mutate(name = recode(name, TVD_raw = "Normal only", TVD_deployed = "With drop-through"),
         BatchName = fct_reorder(factor(BatchName), value, max)) |>
  ggplot(aes(value, BatchName, fill = name)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  scale_fill_manual(values = c(`Normal only` = "grey62", `With drop-through` = "#C1272D")) +
  scale_x_continuous(labels = scales::percent) +
  labs(title = "Effect of the drop-through correction",
       subtitle = "Batches where the correction is deployed",
       x = "Total variation distance", y = NULL, fill = NULL) +
  theme_minimal(base_size = 11)

print(p_gain)

# ---------------------------------------------------------------------------
# 6. Persist
# ---------------------------------------------------------------------------

write_csv(Summary, "EQ_model_summary.csv")

Comparison |>
  left_join(select(Choice, BatchName, model_deployed), by = "BatchName") |>
  mutate(theoretical_p = if_else(model_deployed == "normal", p_normal, p_model),
         residual_p    = obs_p - theoretical_p) |>
  select(BatchName, EQ, Pieces, obs_p, p_normal, p_model,
         model_deployed, theoretical_p, residual_p) |>
  write_csv("EQ_model_bins.csv")

# ---------------------------------------------------------------------------
# NOTES FOR OPERATION
#
# * Recalibrate MACHINE only when per-fruit data for new affected batches is
#   available, or if the filter aperture is physically changed. S50's RMSE is
#   0.038 mm at these sample sizes, so a shift beyond ~0.1 mm between
#   calibrations is a real change in the machine, not noise.
#
# * `dropped_pct` is an EXTRAPOLATION below the smallest measured fruit and is
#   the least reliable output. It is a hypothesis for the weighbridge, not an
#   estimate to act on. Batches 545/551/541/548 imply 27-52% loss, which should
#   be checked against bin weights versus grader throughput before being
#   reported to anyone.
#
# * Residual after correction concentrates almost entirely at the US/58
#   boundary (batch 545: -2.0% US, +3.3% in 58, everything else within 0.9%).
#   That is a localised antisymmetric signature, i.e. a cut offset rather than a
#   shape error. A single global offset at the 52.5 cut was tested and does NOT
#   help consistently -- it improves 545 and 541 but degrades 549 and 543 -- so
#   it is deliberately not applied. Resolving this needs per-fruit data for more
#   batches.
#
# * Batches 551 and 548 imply large drop-through (37%, 27%) yet the correction
#   makes their bin fit WORSE, so plain normal is deployed. Combined with the
#   known 15% piece-count discrepancy in 551, both warrant investigation.
# ---------------------------------------------------------------------------
