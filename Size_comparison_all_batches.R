library(tidyverse)

# ---------------------------------------------------------------------------
# Compare observed grader EQ distributions against a fitted normal, per batch
# Batch 130 excluded: it uses the old EQ cut set.
# ---------------------------------------------------------------------------

#all_batches        <- read_csv("all_batches.csv", show_col_types = FALSE)
#AppleBatchSizeStats <- read_csv("AppleBatchStats.csv", show_col_types = FALSE)

DiscreteEQCuts <- tibble(
  EQ        = c("US", "58", "63 S", "63 M", "63 L", "2PC", "67",
                "67(66.4)", "67(66.6)", "67(66.8)", "72", "OS"),
  EQCutsLow  = c(30, 52.5, 54.4, 56.1, 59.0, 61.5, 62.5, 66.2, 66.4, 66.6, 66.8, 71.9),
  EQCutsHigh = c(52.5, 54.4, 56.1, 59.0, 61.5, 62.5, 66.2, 66.4, 66.6, 66.8, 71.9, 100.0)
) |>
  mutate(EQ = fct_inorder(EQ))

# --- observed --------------------------------------------------------------

Observed <- all_batches |>
  filter(BatchName != 130) |>
  mutate(Elong = str_trim(str_remove(Grade, "^\\([A-Z]\\) ")),
         EQ    = str_trim(str_remove(Size,  "^\\(\\d+\\) "))) |>
  select(-c(Grade, Size)) |>
  summarise(Pieces = sum(Pieces), .by = c(BatchName, EQ)) |>
  # complete() guarantees a row per bin even where a batch packed none
  complete(BatchName, EQ = levels(DiscreteEQCuts$EQ), fill = list(Pieces = 0)) |>
  mutate(PropPieces = Pieces / sum(Pieces), .by = BatchName) |>
  mutate(EQ = factor(EQ, levels = levels(DiscreteEQCuts$EQ)))

# --- theoretical -----------------------------------------------------------

BatchParams <- AppleBatchSizeStats |>
  filter(BATCH_ID %in% unique(Observed$BatchName)) |>
  select(BatchName = BATCH_ID, mean = MEAN_EQ, sd = SD_EQ) |>
  mutate(BatchName = as.character(BatchName))

TheoreticalProbs <- BatchParams |>
  reframe(DiscreteEQCuts |>
            mutate(prob = pnorm(EQCutsHigh, mean, sd) - pnorm(EQCutsLow, mean, sd)),
          .by = c(BatchName, mean, sd))

Comparison <- Observed |>
  left_join(TheoreticalProbs, by = c("BatchName", "EQ")) |>
  mutate(diff = PropPieces - prob) |>
  arrange(BatchName, EQ)

# --- per-batch fit summary -------------------------------------------------
# With n in the 10^5 range a chi-square GOF test rejects on trivial deviations,
# so judge fit on effect size: total variation distance and the worst bin.

FitSummary <- Comparison |>
  summarise(n_pieces     = sum(Pieces),
            mean_EQ      = first(mean),
            sd_EQ        = first(sd),
            TVD          = 0.5 * sum(abs(diff)),
            max_abs_diff = max(abs(diff)),
            worst_bin    = EQ[which.max(abs(diff))],
            chisq_stat   = sum((Pieces - sum(Pieces) * prob)^2 / (sum(Pieces) * prob)),
            .by = BatchName) |>
  arrange(desc(TVD))

print(FitSummary, n = Inf)

# --- is the misfit in the parameters, or in the shape? ---------------------
# Refit mu and sigma by grouped-data maximum likelihood (multinomial over the
# bins). If the refitted parameters match AppleBatchStats and TVD barely moves,
# the stats are right and any residual misfit is genuine non-normality.

fit_binned_normal <- function(dat) {
  nll <- function(par) {
    p <- pnorm(dat$EQCutsHigh, par[1], exp(par[2])) -
         pnorm(dat$EQCutsLow,  par[1], exp(par[2]))
    -sum(dat$Pieces * log(pmax(p, 1e-12)))
  }
  opt <- optim(c(dat$mean[1], log(dat$sd[1])), nll)
  p   <- pnorm(dat$EQCutsHigh, opt$par[1], exp(opt$par[2])) -
         pnorm(dat$EQCutsLow,  opt$par[1], exp(opt$par[2]))
  tibble(mean_mle = opt$par[1],
         sd_mle   = exp(opt$par[2]),
         TVD_mle  = 0.5 * sum(abs(dat$PropPieces - p)))
}

MLECheck <- Comparison |>
  nest(.by = BatchName) |>
  mutate(fit = map(data, fit_binned_normal)) |>
  select(-data) |>
  unnest(fit) |>
  left_join(BatchParams, by = "BatchName") |>
  left_join(select(FitSummary, BatchName, TVD), by = "BatchName") |>
  transmute(BatchName,
            mean, mean_mle, mean_shift = mean_mle - mean,
            sd,   sd_mle,   sd_shift   = sd_mle - sd,
            TVD, TVD_mle,   TVD_gain   = TVD - TVD_mle) |>
  arrange(desc(TVD))

print(MLECheck, n = Inf)

# --- plots -----------------------------------------------------------------

batch_order <- FitSummary |> arrange(TVD) |> pull(BatchName)

p_bins <- Comparison |>
  mutate(BatchName = factor(BatchName, levels = batch_order)) |>
  pivot_longer(c(PropPieces, prob), names_to = "source", values_to = "p") |>
  mutate(source = recode(source, PropPieces = "Observed", prob = "Normal")) |>
  ggplot(aes(EQ, p, fill = source)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  facet_wrap(~BatchName, scales = "free_y") +
  scale_fill_manual(values = c(Observed = "#C1272D", Normal = "grey60")) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Observed vs normal-theoretical EQ distribution by batch",
       subtitle = "Batches ordered by increasing total variation distance",
       x = NULL, y = "Proportion of pieces", fill = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "top")

p_resid <- Comparison |>
  mutate(BatchName = factor(BatchName, levels = batch_order)) |>
  ggplot(aes(EQ, diff)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_col(fill = "#C1272D", width = 0.7) +
  facet_wrap(~BatchName) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Observed minus theoretical, by EQ bin",
       x = NULL, y = "Difference in proportion") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Normal-quantile linearity check: if EQ is normal, qnorm(cumulative observed
# proportion) plotted against the physical cut point is a straight line with
# slope 1/sd and intercept -mean/sd. Curvature = skew or heavy tails.

p_qq <- Comparison |>
  arrange(BatchName, EQ) |>
  mutate(cum = cumsum(PropPieces), .by = BatchName) |>
  filter(cum > 1e-6, cum < 1 - 1e-6) |>
  mutate(z_obs = qnorm(cum), BatchName = factor(BatchName, levels = batch_order)) |>
  ggplot(aes(EQCutsHigh, z_obs)) +
  geom_abline(aes(intercept = -mean / sd, slope = 1 / sd), colour = "grey60") +
  geom_point(colour = "#C1272D", size = 1.2) +
  facet_wrap(~BatchName) +
  labs(title = "Normal-quantile check at the EQ cut points",
       subtitle = "Grey line = normal with AppleBatchStats parameters; curvature indicates departure",
       x = "EQ cut (mm)", y = "qnorm(cumulative observed proportion)") +
  theme_minimal(base_size = 9)

print(p_bins)

write_csv(Comparison, "EQ_normality_comparison.csv")
write_csv(FitSummary, "EQ_normality_fit_summary.csv")
