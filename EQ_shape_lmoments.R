library(tidyverse)
library(lmom)

# ---------------------------------------------------------------------------
# Per-batch EQ shape analysis via L-moments
#
# EQ is recorded at 0.1 mm resolution, so a GROUP BY on rounded EQ is a
# LOSSLESS representation of the batch distribution in ~650 rows. Pull that
# once; everything below is computed from it without further querying.
# ---------------------------------------------------------------------------

sql <- "
SELECT  
    BATCH_ID
    ,ROUND(MINOR, 1) AS EQ
    ,COUNT(*) AS N 
FROM   ROCKIT_DATA_PROD.COMPAC.STG_COMPAC_BATCH
WHERE  MINOR BETWEEN 20 AND 110
AND START_TIME > '2026-01-01 00:00:00.000'
AND NOT (SIZER_GRADE_NAME IN ('Capture','Captures','Rcy','Capture ','Recycle','Doub','Doubles ','Capt','Leaf','Cap') AND NOT SIZE_NAME IN ('US','OS','OOS','OSS'))
AND BATCH_ID <> 0
GROUP  BY 1, 2;
"
con <- DBI::dbConnect(
  odbc::odbc(),
  Driver            = "SnowflakeDSIIDriver",   
  Server            = Sys.getenv("SNOWFLAKE_SERVER"),
  warehouse         = Sys.getenv("SNOWFLAKE_WHOUSE"),
  UID               = Sys.getenv("SNOWFLAKE_UID"),
  Authenticator     = "snowflake_jwt",
  PRIV_KEY_FILE     = Sys.getenv("SNOWFLAKE_KEY_PATH"),
  PRIV_KEY_FILE_PWD = Sys.getenv("SNOWFLAKE_KEY_PWD")
)

EQHist <- DBI::dbGetQuery(con, sql) |> as_tibble() |> rename_with(tolower)
#EQHist <- read_csv("eq_histogram.csv", show_col_types = FALSE)

DBI::dbDisconnect(con)

stopifnot(all(c("batch_id", "eq", "n") %in% names(EQHist)))

# ---------------------------------------------------------------------------
# 1. Shape statistics per batch
#
# L-moments are computed on the expanded vector. This is exact — samlmu()
# takes no weights, and expanding is cheap per batch (200k doubles ~ 1.6 MB).
# Conventional g1/g2 are carried alongside purely for comparison.
# ---------------------------------------------------------------------------

shape_stats <- function(eq, n) {
  x <- rep(eq, n)
  l <- lmom::samlmu(x, nmom = 4)
  m <- mean(x); s <- sd(x); d <- x - m

  # trimmed g2: kurtosis is dominated by extremes, so flag if trimming moves it
  xt <- x[x >= quantile(x, 0.0005) & x <= quantile(x, 0.9995)]
  dt <- xt - mean(xt)

  tibble(n_fruit = length(x),
         mean_eq = m, sd_eq = s, cv = s / m,
         L2      = unname(l["l_2"]),
         t3      = unname(l["t_3"]),        # L-skewness   (normal: 0)
         t4      = unname(l["t_4"]),        # L-kurtosis   (normal: 0.1226)
         g1      = mean(d^3) / s^3,
         g2      = mean(d^4) / s^4 - 3,
         g2_trim = mean(dt^4) / sd(xt)^4 - 3)
}

Shape <- EQHist |>
  summarise(shape_stats(eq, n), .by = batch_id) |>
  mutate(t4_excess = t4 - 0.1226,
         g2_fragile = abs(g2 - g2_trim) > 0.5 * abs(g2))   # outlier-driven g2

# ---------------------------------------------------------------------------
# 2. Null band: how big is t3 under exact normality, at each batch's n?
#
# Rather than rely on asymptotic constants, simulate on a grid of n and fit
# the sd ~ c / sqrt(n) law, then interpolate to each batch. Cheap and exact
# enough; c is estimated to ~2% with 300 reps.
# ---------------------------------------------------------------------------

set.seed(42)
n_grid <- c(2e3, 5e3, 1e4, 2e4, 5e4)
NullSim <- tibble(n = n_grid) |>
  mutate(sd_t3 = map_dbl(n, \(nn) sd(replicate(300, lmom::samlmu(rnorm(nn), 3)[["t_3"]]))),
         sd_t4 = map_dbl(n, \(nn) sd(replicate(300, lmom::samlmu(rnorm(nn), 4)[["t_4"]]))))

c_t3 <- mean(NullSim$sd_t3 * sqrt(NullSim$n))
c_t4 <- mean(NullSim$sd_t4 * sqrt(NullSim$n))

Shape <- Shape |>
  mutate(se_t3 = c_t3 / sqrt(n_fruit),
         se_t4 = c_t4 / sqrt(n_fruit),
         z_t3  = t3 / se_t3)

# NOTE: at n = 100k, se_t3 is around 0.0006. Every batch will be "significant".
# Judge on MAGNITUDE against the reference scale built next, not on z.

# ---------------------------------------------------------------------------
# 3. Reference scale: what t3 do candidate skewed families give?
#
# Computed by simulation rather than from parameterisation formulae, so it is
# correct regardless of how each family is parameterised. t3 for these families
# depends only on the CV, so build a lookup over the observed CV range.
# ---------------------------------------------------------------------------

t3_of <- function(cv, family, nsim = 3e5) {
  x <- switch(family,
    lognormal = { s <- sqrt(log(1 + cv^2)); rlnorm(nsim, -s^2 / 2, s) },
    gamma     = rgamma(nsim, shape = 1 / cv^2, scale = cv^2),
    weibull   = { k <- uniroot(\(k) sqrt(gamma(1 + 2/k) / gamma(1 + 1/k)^2 - 1) - cv,
                               c(0.3, 200))$root
                  rweibull(nsim, shape = k) })
  lmom::samlmu(x, 3)[["t_3"]]
}

RefScale <- expand_grid(cv     = seq(min(Shape$cv), max(Shape$cv), length.out = 12),
                        family = c("lognormal", "gamma", "weibull")) |>
  mutate(t3_ref = map2_dbl(cv, family, t3_of))

# Practical threshold: a lognormal matched to a typical batch CV. Departures
# smaller than this are not worth acting on.
t3_material <- t3_of(median(Shape$cv), "lognormal")

Shape <- Shape |> mutate(material = abs(t3) > 0.5 * t3_material)

message(sprintf("Reference: lognormal at CV %.3f gives t3 = %.4f. Flagging |t3| > %.4f",
                median(Shape$cv), t3_material, 0.5 * t3_material))

# ---------------------------------------------------------------------------
# 4. Does skew scale with fruit size? (the small-fruit finding, on all batches)
# ---------------------------------------------------------------------------

skew_fit <- lm(t3 ~ mean_eq, data = Shape, weights = 1 / se_t3^2)
print(summary(skew_fit))

p_trend <- Shape |>
  ggplot(aes(mean_eq, t3)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_hline(yintercept = c(-1, 1) * t3_material, linetype = 2, colour = "grey65") +
  geom_point(aes(size = n_fruit, colour = material), alpha = 0.6) +
  geom_smooth(method = "lm", colour = "#C1272D", se = TRUE) +
  scale_colour_manual(values = c(`TRUE` = "#C1272D", `FALSE` = "grey55"), guide = "none") +
  scale_size_area(max_size = 4, guide = "none") +
  labs(title = "L-skewness of EQ against batch mean size",
       subtitle = "Dashed lines: half the L-skewness of a matched lognormal",
       x = "Batch mean EQ (mm)", y = expression(tau[3]~"(L-skewness)")) +
  theme_minimal(base_size = 11)

p_trend

# ---------------------------------------------------------------------------
# 5. L-moment ratio diagram — which family does each batch resemble?
#
# Each batch is one point. Family curves are the locus each 3-parameter family
# traces as its shape parameter varies. Whichever curve the cloud sits on is
# the family to fit; the normal is a single point at (0, 0.1226).
# ---------------------------------------------------------------------------

family_curve <- function(family, nsim = 3e5) {
  tibble(cv = seq(0.01, 0.45, length.out = 25)) |>
    mutate(lm = map(cv, \(cv) {
      x <- switch(family,
        lognormal = { s <- sqrt(log(1 + cv^2)); rlnorm(nsim, -s^2 / 2, s) },
        gamma     = rgamma(nsim, shape = 1 / cv^2, scale = cv^2),
        weibull   = { k <- uniroot(\(k) sqrt(gamma(1 + 2/k) / gamma(1 + 1/k)^2 - 1) - cv,
                                   c(0.3, 200))$root
                      rweibull(nsim, shape = k) })
      l <- lmom::samlmu(x, 4); tibble(t3 = l[["t_3"]], t4 = l[["t_4"]])
    })) |>
    unnest(lm) |> mutate(family = family)
}

Curves <- map(c("lognormal", "gamma", "weibull"), family_curve) |> list_rbind()

p_lmrd <- ggplot() +
  geom_line(data = Curves, aes(t3, t4, colour = family), linewidth = 0.7) +
  geom_point(data = Shape, aes(t3, t4), alpha = 0.45, size = 1.4) +
  annotate("point", x = 0, y = 0.1226, shape = 4, size = 4, stroke = 1.2, colour = "#C1272D") +
  annotate("text", x = 0, y = 0.1226, label = "normal", hjust = -0.25, vjust = 1.8,
           colour = "#C1272D", size = 3.5) +
  labs(title = "L-moment ratio diagram",
       subtitle = "One point per batch; curves are the loci of candidate families",
       x = expression(tau[3]~"(L-skewness)"), y = expression(tau[4]~"(L-kurtosis)"),
       colour = NULL) +
  theme_minimal(base_size = 11)

p_lmrd
# lmom's own version, with the standard hydrological family set, for comparison:
# lmrd(Shape$t3, Shape$t4, pch = 16, col = "grey40")

# ---------------------------------------------------------------------------
# 6. QQ panels — normal vs a fitted skewed alternative (Pearson III)
#
# Plotted on the histogram's own grid, so no need to expand to full vectors.
# ---------------------------------------------------------------------------

qq_data <- function(eq, n, l1, l2, t3) {
  p <- (cumsum(n) - 0.5) / sum(n)
  keep <- p > 1e-5 & p < 1 - 1e-5
  tibble(obs    = eq[keep],
         normal = qnorm(p[keep], l1, l2 * sqrt(pi)),
         pe3    = lmom::quape3(p[keep], lmom::pelpe3(c(l1, l2, t3))))
}

QQ <- EQHist |>
  arrange(batch_id, eq) |>
  left_join(select(Shape, batch_id, L2, t3, mean_eq), by = "batch_id") |>
  summarise(qq_data(eq, n, first(mean_eq), first(L2), first(t3)), .by = batch_id)

p_qq <- QQ |>
  filter(batch_id %in% head(arrange(Shape, desc(abs(t3)))$batch_id, 12)) |>
  pivot_longer(c(normal, pe3), names_to = "model", values_to = "theoretical") |>
  ggplot(aes(theoretical, obs, colour = model)) +
  geom_abline(colour = "grey55") +
  geom_line(linewidth = 0.6) +
  facet_wrap(~batch_id, scales = "free") +
  scale_colour_manual(values = c(normal = "grey45", pe3 = "#C1272D")) +
  labs(title = "QQ plots, 12 most skewed batches",
       subtitle = "Pearson III (red) vs normal (grey); departure from the 1:1 line is misfit",
       x = "Theoretical EQ (mm)", y = "Observed EQ (mm)", colour = NULL) +
  theme_minimal(base_size = 9)

print(p_trend); print(p_lmrd)

write_csv(Shape, "EQ_shape_summary.csv")

# ---------------------------------------------------------------------------
# 7. Reconciliation check — RUN THIS FIRST
#
# Confirms the per-fruit query returns the same population that produced
# AppleBatchStats. APPLES currently runs 0.1-6% above sum(Pieces), so there is
# a pre/post-reject distinction somewhere. If mean_eq and sd_eq below do not
# reproduce MEAN_EQ and SD_EQ, the two sources are different populations and
# nothing above is comparable to the earlier binned work.
# ---------------------------------------------------------------------------

Reconcile <- Shape |>
  select(batch_id, n_fruit, mean_eq, sd_eq) |>
  inner_join(AppleBatchSizeStats |> select(batch_id = BATCH_ID, APPLES, MEAN_EQ, SD_EQ),
             by = "batch_id") |>
  mutate(d_n = n_fruit - APPLES, d_mean = mean_eq - MEAN_EQ, d_sd = sd_eq - SD_EQ)

Reconcile |> summarise(across(c(d_n, d_mean, d_sd), \(x) max(abs(x)))) |> print()
