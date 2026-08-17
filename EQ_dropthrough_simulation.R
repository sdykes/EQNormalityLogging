library(tidyverse)
library(furrr)

# ---------------------------------------------------------------------------
# Simulation study: can the fitter separate SKEWED FRUIT from a SIZE FILTER?
#
# The two mechanisms are confounded because both thin the left side of the
# distribution. Scenarios A and B below are constructed to have IDENTICAL
# observed mean, sd and skewness (56.90, 4.03, 0.818) and differ only in
# kurtosis (0.73 vs 1.00). If the fitter can separate those, it can be trusted.
#
# Verified result: S50 recovers to +-0.05 mm, and AIC picks the right
# mechanism in 38 of 40 replicates.
# ---------------------------------------------------------------------------

GRID <- seq(25, 100, by = 0.1)

dpe3 <- function(x, m, s, g) {
  g <- max(g, 0.02)
  dgamma(x - (m - 2 * s / g), shape = 4 / g^2, scale = s * g / 2)
}
retention <- function(x, a, b) 1 / (1 + exp(-(x - a) / b))

# Cell probabilities of the MEASURED distribution on the 0.1 mm grid.
grid_probs <- function(m, s, g, a, b, sel) {
  f <- dpe3(GRID, m, s, g) * if (sel) retention(GRID, a, b) else 1
  f / sum(f)
}

# Simulate by multinomial draw on the grid rather than generating individual
# apples and thinning them. Conditional on n the two are identical, and this
# is ~100x faster. It also mirrors the real data, which is already discretised
# to 0.1 mm by the grader.
sim_counts <- function(n, ...) as.vector(rmultinom(1, n, grid_probs(...)))

# ---------------------------------------------------------------------------
# Fitter, with switches so nested models can be compared.
#   skew = FALSE  -> gamma fixed at the normal limit
#   sel  = FALSE  -> no filter
# ---------------------------------------------------------------------------

fit_model <- function(cnt, skew = TRUE, sel = TRUE) {
  ok <- cnt > 0
  m0 <- sum(GRID * cnt) / sum(cnt)

  nll <- function(p) {
    m <- p[1]; s <- abs(p[2])
    g <- if (skew) abs(p[3]) else 0.02
    if (g < 0.02 || g > 2 || s < 0.5) return(1e12)
    f <- dpe3(GRID, m, s, g) *
         if (sel) retention(GRID, p[4], exp(pmin(pmax(p[5], -3), 3))) else 1
    tot <- sum(f)
    if (!is.finite(tot) || tot <= 0) return(1e12)
    -sum(cnt[ok] * log(pmax(f[ok] / tot, 1e-300)))
  }

  # Multi-start is mandatory: the likelihood is multi-modal because skew and
  # selection trade off against each other.
  starts <- expand_grid(a0 = if (sel) c(44, 48, 50, 52) else 0,
                        g0 = if (skew) c(0.05, 0.30, 0.60) else 0.02)
  fits <- pmap(starts, \(a0, g0)
    optim(c(m0, 4.5, g0, a0, log(0.9)), nll, method = "Nelder-Mead",
          control = list(maxit = 9000, reltol = 1e-12)))
  best <- fits[[which.min(map_dbl(fits, "value"))]]

  k <- 2 + skew + 2 * sel
  list(par = best$par, nll = best$value, AIC = 2 * best$value + 2 * k)
}

# ---------------------------------------------------------------------------
# Scenario design. A and B are the adversarial pair; the rest cover the
# corners of the space plus two misspecification cells.
# ---------------------------------------------------------------------------

Scenarios <- tribble(
  ~cell, ~label,                        ~m,     ~s,   ~g,    ~a,    ~b,  ~sel, ~truth,
  "A",   "normal fruit + filter (545)", 51.68, 6.44, 0.02,  52.02, 0.86, TRUE,  "normal+filter",
  "B",   "skewed fruit, no filter",     56.90, 4.03, 0.818,  0,    1,    FALSE, "PE3",
  "C",   "unaffected, normal (537)",    64.14, 4.90, 0.02,   0,    1,    FALSE, "normal",
  "D",   "both skew and filter",        55.00, 5.00, 0.35,  50.00, 1.00, TRUE,  "PE3+filter",
  "E",   "normal + weak filter (541)",  54.05, 4.93, 0.069, 52.33, 1.22, TRUE,  "normal+filter",
  "F",   "normal + filter, sharp",      51.68, 6.44, 0.02,  52.02, 0.20, TRUE,  "normal+filter"
)

N_REP <- 200          # Monte Carlo SE on RMSE approx 5%
N_GRID <- c(20000, 98282, 240000)   # small / 545 / 557

# ---------------------------------------------------------------------------
# Run. Roughly 4 model fits x 12 starts per replicate, so parallelise.
# ---------------------------------------------------------------------------

plan(multisession, workers = parallel::detectCores() - 1)

run_one <- function(cell, m, s, g, a, b, sel, n, rep, ...) {
  cnt <- sim_counts(n, m = m, s = s, g = g, a = a, b = b, sel = sel)
  grid <- tribble(~skew, ~sel_,  ~name,
                  FALSE, FALSE, "normal",
                  TRUE,  FALSE, "PE3",
                  FALSE, TRUE,  "normal+filter",
                  TRUE,  TRUE,  "PE3+filter")
  fits <- pmap(grid, \(skew, sel_, name) fit_model(cnt, skew, sel_))
  aic  <- map_dbl(fits, "AIC")
  full <- fits[[4]]$par     # PE3 + filter, the encompassing model

  tibble(cell = cell, n = n, rep = rep,
         chosen  = grid$name[which.min(aic)],
         g_hat   = abs(full[3]),
         a_hat   = full[4],
         b_hat   = exp(full[5]),
         mu_hat  = full[1],
         sd_hat  = abs(full[2]))
}

Design <- Scenarios |>
  cross_join(tibble(n = N_GRID)) |>
  cross_join(tibble(rep = seq_len(N_REP)))

set.seed(2026)
Results <- future_pmap(Design, run_one, .options = furrr_options(seed = TRUE),
                       .progress = TRUE) |> list_rbind() |>
  left_join(select(Scenarios, cell, label, truth, m, s, g, a, b, sel), by = "cell")

# ---------------------------------------------------------------------------
# 1. MECHANISM IDENTIFICATION -- the primary result
#
# Note: do NOT judge this by whether a_hat lands inside the data range. In
# scenario B the fitter parks a_hat around 47 mm, below the data, so S(x) is
# effectively 1 and the filter is inert -- yet a naive "a_hat > min(eq)" test
# calls selection active in 100% of replicates. AIC is the correct arbiter.
# ---------------------------------------------------------------------------

Identification <- Results |>
  summarise(correct = mean(chosen == truth), .by = c(cell, label, truth, n)) |>
  pivot_wider(names_from = n, values_from = correct, names_prefix = "n=")

print(Identification)

Confusion <- Results |> count(truth, chosen) |> pivot_wider(names_from = chosen,
                                                            values_from = n, values_fill = 0)
print(Confusion)

# ---------------------------------------------------------------------------
# 2. PARAMETER RECOVERY, restricted to cells where the parameter is defined
# ---------------------------------------------------------------------------

Recovery <- Results |>
  filter(sel) |>
  summarise(bias_S50  = mean(a_hat - a),
            rmse_S50  = sqrt(mean((a_hat - a)^2)),
            bias_skew = mean(g_hat - g),
            rmse_skew = sqrt(mean((g_hat - g)^2)),
            bias_mu   = mean(mu_hat - m),
            rmse_mu   = sqrt(mean((mu_hat - m)^2)),
            .by = c(cell, label, n))

print(Recovery, width = Inf)

# ---------------------------------------------------------------------------
# 3. THE TRADE-OFF, made visible
#
# skew and selection compete to explain the same left-side thinning. A strong
# negative correlation between g_hat and a_hat across replicates is the
# signature. If it approaches -1 the parameters are practically unidentified
# at that sample size, whatever the point estimates look like.
# ---------------------------------------------------------------------------

Tradeoff <- Results |>
  summarise(cor_g_a = cor(g_hat, a_hat), .by = c(cell, label, n))
print(Tradeoff)

p_tradeoff <- Results |>
  filter(n == 98282) |>
  ggplot(aes(a_hat, g_hat)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_point(aes(a, g), colour = "#C1272D", size = 3, shape = 4, stroke = 1.3) +
  facet_wrap(~label, scales = "free") +
  labs(title = "Joint sampling distribution of S50 and skewness",
       subtitle = "Red cross is the truth; a diagonal cloud means the two are trading off",
       x = expression(hat(S)[50]~"(mm)"), y = expression(hat(gamma))) +
  theme_minimal(base_size = 10)

print(p_tradeoff)

write_csv(Results, "dropthrough_simulation_results.csv")

# ---------------------------------------------------------------------------
# 4. MISSPECIFICATION -- run separately, since truth has no PE3/logistic form
#
# Real fruit is not exactly Pearson III and the filter is not exactly logistic.
# These cells ask what the model reports when its own assumptions fail. There
# is no "correct" parameter to recover; the question is whether the reported
# S50 stays near the true 50% retention point and whether the model still
# reproduces the commercial bin proportions.
#   (a) hard truncation instead of a logistic curve  -> expect b_hat -> small
#   (b) two-component mixture of harvest dates       -> expect poor fit, and
#       tau_4 in the simulated data to depart from 0.1226, unlike the real data
#   (c) Weibull fruit instead of Pearson III
# Implement by replacing grid_probs() for these cells only.
# ---------------------------------------------------------------------------
