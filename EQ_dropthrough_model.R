library(tidyverse)

# ---------------------------------------------------------------------------
# Grader drop-through selection model
#
#   presented density   f(x)  = Pearson III(mu, sigma, gamma)   [normal at gamma = 0]
#   retention curve     S(x)  = 1 / (1 + exp(-(x - a) / b))
#   measured density    f*(x) = f(x) S(x) / integral( f S )
#
# gamma is left FREE so the model can explain the data with skew if skew is
# what is there. On batches 537/541/545/557 it does not: fitted gamma comes
# back at 0.02-0.07 once selection is included.
#
# Input: BatchEQ.csv, per-fruit EQ at 0.1 mm resolution (BATCH_ID, EQ).
# ---------------------------------------------------------------------------

BatchEQ <- read_csv("BatchEQ.csv", show_col_types = FALSE) |> rename_with(tolower)

GRID <- seq(25, 100, by = 0.1)

# Pearson III density by (mean, sd, skew). Positive skew only; gamma is
# floored at 0.02 because the gamma parameterisation is singular at zero and
# the normal limit is reached well before then.
dpe3 <- function(x, m, s, g) {
  g <- max(g, 0.02)
  dgamma(x - (m - 2 * s / g), shape = 4 / g^2, scale = s * g / 2)
}

retention <- function(x, a, b) 1 / (1 + exp(-(x - a) / b))

# ---------------------------------------------------------------------------
# Fit one batch. Multi-start is NOT optional: the likelihood is multi-modal
# and trades skew against selection, since both thin the left side. A single
# start from naive values found the wrong optimum for batch 541.
# ---------------------------------------------------------------------------

fit_dropthrough <- function(eq,
                            a_starts = c(35, 45, 48, 50, 52),
                            g_starts = c(0.05, 0.30, 0.60)) {

  idx <- round((eq - GRID[1]) / 0.1) + 1
  cnt <- tabulate(idx, nbins = length(GRID))
  obs <- cnt > 0

  nll <- function(p) {
    m <- p[1]; s <- abs(p[2]); g <- abs(p[3]); a <- p[4]; b <- exp(p[5])
    if (g < 0.02 || g > 2 || s < 0.5 || abs(p[5]) > 3) return(1e12)
    f <- dpe3(GRID, m, s, g) * retention(GRID, a, b)
    tot <- sum(f)
    if (!is.finite(tot) || tot <= 0) return(1e12)
    -sum(cnt[obs] * log(pmax(f[obs] / tot, 1e-300)))
  }

  starts <- expand_grid(a0 = a_starts, g0 = g_starts)
  fits <- pmap(starts, \(a0, g0)
    optim(c(mean(eq), sd(eq), g0, a0, log(0.8)), nll,
          method = "Nelder-Mead",
          control = list(maxit = 20000, reltol = 1e-12)))

  best <- fits[[which.min(map_dbl(fits, "value"))]]
  p <- best$par
  m <- p[1]; s <- abs(p[2]); g <- abs(p[3]); a <- p[4]; b <- exp(p[5])

  f0 <- dpe3(GRID, m, s, g)
  tibble(mu_presented = m, sd_presented = s, skew = g,
         S50 = a, width = b,
         retained = sum(f0 * retention(GRID, a, b)) / sum(f0),
         nll = best$value, n = length(eq))
}

# ---------------------------------------------------------------------------
# Is the filter real? Decide by NESTED MODEL COMPARISON, not by where S50 lands.
#
# A fitted S50 below the data range leaves S(x) ~ 1 and the filter inert, but
# the parameter still has a value -- so testing "S50 > min(eq)" wrongly reports
# selection in 100% of no-filter replicates (see EQ_dropthrough_simulation.R).
# AIC over the four nested models picks the right mechanism in 38 of 40
# simulated replicates.
# ---------------------------------------------------------------------------

fit_nested <- function(eq) {
  idx <- round((eq - GRID[1]) / 0.1) + 1
  cnt <- tabulate(idx, nbins = length(GRID))
  ok  <- cnt > 0
  m0  <- mean(eq)

  one <- function(skew, sel) {
    nll <- function(p) {
      m <- p[1]; s <- abs(p[2]); g <- if (skew) abs(p[3]) else 0.02
      if (g < 0.02 || g > 2 || s < 0.5) return(1e12)
      f <- dpe3(GRID, m, s, g) *
           if (sel) retention(GRID, p[4], exp(pmin(pmax(p[5], -3), 3))) else 1
      tot <- sum(f)
      if (!is.finite(tot) || tot <= 0) return(1e12)
      -sum(cnt[ok] * log(pmax(f[ok] / tot, 1e-300)))
    }
    st <- expand_grid(a0 = if (sel) c(44, 48, 50, 52) else 0,
                      g0 = if (skew) c(0.05, 0.30, 0.60) else 0.02)
    fs <- pmap(st, \(a0, g0) optim(c(m0, sd(eq), g0, a0, log(0.9)), nll,
                                   method = "Nelder-Mead",
                                   control = list(maxit = 9000, reltol = 1e-12)))
    b <- fs[[which.min(map_dbl(fs, "value"))]]
    tibble(nll = b$value, k = 2 + skew + 2 * sel, AIC = 2 * b$value + 2 * (2 + skew + 2 * sel))
  }

  tribble(~model, ~skew, ~sel,
          "normal",        FALSE, FALSE,
          "PE3",           TRUE,  FALSE,
          "normal+filter", FALSE, TRUE,
          "PE3+filter",    TRUE,  TRUE) |>
    mutate(res = map2(skew, sel, one)) |>
    unnest(res) |>
    mutate(dAIC = AIC - min(AIC)) |>
    arrange(AIC)
}

Fits <- BatchEQ |>
  nest(.by = batch_id) |>
  mutate(fit = map(data, \(d) fit_dropthrough(d$eq))) |>
  select(-data) |>
  unnest(fit)

print(Fits, width = Inf)

# Sanity checks worth making explicitly:
#  - S50 should agree ACROSS affected batches. It is a machine constant, not a
#    batch property. Disagreement beyond ~0.5 mm means something is wrong.
#  - Unaffected batches should return selection_active == FALSE, i.e. the model
#    was free to invoke a filter and declined.
#  - `retained` is an extrapolation below the observed minimum and is the least
#    reliable output. Validate against physical throughput, do not rely on it.

Selection <- BatchEQ |>
  nest(.by = batch_id) |>
  mutate(sel = map(data, \(d) fit_nested(d$eq))) |>
  select(-data) |>
  unnest(sel)

print(Selection, n = Inf)   # per batch: which mechanism AIC supports

Fits <- Fits |>
  left_join(Selection |> slice_min(AIC, by = batch_id) |>
              select(batch_id, chosen_model = model), by = "batch_id") |>
  mutate(selection_active = str_detect(chosen_model, "filter"))

Fits |>
  filter(selection_active) |>
  summarise(S50_mean = mean(S50), S50_range = diff(range(S50))) |>
  print()

# ---------------------------------------------------------------------------
# Model comparison at the commercial cuts
# ---------------------------------------------------------------------------

EDGES <- c(0, 52.5, 54.4, 56.1, 59.0, 61.5, 62.5, 66.2, 66.4, 66.6, 66.8, 71.9, 1000)

tvd <- function(o, p) 0.5 * sum(abs(o / sum(o) - p / sum(p)))

grid_probs <- function(dens) {
  cd <- c(0, cumsum(dens / sum(dens)))
  diff(pmax(approx(c(GRID[1] - 0.05, GRID + 0.05), cd,
                   xout = pmin(pmax(EDGES, GRID[1] - 0.05), tail(GRID, 1) + 0.05),
                   rule = 2)$y, 0))
}

Comparison <- BatchEQ |>
  nest(.by = batch_id) |>
  left_join(Fits, by = "batch_id") |>
  mutate(res = pmap(list(data, mu_presented, sd_presented, skew, S50, width),
    \(d, m, s, g, a, b) {
      o <- as.numeric(table(cut(d$eq, EDGES, right = FALSE)))
      tibble(
        TVD_normal = tvd(o, diff(pnorm(EDGES, mean(d$eq), sd(d$eq)))),
        TVD_model  = tvd(o, grid_probs(dpe3(GRID, m, s, g) * retention(GRID, a, b))))
    })) |>
  select(batch_id, res) |>
  unnest(res)

print(Comparison)

# ---------------------------------------------------------------------------
# Diagnostic plot: observed histogram, fitted measured density, retention curve
# ---------------------------------------------------------------------------

Curves <- Fits |>
  mutate(d = pmap(list(mu_presented, sd_presented, skew, S50, width, n),
    \(m, s, g, a, b, n) {
      f <- dpe3(GRID, m, s, g) * retention(GRID, a, b)
      tibble(eq = GRID, fitted = f / sum(f) * n * 0.1 / 0.1,
             S = retention(GRID, a, b))
    })) |>
  select(batch_id, d) |>
  unnest(d)

p_fit <- ggplot() +
  geom_histogram(data = BatchEQ, aes(eq), binwidth = 0.1,
                 fill = "grey70", colour = NA) +
  geom_line(data = Curves, aes(eq, fitted), colour = "#C1272D", linewidth = 0.6) +
  facet_wrap(~batch_id, scales = "free_y") +
  coord_cartesian(xlim = c(40, 85)) +
  labs(title = "Measured EQ with fitted normal x drop-through model",
       x = "EQ (mm)", y = "Count") +
  theme_minimal(base_size = 10)

p_sel <- Curves |>
  semi_join(filter(Fits, selection_active), by = "batch_id") |>
  ggplot(aes(eq, S, colour = factor(batch_id))) +
  geom_line(linewidth = 0.8) +
  coord_cartesian(xlim = c(44, 60)) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Fitted retention curves",
       subtitle = "Independently estimated per batch; should coincide if the filter is a machine constant",
       x = "EQ (mm)", y = "Probability fruit is retained", colour = "Batch") +
  theme_minimal(base_size = 11)

print(p_fit); print(p_sel)

write_csv(Fits, "EQ_dropthrough_fits.csv")
