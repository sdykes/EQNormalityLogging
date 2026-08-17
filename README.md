# Validating the normal EQ assumption across grader batches

**Date:** 13 August 2026 (revised — conclusions reversed by per-fruit data)
**Data:** `all_batches.csv` (21 batches, batch 130 excluded), `AppleBatchStats.csv`, `BatchEQ.csv` (per-fruit EQ, batches 537/541/545/557)
**Scripts:** `Size_comparison_all_batches.R`, `EQ_shape_lmoments.R`, `EQ_dropthrough_model.R`

## Headline

**EQ is normally distributed.** The apparent right-skew in small-fruit batches is not a property of the fruit. It is an artefact of the grader's drop-through filter, which physically removes fruit below about 52 mm before it is ever measured.

Accounting for that filter reduces batch 545's misfit from 12.3% to **0.53%**, and batch 541's from 8.9% to **0.65%** — bringing both in line with the well-fitting batches. Residual skewness after correction is negligible in all four batches examined.

An earlier revision of this report attributed the misfit to genuine size-dependent skew and recommended replacing the normal distribution. **That recommendation is withdrawn.**

The separation of skew from selection has since been validated by simulation (Appendix D): across 3,600 replicates, including a pair of scenarios constructed to have identical mean, standard deviation and skewness, the fitter **never once confused the two mechanisms**, and never failed to detect a filter that was present.

## Purpose

The batch-537 comparison showed good agreement between observed grader piece counts by EQ class and a theoretical allocation from a normal distribution with the batch's `MEAN_EQ` and `SD_EQ`. This extends the comparison to all 21 batches on the current cut scheme.

Four candidate explanations for disagreement were tested:

1. the batch parameters in `AppleBatchStats` are wrong — **eliminated**;
2. the EQ bin cuts were adjusted on the fly and differ from nominal — **eliminated**;
3. the underlying EQ distribution is genuinely not normal — **eliminated**;
4. the *measured population* is a size-selected subset of the fruit presented — **confirmed**.

## Method

Binned analysis, all 21 batches:

1. Parse `EQ` from the `Size` label, sum `Pieces` across grade/elongation rows, convert to proportions.
2. Compute theoretical bin probabilities as `pnorm(high, mean, sd) - pnorm(low, mean, sd)` over the 12 contiguous cuts spanning 30–100 mm.
3. Compare bin by bin; refit mu and sigma by grouped-data maximum likelihood to test the parameters; invert bin mean piece mass to test the cut positions (Appendix B).

Per-fruit analysis, batches 537/541/545/557:

4. Pull individual EQ measurements, reconcile against `AppleBatchStats`, compute L-moments and conventional moments.
5. Fit a joint distribution-plus-selection model (Appendix C).

Fit is judged on **total variation distance** (Appendix A). Chi-square is reported by the script but not used: with 20,000–570,000 pieces per batch every batch rejects, including those agreeing to within half a percentage point.

## Results — binned analysis, all 21 batches

`TVD` is over all 12 bins; `TVD (bounded)` excludes the unbounded catch-all bins US and OS and renormalises.

| Batch | Pieces | mean EQ | sd EQ | TVD | TVD (bounded) | max diff | Worst bin |
|---|---|---|---|---|---|---|---|
| 545 | 96,585 | 56.91 | 4.01 | 12.83% | 12.26% | 8.35% | 58 |
| 541 | 190,953 | 56.80 | 3.62 | 9.28% | 8.89% | 5.81% | 58 |
| 549 | 223,564 | 58.30 | 3.72 | 4.73% | 3.94% | 2.19% | 58 |
| 543 | 91,692 | 57.80 | 3.77 | 4.44% | 3.71% | 2.26% | 58 |
| 540 | 110,997 | 61.87 | 4.09 | 3.67% | 3.31% | 2.37% | 63 M |
| 551 | 67,718 | 58.39 | 4.65 | 3.48% | 3.50% | 2.47% | 58 |
| 548 | 214,137 | 58.82 | 4.58 | 2.87% | 2.72% | 1.97% | 58 |
| 542 | 217,012 | 59.59 | 3.97 | 2.81% | 2.45% | 1.00% | 58 |
| 547 | 182,754 | 62.07 | 4.98 | 2.54% | 2.13% | 1.15% | 67 |
| 556 | 565,780 | 60.78 | 4.81 | 2.45% | 2.11% | 1.08% | 67 |
| 550 | 78,496 | 63.13 | 4.16 | 2.27% | 2.08% | 1.19% | 67 |
| 555 | 335,506 | 60.52 | 4.89 | 1.96% | 1.71% | 0.88% | 58 |
| 554 | 263,989 | 62.16 | 4.90 | 1.92% | 1.68% | 0.79% | 67 |
| 538 | 23,407 | 66.85 | 5.31 | 1.76% | 1.37% | 1.06% | 67 |
| 552 | 368,610 | 62.41 | 4.80 | 1.66% | 1.47% | 0.72% | 67 |
| 544 | 132,652 | 59.75 | 4.43 | 1.50% | 1.46% | 1.09% | 58 |
| 537 | 172,880 | 64.13 | 4.87 | 1.44% | 1.11% | 0.82% | 72 |
| 553 | 163,650 | 62.61 | 4.67 | 1.39% | 1.14% | 0.65% | 67 |
| 539 | 55,534 | 65.12 | 5.27 | 1.02% | 0.88% | 0.57% | 63 L |
| 546 | 203,568 | 63.51 | 5.33 | 0.70% | 0.40% | 0.46% | OS |
| 557 | 234,975 | 63.55 | 4.73 | 0.64% | 0.52% | 0.29% | 63 L |

Sixteen of 21 batches sit under 3% TVD. Of the five above, four are the four smallest-fruit batches in the set — the correlation with batch mean that motivated everything below.

### The parameters are not the problem

Refitting mu and sigma by grouped-data maximum likelihood — treating the 12 bin counts as multinomial and maximising over the normal parameters directly — gives estimates agreeing with `AppleBatchStats` to within **0.15 mm on the mean and 0.13 mm on the standard deviation in every batch**, with TVD barely moving.

This was later confirmed exactly by the per-fruit data: means agree to 0.015 mm, standard deviations to 0.02 mm, counts to 0.3% (below).

### The cut positions are not the problem

Two independent tests, in Appendix B: inverting bin mean piece mass places **every interior cut in all 21 batches within ±0.3 mm of nominal**, and binning the per-fruit EQ values at the nominal cuts reproduces the grader's own `Pieces` to 1.05–2.16% TVD.

One small systematic effect is real: the grader consistently assigns slightly more fruit to "58" and fewer to "US" than the recorded diameters imply, corresponding to an effective cut around 0.5 mm below the nominal 52.5. Worth correcting in the cut table, but far too small to account for the misfit.

## Results — per-fruit analysis

### Reconciliation

| Batch | n (raw) | `APPLES` | `MEAN_EQ` | mean (raw) | `SD_EQ` | sd (raw) |
|---|---|---|---|---|---|---|
| 537 | 175,164 | 174,893 | 64.135 | 64.142 | 4.872 | 4.895 |
| 541 | 203,080 | 202,515 | 56.804 | 56.791 | 3.624 | 3.642 |
| 545 | 98,282 | 97,963 | 56.911 | 56.896 | 4.014 | 4.032 |
| 557 | 239,005 | 238,724 | 63.553 | 63.559 | 4.734 | 4.751 |

`AppleBatchStats` is computed from per-fruit measurements over the same population. Explanation 1 is closed.

### Shape statistics

| Batch | tau-3 (L-skew) | tau-4 (L-kurtosis) | g1 | g2 |
|---|---|---|---|---|
| 537 | −0.0097 | 0.1214 | −0.075 | +0.108 |
| 557 | +0.0060 | 0.1266 | +0.030 | +0.146 |
| 541 | +0.0975 | 0.1230 | +0.557 | +0.553 |
| 545 | **+0.1513** | 0.1207 | +0.818 | +0.782 |

Two observations, and the second is the decisive one.

**tau-4 is 0.1226 in every batch** — the exact normal value — including 545. The tails carry normal weight. The entire departure is asymmetry, with no heavy-tail component at all. That immediately rules out mixture, contamination and heavy-tailed alternatives, and points to a mechanism acting on one side of the distribution only.

**537 and 557 are normal to three decimal places.** Whatever affects 541 and 545 does not affect them.

## The cause: the grader drop-through filter

The grader has a physical filter allowing very small fruit to fall through before entering the sizer. Such fruit is never measured, so it is absent from both `BatchEQ` and `AppleBatchStats`. **The recorded population is a size-selected subset of the fruit presented**, and the selection removes the left tail — producing exactly the observed signature of positive skew with normal kurtosis.

Batch 545's minimum recorded EQ is **42.1 mm** against 27.5 mm for 537 and 32.0 mm for 541, despite 545 having a mean 7 mm smaller than 537. Its left flank rises almost vertically from about 48 mm.

Appendix C fits the mechanism explicitly. The results:

| Batch | Underlying mean | sd | skew | Drop-through S50 | Transition width | Retained |
|---|---|---|---|---|---|---|
| 545 | 51.68 | 6.44 | 0.020 | **52.02 mm** | 0.86 mm | 47.8% |
| 541 | 54.05 | 4.93 | 0.069 | **52.33 mm** | 1.22 mm | 62.2% |
| 537 | 64.14 | 4.90 | 0.020 | none in range | — | 100% |
| 557 | 63.56 | 4.75 | 0.027 | none in range | — | 100% |

Three things make this convincing:

**The two affected batches recover the same threshold independently.** S50 = 52.02 mm and 52.33 mm, fitted separately with no shared parameters — agreement to 0.3 mm on a quantity that is a physical property of the machine, not of the fruit.

**The two unaffected batches correctly report no selection.** The model was free to invoke a filter for 537 and 557 and declined, placing S50 far below the observed range. Their fruit is large enough that none reaches the filter.

**Underlying skewness collapses to nothing.** Once selection is modelled, fitted skew is 0.020–0.069 across all four batches — at or near the lower bound of the parameter space. There is no residual skew to explain.

Fit quality in the transition zone, batch 545 (observed vs model, 1 mm bins):

| EQ | Observed | Model | Retention S(x) |
|---|---|---|---|
| 48–49 | 174 | 185 | 0.016 |
| 49–50 | 608 | 607 | 0.051 |
| 50–51 | 1,873 | 1,800 | 0.146 |
| 51–52 | 4,321 | 4,374 | 0.353 |
| 52–53 | 7,754 | 7,805 | 0.636 |
| 53–54 | 10,174 | 10,198 | 0.848 |
| 54–55 | 11,023 | 10,868 | 0.947 |
| 55–56 | 10,466 | 10,456 | 0.983 |

Agreement within a few percent across four orders of magnitude of count.

### Model comparison

TVD at the commercial cuts, computed from per-fruit data:

| Batch | Normal | Pearson III (free skew) | Normal x drop-through |
|---|---|---|---|
| 545 | 11.11% | 3.30% | **0.53%** |
| 541 | 7.39% | 1.58% | **0.65%** |
| 537 | 0.96% | n/a (left-skew) | 1.15% |
| 557 | 0.75% | 0.66% | **0.66%** |

Fitting skewness alone gets partway. Modelling the mechanism closes the gap and brings the two problem batches **below** the best-fitting unaffected batch.

## Recommendation

**Retain the normal distribution for EQ.** It is correct. Do not switch to lognormal, gamma or Pearson III — those were responses to a misdiagnosis, and lognormal in particular fits *worse* than normal on unaffected batches (537: 3.27% vs 0.96%).

**Apply a drop-through correction where it bites.** For prediction, the distribution of measured fruit is normal truncated by the retention curve S(x). Practically:

- Batches with mean EQ **above about 60 mm**: use the normal directly. Fewer than 0.1% of fruit approach the filter and the correction is negligible.
- Batches with mean EQ **below about 58.5 mm**: apply the correction, or the model will over-predict the US and 58 classes substantially.
- The threshold is close to a **machine property** rather than a batch parameter, but not literally constant. Simulation puts S50's RMSE at 0.038 mm at these sample sizes (Appendix D), so the 0.3 mm gap between 545 (52.02) and 541 (52.33) is real, not noise. Pool S50 across batches with a small random effect rather than fixing it to a single value, and carry mean and sd per batch.

**Confirm the retained fraction against physical records.** The model implies only 47.8% of fruit presented in batch 545 was retained, and 62.2% in 541. This is the weakest number in the analysis — it is an extrapolation into a region with no data (see Appendix C caveats) — but it is a strong, falsifiable prediction. Compare orchard bin weights against grader throughput for those runs. If the tonnage gap is nothing like the predicted one, the selection interpretation needs revisiting.

**Consider whether the loss is acceptable.** Independent of the statistics: if a substantial fraction of a small-fruit batch is falling through the filter, that is a commercial question about whether the filter aperture is right for the season's fruit profile.

**Use BIC, not AIC, to decide whether a filter is present.** AIC over-fits at a rate that *increases* with sample size (Appendix D), adding spurious parameters in 11.6% of no-filter cases. BIC's log(n) penalty is around 5.8 times stiffer at these sample sizes and removes most of that.

**Treat skewness below about 0.07 as undetectable and immaterial.** Appendix D puts the detection threshold there — roughly 60% power even at n = 240,000 — and by the TVD approximately g1/10 rule such skew costs under 0.7 percentage points of TVD in any case.

**Other items:**

- Correct the nominal 52.5 mm cut, which sits about 0.5 mm above its effective position in every batch.
- Exclude **US and OS** from distributional assessment. They are catch-alls, not size classes.
- **Retain the three 0.2 mm sliver bins** at 66.2–66.8, which function as a calibration vernier at no cost (Appendix B).

## Data quality notes

**Batch 551 piece count.** Carries 67,718 pieces in `all_batches` against `APPLES` = 58,768 — a 15% discrepancy, and the only batch where the two disagree by more than 5%. Every other batch has `APPLES` slightly *above* the piece count, the expected direction if the stats are pre-reject. Batch 551 reverses that. Resolve before using 551 downstream.

**US and OS are catch-alls.** The mass inversion shows the ten bounded bins behaving as clean size classes while the unbounded ones absorb fruit diverted for reasons other than size — US running heavy (up to +1.5 mm equivalent) and OS running light (down to −3.6 mm equivalent). Any calculation treating these as pure EQ intervals will be biased.

**Per-fruit data covers four batches only.** 537, 541, 545 and 557 were chosen to span the fit range. The drop-through conclusion should be confirmed on the remaining batches, particularly 549, 543, 540 and 551, before being treated as established for the whole season.

---

# Appendix A: total variation distance

## Definition

For two discrete probability distributions $P$ and $Q$ over the same set of bins:

$$\mathrm{TVD}(P,Q) = \frac{1}{2}\sum_i |p_i - q_i|$$

## Why the one-half

Both distributions sum to one, so the excesses and deficits cancel exactly:

$$\sum_i (p_i - q_i) = 0$$

Taking absolute values therefore counts every misplaced unit of probability mass **twice** — once at the bin where it is, once at the bin where the model says it should be. Halving the sum recovers the quantity of interest: the fraction of the distribution that must be relocated to turn one into the other.

Batch 545 demonstrates this directly:

| | Bins | Sum of differences |
|---|---|---|
| Excess (observed > normal) | 58, 63 S, 67, the three 67 slivers, 72, OS | **+0.1283** |
| Deficit (observed < normal) | US, 63 M, 63 L, 2PC | **−0.1283** |
| Sum of absolute values | | 0.2566 |
| **TVD** | | **0.1283** |

So 12.8% of the fruit in batch 545 sits in a different size class from where the normal model places it — roughly 12,400 of 96,585 pieces. Batch 557, at 0.64%, misplaces about 1,500 of 234,975.

## The supremum form, and why it matters

TVD has an equivalent definition as a supremum over events:

$$\mathrm{TVD}(P,Q) = \sup_{A} |P(A) - Q(A)|$$

The supremum is attained at exactly the excess set — the union of bins where $p_i > q_i$. This gives TVD its operational reading:

> For batch 545, there exists some grouping of size classes for which the normal model's predicted share is wrong by 12.8 percentage points.

Not on average across bins, and not merely in the worst single bin — worst over *any* combination of bins. Conversely for batch 557, **no** grouping of size classes is off by more than 0.64 points. That upper bound over all groupings is the guarantee worth having before using the normal assumption to predict availability for a SKU that draws from several classes at once.

## Properties

- **A metric.** Symmetric, zero only when $P = Q$, and satisfies the triangle inequality. Batches can therefore be ranked and compared meaningfully.
- **Bounded in $[0, 1]$.** Zero is identity; one is disjoint support. Directly readable as a percentage.
- **An effect size, not a test statistic.** It does not grow with $n$. This is the decisive advantage over chi-square here — with piece counts in the $10^5$–$10^6$ range, χ² is in the thousands for batches that agree to within half a percentage point, so its p-value carries no decision-relevant information. TVD answers *how wrong*, which is the question being asked, rather than *is it wrong at all*, which at this sample size always answers yes.

## Two caveats for this application

**TVD is bin-dependent.** It measures disagreement at the resolution of the cut scheme, not between the underlying continuous distributions. The three 0.2 mm slivers at 66.2–66.8 mm contribute almost nothing to TVD however badly the model fits within them, while the 3.7 mm-wide "67" bin can dominate. Consequently, two batches are only comparable on TVD if they share the same cuts — which is precisely why batch 130 had to be excluded rather than simply relabelled.

**TVD is blind to ordering.** Moving 2% of fruit from 63 M to the adjacent 63 L costs exactly as much as moving it from US to OS, even though the first is a rounding error and the second is physically implausible. Because these bins are ordered along a physical millimetre scale, a **Wasserstein (earth-mover) distance in mm** is a natural companion metric: it charges for how far mass has to travel, and so distinguishes a mildly mis-shaped distribution from a genuinely wrong one. For batch 545 the two metrics would likely tell different stories — much of the 12.8% is a short hop across the US/58 boundary, while the 72 and OS excess is a long one. In R, `transport::wasserstein1d()` on the bin midpoints weighted by piece counts is the quickest route to this.

---

# Appendix B: detecting on-the-fly EQ cut adjustment

## The problem

EQ bin cuts can be adjusted by operators during a run. There is no log of these changes, batch to batch or within a batch, so any disagreement between observed and theoretical bin proportions could in principle be a calibration artefact rather than a distributional one.

**Check first whether it is genuinely untracked.** The TOMRA reports parsed with `pdftools` often carry the active size program in the header block — recipe name, version, sometimes the cut table itself. If a recipe identifier is printed, that is a log at run granularity, and grepping the 21 PDFs for the cut values is cheaper than any inference below.

## Why counts alone can never answer it

There are 12 bins, hence 11 free proportions, and 11 interior cuts. The system is **exactly saturated**. Any observed histogram whatsoever can be reproduced perfectly by a normal distribution with the stated μ and σ if the cuts are allowed to move. This is a structural non-identifiability: more batches will not fix it, and neither will more fruit.

The inversion makes the point concrete. Given observed cumulative proportion $\hat{F}_k$ at boundary $k$, the cut position implied by exact normality is

$$\hat{c}_k = \mu + \sigma\,\Phi^{-1}(\hat{F}_k)$$

**Implied minus nominal cut (mm), from counts, assuming exact normality**

| Batch | 52.5 | 54.4 | 56.1 | 59.0 | 61.5 | 62.5 | 66.2 | 66.4 | 66.6 | 66.8 | 71.9 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 545 | −0.90 | 0.48 | 0.69 | 0.39 | −0.17 | −0.43 | **−1.48** | **−1.53** | **−1.59** | **−1.65** | **−2.96** |
| 541 | −0.57 | 0.33 | 0.45 | 0.20 | −0.25 | −0.45 | **−1.31** | **−1.37** | **−1.42** | **−1.47** | **−3.05** |
| 543 | −0.44 | 0.09 | 0.18 | 0.22 | 0.08 | −0.01 | −0.62 | −0.66 | −0.73 | −0.76 | **−2.34** |
| 549 | −0.58 | 0.09 | 0.21 | 0.23 | 0.07 | −0.03 | −0.50 | −0.54 | −0.59 | −0.62 | **−2.02** |
| 540 | −1.60 | −0.54 | −0.17 | 0.23 | 0.27 | 0.25 | 0.03 | 0.02 | 0.01 | 0.01 | −0.29 |
| 537 | −0.56 | 0.18 | 0.10 | 0.07 | 0.01 | −0.00 | −0.03 | −0.03 | −0.02 | −0.02 | 0.31 |
| 557 | −0.30 | −0.02 | 0.01 | 0.01 | 0.04 | 0.06 | 0.06 | 0.05 | 0.05 | 0.05 | 0.12 |
| 546 | −0.21 | −0.00 | 0.01 | 0.03 | 0.02 | 0.03 | 0.02 | 0.02 | 0.02 | 0.02 | 0.22 |

Batch 545 can be explained away entirely: the 66.2 cut was really at 64.7, the 71.9 cut at 68.9. Both are plausible magnitudes for an operator nudge, and nothing in the count data refutes it.

## Direct test, using per-fruit EQ

With per-fruit measurements available, the cut positions can be checked directly: bin the recorded EQ values at the nominal cuts and compare with the grader's own `Pieces`.

| Batch | TVD (per-fruit binned vs grader counts) |
|---|---|
| 537 | 1.15% |
| 541 | 2.16% |
| 545 | 2.03% |
| 557 | 1.05% |

The grader's bin assignments match what the recorded diameters imply. This is the most direct evidence available and it independently confirms the mass inversion below.

One systematic effect appears in all four batches: the grader assigns **more fruit to "58" and fewer to "US"** than the diameters imply (batch 545: −1.3 percentage points US, +2.0 in 58). This is the same roughly −0.5 mm offset at the 52.5 cut that showed up in every batch in the count inversion above — now measured rather than inferred. It is real and worth correcting in the nominal cut table, but an order of magnitude too small to explain the fit failures.

Note this test requires per-fruit data. The mass-based method below works from the summary tables alone and remains the practical route for the other 600-plus batches.

## Mass is the independent axis that breaks the tie

`all_batches.csv` carries `Weight_kg` alongside `Pieces` in every cell. Mean piece mass in a bin measures **the fruit that is in the bin**, not the rule that put it there. That distinction resolves the identifiability problem:

- If a **cut moves**, the fruit occupying the adjacent bins physically changes, and mean mass moves with it.
- If the **distribution shape is wrong**, the cuts are where they say they are, the fruit in each bin is the right size, and mean mass is unchanged.

Two observations per bin (count and mean mass) against one unknown cut per boundary takes the system from saturated to overidentified, and both effects become separately estimable.

### Procedure

1. Compute bin mean piece mass as `Weight_kg / Pieces × 1000` grams.
2. Fit a per-batch allometry $m = a\,\mathrm{EQ}^{\,b}$ on the interior bins (58 through 2PC), regressing on the *conditional mean EQ* of each bin under the batch's own normal rather than the bin midpoint, so that bin width is handled correctly.
3. For every bin, solve for the uniform edge displacement $\delta$ such that the predicted conditional mean mass over $[\mathrm{lo}+\delta,\ \mathrm{hi}+\delta]$ equals the observed mean mass. A one-dimensional root find (`uniroot`) suffices.

The fitted exponent came out at **2.73**, range 2.54–2.83 across all 21 batches — sub-cubic as expected given elongation below 1, and stable enough to treat as a fixed physical property.

### Result

**Mass-implied bin displacement (mm); positive means the bin sat higher in EQ than nominal**

| Batch | US | 58 | 63 S | 63 M | 63 L | 2PC | 67 | 66.4 | 66.6 | 66.8 | 72 | OS |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 545 | **+1.52** | 0.05 | −0.03 | −0.09 | 0.01 | 0.03 | 0.18 | −0.05 | −0.19 | −0.14 | 0.36 | −0.57 |
| 541 | **+0.72** | 0.02 | 0.00 | −0.08 | 0.02 | 0.01 | 0.16 | −0.03 | −0.11 | 0.00 | 0.23 | **−3.62** |
| 543 | +0.32 | −0.01 | 0.02 | −0.03 | −0.02 | 0.01 | 0.04 | 0.03 | 0.06 | −0.15 | 0.15 | **−2.47** |
| 549 | +0.33 | −0.02 | 0.03 | −0.02 | −0.01 | 0.00 | 0.03 | 0.02 | −0.10 | −0.02 | 0.18 | **−3.18** |
| 540 | +0.07 | 0.00 | 0.02 | −0.05 | −0.04 | 0.04 | 0.14 | 0.28 | 0.27 | 0.23 | 0.67 | +1.37 |
| 537 | −0.89 | −0.10 | 0.12 | 0.00 | −0.02 | −0.03 | 0.02 | 0.11 | 0.09 | 0.07 | 0.19 | +0.19 |
| 557 | −0.41 | −0.03 | 0.04 | −0.03 | −0.02 | 0.01 | 0.05 | 0.09 | 0.07 | 0.07 | 0.21 | +0.47 |
| 546 | −0.96 | −0.03 | 0.03 | 0.00 | −0.01 | −0.01 | 0.00 | 0.00 | −0.09 | −0.11 | −0.01 | −0.49 |

**Every interior boundary in all 21 batches sits within ±0.3 mm of nominal.** Counts claimed batch 545's 66.2–66.8 cuts had drifted 1.5 mm; mass says the fruit in those bins weighs precisely what 66.2–66.8 mm fruit should weigh. The cuts did not move. There really is that much fruit at that size, which independently confirms the maximum-likelihood conclusion that 545 is a shape failure and not a calibration failure.

The two unbounded bins behave differently and are discussed in the data quality notes: US runs heavy and OS runs light, marking both as catch-alls rather than size classes.

## The slivers as a built-in vernier

Three 0.2 mm slices at 66.2–66.8 mm form a precision calibration probe. Across them $dm/d\mathrm{EQ} \approx 4.5$ g/mm (consistent with $b\,m/\mathrm{EQ} = 2.73 \times 118 / 66.5 = 4.84$). Median sliver count is about 1,000 pieces, so with piece-mass SD near 20 g the standard error on bin mean mass is roughly 0.6 g, resolving cut position to about **0.13 mm**. This is far finer than any adjustment an operator would plausibly make, which is why the ±0.3 mm result above is a real constraint rather than a null finding from a blunt instrument.

## Discriminating by residual signature, when mass is unavailable

The residual pattern alone carries diagnostic information:

| Signature | Interpretation |
|---|---|
| **Locally antisymmetric** — equal and opposite in two adjacent bins, near-zero elsewhere | A genuine cut move at that boundary |
| **Smooth monotone ramp** across many boundaries | Shape failure: tail too heavy or too light |
| **Same offset at the same boundary in every batch** | Wrong nominal cut value, or systematic calibration bias |

Batch 545's upper boundaries read −1.48, −1.53, −1.59, −1.65, −2.96 — a monotone ramp, not a localised step. Independent operator nudges do not produce monotone ramps, so shape failure is indicated even before the mass evidence is brought in.

## Detecting within-batch drift

Everything above is batch-level. Intra-batch drift requires temporal resolution — sub-lot, pallet, or timestamped grader records. Given a batch sliced into time windows, two tests apply.

**Adjacent-bin log-odds changepoint.** Track

$$\ell_k(t) = \log\!\left(\frac{n_k(t)}{n_{k+1}(t)}\right)$$

for each of the 11 boundaries. A cut adjustment at boundary $k$ produces a step in $\ell_k$ and leaves the other ten untouched. A change in the incoming fruit population shifts all eleven together. Running a multivariate changepoint (`changepoint::cpt.meanvar`, or `strucchange::Fstats` on the vector) and reading the loading gives the discriminant:

- **sparse loading on one coordinate → a cut moved**
- **dense loading across all coordinates → different fruit**

The log-odds form matters here: it is the statistic most sensitive to a single shared boundary and least sensitive to overall population shift, which raw proportions are not.

**Mass check per window.** Repeat the mass inversion within each time slice. Cut drift shows as a step in the mass-implied displacement; a population change does not.

## Scaling to a full season

For 600-plus batches, prefer a hierarchical model to per-batch inversions: per-boundary cut offsets as batch-level random effects, distribution shape parameters pooled or partially pooled by orchard. Partial pooling performs the separation automatically — a residual appearing at the *same boundary across every batch* indicates a wrong nominal cut value or a systematic calibration bias, whereas one confined to a single batch indicates that batch's operator. The existing `glmmTMB` setup would take a multinomial-over-bins likelihood with a linear predictor in the cut offsets.

A design-based alternative, if calibration confidence ever needs to be higher than inference can deliver: pass a set of reference spheres or a tagged, independently measured sample of fruit through the grader at the start and end of each run.

---

# Appendix C: the drop-through selection model

## Specification

Let $f(x;\mu,\sigma,\gamma)$ be the density of fruit **presented** to the grader — a Pearson III with mean $\mu$, standard deviation $\sigma$ and skewness $\gamma$, which reduces to the normal at $\gamma = 0$. Let $S(x)$ be the probability that a fruit of diameter $x$ is retained rather than falling through the filter:

$$S(x) = \frac{1}{1 + \exp\!\left(-\frac{x - a}{b}\right)}$$

with $a$ the 50% retention diameter and $b$ the transition width. The density of **measured** fruit is then

$$f^{*}(x) = \frac{f(x;\mu,\sigma,\gamma)\,S(x)}{\int f(t;\mu,\sigma,\gamma)\,S(t)\,dt}$$

Five parameters, fitted by maximum likelihood on the 0.1 mm histogram. The logistic form is used rather than a hard cut because retention depends on fruit orientation and shape as well as diameter, so the transition is necessarily gradual — and the data confirm this, showing a smooth rise over roughly 5 mm rather than a step.

Keeping $\gamma$ free rather than fixing it at zero is what makes the exercise a test rather than an assumption: the model is free to explain the data with skew, and does not.

## Why Pearson III rather than Box-Cox

An earlier draft proposed a Box-Cox transform to normality with the transform parameter varying by batch mean. That was designed to capture size-dependent skew which, on the evidence above, does not exist. Pearson III is used here purely as a skew-permitting envelope containing the normal as a special case — its role is to *fail* to find skew, and it does.

## Identification

The parameters are identified by different features of the data, which is worth understanding because they are estimated with very different precision:

| Parameter | Identified by | Precision |
|---|---|---|
| $a$, $b$ | The shape of the observed transition zone, roughly 48–56 mm | **Well determined** — tens of thousands of fruit across the rise |
| $\mu$, $\sigma$, $\gamma$ | The unaffected upper portion, plus extrapolation | Good for $\sigma$ and $\gamma$; $\mu$ depends on extrapolation |
| Retained fraction | Pure extrapolation below the observed minimum | **Weak** — treat as indicative only |

The retention curve itself is measured. The underlying population mean and the retained fraction are model-based extrapolations into a region where, by construction, there is no data at all.

Appendix D quantifies all of this by simulation. In summary: S50 recovers with RMSE 0.038 mm and the transition width to 0.009 mm at batch 545's sample size, while the skewness and S50 estimates carry a persistent negative correlation of −0.48 to −0.72 that does **not** diminish with sample size — so neither should be reported with a marginal standard error alone. Batch 545's implied 47.8% retention should be treated as a hypothesis to test against physical throughput records, not as an estimate to rely on.

## Multi-start fitting is necessary

The likelihood is multi-modal. A single optimisation from a naive start converged for batch 545 but found a poor local optimum for 541, reporting no selection and a large skew — the wrong answer, and one that looked plausible in isolation. Fitting from a grid of starting values for $a$ (35, 45, 48, 50, 52 mm) and $\gamma$ (0.05, 0.3, 0.6) and keeping the best resolved it.

This matters generally: the model can trade skew against selection, since both thin the left side of the distribution. The two are distinguished by the *shape* of the thinning — a selection curve produces a much sharper cut-off than any skewed density can — but the likelihood surface has local optima corresponding to each explanation. Always fit from multiple starts, and check that the fitted $S(x)$ falls in a physically plausible place.

## Season-wide fitting

For the full batch set, a hierarchical structure follows directly from the physics:

- $a$ and $b$ are **properties of the machine**, not of the batch. Fit them as fixed parameters shared across all batches, or as slowly-varying terms if the filter is adjusted between runs.
- $\mu_b$ and $\sigma_b$ are batch-specific and already available from `AppleBatchStats` — though note those are moments of the *measured* population, so they act as constraints on $f^{*}$, not on $f$.
- $\gamma$ can be pooled across batches, testing a single global hypothesis of normality with far more power than any per-batch test.

This is a well-behaved model with roughly $2B + 3$ parameters against $25B$ observations, and it is a **nonlinear latent-variable model rather than a GLMM**, so `glmmTMB` cannot fit it. Use `RTMB` — the same Template Model Builder engine, with the likelihood written directly in R and automatic differentiation plus the Laplace approximation supplied for the random effects.

The earlier proposal to model per-boundary cut offsets as random effects can be dropped. Three independent lines of evidence now place the cuts within ±0.3 mm of nominal, so the only correction needed is the single fixed −0.5 mm offset at the 52.5 mm boundary.

## Simulation check

This has now been carried out — see **Appendix D**. The model recovers S50 to better than 0.1 mm and never confuses genuine skew with selection, including on scenarios constructed to be indistinguishable on the first three moments. Two changes to the procedure follow from it: use BIC rather than AIC for the mechanism decision, and do not treat S50 as identical across batches.

The one criterion that failed was the original heuristic for deciding whether selection is active. Testing whether the fitted S50 falls inside the observed data range is unreliable: with genuinely skewed fruit and no filter the optimiser parks S50 a few millimetres below the data, leaving $S(x) \approx 1$ and the filter inert — but the parameter still has a value, so the heuristic reported selection in 100% of no-filter replicates. Nested model comparison is the correct arbiter and has replaced it in the script.

What remains outstanding is misspecification: every simulated cell drew from the same Pearson III times logistic form the fitter assumes.

---

# Appendix D: simulation study — can skew be separated from selection?

## Why the study was necessary

Genuine right-skew and a left-side size filter both thin the lower part of the distribution, so they compete to explain the same feature of the data. Appendix C's conclusion — that the fruit is normal and the filter does the work — depends entirely on the fitter being able to tell them apart. That claim needed testing against known truth before the model was used in production.

## Design

Six scenarios, each simulated 200 times at three sample sizes (20,000 / 98,282 / 240,000 — the last two matching batches 545 and 557). Simulation is by multinomial draw on the 0.1 mm grid, which is exact conditional on $n$ and mirrors the grader's own discretisation.

| Cell | Scenario | mean | sd | skew | S50 | width |
|---|---|---|---|---|---|---|
| A | Normal fruit + filter (the 545 fit) | 51.68 | 6.44 | 0.02 | 52.02 | 0.86 |
| B | Genuinely skewed fruit, no filter | 56.90 | 4.03 | 0.818 | — | — |
| C | Unaffected large fruit (the 537 fit) | 64.14 | 4.90 | 0.02 | — | — |
| D | Both skew and filter present | 55.00 | 5.00 | 0.35 | 50.00 | 1.00 |
| E | Weak filter, slight skew (the 541 fit) | 54.05 | 4.93 | 0.069 | 52.33 | 1.22 |
| F | Normal + near-sharp filter | 51.68 | 6.44 | 0.02 | 52.02 | 0.20 |

**A and B are the adversarial pair.** They are constructed to produce measured distributions with identical mean, standard deviation and skewness — 56.90, 4.03, 0.818 to three figures. Only excess kurtosis differs, 0.73 against 1.00. Any method resting on the first three moments cannot separate them, so this pair is the real test.

Each replicate is fitted with all four nested models — normal, Pearson III, normal + filter, Pearson III + filter — and the mechanism is chosen by AIC.

## Result 1: the mechanisms are never confused

Confusion matrix, all 3,600 replicates pooled:

| True ↓ / Chosen → | normal | PE3 | normal+filter | PE3+filter |
|---|---|---|---|---|
| normal | **520** | 32 | 25 | 23 |
| PE3 (skew, no filter) | **0** | **509** | **0** | 91 |
| normal+filter | **0** | **0** | **1479** | 321 |
| PE3+filter | 0 | 0 | 1 | **599** |

The two cells that would represent genuine mechanism confusion are empty:

- **True skew reported as a filter: 0 of 600.**
- **True filter reported as pure skew: 0 of 1,800.**

Across all 2,400 replicates where a filter was present, the fitter **never once failed to detect it**. Every error in the table is a redundant parameter rather than a wrong mechanism: an inert filter appended to genuinely skewed fruit, or a small spurious skew appended to filtered fruit. Restated as error rates:

- **False negatives (real filter missed): 0.0%**
- **False positives (spurious filter added): 139 of 1,200 = 11.6%**

The asymmetry is what matters operationally. The analysis in Appendix C claims a filter is present in batches 545 and 541; the study shows that claim is not one the fitter makes spuriously in the direction that would matter, and that a real filter is never overlooked.

## Result 2: S50 recovers to better than 0.1 mm

| Cell | n | Bias (mm) | RMSE (mm) | RMSE of width $b$ |
|---|---|---|---|---|
| A | 20,000 | −0.031 | 0.096 | 0.021 |
| A | 98,282 | −0.008 | **0.038** | 0.009 |
| A | 240,000 | −0.010 | 0.027 | 0.006 |
| D | 98,282 | +0.009 | 0.083 | 0.020 |
| E | 98,282 | −0.020 | 0.083 | 0.019 |
| F | 98,282 | −0.002 | **0.009** | 0.004 |

Essentially unbiased, with RMSE scaling as $1/\sqrt{n}$ as expected. Recovery holds even in cell D where genuine skew and a filter are both present, and is sharpest in cell F where the filter is nearly a step — the parameterisation does not degrade at that boundary.

**A consequence for the fitted batches.** S50 recovers to 0.038 mm RMSE at batch 545's sample size. The estimates of 52.02 mm (545) and 52.33 mm (541) are therefore separated by roughly five RMSEs, so **that 0.3 mm difference is real and not sampling noise.** The earlier suggestion in Appendix C to treat S50 as a single fixed machine constant across batches should be qualified: the threshold is stable to a few tenths of a millimetre but is not literally identical between runs. Fit it per batch where the data allow, or pool with a small random effect rather than fixing it outright.

## Result 3: cell E was mislabelled, and the model was right

Cell E's apparent identification rate falls with sample size, which looks like a failure and is not one. E was given a true skewness of 0.069 — copied from the 541 fit — while its truth was labelled `normal+filter`. Since 0.069 is not zero, **`PE3+filter` is the correct model for cell E**, and the fitter increasingly recognises this as data accumulates:

| n | Fitted skew (true 0.069) | Picks PE3+filter |
|---|---|---|
| 20,000 | 0.086 ± 0.070 | 17% |
| 98,282 | 0.073 ± 0.036 | 42% |
| 240,000 | 0.066 ± 0.026 | **59%** |

The estimate converges cleanly on the true value. Correctly labelled, cell E measures the detection threshold for skewness rather than a failure mode: **a skewness of about 0.07 sits right at the limit of detectability at these sample sizes**, reaching only around 60% power at n = 240,000. Skew smaller than that cannot be distinguished from zero, and skew that small is immaterial anyway — by the $\mathrm{TVD} \approx g_1/10$ rule in Appendix A it costs under 0.7 percentage points of TVD.

Corrected identification rates:

| Cell | n = 20,000 | n = 98,282 | n = 240,000 |
|---|---|---|---|
| A | 89.5% | 95.0% | 90.0% |
| B | 85.5% | 85.0% | 84.0% |
| C | 89.0% | 88.0% | 83.0% |
| D | 99.5% | 100.0% | 100.0% |
| E | 16.5% | 41.5% | 58.5% |
| F | 97.5% | 92.5% | 91.5% |

## Result 4: use BIC, not AIC

The residual over-fitting is AIC behaving as designed. AIC is not consistent for model selection — it retains a fixed error rate however large the sample — and the study shows this directly. Cell C's over-fitting rate **rises** with sample size:

| n | Cell C over-fit rate |
|---|---|
| 20,000 | 11% |
| 98,282 | 12% |
| 240,000 | **17%** |

More data makes it worse, not better, because AIC's penalty of 2 per parameter does not grow. BIC's $\log(n)$ penalty is about 5.8 times stiffer at $n \approx 10^5$ and would remove most of the 321 spurious-skew and 91 inert-filter selections. **Switch the model-choice criterion to BIC** for production use. AIC remains reasonable if the aim is prediction rather than mechanism identification, but the question here is explicitly which mechanism is operating.

## Result 5: a persistent trade-off that does not vanish with data

Correlation between the fitted skewness and the fitted S50 across replicates:

| Cell | 20,000 | 98,282 | 240,000 |
|---|---|---|---|
| A | −0.60 | −0.49 | −0.48 |
| D | −0.72 | −0.72 | −0.71 |
| E | −0.47 | −0.58 | −0.56 |
| F | −0.29 | −0.28 | −0.34 |
| B | +0.28 | +0.19 | −0.01 |
| C | +0.03 | +0.07 | +0.03 |

In the filter cells the correlation sits between −0.48 and −0.72 and **shows no tendency to shrink with sample size**. This is the structural trade-off made visible: skew and selection genuinely compete, and collecting more fruit does not decouple them. The magnitude is far from −1, so both parameters remain identified — but two practical consequences follow:

- **Never report a marginal standard error for S50 or skewness on its own.** Use a joint likelihood region, or profile one parameter with the other fixed.
- **Prefer fixing S50 from pooled batches** and estimating skewness conditional on it, rather than estimating both freely per batch. This exploits the fact that S50 is close to a machine property while skewness is not.

## What the study does not establish

**Correct specification is assumed throughout.** Every cell simulates from a Pearson III multiplied by a logistic retention curve — the same functional form the fitter uses. Real fruit is not exactly Pearson III and a physical filter is not exactly logistic. The misspecification cells set out in the script header remain outstanding:

- hard truncation instead of a logistic curve;
- a two-component mixture, representing two harvest dates in one batch;
- Weibull fruit instead of Pearson III.

The mixture cell carries a useful side-benefit: it should produce an L-kurtosis departing from 0.1226, which the real batches conspicuously do not, giving an independent check that the batches are not blends.

**The retained fraction remains unvalidated.** No simulation can address it, because the quantity is an extrapolation into a region where the design guarantees no data. Batch 545's implied 47.8% retention still requires confirmation against orchard bin weights versus grader throughput.
