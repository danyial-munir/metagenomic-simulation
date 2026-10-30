# Phase 4 analysis plan: seasonal scenario power study (FROZEN)

Status: **frozen before any power result was computed.** The notebook `05_power.ipynb` prints the SHA-256 of this file and
stops if it differs from the value recorded in `CLAUDE.md`. Any later change must be added as a dated amendment at the end,
with its reason, never by editing the text above it.

## 1. Question (RQ3)

If a seasonal driver changes the relative abundance of some skin taxa with amplitude *A*, how large must *A* be, and how many
samples *n* are needed, before standard tools detect it? The seasonal signal is **imposed by simulation**; nothing here
estimates real seasonality.

## 2. Data-generating model

- **Primary simulator: M3** (logistic-normal multinomial, μ and Σ from `model_params.npz`, λ = 0.05), D = 22 categories
  (20 named taxa, `Unassigned`, reference `Other`). **Sensitivity simulator: M4** (M3 + structural zeros, `model_params_m4.npz`).
  The choice follows the pre-set rule of Phase 3.
- **Design:** one sample per subject, months balanced (n/12 subjects per month, random order). Cross-sectional, like the real study.
- **Driver:** `s(m) = cos(2π(m − φ)/12)`, φ = 7 (peak in July, trough in January).
- **Effect:** for each simulated study, *K* of the 20 named taxa are chosen at random (so results average over abundant and
  rare taxa); their ALR coordinate gets `z_ij += A · s(m_i) · sign_j`, signs alternating +1, −1, +1, … (some taxa peak in summer,
  some in winter). Peak-to-trough fold change relative to the reference = e^{2A}.
- `Unassigned` and `Other` are never affected. Under M4 presence/absence is not seasonal (simplification).
- **Reads:** Multinomial(depth, softmax(z)).
- simplification: no simulated temperature covariate (it would add a second, correlated driver without changing the question).

## 3. Grid (all values fixed now)

| Factor | Values | Main setting |
|---|---|---|
| Amplitude *A* (log-ratio units) | 0, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5 | all |
| Sample size *n* | 48, 96, 192, 288, 384 (multiples of 12; 288 ≈ the real study, 273) | all |
| Affected taxa *K* | 1, 3, 5 | 3 |
| Depth (reads/sample) | 1000, 10 000 | 1000 |
| Overdispersion (Σ multiplied by) | 0.5, 1, 2 | 1 |
| Simulator | M3, M4 | M3 |

Scenarios:
1. **Main:** M3, K = 3, depth 1000, Σ × 1, full *A* × *n* grid (35 points).
2. **M4 sensitivity:** same full *A* × *n* grid with M4 (35 points).
3. **One-at-a-time sensitivity** (M3, all *A*, *n* ∈ {96, 288}): K = 1 and K = 5; depth 10 000; Σ × 0.5 and Σ × 2 (70 points).
4. **Negative control:** M3, main setting, A = 1.0, every *n*, month labels **shuffled** before testing (5 points).

Total 145 grid points.

## 4. Detection tools (pre-specified)

All on the 21 parts without `Unassigned`; CLR with pseudocount 0.5 (sensitivity checked in Phase 1).

1. **PERMANOVA, Aitchison distance, factor = season** (DJF, MAM, JJA, SON; 4 groups), 199 permutations. *(primary)*
2. **PERMANOVA, Bray–Curtis, factor = season**, 199 permutations. *(secondary)*
3. **Per-taxon cosinor on CLR** for the 20 named taxa: `clr_ij = b0 + b1 cos(2πm/12) + b2 sin(2πm/12) + ε`, F-test of
   b1 = b2 = 0 (2 df), **Benjamini–Hochberg at q = 0.05** over the 20 taxa. *(primary)*

Significance level α = 0.05. 199 permutations make the permutation test exact at α = 0.05 (α·(199+1) = 10 is an integer).

## 5. Outcomes

- **Power** of each test = share of replicates with p ≤ 0.05 (cosinor: at least one BH rejection).
- **Type I error** = the same at A = 0. Accepted if within 0.05 ± 2·MCSE = [0.036, 0.064] (R = 1000); otherwise flagged and reported.
- **Cosinor per-taxon:** sensitivity (true positives / K) and empirical FDR (false rejections / all rejections, over
  replicates with ≥ 1 rejection). Note: shifting K taxa changes the CLR of *all* taxa a little (closure), so some FDR inflation
  above q is expected and will be reported, not hidden.
- **Minimum detectable effect (MDE):** smallest *A* with power ≥ 0.80, by linear interpolation between grid values of *A*;
  reported as "> 1.5" if not reached.
- **Power vs n** at each *A*; power heatmap over (*A*, *n*).
- **Negative control:** rejection rate with shuffled labels must be ≈ 0.05 (same acceptance band as type I error).

## 6. Replicates and Monte Carlo error

R = **1000** replicates per grid point. Monte Carlo standard error of a rejection rate p is √(p(1−p)/R): at most 0.016
(p = 0.5), 0.007 at p = 0.05. Power curves are drawn with 95 % Wilson intervals.
Runtime estimate: about 12 ms per replicate → about 12 minutes on one core for 145 000 replicates (parallelised over grid points).

## 7. Randomness

Each replicate uses its own generator `np.random.default_rng([SEED, 5, point_id, replicate])`, SEED = 20261030, so results
do not depend on execution order or the number of parallel workers.

## 8. Reporting rules

- All grid points are reported, including negative or null results; no grid point, test or replicate is dropped.
- Primary conclusions use the main scenario and the two primary tools; everything else is labelled sensitivity/secondary.
- No change to tools, thresholds or grid after seeing results, except as a dated amendment below with its reason.

## Amendments

Verification: the text above this heading is unchanged since freezing. Its SHA-256 is
`53a474f05a8a813d429dcb5ec30c89fe46aba5aa1f06548b65a9bf63fc49d35d` (the whole file, before this amendment, had
`dfd32c3c0aa12e06cb0c3371384072f52ffe8bc7bd236bf9fb9e36903f330285`). `05_power.ipynb` checks the part above this heading.

**Amendment 1 (2026-10-30, added AFTER seeing the first results).** The FDR in section 5 was defined as the mean false-discovery
proportion over replicates with at least one rejection. That is the *positive* FDR, which Benjamini–Hochberg does not control,
so it cannot be compared with q = 0.05 (at small A almost every rejecting replicate is a chance rejection, so it is close to 1 by
construction). Added: the **standard FDR = E[false rejections / max(rejections, 1)]**, counting replicates without rejections
as 0, which is the quantity BH controls. The pre-specified quantity is still reported unchanged. No test, threshold, grid value
or replicate changes; the simulations are not rerun with different settings (same seeds, identical results).
