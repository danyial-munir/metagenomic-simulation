# Concept sheet (running list, for the oral defence)

Format per concept: plain language, key formula, assumptions, what breaks, likely examiner questions.
Entries are added phase by phase. **Phase 0 entries below.** Empty headings are placeholders for later phases.

## Phase 0 – data audit

### Rarefaction and what the numbers in the file are
- **Plain language:** every sample was cut to the same number of sequencing reads (1000). Each cell is therefore *count / 1000*.
- **Evidence in our data:** every value × 1000 is an integer, smallest non-zero value = 0.001, `Reads = 1000` for all subjects.
- **Why it matters:** with exact counts out of a known total N = 1000 we can use count models (multinomial, Dirichlet-multinomial) instead of treating the numbers as continuous.
- **What breaks / caveat:** rarefying discards reads and information (McMurdie & Holmes 2014). Here it is already done, so we model the data as given and say so.
- **Likely questions:** Why is a column sum of ~1 informative? (Compositional: only ratios carry information.) Why not rescale to percentages? (Loses the count scale, hence the sampling noise.)

### Missing mass and the `Unassigned` category
- 48 of 273 samples sum to < 1000 assigned reads. Decision: add `Unassigned = 1000 − sum` so each sample is a full multinomial of 1000 reads. Sensitivity: drop samples with < 80 % assigned.
- **Likely question:** why not just renormalise? (Renormalising inflates the assigned taxa and invents reads that were never observed.)

### Sparsity and prevalence filter
- 93 % of cells are zero. **Prevalence** = share of samples where a taxon has ≥ 1 read. Taxa below the threshold are merged into `Other`. The sensitivity table shows what changes.
- **Likely question:** why is the threshold arbitrary, and how do you know results are robust to it? (Sensitivity table; compare model conclusions across thresholds.)

### Missing data in covariates
- No imputation. Complete-case analysis for each model, always reporting n. "Unknown" is NaN, but `0 washes/week` is a real value, not missing.
- **Likely question:** why not impute? (Imputation would hide that 42 % of wash frequency is unknown and add made-up information to a small sample.)

## Phase 1 – EDA

### Compositional data, closure, CLR
- **Plain language:** each sample is divided by its total, so only ratios between taxa are informative. If one taxon rises, the others fall in relative terms.
- **Formula:** `clr(x)_i = log x_i − (1/D) Σ_j log x_j`; rows sum to 0; Aitchison distance = Euclidean distance on CLR.
- **Assumptions:** zeros are replaced (pseudocount 0.5 here) so the logs exist.
- **What breaks:** the pseudocount is arbitrary; with 93 % zeros it matters. Checked: Spearman correlation of pairwise distances across pseudocounts 0.1–2 is ≥ 0.89 (min over pairs), so the ordering of samples is fairly stable, but absolute distances change.
- **Questions:** Why not correlations on proportions? (Spurious negative correlations from closure.) Why CLR not ALR? (Symmetric in the parts; ALR needs a reference taxon.) Why can't you take log of 0? (Sampling zeros: 1000 reads miss rare taxa.)

### Alpha diversity
- Richness = number of taxa with ≥ 1 read. Shannon `H = −Σ p ln p`. Gini–Simpson `1 − Σ p²` = probability two random reads are from different taxa.
- **Caveat found:** richness is strongly tied to the number of assigned reads (Spearman −0.69 in our data). Hypothesis, not verified: the source table dropped some rare taxa, so samples with many rare taxa lost more mass. Shannon and Simpson are less sensitive to rare taxa.
- **Questions:** richness vs Shannon? (Richness counts presence; Shannon also weights evenness.) Why is richness not comparable across depths? (Deeper sequencing finds more rare taxa; hence rarefaction.)

### Beta diversity, PCoA, PERMANOVA
- Bray–Curtis `1 − 2Σmin(a,b)/(Σa+Σb)`. PCoA = eigen-decomposition of the double-centred squared distance matrix; Bray–Curtis is not Euclidean so some eigenvalues are negative (we report % explained from positive ones).
- PERMANOVA: pseudo-F = (SS_between/(g−1)) / (SS_within/(n−g)) computed on distances; p-value from permuting group labels (999 here). R² = share of distance variance explained. It tests differences in *centroids and spread*; unequal dispersion can give a false positive.
- Our numbers: covariates explain only ~1–3 % of variance (R² of 0.004–0.03), strongest for age group and `In_NC`.
- **Questions:** Why permute? (No distributional assumption on distances.) What does R² = 0.027 mean? (Small effect even when p < 0.001.) What is the arch/horseshoe in the Aitchison PCoA? (A gradient-induced distortion when the first axis dominates; not a second biological signal.)

### Taylor's law and overdispersion
- Variance vs mean per taxon on log–log axes: `σ² = a μ^b`. A common-probability multinomial gives `σ² = μ(1 − μ/N)`, slope 1.
- Our data: slope b = 1.84 (bootstrap CI 1.82–1.89); every taxon has more than 2× the multinomial variance (median ≈ 22×). So M0 will fail, as predicted, and M1 or M3 are needed.
- **Questions:** What does b≈2 mean? (Constant coefficient of variation: heavy-tailed, lognormal-like abundance across samples.) Does overdispersion prove the DM is right? (No: it only rules out M0; DM implies a specific mean-variance form, to be tested in Phase 3.)

### Multiple testing (first use)
- 12 exploratory alpha tests and 12 PERMANOVAs corrected with Benjamini–Hochberg. BH controls the expected share of false discoveries among the rejections, not the chance of any false positive.

## Phase 2 – models

### Likelihood and maximum likelihood
- **Plain language:** the likelihood is the probability the model gives to the data we actually observed, seen as a function of the parameters. Maximum likelihood (ML) picks the parameters that make the observed data most probable.
- **Formula:** `ℓ(θ) = Σ_i log P(x_i | θ)`; `θ̂ = argmax ℓ(θ)`. We maximise numerically (L-BFGS) with analytic gradients, checked against finite differences in the tests.
- **Assumptions:** samples independent given θ; the model family contains (approximately) the truth.
- **What breaks:** if the model is wrong, ML still gives "the best wrong model" (closest in KL divergence); standard errors from the Hessian may then be off, which is one reason we bootstrap.
- **Questions:** Why the log? (Products → sums; numerically stable; same maximiser.) Why include the multinomial coefficient? (So log-likelihoods of different models are true log-probabilities of the same data and can be compared.)

### Multinomial (M0) vs Dirichlet-multinomial (M1)
- **M0:** every sample draws its N reads with the same probabilities p. `Var x_j = N p_j(1−p_j)`.
- **M1:** each sample first draws its own probabilities `p_i ~ Dirichlet(α)`, then its reads. Mean is unchanged, `p_j = α_j/A` with `A = Σα_j`, but `Var x_j = N p_j(1−p_j)·(N+A)/(1+A) = N p_j(1−p_j)[1+(N−1)ρ]` with `ρ = 1/(1+A)`.
- **Likelihood:** `log P(x|α) = log N! − Σ log x_j! + lnΓ(A) − lnΓ(N+A) + Σ_j [lnΓ(x_j+α_j) − lnΓ(α_j)]`.
- **Checks in the tests:** probabilities sum to 1 over all count vectors; DM → multinomial as A → ∞; simulated mean and variance match the formulas.
- **Limitation:** one ρ for all taxa (one common overdispersion factor) and only negative correlations between taxa. The data show overdispersion that varies a lot between taxa, so M1 is expected to be imperfect.
- **Questions:** What is ρ? (Intra-class correlation: chance that two reads of a sample share their "source"; ρ = 0 is the multinomial.) Why is overdispersion expected biologically? (Real samples differ in composition: people differ.)

### Dirichlet-multinomial regression (M2)
- `log α_ij = b_j + z_i·β_j` for the top-10 taxa. exp(β_j) multiplies α_j, so it multiplies the expected ratio p_j/p_k to unaffected categories k.
- Complete cases only; M1 is refitted on exactly the same samples for comparison.
- **Likelihood-ratio test:** M1 is M2 with all β = 0 (nested), so `LR = 2(ℓ̂_M2 − ℓ̂_M1) ~ χ²_df` under H₀ (Wilks' theorem), df = number of extra parameters. Needs large n and parameters not on the boundary.
- **Questions:** Why can LR be significant while CV says no gain? (LR tests whether *some* effect exists in-sample; CV asks whether the extra parameters improve prediction of new samples. Small real effects with many parameters can fail to pay for their estimation noise.)

### Logistic-normal multinomial (M3), ALR, Laplace approximation, importance sampling
- **Model:** `z_i ~ Normal(μ, Σ)` in ALR coordinates `z_j = log(p_j/p_ref)`; `p_i = softmax([z_i, 0])`; reads ~ Multinomial(N, p_i). Σ allows positive *and* negative correlations between taxa and a separate variance per taxon. No pseudocount is needed: zeros are handled by the multinomial step.
- **Problem:** `p(x) = ∫ Mult(x|softmax z) N(z|μ,Σ) dz` has no closed form.
- **Laplace approximation:** expand log of the integrand to second order around its maximum ẑ (like the saddle-point method): `log p(x) ≈ f(ẑ) − ½ log det Σ − ½ log det H`, H = curvature at ẑ. Exact for a Gaussian integrand; inaccurate when the integrand is skewed (zero counts make it skewed).
- **Importance sampling:** `p(x) = E_q[p(x,z)/q(z)]`, estimated by averaging weights over draws from a proposal q. Unbiased for p(x), reliable only when the effective sample size `ESS = (Σw)²/Σw²` is large.
- **Our finding:** with 114 categories Laplace and IS disagree by several nats and ESS is 2–6 of 20 000 → not trustworthy. At D = 22 they agree within a few tenths of a nat with ESS > 1000. Hence the comparison is done at D = 22 (decision taken together).
- **Fitting:** approximate EM: E-step = Laplace Gaussian for each sample's z; M-step = update μ and Σ; Σ shrunk towards its diagonal by λ (chosen by CV).
- **Questions:** Why shrink Σ? (231 covariance parameters from 273 samples: the raw estimate is noisy; shrinkage trades a little bias for much less variance.) Why does the ALR reference not matter for the likelihood? (Changing reference is a linear map of z; Gaussians stay Gaussian. Only the diagonal shrinkage depends on it.)

### AIC, BIC and when they are not valid
- `AIC = 2k − 2ℓ̂`, `BIC = k ln n − 2ℓ̂`. Lower is better. AIC estimates out-of-sample prediction error; BIC approximates the model evidence and penalises parameters more when n is large.
- Need an exact maximised likelihood and a well-defined parameter count k: not valid for M3 (approximate likelihood, shrinkage makes the effective k smaller than the nominal one), so M3 is judged by CV only.

### Cross-validation
- Split samples into 5 folds; fit on 4, score the log-likelihood of the held-out fold; repeat; every sample is scored exactly once by a model that never saw it.
- Report the mean per sample and the standard error; compare models with **paired** differences on the same samples (removes sample-to-sample variation common to both models).
- **What breaks:** using the same folds to tune λ and to compare models is slightly optimistic; folds must be split by *subject* (here one sample per subject, so fine).
- **Questions:** Why not just compare in-sample likelihood? (A bigger model always fits the training data at least as well.) Why 5 folds? (Trade-off: more folds = less bias, more variance and compute.)

### Bootstrap for parameter uncertainty
- Resample the 273 subjects with replacement, refit, repeat B times; the 2.5 % and 97.5 % percentiles of the refitted values give a 95 % CI.
- **Assumptions:** subjects are independent and the sample is representative; percentile CIs can be off for strongly skewed or boundary parameters.
- **Questions:** Why not the Hessian? (Possible: inverse observed Fisher information gives asymptotic SEs; the bootstrap needs no formula and works for derived quantities like ρ = 1/(1+A).) How many replicates? (Enough that the CI endpoints are stable: 500 for M1/M2, 100 for the slower M3.)

## Phase 3 – validation

### Parametric bootstrap (posterior-predictive-style check)
- **Plain language:** pretend the fitted model is the truth, simulate many datasets with the real design (273 samples × 1000 reads), and compute the same statistic on each. This shows what values the model *can* produce. If the real value lies outside the central 95 % of them, the model fails to reproduce that feature.
- **Formula:** band = [2.5 %, 97.5 %] quantiles of `T(Y_r)`, `Y_r ~ model(θ̂)`, r = 1…R (R = 200).
- **Assumptions/caveat:** uses the plug-in θ̂, so the band ignores the uncertainty of θ̂ (slightly too narrow). With 273 samples and few parameters (M0, M1) the effect is small.
- **Questions:** Difference from the non-parametric bootstrap? (Parametric: draws new data *from the model*, tests the model. Non-parametric: resamples *the real data*, gives the CI of the real statistic, no model needed.) Is "inside the band" proof the model is right? (No: it only fails to reject it for that statistic.)

### Non-parametric bootstrap over samples
- Resample subjects with replacement and recompute the statistic → 95 % CI of the real value.
- **Pitfall handled:** for pairwise statistics (Bray–Curtis) a resampled subject paired with its own copy gives a distance of 0; those pairs are excluded.

### Kolmogorov–Smirnov and Wasserstein distances
- **KS:** `D = max_x |F_real(x) − F_sim(x)|`, the largest vertical gap between the two empirical CDFs. The p-value assumes independent observations (OK for per-sample metrics, *not* for pairwise distances).
- **Wasserstein-1 ("earth mover's"):** `W = ∫ |F_real(x) − F_sim(x)| dx`, the average shift needed to move one distribution onto the other, in the units of the metric (interpretable: "Shannon differs by 0.2 nats on average").
- **Calibration by simulation:** W is compared with its distribution when the model is true (W between other simulated datasets and the reference pool). Monte Carlo p = (1 + #{W_null ≥ W_real})/(1 + R). This works even when observations are dependent.
- **Questions:** Why not only KS? (Its p-value needs independence; KS also focuses on one point of the CDF; W uses the whole shape and has units.) Why does a tiny p with large n not mean a bad model? (With enough data every model is rejected; look at the size of the discrepancy.)

### PERMANOVA real vs simulated, PERMDISP
- Pool real and simulated samples, label them, and test whether the label explains distance variance. A good simulator is *hard to tell apart*: R² ≈ 0.
- PERMANOVA is sensitive to differences in centroid *and* spread; PERMDISP tests equal spread (distance of samples to their group centroid), to tell the two apart.
- **Questions:** Why repeat over 20 simulated datasets? (One simulated dataset is random; the spread of R² shows the Monte Carlo variability.) Why focus on R² not p? (546 samples give power to detect trivial differences.)

### Zero inflation (M4) and the method of moments
- **Plain language:** a zero count can mean "present but not caught in 1000 reads" (sampling zero) or "not there at all" (structural zero). M3 only makes sampling zeros and so produces too few zeros (32 % vs 41 %). M4 adds a probability πⱼ that taxon j is truly absent from a person.
- **Model:** presence bᵢⱼ ~ Bernoulli(1 − πⱼ); z ~ N(μ, Σ); pᵢ = softmax of z over the present taxa only (reference `Other` always present); reads ~ Multinomial(N, pᵢ).
- **Method of moments:** choose parameters so that model expectations equal sample moments. Here 2 per category: prevalence `E[1(x_j>0)] = (1−π_j)·P(x_j>0 | present)` and mean count `E[x_j] = N·E[p_j]`; solved by fixed-point iteration with common random numbers (same random draws every iteration, so the update is smooth). Σ is kept from M3.
- **Assumptions/limits:** absences independent between taxa and of everything else; Σ not re-estimated; not a likelihood fit → no CV log-likelihood or AIC for M4. Prevalence and mean richness are matched *by construction*, so they cannot be used to validate M4 (circular); Shannon, Simpson, dominance, Bray–Curtis, Taylor slope and PERMANOVA remain independent checks.
- **Questions:** Why not maximum likelihood? (Each sample's likelihood sums over all presence patterns of its zero taxa, about 2⁹ on average, each needing a Laplace integral: possible but heavy, and harder to explain; moments give a transparent calibration.) Is matching prevalence "cheating"? (Only for the statistics it was tuned on; that is why those are excluded from M4's validation.) What does πⱼ = 0.5 mean? (Half the people do not carry that taxon at all.)

### Choosing the final model
- Rule fixed **before** looking at the results: the best CV model is the final simulator unless another model reproduces more of the summary statistics. Fixing the rule first avoids picking whichever check happens to favour a preferred model.

## Phase 4 – power study (to do)
p-values, multiple testing (Benjamini–Hochberg FDR), type I error vs power, Monte Carlo error, cosinor regression.
