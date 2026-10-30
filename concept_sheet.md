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

## Phase 2 – models (to do)
Multinomial, Dirichlet-multinomial, logistic-normal, likelihood, AIC/BIC, cross-validation.

## Phase 3 – validation (to do)
Parametric and non-parametric bootstrap, KS/Wasserstein, PERMANOVA.

## Phase 4 – power study (to do)
p-values, multiple testing (Benjamini–Hochberg FDR), type I error vs power, Monte Carlo error, cosinor regression.
