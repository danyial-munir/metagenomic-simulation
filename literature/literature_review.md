# Literature review: simulating compositional microbiome data to study seasonal dynamics

*Project: Statistical Data Analysis for Applied Physics. Topic (proposal 5): simulation of metagenomic data (relative abundances) to study seasonal dynamics of bacterial populations, applied to a human belly-button (navel) skin microbiome dataset.*

All references were checked against an online record (publisher, PubMed/PMC, journal or library catalogue) in October 2026; the verification table at the end lists what was confirmed for each one. Statements about a paper's content are limited to what could be confirmed from its abstract or record.

---

## 1. Scope

The project asks three questions: which generative model best reproduces real skin-microbiome count data (RQ1); how closely data simulated from it match the real data (RQ2); and how large a seasonal change must be, and how many samples are needed, before standard tools detect it (RQ3). Because the dataset has no time dimension, the seasonal signal is imposed by simulation. The review therefore covers six areas: (i) how microbiome sequencing data arise and the dataset used; (ii) compositional data and Aitchison geometry; (iii) count models for compositions (multinomial, Dirichlet-multinomial, logistic-normal, zero inflation); (iv) simulators of microbiome data; (v) diversity metrics, ordination and PERMANOVA; (vi) seasonality in microbial communities, rhythm detection, power analysis and multiple testing.

## 2. Microbiome sequencing data and the dataset

### 2.1 From swab to count table

In amplicon (16S rRNA gene) studies, DNA is extracted from a sample, a variable region of the 16S gene is amplified and sequenced, and reads are clustered or classified into taxa (operational taxonomic units, phylotypes or genera). The result is a table of **read counts** per taxon and sample. The total number of reads per sample, the *library size*, is set by the sequencing run and not by the amount of bacteria on the skin. This has two consequences that shape the whole project. First, only **relative** information is available (Section 3). Second, samples differ in depth, so they are often **rarefied**: randomly subsampled to a common number of reads.

Rarefying is debated. McMurdie and Holmes (2014) argued that it discards information and is statistically inefficient, and recommended model-based normalisation instead. Weiss et al. (2017) compared normalisation and differential-abundance strategies and found that the best choice depends on data characteristics: rarefying clustered samples by biological origin well in ordination, and in their simulations rarefying did not by itself increase the false discovery rate of many differential-abundance methods. In this project the data were already rarefied to 1000 reads per sample before we received them, so we model the data as given (multinomial sampling of 1000 reads) and state rarefaction as a limitation.

### 2.2 The belly-button study

The data come from the *Belly Button Biodiversity* citizen-science project. Its first publication, Hulcr et al. (2012, *PLoS ONE*), sequenced 16S rDNA from navel swabs collected from volunteers in two populations within a nation-wide citizen-science project. After rarefying to 400 reads per sample, belly buttons contained on average 67 bacterial phylotypes. Only 6 phylotypes occurred in more than 80 % of people. These frequent phylotypes were fewer than 0.3 % of all phylotypes but accounted for about one third of the reads, a pattern the authors compared to "oligarchies" in tropical forests. The associated data package on Dryad (Hulcr et al. 2013) contains the raw sequences, aligned representative sequences and a sample mapping file.

**Provenance caveat (unresolved).** The file analysed in this project has 273 subjects, a genus-level taxonomy table rarefied to 1000 reads per sample, and a metadata table with washing frequency and navel type. Popular accounts of Hulcr et al. (2012) describe 60 volunteers, and the paper rarefied to 400 reads. Our table is therefore **not** the published analysis table. It is consistent with a later or extended processing of the same project's samples, but its exact origin could not be verified from public sources. It must be stated as such in the report.

The broader skin literature gives context for what the data should look like. Costello et al. (2009) sampled up to 27 body sites in seven to nine adults on four occasions. Body habitat was the main determinant of community composition, variation *between people* within a habitat was high, and a person's community was fairly stable over time. Grice and Segre (2011) review the skin microbiome as a set of distinct habitats shaped by host and environment. Our exploratory analysis is consistent with this picture: strong dominance by a few taxa (*Corynebacterium*, *Staphylococcaceae*), very sparse tables, and large person-to-person variation.

## 3. Compositional data and Aitchison geometry

### 3.1 The simplex

Because only relative abundances are observed, each sample is a **composition**: a vector of positive parts that carry information only through their ratios. Aitchison (1982, read to the Royal Statistical Society; book 1986) showed that standard multivariate statistics applied to such data can mislead. The data live on the **simplex**, not in ordinary Euclidean space. The constant-sum constraint induces spurious negative correlations, and results can change when a part is added or removed. His solution is to work with **log-ratios**, which map the simplex to real space where ordinary methods apply.

### 3.2 Log-ratio transforms

- **Additive log-ratio (ALR):** log(x_j / x_D) with respect to a reference part. Simple, but not symmetric in the parts. ALR is the natural coordinate system of the logistic-normal distribution (Aitchison 1986) and of the logistic-normal multinomial regression of Xia et al. (2013).
- **Centred log-ratio (CLR):** log(x_j / g(x)), with g(x) the geometric mean. Symmetric. The Euclidean distance between CLR vectors is the **Aitchison distance**. The coordinates sum to zero, so a covariance matrix of CLR values is singular.
- **Isometric log-ratio (ILR):** Egozcue et al. (2003) introduced orthonormal log-ratio coordinates that preserve all metric properties of the simplex and are related to both ALR and CLR. They avoid the singularity of CLR but are harder to interpret, which is why this project does not use them.

### 3.3 Why this matters for microbiome data

Gloor et al. (2017, *Frontiers in Microbiology*) argued that sequencing-based microbiome datasets are compositional at every stage, and that analysing them as if they were absolute abundances leads to errors in correlation, differential abundance and normalisation. They recommend log-ratio methods.

**Zeros** are the practical obstacle: the log of zero is undefined, and microbiome tables are mostly zeros (93 % in our data). The usual remedies are:

- adding a **pseudocount** (simple, but the result depends on its value);
- replacing zeros by small values;
- using a **count model** in which zeros arise naturally from sampling.

This project uses the third option for its generative models. It uses a pseudocount only for descriptive distances and for the per-taxon tests, and checks sensitivity to the pseudocount.

## 4. Count models for compositional data

### 4.1 Multinomial and overdispersion

The simplest model treats each sample's reads as a multinomial draw from a common probability vector. It predicts a binomial variance for each taxon, Var = N·p(1 − p). Ecological counts are almost always more variable than this. Taylor (1961) showed across many species that the variance of population counts grows as a power of the mean, now called **Taylor's law**. In our data the slope is about 1.8 and every taxon is strongly overdispersed relative to the multinomial.

### 4.2 Dirichlet-multinomial (DM)

The DM lets each sample have its own probability vector drawn from a Dirichlet distribution, which inflates the variance by a common factor 1 + (N − 1)ρ.

- **Clustering:** Holmes, Harris and Quince (2012) built a generative framework on mixtures of Dirichlet-multinomials. Each community's taxon probabilities come from one of several Dirichlet components ("metacommunities", linked to enterotypes). They selected the number of components by model evidence, using a Laplace approximation.
- **Testing and power:** La Rosa et al. (2012) used the DM for parameter estimation, multivariate hypothesis tests between groups, and power and sample-size calculations, implemented in the R package *HMP*. They argued that this keeps more of the multivariate structure than testing one taxon at a time.
- **Covariates:** Chen and Li (2013) proposed a **DM regression** in which covariates enter through the Dirichlet parameters, with a sparse group penalty to select covariates. In their simulations, accounting for overdispersion identified microbiome-associated covariates better than ignoring it.

Our M1 (DM) and M2 (DM regression with covariates on the most abundant taxa) follow this line, fitted by plain maximum likelihood.

### 4.3 Logistic-normal multinomial (LNM)

The DM has only one overdispersion parameter and cannot represent positive correlations between taxa. The **logistic-normal** distribution (Aitchison 1986) places a multivariate Gaussian on ALR coordinates. That gives each taxon its own variance and allows arbitrary correlations. Combined with multinomial sampling of reads, it gives the logistic-normal multinomial.

Xia, Chen, Fung and Li (2013, *Biometrics*) developed an LNM regression for microbiome compositions using the ALR transform, with penalised estimation to select covariates associated with taxon composition. Later work extended the framework to clustering and to higher dimensions. A practical difficulty is that the LNM likelihood requires integrating over the latent Gaussian, which has no closed form. In this project we use a Laplace approximation and check it against importance sampling; the check limited the model to the 20 most abundant taxa (Phase 2).

### 4.4 Zero inflation

Real tables often contain more zeros than sampling alone explains: a taxon may be truly absent from a person. Simulators therefore often include **zero inflation**. SparseDOSSA (Ma et al. 2021), for example, models each feature's marginal abundance as a zero-inflated log-normal, with added components for sequencing and for microbe–microbe and microbe–environment structure. Our validation (Phase 3) found exactly this deficit in the LNM: 32 % zeros versus 41 % real. That motivated M4, an LNM with per-taxon structural zeros calibrated by the method of moments.

## 5. Simulators of microbiome data

Simulation is the standard way to evaluate statistical methods for microbiome data, because it provides data with a **known ground truth**. Three examples from the literature:

- **metaSPARSim** (Patuzzi et al. 2019) generates sparse 16S count tables. It models the sequencing step with a multivariate hypergeometric distribution and can estimate its parameters from real data.
- **SparseDOSSA 2** (Ma et al. 2021) fits real community profiles and simulates new ones with known structure, including spiked-in associations, for benchmarking.
- **DM-based simulation** for power and sample size is available through the *HMP* package (La Rosa et al. 2012).
- **Distance-based simulation:** Kelly et al. (2015) simulate distance matrices to estimate PERMANOVA power (*micropower*). The pseudo-F statistic does not follow Fisher's F distribution, so parametric ANOVA power formulas do not apply.

This project follows the same logic with deliberately simple, explainable models (M0–M4). The emphasis is on validating the simulator against the real data **before** using it for a power study. That validation step is the main difference from simply adopting an existing simulator.

## 6. Diversity, ordination and PERMANOVA

- **Alpha diversity** summarises one sample. **Richness** counts the taxa present. **Shannon** entropy, H = −Σ p ln p, comes from Shannon's (1948) theory of communication and weights both richness and evenness. **Simpson's** (1949) index, Σ p², is the probability that two individuals drawn at random belong to the same type; we report 1 − Σ p².
- **Beta diversity** compares samples. The **Bray–Curtis** dissimilarity (Bray and Curtis 1957) is based on shared abundance and ignores joint absences. The **Aitchison distance** (Section 3) is its compositional alternative.
- **Principal coordinates analysis** (PCoA; Gower 1966) recovers coordinates of points from their pairwise distances and plots samples in two dimensions. Non-Euclidean distances such as Bray–Curtis give some negative eigenvalues.
- **PERMANOVA** (Anderson 2001) partitions the variation of a distance matrix among groups and tests it by permuting group labels. It needs few distributional assumptions and handles complex designs. It reacts to differences in group centroids *and* in dispersion.
- **PERMDISP** (Anderson 2006), a distance-based extension of Levene's test, checks the dispersion part. We use both to test whether simulated samples can be told apart from real ones (Phase 3) and whether seasonal groups differ (Phase 4).

## 7. Seasonality in microbial communities

Seasonal dynamics are well documented in environmental microbiomes:

- **Marine:** Gilbert et al. (2012) analysed a six-year monthly 16S time series from a coastal site off Plymouth (UK). Surface bacterial communities showed repeatable seasonal patterns, with winter peaks in diversity, and day length alone explained more than 65 % of the variance in diversity.
- **Wastewater** (the setting of the original project proposal): LaMartina, Mohaimani and Newton (2021) followed sewage influent monthly for five years. The bacterial community assembled into two distinct states each year, driven by the extremes of wastewater temperature, while the human-derived fraction varied more randomly.

For **human-associated** communities, the evidence for seasonality is weaker and comes from small studies:

- **Gut:** Davenport et al. (2014) followed the gut microbiome of 60 Hutterites, whose diet varies with the season. Individuals stayed fairly stable, but population-level composition and diversity differed between winter and summer.
- **Skin, long-term stability:** Oh et al. (2016) found skin bacterial, fungal and viral communities largely stable over time. Stability depended on body site and individual, and was driven by persistence of strains.
- **Skin, seasonal shifts:** Seo et al. (2023) sampled the facial skin of 10 adults weekly for a year. *Cutibacterium* was more abundant in winter and *Corynebacterium*, *Staphylococcus* and *Streptococcus* in summer, together with changes in skin-barrier measures. The authors note that mask-wearing during the COVID-19 pandemic may confound this.

In short, seasonal shifts in human skin communities are plausible but probably modest compared with differences between people. This is the regime the power study explores.

**Rhythm detection.** The standard tool for a periodic signal of known period is **cosinor** regression: fitting y = b₀ + b₁ cos(2πt/T) + b₂ sin(2πt/T) by least squares. Cornelissen (2014) reviews cosinor-based rhythmometry. Framed as a regression, it handles unequally spaced data and gives confidence intervals for the amplitude and phase. Because cos and sin together represent any phase, the peak month need not be known in advance.

## 8. Power analysis, multiple testing and resampling

**Power** is the probability that a test detects an effect of a given size. For microbiome data it is usually estimated by simulation:

- La Rosa et al. (2012) derive power and sample size under the DM.
- Kelly et al. (2015) simulate distance matrices for PERMANOVA.

In both cases power depends strongly on the assumed between-sample variability. This is why our simulator is first validated against the real data.

**Multiple testing.** When many taxa are tested, Benjamini and Hochberg (1995) control the **false discovery rate**, the expected proportion of false rejections among all rejections. Their step-up procedure controls the FDR for independent test statistics and is more powerful than family-wise control.

**Resampling.** The **bootstrap** (Efron 1979) estimates uncertainty by resampling the data (non-parametric) or by simulating from a fitted model (parametric). We use both, for parameter uncertainty and for comparing real and simulated data. Rejection rates estimated by simulation carry a Monte Carlo standard error √(p(1 − p)/R), which sets the number of replicates R.

## 9. Gap addressed by this project

The literature provides count models for compositional data (DM, LNM, zero-inflated variants), simulators, and power frameworks for group comparisons. Three gaps remain:

1. Few studies quantify how detectable a **seasonal** signal would be in a **cross-sectional** skin-microbiome design, where each person is sampled once and seasonal variation competes with large differences between people.
2. Existing simulators are rarely **validated** against the specific dataset before being used for design questions.
3. The compositional side-effects of a seasonal shift are seldom reported. Changing a few taxa moves the CLR of all of them, which affects false discovery rates.

This project addresses these points with a sequence of simple, explainable models fitted to the belly-button data, an explicit validation step, and a pre-registered (frozen) simulation plan. The answer is framed honestly as *"what could be detected if"*, not as a measurement of real seasonality.

---

## References

Aitchison, J. (1982). The statistical analysis of compositional data (with discussion). *Journal of the Royal Statistical Society: Series B*, 44(2), 139–177.

Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*. Monographs on Statistics and Applied Probability 25. London: Chapman and Hall. (Reprinted by Blackburn Press.)

Anderson, M. J. (2001). A new method for non-parametric multivariate analysis of variance. *Austral Ecology*, 26, 32–46.

Anderson, M. J. (2006). Distance-based tests for homogeneity of multivariate dispersions. *Biometrics*, 62(1), 245–253.

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B*, 57(1), 289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

Bray, J. R., & Curtis, J. T. (1957). An ordination of the upland forest communities of southern Wisconsin. *Ecological Monographs*, 27(4), 325–349. https://doi.org/10.2307/1942268

Chen, J., & Li, H. (2013). Variable selection for sparse Dirichlet-multinomial regression with an application to microbiome data analysis. *Annals of Applied Statistics*, 7(1), 418–442. https://doi.org/10.1214/12-AOAS592

Cornelissen, G. (2014). Cosinor-based rhythmometry. *Theoretical Biology and Medical Modelling*, 11, 16. https://doi.org/10.1186/1742-4682-11-16

Costello, E. K., Lauber, C. L., Hamady, M., Fierer, N., Gordon, J. I., & Knight, R. (2009). Bacterial community variation in human body habitats across space and time. *Science*, 326(5960), 1694–1697. https://doi.org/10.1126/science.1177486

Davenport, E. R., et al. (2014). Seasonal variation in human gut microbiome composition. *PLoS ONE*, 9(3), e90731. https://doi.org/10.1371/journal.pone.0090731

Efron, B. (1979). Bootstrap methods: another look at the jackknife. *Annals of Statistics*, 7(1), 1–26. https://doi.org/10.1214/aos/1176344552

Egozcue, J. J., Pawlowsky-Glahn, V., Mateu-Figueras, G., & Barceló-Vidal, C. (2003). Isometric logratio transformations for compositional data analysis. *Mathematical Geology*, 35(3), 279–300.

Gilbert, J. A., et al. (2012). Defining seasonal marine microbial community dynamics. *The ISME Journal*, 6(2), 298–308. https://doi.org/10.1038/ismej.2011.107

Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V., & Egozcue, J. J. (2017). Microbiome datasets are compositional: and this is not optional. *Frontiers in Microbiology*, 8, 2224. https://doi.org/10.3389/fmicb.2017.02224

Gower, J. C. (1966). Some distance properties of latent root and vector methods used in multivariate analysis. *Biometrika*, 53(3–4), 325–338. https://doi.org/10.1093/biomet/53.3-4.325

Grice, E. A., & Segre, J. A. (2011). The skin microbiome. *Nature Reviews Microbiology*, 9(4), 244–253. https://doi.org/10.1038/nrmicro2537

Holmes, I., Harris, K., & Quince, C. (2012). Dirichlet multinomial mixtures: generative models for microbial metagenomics. *PLoS ONE*, 7(2), e30126. https://doi.org/10.1371/journal.pone.0030126

Hulcr, J., Latimer, A. M., Henley, J. B., Rountree, N. R., Fierer, N., Lucky, A., Lowman, M. D., & Dunn, R. R. (2012). A jungle in there: bacteria in belly buttons are highly diverse, but predictable. *PLoS ONE*, 7(11), e47712. https://doi.org/10.1371/journal.pone.0047712

Hulcr, J., Latimer, A. M., Henley, J. B., Rountree, N. R., Fierer, N., Lucky, A., Lowman, M. D., & Dunn, R. R. (2013). Data from: A jungle in there: bacteria in belly buttons are highly diverse, but predictable [Dataset]. Dryad. https://doi.org/10.5061/dryad.q8f74

Kelly, B. J., et al. (2015). Power and sample-size estimation for microbiome studies using pairwise distances and PERMANOVA. *Bioinformatics*, 31(15), 2461–2468. https://doi.org/10.1093/bioinformatics/btv183

La Rosa, P. S., et al. (2012). Hypothesis testing and power calculations for taxonomic-based human microbiome data. *PLoS ONE*, 7(12), e52078. https://doi.org/10.1371/journal.pone.0052078

LaMartina, E. L., Mohaimani, A. A., & Newton, R. J. (2021). Urban wastewater bacterial communities assemble into seasonal steady states. *Microbiome*, 9. https://doi.org/10.1186/s40168-021-01038-5

Ma, S., et al. (2021). A statistical model for describing and simulating microbial community profiles. *PLoS Computational Biology*, 17(9), e1008913. https://doi.org/10.1371/journal.pcbi.1008913

McMurdie, P. J., & Holmes, S. (2014). Waste not, want not: why rarefying microbiome data is inadmissible. *PLoS Computational Biology*, 10(4), e1003531. https://doi.org/10.1371/journal.pcbi.1003531

Oh, J., Byrd, A. L., Park, M., Kong, H. H., & Segre, J. A. (2016). Temporal stability of the human skin microbiome. *Cell*, 165(4), 854–866. https://doi.org/10.1016/j.cell.2016.04.008

Patuzzi, I., Baruzzo, G., Losasso, C., Ricci, A., & Di Camillo, B. (2019). metaSPARSim: a 16S rRNA gene sequencing count data simulator. *BMC Bioinformatics*, 20(Suppl 9), 416. https://doi.org/10.1186/s12859-019-2882-6

Seo, J. Y., You, S. W., Gu, K.-N., Kim, H., Shin, J.-G., Leem, S., Hwang, B. K., Kim, Y., & Kang, N. G. (2023). Longitudinal study of the interplay between the skin barrier and facial microbiome over 1 year. *Frontiers in Microbiology*, 14, 1298632. https://doi.org/10.3389/fmicb.2023.1298632

Shannon, C. E. (1948). A mathematical theory of communication. *Bell System Technical Journal*, 27, 379–423 and 623–656.

Simpson, E. H. (1949). Measurement of diversity. *Nature*, 163, 688.

Taylor, L. R. (1961). Aggregation, variance and the mean. *Nature*, 189(4766), 732–735. https://doi.org/10.1038/189732a0

Weiss, S., et al. (2017). Normalization and microbial differential abundance strategies depend upon data characteristics. *Microbiome*, 5, 27. https://doi.org/10.1186/s40168-017-0237-y

Xia, F., Chen, J., Fung, W. K., & Li, H. (2013). A logistic normal multinomial regression model for microbiome compositional data analysis. *Biometrics*, 69(4), 1053–1063.

---

## Reference verification (October 2026)

Every reference was found in an online record. "Partial" means some fields could not be confirmed and are flagged.

| Reference | Confirmed via | Status / flags |
|---|---|---|
| Aitchison 1982 | biographical and reappraisal records citing JRSS B 44, 139–177 | verified; issue number (2) from standard citation, not seen in a record |
| Aitchison 1986 | library catalogues (Chapman and Hall, series vol. 25) | verified |
| Anderson 2001 | full-text PDF copy, citation records | verified |
| Anderson 2006 | IDEAS/RePEc record, scikit-bio reference | verified |
| Benjamini & Hochberg 1995 | DOI record, citation databases | verified |
| Bray & Curtis 1957 | JSTOR table of contents, DOI | verified |
| Chen & Li 2013 | Project Euclid, DOI, arXiv | verified |
| Cornelissen 2014 | PMC, publisher page | verified |
| Costello et al. 2009 | institutional records, DOI | verified |
| Davenport et al. 2014 | PMC, PLoS | **partial**: co-author list not checked (cited as "et al.") |
| Efron 1979 | Project Euclid, DOI | verified |
| Egozcue et al. 2003 | repository record | verified; no DOI recorded |
| Gilbert et al. 2012 | PMC, ISME issue page, DOI | **partial**: 12 co-authors not listed individually (cited as "et al.") |
| Gloor et al. 2017 | PubMed, PMC, Frontiers | verified |
| Gower 1966 | Rothamsted repository, DOI | verified |
| Grice & Segre 2011 | author manuscript record, DOI | verified |
| Holmes, Harris & Quince 2012 | PMC, DOAJ, DOI | verified |
| Hulcr et al. 2012 | PLoS/PMC record, Dryad landing page | verified (full text not read; content from the abstract) |
| Hulcr et al. 2013 (Dryad data) | Dryad landing page | verified; files could not be downloaded (HTTP 403) |
| Kelly et al. 2015 | PMC, DOI | **partial**: co-author list not checked |
| La Rosa et al. 2012 | PMC, DOI | **partial**: co-author list not checked |
| LaMartina et al. 2021 | PMC, DOAJ, DOI | **partial**: volume 9 confirmed, article number not confirmed |
| Ma et al. 2021 | PMC, DOAJ, DOI | **partial**: co-author list not checked |
| McMurdie & Holmes 2014 | PMC, DOI | verified |
| Oh et al. 2016 | Duke scholars record, DOI | verified |
| Patuzzi et al. 2019 | PMC, publisher | verified |
| Seo et al. 2023 | PMC, Frontiers | verified |
| Shannon 1948 | multiple citation records | verified |
| Simpson 1949 | reprint page, citation records | verified |
| Taylor 1961 | Rothamsted repository, DOI | verified |
| Weiss et al. 2017 | PMC, publisher | **partial**: co-author list not checked |
| Xia, Chen, Fung & Li 2013 | HKU repository record | verified; DOI not recorded; first names of authors abbreviated |
