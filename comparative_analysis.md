# Comparative Analysis: NLS vs. Log-Linear SAR Fitting

A companion to [`analysis_updates.ipynb`](analysis_updates.ipynb), comparing the nonlinear least squares (NLS) results from the [`sars`](https://pypi.org/project/sars/) Python library against the original log-linear OLS calculations published in Jessee, Stout & McMeen (2022).

---

## Fitting Method Differences

The original calculations used log₁₀-transformed OLS regression:

> log₁₀(S) = z · log₁₀(A) + log₁₀(c)

This minimizes squared residuals in log-space and implicitly assumes multiplicative error — i.e., error proportional to the predicted value. The c and z parameters are recovered as c = 10^(intercept) and z = slope.

The `sars` library fits the power model directly in arithmetic space via nonlinear least squares (`scipy.optimize.least_squares`):

> S = cA^z

This minimizes squared residuals in the original species-count scale and assumes additive error. The NLS approach is generally preferred in modern SAR literature because it does not impose the variance structure that log-transformation assumes (Tjørve & Tjørve 2021).

## Parameter Estimates

| Group | Method | c | z | R² |
|---|---|---|---|---|
| Total Herpetofauna | Original (log-linear) | 29.54 | 0.1137 | 0.86 (adj.) |
| Total Herpetofauna | `sars` NLS | 22.49 | 0.1478 | 0.9162 |
| Amphibians | Original (log-linear) | 14.09 | 0.1449 | 0.98 (adj.) |
| Amphibians | `sars` NLS | 12.28 | 0.1616 | 0.9849 |
| Reptiles | Original (log-linear) | 15.45 | 0.0747 | 0.33 (adj.) |
| Reptiles | `sars` NLS | 10.96 | 0.1205 | 0.6619 |

> **R² values are not directly comparable between rows.** The original R² is adjusted R² computed in log-space (measuring fit of log S vs. log A), while the `sars` R² is unadjusted and computed in arithmetic space (measuring fit of S vs. A). Both metrics are retained here to match their respective source methods, but side-by-side numerical comparison is not meaningful.

The c parameter (baseline richness) is systematically lower under NLS, while z (the scaling exponent) is higher. This shift is expected: log-linear regression overweights small-area sites relative to their absolute residuals, pulling c upward and z downward. The discrepancy is most dramatic for reptiles, where z increases by 61% (0.0747 → 0.1205).

The original R² values are adjusted R² in log-space, while the `sars` values are unadjusted R² in arithmetic space. These are not directly comparable, but both tell the same story: herpetofauna and amphibians fit the power model well; reptiles do not.

## z-Value Interpretation

All three NLS z-values fall within the canonical mainland/nested range of 0.10–0.20 (Rosenzweig 1995, Drakare et al. 2006), consistent with the study's nested sampling design (Steele Creek Park ⊂ Sullivan County ⊂ NE Tennessee ⊂ Eastern Tennessee).

| Group | z (Original) | z (NLS) | Range Context |
|---|---|---|---|
| Total Herpetofauna | 0.1137 | 0.1478 | Mainland/nested (0.10–0.20) |
| Amphibians | 0.1449 | 0.1616 | Mainland/nested |
| Reptiles | 0.0747 | 0.1205 | Below range (original); within range (NLS) |

The original reptile z of 0.0747 fell below the expected mainland range, suggesting almost no area effect. Under NLS the estimate shifts to 0.12, placing it within the expected range — though the poor overall fit (R² = 0.66) means neither estimate is reliable.

## The Reptile Anomaly

The most analytically significant finding is the reptile SAR's departure from the power model. Examining the raw data reveals why:

| Site | Area (km²) | Reptile Species |
|---|---|---|
| Steele Creek Park | 9.3 | 21 |
| Sullivan County | 1,114 | 21 |
| NE Tennessee | 4,137 | 24 |
| Eastern Tennessee | 37,438 | 44 |

Reptile richness is flat from 9.3 km² to 1,114 km² (21 species at both scales), increases modestly to 24 at 4,137 km², then jumps to 44 at 37,438 km². This step-function pattern is fundamentally non-power-law. The residuals confirm the systematic misfit:

| Site | Power Residual | Linear Residual |
|---|---|---|
| Steele Creek Park | +6.66 | +0.09 |
| Sullivan County | −4.52 | −0.59 |
| NE Tennessee | −5.89 | +0.54 |
| Eastern Tennessee | +5.02 | −0.04 |

The linear model (S = 20.90 + 0.000618·A) produces dramatically smaller residuals and achieves R² = 0.998 vs. 0.662 for the power model. The multi-model comparison confirms this: the linear model ranks first among 2-parameter models for reptiles by both AIC and BIC.

This suggests that reptile richness in this region may be driven by a threshold effect at the largest spatial scale rather than by smooth area-dependent scaling — potentially reflecting the geographic range boundaries of southern Appalachian reptile species, which may require the full eastern Tennessee extent to capture meaningful turnover.

## Multi-Model Rankings

### The AICc Limitation

With n = 4 observations and k ≥ 3 estimated parameters (the `sars` library counts σ), AICc is undefined for every model:

> AICc = AIC + 2k(k+1) / (n − k − 1)

When n ≤ k + 1, the denominator is ≤ 0, producing AICc = ∞. This means Akaike weights are all zero and model averaging via `sars.sar_average()` cannot produce meaningful results.

This is not a limitation of the library but a fundamental constraint of the sample size. AIC and BIC remain finite and are used here for cautious ranking.

The same constraint affects `sars.bootstrap_ci()`, which generates confidence bands by resampling rows with replacement and computing model-averaged predictions across bootstrap replicates. Because the underlying model averaging relies on Akaike weights, the bootstrap intervals are not meaningful when those weights are uniformly zero. With n ≥ 6–7, both tools become available and would substantially strengthen uncertainty quantification.

### Rankings by BIC

**Herpetofauna** — The top models are 3-parameter models (epm2, powerR, epm1) with R² > 0.99, but these have only 1 residual degree of freedom. The 2-parameter power model (R² = 0.92) is more trustworthy at this sample size.

**Amphibians** — Similar ranking with epm2 > powerR > epm1, but the standard power model already provides an excellent fit (R² = 0.98) that is only marginally improved by adding a third parameter.

**Reptiles** — The linear model ranks first or near-first among 2-parameter models (BIC = 8.25), outperforming the power model (BIC = 29.27) by a wide margin. The top 3-parameter models (asymp, ratio, powerR) achieve R² > 0.998 but with effectively linear behavior — the asymptotic model converges to d = 78.5 with z = 1.38 × 10⁻⁵, producing a curve indistinguishable from a straight line.

### Curve Shape

For herpetofauna and amphibians, all well-fitting models produce **convex** shapes (upward curvature without asymptote), consistent with the original finding and with typical nested SAR data (Tjørve & Tjørve 2021).

For reptiles, although fitted models are technically classified as convex, the best-fitting asymptotic models have approach rates so low that they are functionally linear across the observed range. The data do not support distinguishing between convex, asymptotic, and linear forms for this group.

## Predictive Validation

Testing both fitting methods against independent sites from the original publication:

| Site | Area (km²) | Amp. Reported | Amp. NLS | Amp. Original | Rep. Reported | Rep. NLS | Rep. Original |
|---|---|---|---|---|---|---|---|
| GSMNP | 2,114 | 43 | 42.3 | 42.7 | 38 | 27.6 | 27.4 |
| Cumberland Gap NHP | 82.74 | 28 | 25.1 | 26.7 | 20 | 18.7 | 21.5 |
| Radford AAP | 32.7 | 19 | 21.6 | 23.3 | 14 | 16.7 | 20.1 |
| Henderson Wetland | 0.101 | 12 | 8.5 | 10.1 | —* | 8.3 | 13.0 |

\* Henderson Wetland reptile count is unavailable because the site was not surveyed for reptiles in the original study. Reptile predictions are shown for completeness but cannot be validated.

For amphibians, both methods produce similar predictions. The NLS model gives slightly lower estimates that happen to be closer to reported values at GSMNP (42.3 vs. 43) and Radford AAP, but slightly farther from Cumberland Gap.

For reptiles, both methods substantially underpredict GSMNP (reported 38 vs. predicted ~27–28), confirming the weak predictive power of the reptile SAR regardless of fitting method. This underprediction was noted in the original publication and attributed to the non-significant regression relationship.

At small areas (Henderson Wetland, 0.101 km²), the methods diverge more noticeably: the NLS model predicts 8.5 amphibians vs. 10.1 from the original. The reported count of 12 suggests that both models underpredict at very small, hyper-diverse sites — consistent with the original discussion.

## Threshold Analysis

`sars.sar_threshold()` tests for breakpoints (small-island effect) by comparing continuous two-slope, left-horizontal + right slope, and simple linear models. For all three taxonomic groups, the analysis selects the simple linear model with no breakpoint detected.

This is expected at n = 4 — insufficient data to support a piecewise structure — but documents the absence of statistical support for a small-island effect at this sample size.

## Limitations

All conclusions in this analysis are based on n = 4 data points, which is extremely small for any regression methodology. With only 4 observations:

- **Degrees of freedom are minimal.** Two-parameter models have only 2 residual degrees of freedom; 3-parameter models have just 1. Any model with 3+ parameters can fit the data nearly perfectly regardless of its ecological plausibility.
- **AICc is undefined** for all models (see "The AICc Limitation" above), precluding standard model averaging and Akaike weight comparisons.
- **Statistical power is negligible.** Confidence intervals on parameters are wide, and distinguishing between competing model forms (power vs. logarithmic vs. linear) is not statistically meaningful at this sample size.
- **The NLS vs. OLS comparison is illustrative, not definitive.** Both methods fit the amphibian data well and the reptile data poorly. The parameter shifts between methods, while consistent with theoretical expectations, should not be over-interpreted given the sample size.

These limitations apply equally to the original publication's log-linear analysis and to the NLS reanalysis presented here. The results are best viewed as a baseline that motivates expanded sampling at intermediate spatial scales.

## Recommendations for Further Work

**Expand the dataset.** The single most impactful improvement would be adding intermediate spatial scales (individual watersheds, additional well-surveyed parks, or state-level data). Even 2–3 additional data points would make AICc finite for 2-parameter models and unlock model averaging and bootstrap confidence intervals.

**Taxonomic decomposition.** Disaggregating reptiles into lizards, snakes, and turtles could reveal whether the weak SAR is driven by a specific subgroup — particularly turtles, whose richness may be more habitat-specific (aquatic habitat availability) than area-dependent.

**Integration with iNaturalist.** The `sars` library supports `sars.from_pyinaturalist()` for direct ingestion of citizen-science observation data, which could update species counts since the 2022 publication.

**GIS integration.** If park boundary and county shapefiles are available, `sars.from_geodataframe()` can ingest GeoDataFrames directly with automatic area calculation, streamlining the addition of new spatial scales.

## Summary of Key Findings

1. The NLS fitting method in `sars` produces systematically lower c and higher z values than the original log-linear regression, as expected from the different error assumptions.

2. The qualitative conclusions of the original publication are fully supported: amphibians show a strong power-law SAR (R² ≈ 0.98), total herpetofauna are moderately well-fit (R² ≈ 0.92), and reptiles are poorly fit by the power model (R² ≈ 0.66).

3. Multi-model comparison reveals that reptile richness in this region is better described by a linear relationship than by any curvilinear SAR model, driven by the identical species count at 9.3 km² and 1,114 km² (both 21 species).

4. AICc-based model selection, model averaging, and bootstrap confidence intervals are fundamentally inapplicable at n = 4. BIC and AIC remain usable for cautious ranking.

5. All z-values fall within the expected mainland/nested range (0.10–0.20), consistent with the study's spatial design.

6. Predictive accuracy at independent test sites is similar between the two methods for amphibians. Both methods fail for reptiles, confirming the original conclusion that reptile richness in eastern Tennessee is not reliably predicted by a simple area-based model.

## References

- Drakare, S., Lennon, J.J. & Hillebrand, H. (2006). The imprint of the geographical, evolutionary and ecological context on species–area relationships. *Ecology Letters* 9:215–227.
- Jessee, L.D., Stout, J.B. & McMeen, J.N. (2022). Herpetofauna of Steele Creek Park (Sullivan County, TN), with Comments on Species–Area Relationships of Amphibians and Reptiles in Eastern Tennessee. *Southeastern Naturalist* 21(1):63–73.
- Rosenzweig, M.L. (1995). *Species Diversity in Space and Time*. Cambridge University Press.
- Tjørve, E. & Tjørve, K.M.C. (2021). Mathematical expressions for the species–area relationship and the assumptions behind the models. In: Matthews, T.J., Triantis, K.A. & Whittaker, R.J. (Eds.), *The Species–Area Relationship: Theory and Application*. Cambridge University Press, pp. 157–184.
