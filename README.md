# Weighted Anchors: An Uncertainty-Weighted Multiple Comparisons Method

A new post-hoc multiple comparisons procedure for one-way ANOVA that allocates significance budget based on how "uncertain" each pairwise comparison is.

---

## Method Description

Standard multiple comparison methods (like Tukey's HSD) apply a uniform significance threshold to all pairwise comparisons. This method takes a different approach: comparisons whose p-values fall near the significance boundary are the most "uncertain," and so they receive a larger share of the significance budget. However, so that we can still maintain our desired family-wise error rate, we take this budget away from the "certain" p-values, those near 0 and 1. This adaptive allocation is controlled by a tuning parameter `eta` (η).

---

## Mathematical Formula

### Setup

Let there be $I$ groups each with $J$ observations. For each pair $(i, j)$, compute the test statistic:

$$t_{ij} = \frac{|\bar{x}_i - \bar{x}_j|}{\sqrt{2 \cdot \text{MSE} / J}}$$

with degrees of freedom $df = I(J-1)$, giving a two-sided p-value:

$$p_{ij} = 2\left(1 - pt(|t_{ij}|, df)\right)$$

### Šidák Anchor

The default significance threshold (anchor) is set using the Šidák correction for $m = \binom{I}{2}$ comparisons:

$$\text{anchor} = 1 - (1 - \alpha)^{1/m}$$

### Uncertainty
Compute the uncertainty of each p-value $p_i$ as:

$$\mathrm{uncertainty} = \frac{1}{\mathrm{distance}} = \frac{1}{\left|\log \frac{\alpha_1}{p_i}\right|}$$

Distance here is uniquely defined as the absolute log-ratio of the p-value to the anchor. Taking the logarithm "tames" exponential differences, and the absolute value ensures distances are the same in either direction. Uncertainty is the reciprocal of distance: the closer a p-value is to the anchor, the more uncertain we are.

### Weights
Each p-value's weight is its proportion of total uncertainty:

$$w_i = \frac{\mathrm{uncertainty}_i}{\sum_k \mathrm{uncertainty}_k}$$

### Threshold for $p_i > \alpha_1$ (non-rejected comparisons)
Add budget proportional to weight, scaled by $\eta$:

$$\mathrm{new\ threshold} = \alpha_1 + w_i \eta$$

### Threshold for $p_i \leq \alpha_1$ (rejected comparisons)
Subtract budget proportional to weight, scaled by $\eta$, with a floor to ensure the threshold stays positive:

$$\mathrm{new\ threshold} = \alpha_1 - \min\left(\frac{1}{\mathrm{weight}},\ \left(\frac{\alpha_1}{\eta} - \epsilon\right) * \eta\right)$$

The $\left(\frac{\alpha_1}{\eta} - \epsilon\right)$ term guarantees the threshold stays above 0, since:

$$\alpha_1 - \left(\frac{\alpha_1}{\eta} - \epsilon\right) \cdot \eta = \alpha_1 - \alpha_1 + \epsilon\eta = \epsilon\eta > 0$$

Note that $\frac{1}{w_i}$ is used here (rather than $w_i$ as in step (e)) because we want to de-allocate more budget from comparisons with smaller weights.

### Significance
A comparison is declared significant if $p_i < \mathrm{new\ threshold}$.


---

## Installation / Usage

### Requirements

- R (Latest version recommended)
- The following R packages:

```r
install.packages("progress")
install.packages("knitr")
install.packages("tidyr")
install.packages("ggplot2")
```

### Repository Structure

| File | Description |
|------|-------------|
| `mc_project budget (main).R` | All function definitions (`my_comparison`, `calibrate_eta`, `prop_dist`, `simulate_fwer`, `simulate_power`, `tukey_hsd`) |
| `simulation.R` | Executable code — runs the oil filter demo, FWER simulations, and power analysis |
| `example_data_demo.R` | Standalone example using teaching method data |

### Setup

Clone this repository, then source the functions file to load all functions into your environment:

```r
source("mc_project budget (main).R")
```

To reproduce the full simulation results and power analysis:

```r
source("simulation.R")
```


### Key Functions

#### `my_comparison(means, J, MSE, alpha = 0.05, eta = 1)`

The main comparison function.

| Parameter | Description |
|-----------|-------------|
| `means`   | Numeric vector of group sample means |
| `J`       | Number of observations per group |
| `MSE`     | Pooled mean squared error from ANOVA |
| `alpha`   | Desired FWER level (default: 0.05) |
| `eta`     | Budget redistribution tuning parameter (default: 1) |

Returns a data frame with one row per pair containing: `group1`, `group2`, `diff`, `pval`, `dist`, `uncertainty`, `weights`, `threshold`, and `significant`.

---

#### `calibrate_eta(I, J, target_fwer = 0.05, n_sim = 10000, tolerance = 0.001, precision = 0.001)`

Uses binary search to find the optimal `eta` for a given group structure that achieves the target FWER. You only need to calibrate once per unique combination of `I` and `J`.

| Parameter | Description |
|-----------|-------------|
| `I`       | Number of groups |
| `J`       | Number of observations per group |
| `target_fwer` | Desired FWER (default: 0.05) |
| `n_sim`   | Number of simulations for estimation (default: 10,000) |
| `tolerance` | Acceptable distance below target FWER (default: 0.001) |
| `precision` | Binary search stopping criterion (default: 0.001) |

```r
eta <- calibrate_eta(I = 5, J = 10)
results <- my_comparison(means, J, MSE, eta = eta)
```

#### `prop_dist(x, y)`

Helper function used internally. Returns the log-ratio distance between a p-value and the anchor:

$$\mathrm{prop\_dist}(\alpha_1, p_i) = \left|\log \frac{\alpha_1}{p_i}\right|$$

---

## Example Output

The following example uses the oil filter dataset from "Modern Mathematical Statistics with Applications" by Jay L. Devore, et al., Example 11.5, with 5 groups, 9 observations per group, and a pooled MSE of 0.088.

```r
oil_means <- c(14.5, 13.8, 13.3, 14.3, 13.1)
J <- 9
MSE <- 0.088
global_eta <- 0.001

demonstration <- my_comparison(oil_means, J, MSE, eta = global_eta)
kable(demonstration[order(demonstration$diff), ])
```

```
| group1 | group2 | diff |      pval | dist      | uncertainty | weights   | threshold | significant |
|-------:|-------:|-----:|----------:|----------:|------------:|----------:|----------:|:------------|
|      1 |      4 |  0.2 | 0.1604268 |  3.445426 |   0.2902398 | 0.1250652 | 0.0052413 | FALSE       |
|      3 |      5 |  0.2 | 0.1604268 |  3.445426 |   0.2902398 | 0.1250652 | 0.0052413 | FALSE       |
|      2 |      5 |  0.5 | 0.0009317 |  1.703169 |   0.5871409 | 0.2530007 | 0.0011636 | TRUE        |
|      2 |      4 |  0.5 | 0.0009317 |  1.703169 |   0.5871409 | 0.2530007 | 0.0011636 | TRUE        |
|      1 |      2 |  0.7 | 0.0000116 |  6.086777 |   0.1642906 | 0.0707933 | 0.0011636 | TRUE        |
|      2 |      5 |  0.7 | 0.0000116 |  6.086777 |   0.1642906 | 0.0707933 | 0.0011636 | TRUE        |
|      3 |      4 |  1.0 | 0.0000000 | 12.999620 |   0.0769253 | 0.0331473 | 0.0011636 | TRUE        |
|      1 |      3 |  1.2 | 0.0000000 | 17.479246 |   0.0572107 | 0.0246523 | 0.0011636 | TRUE        |
|      4 |      5 |  1.2 | 0.0000000 | 17.479246 |   0.0572107 | 0.0246523 | 0.0011636 | TRUE        |
|      1 |      5 |  1.4 | 0.0000000 | 21.730182 |   0.0460189 | 0.0198297 | 0.0011636 | TRUE        |
```


---

## When to Use This Method vs. Tukey's HSD

### Use this method when:
- You want to detect real, but small differences
- You have a **moderate number of groups** where some comparisons are near the boundary of significance and careful budget allocation matters
- You are willing to accept **slightly variable FWER** in exchange for improved power on borderline comparisons
- You want a **tunable procedure** — `eta` allows you to control how aggressively budget is redistributed

### Use Tukey's HSD when:
- You need to detect real differences with higher effect sizes
- You have **many groups** where the behavior of the uncertainty weights becomes harder to interpret
- You need a method that is **immediately recognized** by reviewers and journals
- **Simplicity and reproducibility** are priorities over adaptive behavior

---

## License

MIT License. See `LICENSE` for details.

---

## Acknowledgments
Simulation code received from Dr. Sybil Prince-Nelson; Department of Mathematics, Washington and Lee University.

Feedback and debugging assistance provided with help from Claude (Anthropic).
