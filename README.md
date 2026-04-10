# Uncertainty-Weighted Multiple Comparisons Method

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

**(c) Uncertainty.** Compute the uncertainty of each p-value $p_i$ as:

$$\mathrm{uncertainty} = \frac{1}{\mathrm{distance}} = \frac{1}{\left|\log \frac{\alpha_1}{p_i}\right|}$$

Distance here is uniquely defined as the absolute log-ratio of the p-value to the anchor. Taking the logarithm "tames" exponential differences, and the absolute value ensures distances are the same in either direction. Uncertainty is the reciprocal of distance — the closer a p-value is to the anchor, the more uncertain we are.

**(d) Weights.** Each p-value's weight is its proportion of total uncertainty:

$$w_i = \frac{\mathrm{uncertainty}_i}{\sum_k \mathrm{uncertainty}_k}$$

**(e) Threshold for $p_i > \alpha_1$ (non-rejected comparisons).** Add budget proportional to weight, scaled by $\eta$:

$$\mathrm{new\ threshold} = \alpha_1 + w_i \cdot \eta$$

**(f) Threshold for $p_i \leq \alpha_1$ (rejected comparisons).** Subtract budget proportional to weight, scaled by $\eta$, with a floor to ensure the threshold stays positive:

$$\mathrm{new\ threshold} = \alpha_1 - \min\!\left(\frac{1}{w_i},\ \frac{\alpha_1}{\eta} - \epsilon\right) \cdot \eta$$

The $\left(\frac{\alpha_1}{\eta} - \epsilon\right)$ term guarantees the threshold stays above 0, since:

$$\alpha_1 - \left(\frac{\alpha_1}{\eta} - \epsilon\right) \cdot \eta = \alpha_1 - \alpha_1 + \epsilon\eta = \epsilon\eta > 0$$

Note that $\frac{1}{w_i}$ is used here (rather than $w_i$ as in step (e)) because we want to de-allocate more budget from comparisons with smaller weights.

**(g) Significance.** A comparison is declared significant if $p_i < \mathrm{new\ threshold}$.


---

## Installation / Usage

### Requirements

- R (version 4.0 or later)
- No additional packages required (uses base R only)

### Setup

Clone this repository or download `my_comparison.R` directly:

```r
source("my_comparison.R")
```

### Function Signature

```r
my_comparison(means, J, MSE, alpha = 0.05, eta = 1)
```

| Parameter | Description |
|-----------|-------------|
| `means`   | Numeric vector of group sample means |
| `J`       | Number of observations per group |
| `MSE`     | Mean squared error (pooled) from ANOVA |
| `alpha`   | Familywise error rate target (default: 0.05) |
| `eta`     | Tuning parameter controlling budget redistribution (default: 1) |

### Returns

A data frame with one row per pair containing:
- `group1`, `group2` — group indices
- `diff` — absolute difference in means
- `pval` — two-sided p-value
- `weights` — uncertainty weight
- `threshold` — adaptive significance threshold
- `significant` — logical: is the comparison significant?

---

## Example Output

```r
# Simulated example with 3 groups, 10 obs each
set.seed(42)
J <- 10
data <- matrix(rnorm(3 * J, mean = c(0, 0.5, 1.5)[rep(1:3, each = J)]), nrow = 3, byrow = TRUE)
means <- rowMeans(data)
MSE <- sum(apply(data, 1, function(row) sum((row - mean(row))^2))) / (3 * (J - 1))

results <- my_comparison(means, J, MSE, eta = 1)
print(results[, c("group1", "group2", "diff", "pval", "threshold", "significant")])
```

```
  group1 group2      diff      pval threshold significant
1      1      2 0.4821430 0.1843275 0.0177132       FALSE
2      1      3 1.5203847 0.0000412 0.0162891        TRUE
3      2      3 1.0382417 0.0063201 0.0191045        TRUE
```

---

## When to Use This Method vs. Tukey's HSD

### Use this method when:

- You have a **moderate number of groups** (3–6) where some comparisons are near the boundary of significance and careful budget allocation matters
- You are willing to accept **slightly variable FWER** in exchange for improved power on borderline comparisons
- You want a **tunable procedure** — `eta` allows you to control how aggressively budget is redistributed
- Comparisons are of **unequal scientific importance** and you want the procedure to naturally focus sensitivity where it is most needed

### Use Tukey's HSD when:

- You want a **well-established, peer-reviewed** procedure with guaranteed FWER control
- You have **unequal group sizes** (Tukey's HSD handles this via the Tukey-Kramer adjustment; this method currently assumes equal $J$)
- You need a method that is **immediately recognized** by reviewers and journals
- You have **many groups** (6+), where the behavior of the uncertainty weights becomes harder to interpret
- **Simplicity and reproducibility** are priorities over adaptive behavior

---

## FWER Simulation Results

Empirical FWER estimated over 10,000 simulations under the global null (all group means equal):

| Scenario | I | J | FWER | SE |
|----------|---|---|------|----|
| Small    | 3 | 5 | ...  | ... |
| Medium   | 4 | 10| ...  | ... |
| Large    | 5 | 20| ...  | ... |
| Many groups | 6 | 10 | ... | ... |
| Large J  | 3 | 30| ...  | ... |

*(Fill in with your simulation results)*

---

## License

MIT License. See `LICENSE` for details.

---

## Acknowledgments

Feedback and debugging assistance provided with help from Claude (Anthropic).
