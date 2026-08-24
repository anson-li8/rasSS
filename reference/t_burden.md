# Summary-Statistic Weighted Burden Statistic

Calculates the 1-degree-of-freedom weighted burden statistic \\T = w'Z /
\sqrt{w'Rw}\\ of a given window of SNPs. Under the null hypothesis of no
association in the window, \\T\\ asymptotically follows a standard
normal distribution \\N(0,1)\\.

## Usage

``` r
t_burden(w, Z, R)
```

## Arguments

- w:

  Numeric vector. Discovery weights (e.g., marginal effect sizes from
  the discovery cohort) for the SNPs in the window.

- Z:

  Numeric vector. Marginal Z-scores from the target cohort for the SNPs
  in the window. Must be the same length as `w`.

- R:

  Numeric matrix. Linkage Disequilibrium (LD) correlation matrix for the
  SNPs in the window. Dimensions must be `length(w)` x `length(w)`.

## Value

A single numeric value representing the burden test statistic. Returns
`0` if the variance denominator is non-positive or non-finite.
