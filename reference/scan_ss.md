# Summary-Statistic RAS Profile Generation

Generates a Regional Association Score (RAS) profile using only summary
statistics. For each pivotal SNP, it runs an adaptive window scan,
calculates the burden statistic
[`t_burden`](https://anson-li8.github.io/rasSS/reference/t_burden.md),
and keeps the minimum p-value (highest RAS) across the grid.

## Usage

``` r
scan_ss(
  b_disc,
  z_targ,
  R,
  mask,
  skip1 = 10,
  skip2 = 3,
  min_window_size = 3,
  max_window_size = 30
)
```

## Arguments

- b_disc:

  Numeric vector. Genome-wide discovery marginal effect sizes.

- z_targ:

  Numeric vector. Genome-wide target marginal Z-scores.

- R:

  Numeric matrix. Genome-wide LD correlation matrix (e.g., from an
  external reference panel).

- mask:

  Logical vector. Indicates which SNPs are retained after LD pruning.

- skip1:

  Integer. Step size for pivotal SNPs (e.g., evaluate every 10th SNP).

- skip2:

  Integer. Step size for the adaptive window scan at each SNP.

- min_window_size:

  Integer. Minimum half-window size (in SNPs).

- max_window_size:

  Integer. Maximum half-window size (in SNPs).

## Value

A list with two elements:

- `x`: Integer vector of the genomic indices of the pivotal SNPs
  evaluated.

- `y`: Numeric vector of the \\-\log\_{10}(p)\\-values from the most
  significant adaptive window at each pivotal SNP.
