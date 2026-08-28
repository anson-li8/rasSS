# Summary-Statistic RAS Profile Generation.

Generates a Regional Association Score (RAS) profile using only summary
statistics. For each pivotal SNP, the function scans windows at adaptive
sizes, calculates the burden statistic
[`t_burden`](https://anson-li8.github.io/rasSS/reference/t_burden.md),
and retains the minimum p-value (highest RAS) across window sizes to
build the profile.

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

  Numeric vector. Genome-wide discovery marginal effect sizes (used as
  weights \\w\\).

- z_targ:

  Numeric vector. Genome-wide target marginal Z-scores (vector \\Z\\).

- R:

  Numeric matrix. Genome-wide LD correlation matrix (e.g., from an
  external reference panel).

- mask:

  Logical vector. Indicates which SNPs are retained after LD pruning
  (e.g., from
  [`prune_ld()`](https://anson-li8.github.io/rasSS/reference/prune_ld.md)).

- skip1:

  Integer. Step size for pivotal SNPs (default: 10).

- skip2:

  Integer. Step size for adaptive window size increment (default: 3).

- min_window_size:

  Integer. Minimum half-window size in SNPs (default: 3).

- max_window_size:

  Integer. Maximum half-window size in SNPs (default: 30).

## Value

A list with two elements:

- `x`: Integer vector of the genomic indices of the evaluated pivotal
  SNPs .

- `y`: Numeric vector of the \\-\log\_{10}(p)\\-values from the most
  significant window at each pivotal SNP.
