# Summary-statistic RAS scan

Summary-statistic RAS scan

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

  Discovery marginal effect-size weights.

- z_targ:

  Target marginal Z-scores.

- R:

  LD correlation matrix.

- mask:

  Logical SNP mask.

- skip1:

  Pivotal SNP step size.

- skip2:

  Window-size step size.

- min_window_size:

  Minimum half-window size.

- max_window_size:

  Maximum half-window size.

## Value

List with `x` and `y`.
