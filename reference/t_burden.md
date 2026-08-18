# Summary-statistic weighted burden statistic

Summary-statistic weighted burden statistic

## Usage

``` r
t_burden(w, Z, R)
```

## Arguments

- w:

  Weight vector.

- Z:

  Marginal Z-score vector.

- R:

  LD correlation matrix.

## Value

Numeric burden statistic.

## Examples

``` r
w <- c(1, 1); Z <- c(1, 1); R <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
t_burden(w, Z, R)
#> [1] 1.154701
```
