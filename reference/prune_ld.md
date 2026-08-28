# Greedy LD pruning

Removes redundant SNPs based on the criterion of squared correlation
(\\r^2\\) exceeding a specified threshold. Traverses SNPs in order and
drops any following SNPs that are redundant with a retained one.

## Usage

``` r
prune_ld(R, thresh)
```

## Arguments

- R:

  Numeric matrix. LD correlation matrix.

- thresh:

  Numeric. Squared correlation threshold (e.g., 0.2). SNPs with \\r^2 \>
  thresh\\ compared to a retained SNP are eliminated.

## Value

Logical vector of length `nrow(R)`; retained SNPs are indicated by
`TRUE`.
