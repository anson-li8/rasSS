# Simulate genotypes from an LD matrix

Using a multivariate normal distribution constructed by the provided LD
correlation matrix, this function simulates genotypes in a population,
then rounds and clips to integer genotypic dosages (0, 1, 2).

## Usage

``` r
sim_genotypes(n, R)
```

## Arguments

- n:

  Integer. Number of individuals to simulate.

- R:

  Numeric matrix. LD correlation matrix.

## Value

An `n` x `ncol(R)` integer matrix of simulated genotype dosages (values
0, 1, or 2).
