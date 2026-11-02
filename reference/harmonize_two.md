# Harmonize two independent sumstats against an LD reference

Used for valid inference where discovery weights (`w`) and target
Z-scores (`Z`) must come from independent cohorts.

## Usage

``` r
harmonize_two(disc, targ, bim, ld, freq_tol = 0.15)
```

## Arguments

- disc:

  Discovery sumstats (provides weights `w`).

- targ:

  Target sumstats (provides Z-scores `Z`).

- bim:

  Output of
  [`read_bim`](https://anson-li8.github.io/rasSS/reference/read_bim.md).

- ld:

  Output of
  [`ld_matrix`](https://anson-li8.github.io/rasSS/reference/ld_matrix.md).

- freq_tol:

  Maximum \|EAF - reference frequency\| tolerance.

## Value

List with w, Z, bp, R, ids, n_kept, all aligned to each other.
