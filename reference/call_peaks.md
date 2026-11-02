# Run RAS change-point detection and map peaks to base pairs

Wraps [`RAS::ras_detect`](https://rdrr.io/pkg/RAS/man/ras_detect.html)
and
[`RAS::ras_validate`](https://rdrr.io/pkg/RAS/man/ras_validate.html).
Correctly maps the validated `tau_hats` (which are full-SNP-scale
indices) back to base-pair positions.

## Usage

``` r
call_peaks(
  scan,
  bp_per_snp,
  window_size = 100,
  slope_check_window_size = 40,
  slope_threshold = 1e-04,
  davies_threshold = 1e-05,
  min_signal = 2.5,
  skip = 10
)
```

## Arguments

- scan:

  Output of
  [`scan_ss`](https://shan-team.github.io/RAS_GWAS_Summary/reference/scan_ss.md).

- bp_per_snp:

  Base-pair positions aligned to the SNPs given to scan_ss.

- window_size:

  CPD window size.

- slope_check_window_size:

  CPD slope check window.

- slope_threshold:

  CPD slope p-value threshold.

- davies_threshold:

  CPD Davies p-value threshold.

- min_signal:

  Minimum RAS height to validate a peak.

- skip:

  The `skip1` parameter used in `scan_ss`.

## Value

List with snp_index, grid_index, and bp.
