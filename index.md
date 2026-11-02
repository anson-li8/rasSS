# rasSS

rasSS (*R*egional *A*ssociation *S*core *S*ummary *S*tatistic extension)
is an R-package that extends the functionality of RAS to accept summary
statistics with a comparable power and Type I error rate as the original
method. The original method is described in [Jiang et
al. (2025)](https://doi.org/10.1073/pnas.2419721122) and found to have
significant advantages over other GWAS algorithms, such as SKAT,
CauchyGM, and Burden. In search of a faithful summary-statistic analog
to the regional association scores derived from regressions on the
individual genotypic and phenotypic data, rasSS implements a one
degree-of-freedom weighted burden test statistic, evaluating the
association of the adaptive genomic windows from typical GWAS summary
statistics. The algorithm takes these evaluated association scores and
runs the same workflow of Change Point Detection (CPD) to find the
significant genomic regions. The package’s main contribution is access,
allowing the RAS algorithm to be used on studies that only publish
summary statistics, and avoiding the restrictions of privacy and data
storage.

## Installation

You can install the development version of rasSS from Github with:

``` r

# install.packages("pak")
pak::pak("shan-team/RAS_GWAS_Summary")
```

## Example

The core of `rasSS` is the
[`scan_block()`](https://anson-li8.github.io/rasSS/reference/scan_block.md)
function, which ingests discovery and target GWAS summary statistics,
extracts the corresponding LD matrix from a PLINK reference panel,
harmonizes the alleles, and runs the regional association scan.

``` r

library(rasSS)

# Run a full end-to-end scan on a genomic block
res <- scan_block(
  disc_path   = "data/sumstats/discovery_gwas.tsv.gz",
  targ_path   = "data/sumstats/target_gwas.tsv.gz",
  chr         = 6,
  start_bp    = 29000000,
  end_bp      = 34000000,
  ref_prefix  = "data/1000G/chr6_eur",
  out_rds     = "data/results/chr6_mhc.rds"
)

# View validated peaks mapped to base-pair positions
print(res$peaks)
```

For a complete walkthrough using real Schizophrenia GWAS data and the
1000 Genomes LD reference, see the **[Real-Data
Analysis](https://anson-li8.github.io/rasSS/articles/real-data.md)**
vignette.
