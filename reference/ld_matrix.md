# Build LD matrix from a PLINK dosage export

Reads a PLINK2 `--export A` .raw file. Extracts the "counted" allele
from the column name suffix (e.g., `rs123_A` -\> `A`) to provide exact
allele orientation without relying on EAF.

## Usage

``` r
ld_matrix(raw_path)
```

## Arguments

- raw_path:

  Path to the .raw file.

## Value

List with G (dosages), R (correlation LD), ids, freq, and counted
allele.
