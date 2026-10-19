# qountstat

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/qountstat)](https://CRAN.R-project.org/package=qountstat)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- DOI-BADGE -->
<!-- badges: end -->

Statistical analysis of count data and quantal data in ecotoxicology.

The package determines no/lowest observed (adverse) effect concentrations
(NOEC/LOEC) from the data of laboratory and field biotests and provides
diagnostics for choosing an appropriate method.

## Methods

| Function | Purpose |
|---|---|
| `CPCAT()` | Closure Principle Computational Approach Test for count data |
| `Dunnett.GLM()` | Dunnett contrasts on a GLM, with a quasi-Poisson, negative binomial or Poisson distribution |
| `CPFISH()` | Closure Principle Fisher-Freeman-Halton test for quantal data |
| `Overdispersion.test()` | Simulation-calibrated diagnostic for overdispersion |

Each method is accompanied by functions for the statistical power and for the
bootstrapped minimum detectable difference, for example `CPCAT.power()` and
`CPCAT.bMDD()`.

## Installation

From CRAN:

```r
install.packages("qountstat")
```

The development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("YOUR-GITHUB-USERNAME/qountstat")
```

## Example

```r
library(qountstat)

# Check whether the counts are compatible with a Poisson distribution
Overdispersion.test(groups = Daphnia.counts$Concentration,
                    counts = Daphnia.counts$Number_Young)

# Determine the NOEC with the CPCAT
CPCAT(groups = Daphnia.counts$Concentration,
      counts  = Daphnia.counts$Number_Young,
      use.fixed.random.seed = 123)
```

## Citation

```r
citation("qountstat")
```

## References

Lehmann, R.; Bachmann, J.; Maletzki, D.; Polleichtner, C.; Ratte, H.; Ratte, M. (2016):
A new approach to overcome shortcomings with multiple testing of reproduction data in
ecotoxicology. *Stochastic Environmental Research and Risk Assessment* 30(3), 871-882.
<https://doi.org/10.1007/s00477-015-1079-4>

Hothorn, L.; Kluxen, F. (2020): Statistical analysis of no observed effect concentrations
or levels in eco-toxicological assays with overdispersed count endpoints. *bioRxiv*.
<https://doi.org/10.1101/2020.01.15.907881>

Daniels, B.; Kotschik, P.; Walter-Rohde, S.; Graeff, T. (2026): When Counts Deviate:
Robust Overdispersion Diagnostics for Ecotoxicity Data. Poster, 36th SETAC Europe
Annual Meeting, Maastricht.

## License

MIT, see [LICENSE](LICENSE).
