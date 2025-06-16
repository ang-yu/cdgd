
# cdgd

[![R-CMD-check](https://github.com/ang-yu/cdgd/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ang-yu/cdgd/actions/workflows/R-CMD-check.yaml)
![CRAN Downloads
overall](https://cranlogs.r-pkg.org/badges/grand-total/cdgd)

The package cdgd implements the causal decompositions of group
disparities in [Yu and Elwert
(2025)](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-19/issue-1/Nonparametric-causal-decomposition-of-group-disparities/10.1214/24-AOAS1990.full).

## Installation

The latest release of the package can be installed through CRAN.

``` r
install.packages("cdgd")
```

The current development version can be installed from source using
devtools.

``` r
devtools::install_github("ang-yu/cdgd")
```

## Examples

``` r
library(cdgd)  

# load the simulated example data
data(exp_data)
head(exp_data)
#>       outcome treatment  confounder          Q group_a
#> 748 1.4608165         1  0.26306864  0.6748330       0
#> 221 0.4777308         0  1.30296394  0.5920512       1
#> 24  0.8760129         1 -1.49971226  1.6294327       1
#> 497 0.4131192         1 -1.17219619 -0.8391873       1
#> 249 2.0483222         1  1.71790879  2.9546966       1
#> 547 0.1912013         0 -0.02438458 -0.3704544       0
```

### Use cdgd0_ml, cdgd0_pa, or cdgd0_manual for unconditional decomposition

``` r
results0 <- cdgd0_pa(Y="outcome",D="treatment",G="group_a",X=c("confounder","Q"),data=exp_data,alpha=0.05)

round(results0$results, 4)
#>              point     se p_value CI_lower CI_upper
#> total       0.2675 0.0390  0.0000   0.1911   0.3439
#> baseline    0.0421 0.0131  0.0013   0.0164   0.0678
#> prevalence  0.2579 0.0337  0.0000   0.1919   0.3240
#> effect     -0.1372 0.0209  0.0000  -0.1781  -0.0963
#> selection   0.1047 0.0150  0.0000   0.0754   0.1340
```

### Use cdgd1_ml, cdgd1_pa, or cdgd1_manual for conditional decomposition

``` r
results1 <- cdgd1_pa(Y="outcome",D="treatment",G="group_a",X="confounder",Q="Q",data=exp_data,alpha=0.05)

round(results1, 4)
#>                                 point     se p_value CI_lower CI_upper
#> total                          0.2675 0.0390  0.0000   0.1911   0.3439
#> baseline                       0.0421 0.0131  0.0013   0.0164   0.0678
#> conditional prevalence         0.2032 0.0371  0.0000   0.1305   0.2760
#> conditional effect            -0.1644 0.0220  0.0000  -0.2076  -0.1212
#> conditional selection          0.0875 0.0143  0.0000   0.0595   0.1156
#> Q distribution                 0.0990 0.0188  0.0000   0.0621   0.1359
#> conditional Jackson reduction  0.2362 0.0378  0.0000   0.1621   0.3103
```
