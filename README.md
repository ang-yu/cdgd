
# cdgd

<!-- badges: start -->
<!-- badges: end -->

The goal of cdgd is to implement the causal decomposition of group
disparities of Yu and Elwert (2023).

## Installation

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

### Use cdgd0_ml, cdgd0_pa, or cdgd0_manual to get unconditional decomposition

``` r
results0 <- cdgd0_pa(Y="outcome",D="treatment",G="group_a",X=c("confounder","Q"),data=exp_data,alpha=0.05)

round(results0$results, 4)
#>              point     se p_value CI_lower CI_upper
#> total       0.2675 0.0390   0.000   0.1911   0.3439
#> baseline    0.0400 0.0129   0.002   0.0146   0.0653
#> prevalence  0.2565 0.0335   0.000   0.1908   0.3222
#> effect     -0.1365 0.0208   0.000  -0.1773  -0.0957
#> selection   0.1075 0.0139   0.000   0.0802   0.1348
```

### Use cdgd1_ml, cdgd1_pa, or cdgd1_manual to get conditional decomposition

``` r
results1 <- cdgd1_pa(Y="outcome",D="treatment",G="group_a",X="confounder",Q="Q",data=exp_data,alpha=0.05)

round(results1, 4)
#>                                 point     se p_value CI_lower CI_upper
#> total                          0.2675 0.0390   0.000   0.1911   0.3439
#> baseline                       0.0400 0.0129   0.002   0.0146   0.0653
#> conditional prevalence         0.2086 0.0351   0.000   0.1398   0.2774
#> conditional effect            -0.1529 0.0205   0.000  -0.1930  -0.1128
#> conditional selection          0.0888 0.0135   0.000   0.0624   0.1151
#> Q distribution                 0.0830 0.0110   0.000   0.0615   0.1045
#> conditional Jackson reduction  0.2393 0.0355   0.000   0.1697   0.3089
```
