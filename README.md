
# cdgd

<!-- badges: start -->
<!-- badges: end -->

The goal of cdgd is to implement the causal decomposition of group
disparities of Yu and Elwert (2022).

## Installation

``` r
devtools::install_github("ang-yu/cdgd")
```

## Example

``` r
library(cdgd)

# load the simulated example data
data(exp_data)
head(exp_data)
#>       outcome treatment  confounder          Q group_a group_b
#> 748 0.2723801         0  0.17537909  0.4498887       0       1
#> 221 1.9326514         1  2.20197596  1.7280341       1       0
#> 24  1.3544340         1  0.33352516  2.4196218       1       0
#> 497 0.4056452         0  0.55186921  0.7738751       1       0
#> 249 2.4310022         1  2.47860586  3.3031311       1       0
#> 547 0.1561152         0 -0.01625638 -0.2469696       0       1
```

### Use cdgd0 to get point estimates

``` r
set.seed(1)
results0 <- cdgd0(Y="outcome",D="treatment",G1="group_a",G2="group_b",X=c("confounder","Q"),Q="Q",data=exp_data,t=0.05,algorithm="nnet")

results0
#>               item     point
#> 1        mean_Y_G1 1.8074994
#> 2        mean_Y_G2 0.5232842
#> 3       mean_Y0_G1 0.6251887
#> 4       mean_Y0_G2 0.4131599
#> 5        mean_D_G1 0.8440000
#> 6        mean_D_G2 0.1580000
#> 7           ATE_G1 1.3009411
#> 8           ATE_G2 1.0817598
#> 9           ATT_G1 1.4008421
#> 10          ATT_G2 0.6969896
#> 11     diff_mean_Y 1.2842152
#> 12    diff_mean_Y0 0.2120289
#> 13     diff_mean_D 0.6860000
#> 14        diff_ATE 0.2191813
#> 15         dff_ATT 0.7038524
#> 16           total 1.2842152
#> 17        baseline 0.2120289
#> 18      prevalence 0.7420872
#> 19          effect 0.1849890
#> 20       selection 0.1451101
#> 21 cond_prevalence 0.5328118
#> 22     cond_effect 0.1989146
#> 23  cond_selection 0.1683844
#> 24          Q_dist 0.1720756
```

### Use cdgd to get point estimates and confidence intervals.

``` r
# This may take a minute or so.

set.seed(1)
results <- cdgd(Y="outcome",D="treatment",G1="group_a",G2="group_b",X=c("confounder","Q"),Q="Q",data=exp_data,t=0.05,algorithm="nnet",alpha=0.05,k=20)

results
#>               item     point          se       lower     upper
#> 1        mean_Y_G1 1.8074994 0.036099772  1.74561513 1.8564409
#> 2        mean_Y_G2 0.5232842 0.011003712  0.49540941 0.5279628
#> 3       mean_Y0_G1 0.6251887 0.097155561  0.45240846 0.7996548
#> 4       mean_Y0_G2 0.4131599 0.007792466  0.38797411 0.4143755
#> 5        mean_D_G1 0.8440000 0.016761562  0.80830040 0.8726514
#> 6        mean_D_G2 0.1580000 0.012218882  0.13552361 0.1670061
#> 7           ATE_G1 1.3009411 0.091210580  1.07692064 1.4121789
#> 8           ATE_G2 1.0817598 0.116334014  0.87899330 1.0865605
#> 9           ATT_G1 1.4008421 0.112407219  1.14400179 1.5560464
#> 10          ATT_G2 0.6969896 0.043326456  0.66336890 0.8107208
#> 11     diff_mean_Y 1.2842152 0.038215434  1.21797808 1.3442409
#> 12    diff_mean_Y0 0.2120289 0.098689084  0.04802229 0.3939066
#> 13     diff_mean_D 0.6860000 0.020923182  0.65040566 0.7260189
#> 14        diff_ATE 0.2191813 0.122946768 -0.03054613 0.3307919
#> 15         dff_ATT 0.7038524 0.119996753  0.35533672 0.7601443
#> 16           total 1.2842152 0.038215434  1.21797808 1.3442409
#> 17        baseline 0.2120289 0.098689084  0.04802229 0.3939066
#> 18      prevalence 0.7420872 0.075900188  0.61179383 0.7451073
#> 19          effect 0.1849890 0.103947993 -0.02500853 0.2768375
#> 20       selection 0.1451101 0.030839840  0.06897247 0.1373094
#> 21 cond_prevalence 0.5328118 0.106490435  0.31247713 0.6450340
#> 22     cond_effect 0.1989146 0.308732662 -0.73130452 0.4469957
#> 23  cond_selection 0.1683844 0.028872747  0.05518119 0.1616613
#> 24          Q_dist 0.1720756 0.299023721  0.05323203 0.8807518
```
