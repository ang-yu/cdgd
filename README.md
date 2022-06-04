
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

# simulated example data
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
#>               item      point
#> 1        mean_Y_G1 1.80749942
#> 2        mean_Y_G2 0.52328423
#> 3       mean_Y0_G1 0.78431665
#> 4       mean_Y0_G2 0.41323307
#> 5        mean_D_G1 0.84400000
#> 6        mean_D_G2 0.15800000
#> 7           ATE_G1 1.14177258
#> 8           ATE_G2 0.90254558
#> 9           ATT_G1 1.21230186
#> 10          ATT_G2 0.69652631
#> 11     diff_mean_Y 1.28421519
#> 12    diff_mean_Y0 0.37108358
#> 13     diff_mean_D 0.68600000
#> 14        diff_ATE 0.23922700
#> 15         dff_ATT 0.51577555
#> 16           total 1.28421519
#> 17        baseline 0.37108358
#> 18      prevalence 0.61914627
#> 19          effect 0.20190759
#> 20       selection 0.09207776
#> 21 cond_prevalence 0.43522769
#> 22     cond_effect 0.22579569
#> 23  cond_selection 0.08103174
#> 24          Q_dist 0.17107650
```

### Use cdgd to get point estimates and confidence intervals.

``` r
# This may take a minute or so.

set.seed(1)
results <- cdgd(Y="outcome",D="treatment",G1="group_a",G2="group_b",X=c("confounder","Q"),Q="Q",data=exp_data,t=0.05,algorithm="nnet",alpha=0.05,k=20)

results
#>               item      point          se       lower      upper
#> 1        mean_Y_G1 1.80749942 0.036099772  1.74561513 1.85644092
#> 2        mean_Y_G2 0.52328423 0.011003712  0.49540941 0.52796277
#> 3       mean_Y0_G1 0.78431665 0.027520347  0.74979367 0.83108238
#> 4       mean_Y0_G2 0.41323307 0.007348987  0.38790382 0.41203214
#> 5        mean_D_G1 0.84400000 0.016761562  0.80830040 0.87265136
#> 6        mean_D_G2 0.15800000 0.012218882  0.13552361 0.16700611
#> 7           ATE_G1 1.14177258 0.027771683  1.06941747 1.17170515
#> 8           ATE_G2 0.90254558 0.041319786  0.83097975 0.97990939
#> 9           ATT_G1 1.21230186 0.028229266  1.13432462 1.23330626
#> 10          ATT_G2 0.69652631 0.039727746  0.68841183 0.81120702
#> 11     diff_mean_Y 1.28421519 0.038215434  1.21797808 1.34424089
#> 12    diff_mean_Y0 0.37108358 0.027131041  0.34559965 0.42606636
#> 13     diff_mean_D 0.68600000 0.020923182  0.65040566 0.72601890
#> 14        diff_ATE 0.23922700 0.049704942  0.13963678 0.29331293
#> 15         dff_ATT 0.51577555 0.041181987  0.34582055 0.49618313
#> 16           total 1.28421519 0.038215434  1.21797808 1.34424089
#> 17        baseline 0.37108358 0.027131041  0.34559965 0.42606636
#> 18      prevalence 0.61914627 0.032832889  0.56583509 0.66375856
#> 19          effect 0.20190759 0.043695607  0.11714127 0.25447088
#> 20       selection 0.09207776 0.009267851  0.06034112 0.09313772
#> 21 cond_prevalence 0.43522769 0.067125567  0.29041768 0.48077095
#> 22     cond_effect 0.22579569 0.116696517 -0.07777706 0.40401287
#> 23  cond_selection 0.08103174 0.008604817  0.04970048 0.07713673
#> 24          Q_dist 0.17107650 0.130527796 -0.06887207 0.33492071
```
