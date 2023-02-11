
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

results0$results
#>        names              point                 se           CI_lower
#> 1      total  0.267479354808872 0.0389768846815615  0.191086064603441
#> 2   baseline  0.039997497539402 0.0129327144177851 0.0146498430582013
#> 3 prevalence  0.256512539568984  0.033532390580103  0.190790261716452
#> 4     effect -0.136516278368955  0.020809641461047 -0.177302426163798
#> 5  selection  0.107485596069441 0.0139415961675306  0.080160569694079
#>              CI_upper
#> 1   0.343872645014304
#> 2  0.0653451520206026
#> 3   0.322234817421516
#> 4 -0.0957301305741109
#> 5   0.134810622444802
```

### Use cdgd1_ml, cdgd1_pa, or cdgd1_manual to get conditional decomposition

``` r
results1 <- cdgd1_pa(Y="outcome",D="treatment",G="group_a",X="confounder",Q="Q",data=exp_data,alpha=0.05)

results1
#>                           names              point                 se
#> 1                         total  0.267479354808872 0.0389768846815615
#> 2                      baseline  0.039997497539402 0.0129327144177851
#> 3        conditional prevalence  0.208604685983608 0.0351067274912385
#> 4            conditional effect -0.152921251811951 0.0204569657882307
#> 5         conditional selection 0.0887548333345764 0.0134667561793527
#> 6                Q distribution 0.0830435897632368 0.0109725233258662
#> 7 conditional Jackson reduction  0.239315359423403 0.0355208384410009
#>             CI_lower           CI_upper
#> 1  0.191086064603441  0.343872645014304
#> 2 0.0146498430582014 0.0653451520206026
#> 3  0.139796764485719  0.277412607481498
#> 4 -0.193016167989852 -0.112826335634051
#> 5 0.0623604762344629   0.11514919043469
#> 6 0.0615378392250134   0.10454934030146
#> 7  0.169695795378375   0.30893492346843
```
