
# cdgd

<!-- badges: start -->
<!-- badges: end -->

The goal of cdgd is to implement the causal decomposition of group
disparities of Yu and Elwert (2022).

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
#> 3 prevalence  0.256512539568985 0.0327325145064224  0.192357990012962
#> 4     effect -0.136516278368955 0.0207371744191902 -0.177160393371693
#> 5  selection  0.107485596069441 0.0138872229169067 0.0802671393070242
#>             CI_upper
#> 1  0.343872645014304
#> 2 0.0653451520206026
#> 3  0.320667089125008
#> 4 -0.095872163366217
#> 5  0.134704052831857
```

### Use cdgd1_ml, cdgd1_pa, or cdgd1_manual to get conditional decomposition

``` r
results1 <- cdgd1_pa(Y="outcome",D="treatment",G="group_a",X="confounder",Q="Q",data=exp_data,alpha=0.05)

results1
#>                           names              point                 se
#> 1                         total  0.267479354808872 0.0389768846815615
#> 2                      baseline  0.039997497539402 0.0129327144177851
#> 3        conditional prevalence  0.209003240763591 0.0338235692692773
#> 4            conditional effect 0.0661663537818592 0.0778173202907227
#> 5         conditional selection 0.0887548333346887 0.0588644397620432
#> 6                Q distribution -0.136442570610669 0.0726359768772961
#> 7 conditional Jackson reduction  0.239713914203476 0.0352865635907051
#>              CI_lower            CI_upper
#> 1   0.191086064603441   0.343872645014304
#> 2  0.0146498430582014  0.0653451520206026
#> 3   0.142710263167212    0.27529621835997
#> 4 -0.0863527913613753   0.218685498925094
#> 5 -0.0266173485690436   0.204127015238421
#> 6  -0.278806469272053 0.00592132805071557
#> 7   0.170553520427511    0.30887430797944
```
