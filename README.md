
# cdgd

<!-- badges: start -->
<!-- badges: end -->

The goal of cdgd is to implement the causal decomposition of group
disparities of Yu and Elwert (2022).

## Installation

You can install the development version of cdgd from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ang-yu/cdgd")
#> Downloading GitHub repo ang-yu/cdgd@HEAD
#> 
#> * checking for file ‘/private/var/folders/0q/j5fgjmvn167f9h6g2mrdv91c0000gn/T/Rtmp5k0MA5/remotes16bc7652b5cdb/ang-yu-cdgd-59f8f80/DESCRIPTION’ ... OK
#> * preparing ‘cdgd’:
#> * checking DESCRIPTION meta-information ... OK
#> * checking for LF line-endings in source and make files and shell scripts
#> * checking for empty or unneeded directories
#> * building ‘cdgd_0.0.0.9000.tar.gz’
#> Installing package into '/private/var/folders/0q/j5fgjmvn167f9h6g2mrdv91c0000gn/T/Rtmppqng7D/temp_libpatha6b06cadd240'
#> (as 'lib' is unspecified)
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(cdgd)

head(exp_data)
#>    outcome treatment confounder         Q group_a group_b
#> 1 1.318643         1   1.497808 0.5537125       1       0
#> 2 2.054619         1   2.131531 2.3158558       1       0
#> 3 1.803570         1   1.921083 1.6572525       1       0
#> 4 1.668810         1   2.886785 0.0686469       1       0
#> 5 2.030661         1   2.116971 2.2428210       1       0
#> 6 1.955015         1   2.318630 1.6372321       1       0

set.seed(1)
results0 <- cdgd0(Y="outcome",D="treatment",G1="group_a",G2="group_b",X=c("confounder","Q"),Q="Q",data=exp_data,t=0.05,algorithm="nnet")

results0
#>               item      point
#> 1        mean_Y_G1 1.82626302
#> 2        mean_Y_G2 0.51487709
#> 3       mean_Y0_G1 0.80204650
#> 4       mean_Y0_G2 0.41296185
#> 5        mean_D_G1 0.86800000
#> 6        mean_D_G2 0.14600000
#> 7           ATE_G1 1.12179482
#> 8           ATE_G2 0.86344271
#> 9           ATT_G1 1.17997295
#> 10          ATT_G2 0.69804959
#> 11     diff_mean_Y 1.31138593
#> 12    diff_mean_Y0 0.38908465
#> 13     diff_mean_D 0.72200000
#> 14        diff_ATE 0.25835211
#> 15         dff_ATT 0.48192337
#> 16           total 1.31138593
#> 17        baseline 0.38908465
#> 18      prevalence 0.62340564
#> 19          effect 0.22424963
#> 20       selection 0.07464602
#> 21 cond_prevalence 0.44547350
#> 22     cond_effect 0.36366893
#> 23  cond_selection 0.06141283
#> 24          Q_dist 0.05174603
```
