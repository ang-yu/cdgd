
# cdgd

<!-- badges: start -->
<!-- badges: end -->

The goal of cdgd is to implement the causal decomposition of group
disparities of Yu and Elwert (2022).

## Installation

You can install the development version of cdgd from
[GitHub](https://github.com/) with:

``` r
devtools::install_github("ang-yu/cdgd")
#> Downloading GitHub repo ang-yu/cdgd@HEAD
#> 
#> * checking for file ‘/private/var/folders/0q/j5fgjmvn167f9h6g2mrdv91c0000gn/T/RtmppvmStN/remotes178b179a94480/ang-yu-cdgd-27b425b/DESCRIPTION’ ... OK
#> * preparing ‘cdgd’:
#> * checking DESCRIPTION meta-information ... OK
#> * checking for LF line-endings in source and make files and shell scripts
#> * checking for empty or unneeded directories
#> * building ‘cdgd_0.0.0.9000.tar.gz’
#> Installing package into '/private/var/folders/0q/j5fgjmvn167f9h6g2mrdv91c0000gn/T/Rtmppqng7D/temp_libpatha6b06166a1cd'
#> (as 'lib' is unspecified)
```

## Example

``` r
library(cdgd)

# simulated example data
head(exp_data)
#>    outcome treatment confounder         Q group_a group_b
#> 1 1.318643         1   1.497808 0.5537125       1       0
#> 2 2.054619         1   2.131531 2.3158558       1       0
#> 3 1.803570         1   1.921083 1.6572525       1       0
#> 4 1.668810         1   2.886785 0.0686469       1       0
#> 5 2.030661         1   2.116971 2.2428210       1       0
#> 6 1.955015         1   2.318630 1.6372321       1       0
```

# Use cdgd0 to get point estimates

``` r
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

# Use cdgd to get point estimates and confidence intervals. Will take a minute or so.

``` r
set.seed(1)
results <- cdgd(Y="outcome",D="treatment",G1="group_a",G2="group_b",X=c("confounder","Q"),Q="Q",data=exp_data,t=0.05,algorithm="nnet",alpha=0.05,k=20)

results
#>               item      point          se       lower      upper
#> 1        mean_Y_G1 1.82626302 0.028601434  1.78001565 1.86099860
#> 2        mean_Y_G2 0.51487709 0.009368193  0.49267817 0.52674606
#> 3       mean_Y0_G1 0.80204650 0.036528592  0.76695783 0.86922163
#> 4       mean_Y0_G2 0.41296185 0.006388473  0.40108986 0.41989096
#> 5        mean_D_G1 0.86800000 0.014778835  0.84265010 0.88932806
#> 6        mean_D_G2 0.14600000 0.011798501  0.11627907 0.16194332
#> 7           ATE_G1 1.12179482 0.035175873  1.03549622 1.16847126
#> 8           ATE_G2 0.86344271 0.054719465  0.71362095 0.91581147
#> 9           ATT_G1 1.17997295 0.035738735  1.09401099 1.22648759
#> 10          ATT_G2 0.69804959 0.047109306  0.53556758 0.74062652
#> 11     diff_mean_Y 1.31138593 0.029575160  1.25205620 1.35656478
#> 12    diff_mean_Y0 0.38908465 0.034230003  0.35670965 0.45820574
#> 13     diff_mean_D 0.72200000 0.018343084  0.68801277 0.74851162
#> 14        diff_ATE 0.25835211 0.058365332  0.16195366 0.36909519
#> 15         dff_ATT 0.48192337 0.039894640  0.41861227 0.55704602
#> 16           total 1.31138593 0.029575160  1.25205620 1.35656478
#> 17        baseline 0.38908465 0.034230003  0.35670965 0.45820574
#> 18      prevalence 0.62340564 0.043007047  0.50534449 0.66606421
#> 19          effect 0.22424963 0.051233810  0.14033833 0.31892798
#> 20       selection 0.07464602 0.008519691  0.05276401 0.07995860
#> 21 cond_prevalence 0.44547350 0.085080809  0.20947097 0.54425993
#> 22     cond_effect 0.36366893 0.085220461  0.18388420 0.50795542
#> 23  cond_selection 0.06141283 0.009109305  0.04066190 0.07116794
#> 24          Q_dist 0.05174603 0.086608148 -0.14695387 0.17040513
```
