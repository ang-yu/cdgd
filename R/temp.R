

d=1
g1=0
g2=1
g3=0

if (d==0 & g1==0) {
  IPO_arg <- IPO_D0G0
  YdgivenQ.Pred_arg <- Y0givenQ.Pred_G0
  g1givenQ.Pred_arg <- 1-GgivenQ.Pred}
if (d==1 & g1==0) {
  IPO_arg <- IPO_D1G0
  YdgivenQ.Pred_arg <- Y1givenQ.Pred_G0
  g1givenQ.Pred_arg <- 1-GgivenQ.Pred}
if (d==0 & g1==1) {
  IPO_arg <- IPO_D0G1
  YdgivenQ.Pred_arg <- Y0givenQ.Pred_G1
  g1givenQ.Pred_arg <- GgivenQ.Pred}
if (d==1 & g1==1) {
  IPO_arg <- IPO_D1G1
  YdgivenQ.Pred_arg <- Y1givenQ.Pred_G1
  g1givenQ.Pred_arg <- GgivenQ.Pred}

if (g2==0) {
  DgivenQ.Pred_arg <- DgivenQ.Pred_G0
  g2givenQ.Pred_arg <- 1-GgivenQ.Pred
}
if (g2==1) {
  DgivenQ.Pred_arg <- DgivenQ.Pred_G1
  g2givenQ.Pred_arg <- GgivenQ.Pred
}

if (g3==0) {
  g3givenQ.Pred_arg <- 1-GgivenQ.Pred
}
if (g3==1) {
  g3givenQ.Pred_arg <- GgivenQ.Pred
}

stab1 <- mean(as.numeric(data[,G]==g1)/mean(data[,G]==g3)*g3givenQ.Pred_arg/g1givenQ.Pred_arg)
stab2 <- mean(as.numeric(data[,G]==g2)/mean(data[,G]==g3)*g3givenQ.Pred_arg/g2givenQ.Pred_arg)

mean(as.numeric(data[,G]==g1)/mean(data[,G]==g3)*g3givenQ.Pred_arg/g1givenQ.Pred_arg/stab1*(IPO_arg)*DgivenQ.Pred_arg +
      as.numeric(data[,G]==g2)/mean(data[,G]==g3)*g3givenQ.Pred_arg/g2givenQ.Pred_arg/stab2*(data[,D]-DgivenQ.Pred_arg)*YdgivenQ.Pred_arg )

mean(as.numeric(data[,G]==g1)/mean(data[,G]==g3)*g3givenQ.Pred_arg/g1givenQ.Pred_arg/stab1*(IPO_arg)*DgivenQ.Pred_arg)
mean(as.numeric(data[,G]==g2)/mean(data[,G]==g3)*g3givenQ.Pred_arg/g2givenQ.Pred_arg/stab2*(data[,D]-DgivenQ.Pred_arg)*YdgivenQ.Pred_arg)


d=0
g1=0
g2=1
g3=0

0.2222894
0.1702883

0.07992304
0.1250708


# conditional (gbm):
0.1996398
0.08762225

# unconditional (gbm):
0.2267091
0.1008862

# cor.test(IPO_D1G0[data[G]==0]-IPO_D0G0[data[G]==0], data[data[G]==0,Q])


## the conditional prevalence is basically the same as the unconditional prevalence. This is because although \xi_{1bab}-\xi_{0bab}
# is unsurprisingly smaller than \xi_{1ba}-\xi_{0ba}, \xi_{1bbb}-\xi_{0bbb} is also smaller than \xi_{1bb}-\xi_{obb}. The latter is
# because (in group b) within Q, the positive selection is smaller than the unconditional selection. This is in turn because Q is
# negatively correlated with the treatment effect
# cor.test(IPO_D1G0[data[G]==0]-IPO_D0G0[data[G]==0], data[data[G]==0,Q])

0.3863603-0.2674638


mean(as.numeric(data[,G]==g1)/mean(data[,G]==g3)*g3givenQ.Pred_arg/g1givenQ.Pred_arg*(IPO_arg)*DgivenQ.Pred_arg +
        as.numeric(data[,G]==g2)/mean(data[,G]==g3)*g3givenQ.Pred_arg/g2givenQ.Pred_arg*(data[,D]-DgivenQ.Pred_arg)*YdgivenQ.Pred_arg )

(1-GgivenQ.Pred)/GgivenQ.Pred

plot(lowess(data[,Q],GgivenQ.Pred))
plot(data[,Q],GgivenQ.Pred)

psi_dgg(1,0,1)
[1] 0.2439536
> psi_dgg(0,0,1)
[1] 0.1160645
> psi_dgg(1,0,0)
[1] 0.02379137
> psi_dgg(0,0,0)
[1] 0.0113191


> psi_dggg(1,0,1,0)
[1] 2.542514
> psi_dggg(0,0,1,0)
[1] 2.029086
> psi_dggg(1,0,0,0)
[1] 0.02574895
> psi_dggg(0,0,0,0)
[1] -0.04823942



Y="adult_income_rank"
D="completion"
G="pincome_1"
X=c("AFQT","gender","medu","parental_presence",
    "n_sib","urban","edu_exp","age","friend_edu_exp","rotter_score","rosenberg_irt_score",
    "sig_other_exp2","sig_other_exp3","sig_other_exp4","foreign_lang",
    "SMSA2","SMSA3","SMSA4","mother_seperate",
    "school_satis2","school_satis3","fm_foreign_born",
    "region2","region3","region4","m_work","race2","race3")
data=data
algorithm="gbm"
alpha=0.05
set.seed(1)




Y="adult_income_rank"
D="completion"
G="pincome_1"
Q="AFQT"
X=c("AFQT","gender","medu","parental_presence",
    "n_sib","urban","edu_exp","age","friend_edu_exp","rotter_score","rosenberg_irt_score",
    "sig_other_exp2","sig_other_exp3","sig_other_exp4","foreign_lang",
    "SMSA2","SMSA3","SMSA4","mother_seperate",
    "school_satis2","school_satis3","fm_foreign_born",
    "region2","region3","region4","m_work","race2","race3")
Q="AFQT"
data=data
algorithm="gbm"
alpha=0.05
set.seed(1)



> mean(DgivenQ.Pred_G1-DgivenQ.Pred_G0)
[1] 0.1703779
> mean(data[data[,G]==1,D])-mean(data[data[,G]==0,D])
[1] 0.3771958


mean( ( (DgivenQ.Pred_G1-DgivenQ.Pred_G0)*(Y1givenQ.Pred_G0-Y0givenQ.Pred_G0) )[data[,G]==0] )




> cond_result_nnet
names               point                  se            CI_lower            CI_upper
1                  total    0.28558054882436   0.015065439248745   0.256052830485544   0.315108267163177
2               baseline   0.261498447987045  0.0220916356632679   0.218199637727459   0.304797258246631
3 conditional prevalence    0.15663724165777  0.0126521959055161   0.131839393357614   0.181435089957927
4     conditional effect  -0.259959754868767  0.0567473802044695  -0.371182576286528  -0.148736933451006
5  conditional selection -0.0333240218690801 0.00988756011173774 -0.0527032835830609 -0.0139447601550993
6         Q distribution   0.160728635917392  0.0469400019704429  0.0687279226210847   0.252729349213699

> cond_result_ranger
names               point                  se           CI_lower            CI_upper
1                  total    0.28558054882436   0.015065439248745  0.256052830485544   0.315108267163177
2               baseline   0.258649885458308  0.0202481504956085  0.218964239733368   0.298335531183247
3 conditional prevalence   0.126345931013473 0.00970935770593622  0.107315939596822   0.145375922430125
4     conditional effect  -0.112002632354042  0.0375085218146465 -0.185517984224084     -0.038487280484
5  conditional selection -0.0722358707884052  0.0182694242102866 -0.108043284258851 -0.0364284573179594
6         Q distribution  0.0848232354950263  0.0258108685006706 0.0342348628240127    0.13541160816604

> cond_result_gbm
names              point                  se            CI_lower               CI_upper
1                  total   0.28558054882436   0.015065439248745   0.256052830485544      0.315108267163177
2               baseline  0.256966503163197  0.0227143732752011   0.212447149612404       0.30148585671399
3 conditional prevalence  0.105884324456687 0.00779619528823791  0.0906040624752999      0.121164586438074
4     conditional effect -0.304992115535094  0.0735272312374309  -0.449102840643407     -0.160881390426781
5  conditional selection  -0.01459029253244 0.00743159659789168 -0.0291559542119381 -0.0000246308529419283
6         Q distribution   0.24231212927201  0.0621265869573617   0.120546256353185      0.364078002190835




