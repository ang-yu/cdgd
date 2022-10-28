

d=1
g1=0
g2=0
if (d==0 & g1==0) {
  IPO_arg <- IPO_D0G0
  YgivenX.Pred_arg <- YgivenX.Pred_D0G0}
if (d==1 & g1==0) {
  IPO_arg <- IPO_D1G0
  YgivenX.Pred_arg <- YgivenX.Pred_D1G0}
if (d==0 & g1==1) {
  IPO_arg <- IPO_D0G1
  YgivenX.Pred_arg <- YgivenX.Pred_D0G1}
if (d==1 & g1==1) {
  IPO_arg <- IPO_D1G1
  YgivenX.Pred_arg <- YgivenX.Pred_D1G1}
as.numeric(data[,G]==g1)/mean(data[,G]==g1)*IPO_arg*mean(data[,G]==g2/mean(data[,G]==g2)*data[,D]) +
  as.numeric(data[,G]==g2)/mean(data[,G]==g2)*mean(data[,G]==g1/mean(data[,G]==g1)*YgivenX.Pred_arg)*(data[,D]-mean(as.numeric(data[,G]==g2)/mean(data[,G]==g2)*data[,D]))









