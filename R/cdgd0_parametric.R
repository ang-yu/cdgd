


#' Perform unconditional decomposition
#'
#' @param Y Outcome. The name of a continuous variable.
#' @param D Treatment status. The name of a binary numeric variable taking values of 0 and 1.
#' @param G1 Group 1 membership. The name of a binary factor variable taking values of 0 and 1.
#' @param G2 Group 2 membership. The name of a binary factor variable taking values of 0 and 1.
#' @param X Confounders. The vector of the names of numeric variables.
#' @param data A data frame.
#' @param alpha 1-alpha confidence interval.
#'
#' @return A data frame of point estimates.
#'
#' @export
#'
#' @examples
#' data(exp_data)
#'
#' set.seed(1)
#'
#' results0 <- cdgd0(
#' Y="outcome",
#' D="treatment",
#' G1="group_a",
#' G2="group_b",
#' X=c("confounder"),
#' data=exp_data,
#' algorithm="nnet")
#'
#' results0




cdgd0_parametric <- function(Y,D,G,X,data,algorithm,alpha=0.05) {

  data <- as.data.frame(data)

  ### outcome regression model
  YgivenDGX.Model <- lm(stats::as.formula(paste(Y, paste(D,G,paste(X,collapse="+"),sep="+"), sep="~")), data=data)

  ### propensity score model
  DgivenGX.Model <- glm(stats::as.formula(paste(D, paste(G,paste(X,collapse="+"),sep="+"), sep="~")), data=data, family=binomial(link="logit"))

  ### predictions
  YgivenX.Pred_D0G0 <- YgivenX.Pred_D1G0 <- YgivenX.Pred_D0G1 <- YgivenX.Pred_D1G1 <- DgivenX.Pred_G0 <- DgivenX.Pred_G1 <- rep(NA, nrow(data))

  pred_data <- data
  pred_data[,D] <- 0
  pred_data[,G] <- 0
  YgivenX.Pred_D0G0 <- stats::predict(YgivenDGX.Model, newdata = pred_data)

  pred_data <- data
  pred_data[,D] <- 1
  pred_data[,G] <- 0
  YgivenX.Pred_D1G0 <- stats::predict(YgivenDGX.Model, newdata = pred_data)

  pred_data <- data
  pred_data[,D] <- 0
  pred_data[,G] <- 1
  YgivenX.Pred_D0G1 <- stats::predict(YgivenDGX.Model, newdata = pred_data)

  pred_data <- data
  pred_data[,D] <- 1
  pred_data[,G] <- 1
  YgivenX.Pred_D1G1 <- stats::predict(YgivenDGX.Model, newdata = pred_data)

  pred_data <- data
  pred_data[,G] <- 0
  DgivenX.Pred_G0 <- stats::predict(DgivenGX.Model, newdata = pred_data, type="response")

  pred_data <- data
  pred_data[,G] <- 1
  DgivenX.Pred_G1 <- stats::predict(DgivenGX.Model, newdata = pred_data, type="response")

  ### The "IPO" (individual potential outcome) function
  # For each d and g value, we have IE(d,g)=\frac{\one(D=d)}{\pi(d,X,g)}[Y-\mu(d,X,g)]+\mu(d,X,g)
  # We stablize the weight by dividing the sample average of estimated weights

  IPO_D0G0 <- (1-data[,D])/(1-DgivenX.Pred_G0)/mean((1-data[,D])/(1-DgivenX.Pred_G0))*(data[,Y]-YgivenX.Pred_D0G0) + YgivenX.Pred_D0G0
  IPO_D1G0 <- data[,D]/DgivenX.Pred_G0/(mean(data[,D]/DgivenX.Pred_G0))*(data[,Y]-YgivenX.Pred_D1G0) + YgivenX.Pred_D1G0
  IPO_D0G1 <- (1-data[,D])/(1-DgivenX.Pred_G0)/(mean((1-data[,D])/(1-DgivenX.Pred_G0)))*(data[,Y]-YgivenX.Pred_D0G1) + YgivenX.Pred_D0G1
  IPO_D1G1 <- data[,D]/DgivenX.Pred_G0/mean(data[,D]/DgivenX.Pred_G0)*(data[,Y]-YgivenX.Pred_D1G1) + YgivenX.Pred_D1G1

  ### The one-step estimate of \xi_{dg} and \xi_{dgg'}
  psi_00 <- mean( (1-data[,G])/(1-mean(data[,G]))*IPO_D0G0 )
  psi_01 <- mean( data[,G]/mean(data[,G])*IPO_D0G1 )

  # There are 8 dgg' combinations, so we define a function first
  psi_dgg <- function(d,g1,g2) {
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

    psi_dgg <- mean( as.numeric(data[,G]==g1)/mean(data[,G]==g1)*IPO_arg*mean(as.numeric(data[,G]==g2)/mean(data[,G]==g2)*data[,D]) +
                       as.numeric(data[,G]==g2)/mean(data[,G]==g2)*mean(as.numeric(data[,G]==g1)/mean(data[,G]==g1)*YgivenX.Pred_arg)*(data[,D]-mean(as.numeric(data[,G]==g2)/mean(data[,G]==g2)*data[,D])) )

    return(psi_dgg)
  }

  ### point estimates
  Y_G0 <- mean((1-data[,G])/(1-mean(data[,G]))*data[,Y])       # mean outcome estimate for group 0
  Y_G1 <- mean(data[,G]/mean(data[,G])*data[,Y])               # mean outcome estimate for group 1
  total <- Y_G1-Y_G0

  baseline <- psi_01-psi_00
  prevalence <- psi_dgg(1,0,1)-psi_dgg(1,0,0)-psi_dgg(0,0,1)+psi_dgg(0,0,0)
  effect <- psi_dgg(1,1,1)-psi_dgg(0,1,1)-psi_dgg(1,0,1)+psi_dgg(0,0,1)
  selection <- total-baseline-prevalence-effect

  ### standard error estimates
  se <- function(x) {sqrt( mean(x^2)/nrow(data) )}
  total_se <- se( data[,G]/mean(data[,G])*(data[,Y]-Y_G1) - (1-data[,G])/(1-mean(data[,G]))*(data[,Y]-Y_G0) )
  baseline_se <- se( data[,G]/mean(data[,G])*(IPO_D0G1-psi_01) - (1-data[,G])/(1-mean(data[,G]))*(IPO_D0G0-psi_00) )

  EIF_dgg <- function(d,g1,g2) {
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

    return(
      as.numeric(data[,G]==g1)/mean(data[,G]==g1)*IPO_arg*mean(as.numeric(data[,G]==g2)/mean(data[,G]==g2)*data[,D]) +
        as.numeric(data[,G]==g2)/mean(data[,G]==g2)*mean(as.numeric(data[,G]==g1)/mean(data[,G]==g1)*YgivenX.Pred_arg)*(data[,D]-mean(as.numeric(data[,G]==g2)/mean(data[,G]==g2)*data[,D])) -
        as.numeric(data[,G]==g1)/mean(data[,G]==g1)*psi_dgg(d,g1,g2)
    )
  }

  prevalence_se <- se( EIF_dgg(1,0,1)-EIF_dgg(1,0,0)-EIF_dgg(0,0,1)+EIF_dgg(0,0,0) )
  effect_se <- se( EIF_dgg(1,1,1)-EIF_dgg(0,1,1)-EIF_dgg(1,0,1)+EIF_dgg(0,0,1) )
  selection_se <- se( data[,G]/mean(data[,G])*(data[,Y]-Y_G1) - (1-data[,G])/(1-mean(data[,G]))*(data[,Y]-Y_G0) -
                        ( data[,G]/mean(data[,G])*(IPO_D0G1-psi_01) - (1-data[,G])/(1-mean(data[,G]))*(IPO_D0G0-psi_00) ) -
                        ( EIF_dgg(1,0,1)-EIF_dgg(1,0,0)-EIF_dgg(0,0,1)+EIF_dgg(0,0,0) ) -
                        ( EIF_dgg(1,1,1)-EIF_dgg(0,1,1)-EIF_dgg(1,0,1)+EIF_dgg(0,0,1) ) )

  ### output results
  point <- c(total,
             baseline,
             prevalence,
             effect,
             selection)
  se <- c(total_se,
          baseline_se,
          prevalence_se,
          effect_se,
          selection_se)
  CI_lower <- point - qnorm(1-alpha/2)*se
  CI_upper <- point + qnorm(1-alpha/2)*se
  names <- c("total",
             "baseline",
             "prevalence",
             "effect",
             "selection")

  output <- as.data.frame(cbind(names, point,se,CI_lower,CI_upper))

  return(output)
}

