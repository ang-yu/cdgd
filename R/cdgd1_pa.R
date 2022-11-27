


#' Perform conditional decomposition via parametric models
#'
#' @param Y Outcome. The name of a continuous variable.
#' @param D Treatment status. The name of a binary numeric variable taking values of 0 and 1.
#' @param G Advantaged group membership. The name of a binary numeric variable taking values of 0 and 1.
#' @param Q Conditional set. The vector of the names of numeric variables.
#' @param X Confounders. The vector of the names of numeric variables.
#' @param data A data frame.
#' @param algorithm The ML algorithm for modelling. "nnet" for neural network, "ranger" for random forests, "gbm" for generalized boosted models.
#' @param alpha 1-alpha confidence interval.
#'
#' @return A data frame of point estimates.
#'
#' @export
#'
#' @examples
#' data(exp_data)
#'
#' results0 <- cdgd0(
#' Y="outcome",
#' D="treatment",
#' G="group_a",
#' #' X=c("confounder","Q"),
#' Q="Q",
#' data=exp_data,
#' t=0.05,
#' algorithm="nnet")
#'
#' results0



cdgd1_pa <- function(Y,D,G,X,Q,data,algorithm,alpha=0.05) {

  data <- as.data.frame(data)

  ### outcome regression model
  YgivenDGXQ.Model <- lm(stats::as.formula(paste(Y, paste(paste(D,c(G,Q,X),sep="*"),collapse="+"), sep="~")), data=data)

  ### propensity score model
  DgivenGXQ.Model <- glm(stats::as.formula(paste(D, paste(G,Q,paste(X,collapse="+"),sep="+"), sep="~")), data=data, family=binomial(link="logit"))

  ### predictions
  YgivenXQ.Pred_D0G0 <- YgivenXQ.Pred_D1G0 <- YgivenXQ.Pred_D0G1 <- YgivenXQ.Pred_D1G1 <- DgivenXQ.Pred_G0 <- DgivenXQ.Pred_G1 <- rep(NA, nrow(data))

  pred_data <- data
  pred_data[,D] <- 0
  pred_data[,G] <- 0
  YgivenXQ.Pred_D0G0 <- stats::predict(YgivenDGX.Model, newdata = pred_data)

  pred_data <- data
  pred_data[,D] <- 1
  pred_data[,G] <- 0
  YgivenXQ.Pred_D1G0 <- stats::predict(YgivenDGX.Model, newdata = pred_data)

  pred_data <- data
  pred_data[,D] <- 0
  pred_data[,G] <- 1
  YgivenXQ.Pred_D0G1 <- stats::predict(YgivenDGX.Model, newdata = pred_data)

  pred_data <- data
  pred_data[,D] <- 1
  pred_data[,G] <- 1
  YgivenXQ.Pred_D1G1 <- stats::predict(YgivenDGX.Model, newdata = pred_data)

  pred_data <- data
  pred_data[,G] <- 0
  DgivenXQ.Pred_G0 <- stats::predict(DgivenGXQ.Model, newdata = pred_data, type="response")

  pred_data <- data
  pred_data[,G] <- 1
  DgivenXQ.Pred_G1 <- stats::predict(DgivenGXQ.Model, newdata = pred_data, type="response")

  ### Estimate E(Y_d | Q,g)
  data_temp <- data[,c(G,Q)]

  pred_data <- data
  pred_data[,D] <- 1
  data_temp$YgivenGXQ.Pred_D0_ncf <- stats::predict(YgivenDGXQ.Model, newdata = pred_data)

  pred_data <- data
  pred_data[,D] <- 0
  data_temp$YgivenGXQ.Pred_D1_ncf <- stats::predict(YgivenDGXQ.Model, newdata = pred_data)

  Y0givenGQ.Model <- lm(stats::as.formula(paste("YgivenGXQ.Pred_D0_ncf", paste(G,Q,sep="*"), sep="~")), data=data_temp)
  Y1givenGQ.Model <- lm(stats::as.formula(paste("YgivenGXQ.Pred_D1_ncf", paste(G,Q,sep="*"), sep="~")), data=data_temp)

  Y0givenQ.Pred_G0 <- Y0givenQ.Pred_G1 <- Y1givenQ.Pred_G0 <- Y1givenQ.Pred_G1 <- rep(NA, nrow(data))

  pred_data <- data
  pred_data[,G] <- 1
  Y0givenQ.Pred_G1 <- stats::predict(Y0givenGQ.Model, newdata = pred_data)
  Y1givenQ.Pred_G1 <- stats::predict(Y1givenGQ.Model, newdata = pred_data)

  pred_data <- data
  pred_data[,G] <- 0
  Y0givenQ.Pred_G0 <- stats::predict(Y0givenGQ.Model, newdata = pred_data)
  Y1givenQ.Pred_G0 <- stats::predict(Y1givenGQ.Model, newdata = pred_data)

  ### The "IPO" (individual potential outcome) function
  # For each d and g value, we have IE(d,g)=\frac{\one(D=d)}{\pi(d,X,g)}[Y-\mu(d,X,g)]+\mu(d,X,g)
  # We stablize the weight by dividing the sample average of estimated weights

  IPO_D0G0 <- (1-data[,D])/(1-DgivenXQ.Pred_G0)/mean((1-data[,D])/(1-DgivenXQ.Pred_G0))*(data[,Y]-YgivenXQ.Pred_D0G0) + YgivenXQ.Pred_D0G0
  IPO_D1G0 <- data[,D]/DgivenXQ.Pred_G0/(mean(data[,D]/DgivenXQ.Pred_G0))*(data[,Y]-YgivenXQ.Pred_D1G0) + YgivenXQ.Pred_D1G0
  IPO_D0G1 <- (1-data[,D])/(1-DgivenXQ.Pred_G0)/(mean((1-data[,D])/(1-DgivenXQ.Pred_G0)))*(data[,Y]-YgivenXQ.Pred_D0G1) + YgivenXQ.Pred_D0G1
  IPO_D1G1 <- data[,D]/DgivenXQ.Pred_G0/mean(data[,D]/DgivenXQ.Pred_G0)*(data[,Y]-YgivenXQ.Pred_D1G1) + YgivenXQ.Pred_D1G1

  ### Estimate E(D | Q,g')
  DgivenGQ.Model <- glm(stats::as.formula(paste(D, paste(G,Q,sep="*"), sep="~")), data=data, family=binomial(link="logit"))

  DgivenQ.Pred_G0 <- DgivenQ.Pred_G1 <- rep(NA, nrow(data))

  pred_data <- data
  pred_data[,G] <- 0
  DgivenQ.Pred_G0 <- stats::predict(DgivenGQ.Model, newdata = pred_data, type="response")

  pred_data <- data
  pred_data[,G] <- 1
  DgivenQ.Pred_G1 <- stats::predict(DgivenGQ.Model, newdata = pred_data, type="response")

  ### Estimate p_g(Q)=Pr(G=g | Q)
  GgivenQ.Model <- glm(stats::as.formula(paste(G, paste(Q,sep="+"), sep="~")), data=data, family=binomial(link="logit"))

  GgivenQ.Pred <- rep(NA, nrow(data))
  GgivenQ.Pred <- stats::predict(DgivenGQ.Model, newdata = data, type="response")

  ### The one-step estimate of \xi_{dg}
  psi_00 <- mean( (1-data[,G])/(1-mean(data[,G]))*IPO_D0G0 )
  psi_01 <- mean( data[,G]/mean(data[,G])*IPO_D0G1 )

  ### The one-step estimate of \xi_{dgg'g''}
  # There are 8 possible dgg'g'' combinations, so we define a function first
  psi_dggg <- function(d,g1,g2,g3) {
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

    # denominators for weight stabilization using the fact that E( \frac{\one(G=g)p_{g''}(Q)}{p_g(Q)p_{g''}} ) and E( \frac{\one(G=g')p_{g''}(Q)}{p_{g'}(Q)p_{g''}} ) are both 1.
    stab1 <- mean(as.numeric(data[,G]==g1)/mean(data[,G]==g3)*g3givenQ.Pred_arg/g1givenQ.Pred_arg)
    stab2 <- mean(as.numeric(data[,G]==g2)/mean(data[,G]==g3)*g3givenQ.Pred_arg/g2givenQ.Pred_arg)

    psi_dggg <- mean( as.numeric(data[,G]==g3)/mean(data[,G]==g3)*YdgivenQ.Pred_arg*DgivenQ.Pred_arg +
                        as.numeric(data[,G]==g1)/mean(data[,G]==g3)*g3givenQ.Pred_arg/g1givenQ.Pred_arg/stab1*(IPO_arg-YdgivenQ.Pred_arg)*DgivenQ.Pred_arg +
                        as.numeric(data[,G]==g2)/mean(data[,G]==g3)*g3givenQ.Pred_arg/g2givenQ.Pred_arg/stab2*(data[,D]-DgivenQ.Pred_arg)*YdgivenQ.Pred_arg )

    return(psi_dggg)
  }

  ### point estimates
  Y_G0 <- mean((1-data[,G])/(1-mean(data[,G]))*data[,Y])       # mean outcome estimate for group 0
  Y_G1 <- mean(data[,G]/mean(data[,G])*data[,Y])               # mean outcome estimate for group 1
  total <- Y_G1-Y_G0

  baseline <- psi_01-psi_00
  cond_prevalence <- psi_dggg(1,0,1,0)-psi_dggg(0,0,1,0)-psi_dggg(1,0,0,0)+psi_dggg(0,0,0,0)
  cond_effect <- psi_dggg(1,1,1,1)-psi_dggg(0,1,1,1)-psi_dggg(1,0,1,1)+psi_dggg(0,0,1,1)
  Q_dist <- psi_dggg(1,0,1,1)-psi_dggg(0,0,1,1)-psi_dggg(1,0,1,0)+psi_dggg(0,0,1,0)
  cond_selection <- total-baseline-cond_prevalence-cond_effect-Q_dist

  ### standard error estimates
  se <- function(x) {sqrt( mean(x^2)/nrow(data) )}
  total_se <- se( data[,G]/mean(data[,G])*(data[,Y]-Y_G1) - (1-data[,G])/(1-mean(data[,G]))*(data[,Y]-Y_G0) )
  baseline_se <- se( data[,G]/mean(data[,G])*(IPO_D0G1-psi_01) - (1-data[,G])/(1-mean(data[,G]))*(IPO_D0G0-psi_00) )

  EIF_dggg <- function(d,g1,g2,g3) {
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

    # denominators for weight stabilization using the fact that E( \frac{\one(G=g)p_{g''}(Q)}{p_g(Q)p_{g''}} ) and E( \frac{\one(G=g')p_{g''}(Q)}{p_{g'}(Q)p_{g''}} ) are both 1.
    stab1 <- mean(as.numeric(data[,G]==g1)/mean(data[,G]==g3)*g3givenQ.Pred_arg/g1givenQ.Pred_arg)
    stab2 <- mean(as.numeric(data[,G]==g2)/mean(data[,G]==g3)*g3givenQ.Pred_arg/g2givenQ.Pred_arg)

    return(
      as.numeric(data[,G]==g3)/mean(data[,G]==g3)*(YdgivenQ.Pred_arg*DgivenQ.Pred_arg-psi_dggg(d,g1,g2,g3)) +
        as.numeric(data[,G]==g1)/mean(data[,G]==g3)*g3givenQ.Pred_arg/g1givenQ.Pred_arg/stab1*(IPO_arg-YdgivenQ.Pred_arg)*DgivenQ.Pred_arg +
        as.numeric(data[,G]==g2)/mean(data[,G]==g3)*g3givenQ.Pred_arg/g2givenQ.Pred_arg/stab2*(data[,D]-DgivenQ.Pred_arg)*YdgivenQ.Pred_arg
    )
  }

  cond_prevalence_se <- se( EIF_dggg(1,0,1,0)-EIF_dggg(0,0,1,0)-EIF_dggg(1,0,0,0)+EIF_dggg(0,0,0,0) )
  cond_effect_se <- se( EIF_dggg(1,1,1,1)-EIF_dggg(0,1,1,1)-EIF_dggg(1,0,1,1)+EIF_dggg(0,0,1,1) )
  Q_dist_se <- se( EIF_dggg(1,0,1,1)-EIF_dggg(0,0,1,1)-EIF_dggg(1,0,1,0)+EIF_dggg(0,0,1,0) )
  cond_selection_se <- se( data[,G]/mean(data[,G])*(data[,Y]-Y_G1) - (1-data[,G])/(1-mean(data[,G]))*(data[,Y]-Y_G0) -
                             ( data[,G]/mean(data[,G])*(IPO_D0G1-psi_01) - (1-data[,G])/(1-mean(data[,G]))*(IPO_D0G0-psi_00) ) -
                             ( EIF_dggg(1,0,1,0)-EIF_dggg(0,0,1,0)-EIF_dggg(1,0,0,0)+EIF_dggg(0,0,0,0) ) -
                             ( EIF_dggg(1,1,1,1)-EIF_dggg(0,1,1,1)-EIF_dggg(1,0,1,1)+EIF_dggg(0,0,1,1) ) -
                             ( EIF_dggg(1,0,1,1)-EIF_dggg(0,0,1,1)-EIF_dggg(1,0,1,0)+EIF_dggg(0,0,1,0) ))

  ### output results
  point <- c(total,
             baseline,
             cond_prevalence,
             cond_effect,
             cond_selection,
             Q_dist)
  se <- c(total_se,
          baseline_se,
          cond_prevalence_se,
          cond_effect_se,
          cond_selection_se,
          Q_dist_se)
  CI_lower <- point - qnorm(1-alpha/2)*se
  CI_upper <- point + qnorm(1-alpha/2)*se
  names <- c("total",
             "baseline",
             "conditional prevalence",
             "conditional effect",
             "conditional selection",
             "Q distribution")

  output <- as.data.frame(cbind(names, point,se,CI_lower,CI_upper))

  return(output)
}
