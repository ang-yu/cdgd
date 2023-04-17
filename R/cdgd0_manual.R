
#' Perform unconditional decomposition with nuisance functions estimated beforehand
#'
#' This function gives the user full control over the estimation of the nuisance functions.
#' For the unconditional decomposition, three nuisance functions (YgivenGX.Pred_D0, YgivenGX.Pred_D1, and DgivenGX.Pred) need to be estimated.
#' The nuisance functions should be estimated using cross-fitting if Donsker class is not assumed.
#'
#' @param Y Outcome. The name of a continuous variable.
#' @param D Treatment status. The name of a binary numeric variable taking values of 0 and 1.
#' @param G Advantaged group membership. The name of a binary numeric variable taking values of 0 and 1.
#' @param YgivenGX.Pred_D0 A numeric vector of predicted Y values given X, G, and D=0. Vector length=nrow(data).
#' @param YgivenGX.Pred_D1 A numeric vector of predicted Y values given X, G, and D=1. Vector length=nrow(data).
#' @param DgivenGX.Pred A numeric vector of predicted D values given X and G. Vector length=nrow(data).
#' @param data A data frame.
#' @param alpha 1-alpha confidence interval.
#'
#' @return A list of estimates.
#'
#' @export
#'
#' @examples
#' # This example will take a minute to run.
#' \donttest{
#' data(exp_data)
#'
#' Y="outcome"
#' D="treatment"
#' G="group_a"
#' X=c("Q","confounder")
#' data=exp_data
#'
#' set.seed(1)
#'
#' ### estimate the nuisance functions with cross-fitting
#' sample1 <- sample(nrow(data), floor(nrow(data)/2), replace=FALSE)
#' sample2 <- setdiff(1:nrow(data), sample1)
#'
#' ### outcome regression model
#'
#' message <- utils::capture.output( YgivenDGX.Model.sample1 <-
#'     caret::train(stats::as.formula(paste(Y, paste(D,G,paste(X,collapse="+"),sep="+"), sep="~")),
#'              data=data[sample1,], method="ranger", trControl=caret::trainControl(method="cv"),
#'              tuneGrid=expand.grid(mtry=c(2,4),splitrule=c("variance"),min.node.size=c(50,100))) )
#' message <- utils::capture.output( YgivenDGX.Model.sample2 <-
#'     caret::train(stats::as.formula(paste(Y, paste(D,G,paste(X,collapse="+"),sep="+"), sep="~")),
#'              data=data[sample2,], method="ranger", trControl=caret::trainControl(method="cv"),
#'              tuneGrid=expand.grid(mtry=c(2,4),splitrule=c("variance"),min.node.size=c(50,100))) )
#'
#' ### propensity score model
#' data[,D] <- as.factor(data[,D])
#' levels(data[,D]) <- c("D0","D1")  # necessary for caret implementation of ranger
#'
#' message <- utils::capture.output( DgivenGX.Model.sample1 <-
#'     caret::train(stats::as.formula(paste(D, paste(G,paste(X,collapse="+"),sep="+"), sep="~")),
#'              data=data[sample1,], method="ranger",
#'              trControl=caret::trainControl(method="cv", classProbs=TRUE),
#'              tuneGrid=expand.grid(mtry=c(1,2),splitrule=c("gini"),min.node.size=c(50,100))) )
#' message <- utils::capture.output( DgivenGX.Model.sample2 <-
#'     caret::train(stats::as.formula(paste(D, paste(G,paste(X,collapse="+"),sep="+"), sep="~")),
#'              data=data[sample2,], method="ranger",
#'              trControl=caret::trainControl(method="cv", classProbs=TRUE),
#'              tuneGrid=expand.grid(mtry=c(1,2),splitrule=c("gini"),min.node.size=c(50,100))) )
#'
#' data[,D] <- as.numeric(data[,D])-1
#'
#' ### cross-fitted predictions
#' YgivenGX.Pred_D0 <- YgivenGX.Pred_D1 <- DgivenGX.Pred <- rep(NA, nrow(data))
#'
#' pred_data <- data
#' pred_data[,D] <- 0
#' YgivenGX.Pred_D0[sample2] <- stats::predict(YgivenDGX.Model.sample1, newdata = pred_data[sample2,])
#' YgivenGX.Pred_D0[sample1] <- stats::predict(YgivenDGX.Model.sample2, newdata = pred_data[sample1,])
#'
#' pred_data <- data
#' pred_data[,D] <- 1
#' YgivenGX.Pred_D1[sample2] <- stats::predict(YgivenDGX.Model.sample1, newdata = pred_data[sample2,])
#' YgivenGX.Pred_D1[sample1] <- stats::predict(YgivenDGX.Model.sample2, newdata = pred_data[sample1,])
#'
#' pred_data <- data
#' DgivenGX.Pred[sample2] <- stats::predict(DgivenGX.Model.sample1,
#'     newdata = pred_data[sample2,], type="prob")[,2]
#' DgivenGX.Pred[sample1] <- stats::predict(DgivenGX.Model.sample2,
#'     newdata = pred_data[sample1,], type="prob")[,2]
#'
#' results <- cdgd0_manual(Y=Y,D=D,G=G,
#'                        YgivenGX.Pred_D0=YgivenGX.Pred_D0,
#'                        YgivenGX.Pred_D1=YgivenGX.Pred_D1,
#'                        DgivenGX.Pred=DgivenGX.Pred,
#'                        data=data)
#'
#' results}



cdgd0_manual <- function(Y,D,G,YgivenGX.Pred_D1,YgivenGX.Pred_D0,DgivenGX.Pred,data,alpha=0.05) {

  data <- as.data.frame(data)

  zero_one <- sum(DgivenGX.Pred==0)+sum(DgivenGX.Pred==1)
  if ( zero_one>0 ) {
    stop(
      paste("D given G and X and are exact 0 or 1 in", zero_one, "cases.", sep=" "),
      call. = FALSE
    )
  }

  ### The "IPO" (individual potential outcome) function
  # For each d and g value, we have IE(d,g)=\frac{\one(D=d)}{\pi(d,X,g)}[Y-\mu(d,X,g)]+\mu(d,X,g)
  # We stabilize the weight by dividing the sample average of estimated weights

  IPO_D0 <- (1-data[,D])/(1-DgivenGX.Pred)/mean((1-data[,D])/(1-DgivenGX.Pred))*(data[,Y]-YgivenGX.Pred_D0) + YgivenGX.Pred_D0
  IPO_D1 <- data[,D]/DgivenGX.Pred/mean(data[,D]/DgivenGX.Pred)*(data[,Y]-YgivenGX.Pred_D1) + YgivenGX.Pred_D1

  ### The one-step estimate of \xi_{dg} and \xi_{dgg'}
  psi_00 <- mean( (1-data[,G])/(1-mean(data[,G]))*IPO_D0 )
  psi_01 <- mean( data[,G]/mean(data[,G])*IPO_D0 )
  psi_10 <- mean( (1-data[,G])/(1-mean(data[,G]))*IPO_D1 )
  psi_11 <- mean( data[,G]/mean(data[,G])*IPO_D1 )

  # Note that this is basically DML2. We could also use DML1:
  #psi_00_S1 <- mean( (1-data[sample1,G])/(1-mean(data[sample1,G]))*IPO_D0[sample1] )     # sample 1 estimate
  #psi_00_S2 <- mean( (1-data[sample2,G])/(1-mean(data[sample2,G]))*IPO_D0[sample2] )     # sample 2 estimate
  #psi_00 <- (1/2)*(psi_00_S1+psi_00_S2)
  #psi_01_S1 <- mean( data[sample1,G]/mean(data[sample1,G])*IPO_D0[sample1] )     # sample 1 estimate
  #psi_01_S2 <- mean( data[sample1,G]/mean(data[sample1,G])*IPO_D0[sample2] )     # sample 2 estimate
  #psi_01 <- (1/2)*(psi_01_S1+psi_01_S2)

  # There are 8 dgg' combinations, so we define a function first
  psi_dgg <- function(d,g1,g2) {
    if (d==0 & g1==0) {
      IPO_arg <- IPO_D0
      YgivenX.Pred_arg <- YgivenGX.Pred_D0}
    if (d==1 & g1==0) {
      IPO_arg <- IPO_D1
      YgivenX.Pred_arg <- YgivenGX.Pred_D1}
    if (d==0 & g1==1) {
      IPO_arg <- IPO_D0
      YgivenX.Pred_arg <- YgivenGX.Pred_D0}
    if (d==1 & g1==1) {
      IPO_arg <- IPO_D1
      YgivenX.Pred_arg <- YgivenGX.Pred_D1}

    psi_dgg <- mean( as.numeric(data[,G]==g1)/mean(data[,G]==g1)*IPO_arg*mean(as.numeric(data[,G]==g2)/mean(data[,G]==g2)*data[,D]) )
    # Note that this is basically DML2. We could also use DML1:
    #psi_dgg_S1 <- mean( as.numeric(data[sample1,G]==g1)/mean(data[sample1,G]==g1)*IPO_arg[sample1]*mean(as.numeric(data[sample1,G]==g2)/mean(data[sample1,G]==g2)*data[sample1,D]) +
    #                      as.numeric(data[sample1,G]==g2)/mean(data[sample1,G]==g2)*mean(as.numeric(data[sample1,G]==g1)/mean(data[sample1,G]==g1)*YgivenX.Pred_arg)*(data[sample1,D]-mean(as.numeric(data[sample1,G]==g2)/mean(data[sample1,G]==g2)*data[sample1,D])) )
    #psi_dgg_S2 <- mean( as.numeric(data[sample2,G]==g1)/mean(data[sample2,G]==g1)*IPO_arg[sample2]*mean(as.numeric(data[sample2,G]==g2)/mean(data[sample2,G]==g2)*data[sample2,D]) +
    #                      as.numeric(data[sample2,G]==g2)/mean(data[sample2,G]==g2)*mean(as.numeric(data[sample2,G]==g1)/mean(data[sample2,G]==g1)*YgivenX.Pred_arg)*(data[sample2,D]-mean(as.numeric(data[sample2,G]==g2)/mean(data[sample2,G]==g2)*data[sample2,D])) )
    #psi_dgg <- (1/2)*(psi_dgg_S1+psi_dgg_S2)

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

  Jackson_reduction <- psi_00+psi_dgg(1,0,1)-psi_dgg(0,0,1)-mean((1-data[,G])/(1-mean(data[,G]))*data[,Y])

  ### standard error estimates
  se <- function(x) {sqrt( mean(x^2)/nrow(data) )}
  total_se <- se( data[,G]/mean(data[,G])*(data[,Y]-Y_G1) - (1-data[,G])/(1-mean(data[,G]))*(data[,Y]-Y_G0) )
  baseline_se <- se( data[,G]/mean(data[,G])*(IPO_D0-psi_01) - (1-data[,G])/(1-mean(data[,G]))*(IPO_D0-psi_00) )
  # Alternatively, we could use
  # se( c( data[sample1,G]/mean(data[sample1,G])*(IPO_D0[sample1]-psi_01) - (1-data[sample1,G])/(1-mean(data[sample1,G]))*(IPO_D0[sample1]-psi_00),
  #         data[sample2,G]/mean(data[sample2,G])*(IPO_D0[sample2]-psi_01) - (1-data[sample2,G])/(1-mean(data[sample2,G]))*(IPO_D0[sample2]-psi_00) ) )
  # But there isn't a theoretically strong reason to prefer one over the other.

  EIF_dgg <- function(d,g1,g2) {
    if (d==0 & g1==0) {
      IPO_arg <- IPO_D0
      YgivenX.Pred_arg <- YgivenGX.Pred_D0
      psi_arg <- psi_00}
    if (d==1 & g1==0) {
      IPO_arg <- IPO_D1
      YgivenX.Pred_arg <- YgivenGX.Pred_D1
      psi_arg <- psi_10}
    if (d==0 & g1==1) {
      IPO_arg <- IPO_D0
      YgivenX.Pred_arg <- YgivenGX.Pred_D0
      psi_arg <- psi_01}
    if (d==1 & g1==1) {
      IPO_arg <- IPO_D1
      YgivenX.Pred_arg <- YgivenGX.Pred_D1
      psi_arg <- psi_11}

    return(
      as.numeric(data[,G]==g1)/mean(data[,G]==g1)*IPO_arg*mean(as.numeric(data[,G]==g2)/mean(data[,G]==g2)*data[,D]) +
        as.numeric(data[,G]==g2)/mean(data[,G]==g2)*psi_arg*(data[,D]-mean(as.numeric(data[,G]==g2)/mean(data[,G]==g2)*data[,D])) -
        as.numeric(data[,G]==g1)/mean(data[,G]==g1)*psi_dgg(d,g1,g2)
    )
  }
  # Alternatively, we could use
  # return(
  #   c(as.numeric(data[sample1,G]==g1)/mean(data[sample1,G]==g1)*IPO_arg[sample1]*mean(as.numeric(data[sample1,G]==g2)/mean(data[sample1,G]==g2)*data[sample1,D]) +
  #       as.numeric(data[sample1,G]==g2)/mean(data[sample1,G]==g2)*mean(as.numeric(data[sample1,G]==g1)/mean(data[sample1,G]==g1)*YgivenX.Pred_arg)*(data[sample1,D]-mean(as.numeric(data[sample1,G]==g2)/mean(data[sample1,G]==g2)*data[sample1,D])) -
  #       as.numeric(data[sample1,G]==g1)/mean(data[sample1,G]==g1)*psi_dgg(d,g1,g2)
  #     , as.numeric(data[sample2,G]==g1)/mean(data[sample2,G]==g1)*IPO_arg[sample2]*mean(as.numeric(data[sample2,G]==g2)/mean(data[sample2,G]==g2)*data[sample2,D]) +
  #       as.numeric(data[sample2,G]==g2)/mean(data[sample2,G]==g2)*mean(as.numeric(data[sample2,G]==g1)/mean(data[sample2,G]==g1)*YgivenX.Pred_arg)*(data[sample2,D]-mean(as.numeric(data[sample2,G]==g2)/mean(data[sample2,G]==g2)*data[sample2,D])) -
  #       as.numeric(data[sample2,G]==g1)/mean(data[sample2,G]==g1)*psi_dgg(d,g1,g2))
  # )
  # But there isn't a theoretically strong reason to prefer one over the other.

  prevalence_se <- se( EIF_dgg(1,0,1)-EIF_dgg(1,0,0)-EIF_dgg(0,0,1)+EIF_dgg(0,0,0) )
  effect_se <- se( EIF_dgg(1,1,1)-EIF_dgg(0,1,1)-EIF_dgg(1,0,1)+EIF_dgg(0,0,1) )
  selection_se <- se( data[,G]/mean(data[,G])*(data[,Y]-Y_G1) - (1-data[,G])/(1-mean(data[,G]))*(data[,Y]-Y_G0) -
                        ( data[,G]/mean(data[,G])*(IPO_D0-psi_01) - (1-data[,G])/(1-mean(data[,G]))*(IPO_D0-psi_00) ) -
                        ( EIF_dgg(1,0,1)-EIF_dgg(1,0,0)-EIF_dgg(0,0,1)+EIF_dgg(0,0,0) ) -
                        ( EIF_dgg(1,1,1)-EIF_dgg(0,1,1)-EIF_dgg(1,0,1)+EIF_dgg(0,0,1) ) )

  Jackson_reduction_se <- se( (1-data[,G])/(1-mean(data[,G]))*(IPO_D0-psi_00)+EIF_dgg(1,0,1)-EIF_dgg(0,0,1)-(1-data[,G])/(1-mean(data[,G]))*(data[,Y]-Y_G0) )

  ### output results
  point <- c(total,
             baseline,
             prevalence,
             effect,
             selection)

  point_specific <- c(Y_G1,
                      Y_G0,
                      psi_01,
                      psi_00,
                      mean(data[,G]/mean(data[,G])*data[,D]),
                      mean((1-data[,G])/(1-mean(data[,G]))*data[,D]),
                      mean(data[,G]/mean(data[,G])*data[,D])-mean((1-data[,G])/(1-mean(data[,G]))*data[,D]),
                      mean(data[,G]/mean(data[,G])*(IPO_D1-IPO_D0)),
                      mean((1-data[,G])/(1-mean(data[,G]))*(IPO_D1-IPO_D0)),
                      mean(data[,G]/mean(data[,G])*(IPO_D1-IPO_D0)) - mean((1-data[,G])/(1-mean(data[,G]))*(IPO_D1-IPO_D0)),
                      Y_G1-psi_01-psi_dgg(1,1,1)+psi_dgg(0,1,1),
                      Y_G0-psi_00-psi_dgg(1,0,0)+psi_dgg(0,0,0),
                      Jackson_reduction)

  se_est <- c(total_se,
              baseline_se,
              prevalence_se,
              effect_se,
              selection_se)

  se_est_specific <- c(se( data[,G]/mean(data[,G])*(data[,Y]-Y_G1) ),
                       se( (1-data[,G])/(1-mean(data[,G]))*(data[,Y]-Y_G0) ),
                       se( data[,G]/mean(data[,G])*(IPO_D0-psi_01)),
                       se( (1-data[,G])/(1-mean(data[,G]))*(IPO_D0-psi_00)),
                       se( data[,G]/mean(data[,G])*(data[,D]-mean(data[,G]/mean(data[,G])*data[,D])) ),
                       se( (1-data[,G])/(1-mean(data[,G]))*(data[,D]-mean((1-data[,G])/(1-mean(data[,G]))*data[,D])) ),
                       se( data[,G]/mean(data[,G])*(data[,D]-mean(data[,G]/mean(data[,G])*data[,D])) - (1-data[,G])/(1-mean(data[,G]))*(data[,D]-mean((1-data[,G])/(1-mean(data[,G]))*data[,D])) ),
                       se( data[,G]/mean(data[,G])*(IPO_D1-IPO_D0-mean(data[,G]/mean(data[,G])*(IPO_D1-IPO_D0))) ),
                       se( (1-data[,G])/(1-mean(data[,G]))*(IPO_D1-IPO_D0-mean((1-data[,G])/(1-mean(data[,G]))*(IPO_D1-IPO_D0))) ),
                       se( data[,G]/mean(data[,G])*(IPO_D1-IPO_D0-mean(data[,G]/mean(data[,G])*(IPO_D1-IPO_D0))) - (1-data[,G])/(1-mean(data[,G]))*(IPO_D1-IPO_D0-mean((1-data[,G])/(1-mean(data[,G]))*(IPO_D1-IPO_D0))) ),
                       se( data[,G]/mean(data[,G])*(data[,Y]-Y_G1)-data[,G]/mean(data[,G])*(IPO_D0-psi_01)-EIF_dgg(1,1,1)+EIF_dgg(0,1,1) ),
                       se( (1-data[,G])/(1-mean(data[,G]))*(data[,Y]-Y_G0)-(1-data[,G])/(1-mean(data[,G]))*(IPO_D0-psi_00)-EIF_dgg(1,0,0)+EIF_dgg(0,0,0) ),
                       Jackson_reduction_se)

  p_value <- (1-stats::pnorm(abs(point/se_est)))*2
  CI_lower <- point - stats::qnorm(1-alpha/2)*se_est
  CI_upper <- point + stats::qnorm(1-alpha/2)*se_est

  p_value_specific <- (1-stats::pnorm(abs(point_specific/se_est_specific)))*2
  CI_lower_specific <- point_specific - stats::qnorm(1-alpha/2)*se_est_specific
  CI_upper_specific <- point_specific + stats::qnorm(1-alpha/2)*se_est_specific

  names <- c("total",
             "baseline",
             "prevalence",
             "effect",
             "selection")

  names_specific <- c("Y_G1",
                      "Y_G0",
                      "Y0_G1",
                      "Y0_G0",
                      "D_G1",
                      "D_G0",
                      "D_G1-D_G0",
                      "ATE_G1",
                      "ATE_G0",
                      "ATE_G1-ATE_G0",
                      "Cov_G1",
                      "Cov_G0",
                      "Jackson reduction")

  results <- as.data.frame(cbind(point,se_est,p_value,CI_lower,CI_upper))
  results_specific <- as.data.frame(cbind(point_specific,se_est_specific,p_value_specific,CI_lower_specific,CI_upper_specific))
  rownames(results) <- names
  rownames(results_specific) <- names_specific
  colnames(results) <- colnames(results_specific) <- c("point","se","p_value","CI_lower","CI_upper")

  output <- list(results=results, results_specific=results_specific)

  return(output)
}

