


#' Perform unconditional decomposition via parametric models
#'
#' @param Y Outcome. The name of a numeric variable (can be binary and take values of 0 and 1).
#' @param D Treatment status. The name of a binary numeric variable taking values of 0 and 1.
#' @param G Advantaged group membership. The name of a binary numeric variable taking values of 0 and 1.
#' @param X Confounders. A vector of variable names.
#' @param data A data frame.
#' @param alpha 1-alpha confidence interval.
#' @param trim Threshold for trimming the propensity score. When trim=a, individuals with propensity scores lower than a or higher than 1-a will be dropped.
#' @param weight Sampling weights. The name of a numeric variable. If unspecified, equal weights are used. Technically, the weight should be a deterministic function of X and G.
#'
#' @return A list of estimates.
#'
#' @export
#'
#' @examples
#' data(exp_data)
#'
#' results <- cdgd0_pa(
#' Y="outcome",
#' D="treatment",
#' G="group_a",
#' X=c("Q","confounder"),
#' data=exp_data)
#'
#' results[[1]]




cdgd0_pa <- function(Y,D,G,X,data,alpha=0.05,trim=0,weight=NULL) {

  data <- as.data.frame(data)

  if ( sum(is.na(data[,c(Y,D,G,X)]))>0 ) {
    stop(
      "There are missing values in key variables.",
      call. = FALSE
    )
  }

  # treatment model
  DgivenGX.Model <- stats::glm(stats::as.formula(paste(D, paste(G,paste(X,collapse="+"),sep="+"), sep="~")), data=data, family=stats::binomial(link="logit"))

  # treatment predictions
  DgivenGX.Pred <- rep(NA, nrow(data))
  DgivenGX.Pred <- stats::predict(DgivenGX.Model, newdata = data, type="response")

  # trim the sample based on the propensity score
  dropped <- sum(DgivenGX.Pred<trim | DgivenGX.Pred>1-trim)  # the number of dropped obs
  data <- data[DgivenGX.Pred>=trim & DgivenGX.Pred<=1-trim, ]
  DgivenGX.Pred <- DgivenGX.Pred[DgivenGX.Pred>=trim & DgivenGX.Pred<=1-trim]

  zero_one <- sum(DgivenGX.Pred==0)+sum(DgivenGX.Pred==1)
  if ( zero_one>0 ) {
    stop(
      paste("D given X and G are exact 0 or 1 in", zero_one, "cases.", sep=" "),
      call. = FALSE
    )
  }

  ### outcome regression model
  YgivenDGX.Model <- stats::lm(stats::as.formula(paste(Y, paste(paste(D,c(G,X),sep="*"),collapse="+"), sep="~")), data=data)

  # outcome predictions
  YgivenGX.Pred_D0 <- YgivenGX.Pred_D1 <- rep(NA, nrow(data))

  pred_data <- data
  pred_data[,D] <- 0
  YgivenGX.Pred_D0 <- stats::predict(YgivenDGX.Model, newdata = pred_data)

  pred_data <- data
  pred_data[,D] <- 1
  YgivenGX.Pred_D1 <- stats::predict(YgivenDGX.Model, newdata = pred_data)

  ### The "IPO" (individual potential outcome) function
  # For each d and g value, we have IE(d,g)=\frac{\one(D=d)}{\pi(d,X,g)}[Y-\mu(d,X,g)]+\mu(d,X,g)
  # We stabilize the weight by dividing the sample average of estimated weights

  IPO_D0 <- (1-data[,D])/(1-DgivenGX.Pred)/mean((1-data[,D])/(1-DgivenGX.Pred))*(data[,Y]-YgivenGX.Pred_D0) + YgivenGX.Pred_D0
  IPO_D1 <- data[,D]/DgivenGX.Pred/mean(data[,D]/DgivenGX.Pred)*(data[,Y]-YgivenGX.Pred_D1) + YgivenGX.Pred_D1

  if (is.null(weight)) {
    weight <- rep(1, nrow(data))
  } else {
    weight <- data[,weight]
  }
  weight0 <- (1-data[,G])/(1-mean(data[,G]))*weight/mean((1-data[,G])/(1-mean(data[,G]))*weight)
  weight1 <- data[,G]/mean(data[,G])*weight/mean(data[,G]/mean(data[,G])*weight)

  ### The one-step estimate of \xi_{dg} and \xi_{dgg'}
  psi_00 <- mean( weight0*IPO_D0 )
  psi_01 <- mean( weight1*IPO_D0 )
  psi_10 <- mean( weight0*IPO_D1 )
  psi_11 <- mean( weight1*IPO_D1 )

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

    weight_g1 <- as.numeric(data[,G]==g1)/mean(data[,G]==g1)*weight/mean(as.numeric(data[,G]==g1)/mean(data[,G]==g1)*weight)
    weight_g2 <- as.numeric(data[,G]==g2)/mean(data[,G]==g2)*weight/mean(as.numeric(data[,G]==g2)/mean(data[,G]==g2)*weight)

    psi_dgg <- mean( weight_g1*IPO_arg*mean(weight_g2*data[,D]) )

    return(psi_dgg)
  }

  ### point estimates
  Y_G0 <- mean(weight0*data[,Y])       # mean outcome estimate for group 0
  Y_G1 <- mean(weight1*data[,Y])       # mean outcome estimate for group 1
  total <- Y_G1-Y_G0

  baseline <- psi_01-psi_00
  prevalence <- psi_dgg(1,0,1)-psi_dgg(1,0,0)-psi_dgg(0,0,1)+psi_dgg(0,0,0)
  effect <- psi_dgg(1,1,1)-psi_dgg(0,1,1)-psi_dgg(1,0,1)+psi_dgg(0,0,1)
  selection <- total-baseline-prevalence-effect

  Jackson_reduction <- psi_00+psi_dgg(1,0,1)-psi_dgg(0,0,1)-Y_G0

  ### standard error estimates
  se <- function(x) {sqrt( mean(x^2)/nrow(data) )}
  total_se <- se( weight1*(data[,Y]-Y_G1) - weight0*(data[,Y]-Y_G0) )
  baseline_se <- se( weight1*(IPO_D0-psi_01) - weight0*(IPO_D0-psi_00) )

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

    weight_g1 <- as.numeric(data[,G]==g1)/mean(data[,G]==g1)*weight/mean(as.numeric(data[,G]==g1)/mean(data[,G]==g1)*weight)
    weight_g2 <- as.numeric(data[,G]==g2)/mean(data[,G]==g2)*weight/mean(as.numeric(data[,G]==g2)/mean(data[,G]==g2)*weight)

    return(
      weight_g1*IPO_arg*mean(weight_g2*data[,D]) +
        weight_g2*psi_arg*(data[,D]-mean(weight_g2*data[,D])) -
        weight_g1*psi_dgg(d,g1,g2)
    )
  }

  prevalence_se <- se( EIF_dgg(1,0,1)-EIF_dgg(1,0,0)-EIF_dgg(0,0,1)+EIF_dgg(0,0,0) )
  effect_se <- se( EIF_dgg(1,1,1)-EIF_dgg(0,1,1)-EIF_dgg(1,0,1)+EIF_dgg(0,0,1) )
  selection_se <- se( weight1*(data[,Y]-Y_G1) - weight0*(data[,Y]-Y_G0) -
                        ( weight1*(IPO_D0-psi_01) - weight0*(IPO_D0-psi_00) ) -
                        ( EIF_dgg(1,0,1)-EIF_dgg(1,0,0)-EIF_dgg(0,0,1)+EIF_dgg(0,0,0) ) -
                        ( EIF_dgg(1,1,1)-EIF_dgg(0,1,1)-EIF_dgg(1,0,1)+EIF_dgg(0,0,1) ) )

  Jackson_reduction_se <- se( weight0*(IPO_D0-psi_00)+EIF_dgg(1,0,1)-EIF_dgg(0,0,1)-weight0*(data[,Y]-Y_G0) )

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
                      mean(weight1*data[,D]),
                      mean(weight0*data[,D]),
                      mean(weight1*data[,D])-mean(weight0*data[,D]),
                      mean(weight1*(IPO_D1-IPO_D0)),
                      mean(weight0*(IPO_D1-IPO_D0)),
                      mean(weight1*(IPO_D1-IPO_D0)) - mean(weight0*(IPO_D1-IPO_D0)),
                      Y_G1-psi_01-psi_dgg(1,1,1)+psi_dgg(0,1,1),
                      Y_G0-psi_00-psi_dgg(1,0,0)+psi_dgg(0,0,0),
                      Jackson_reduction)

  se_est <- c(total_se,
              baseline_se,
              prevalence_se,
              effect_se,
              selection_se)

  se_est_specific <- c(se( weight1*(data[,Y]-Y_G1) ),
                       se( weight0*(data[,Y]-Y_G0) ),
                       se( weight1*(IPO_D0-psi_01)),
                       se( weight0*(IPO_D0-psi_00)),
                       se( weight1*(data[,D]-mean(weight1*data[,D])) ),
                       se( weight0*(data[,D]-mean(weight0*data[,D])) ),
                       se( weight1*(data[,D]-mean(weight1*data[,D])) - weight0*(data[,D]-mean(weight0*data[,D])) ),
                       se( weight1*(IPO_D1-IPO_D0-mean(weight1*(IPO_D1-IPO_D0))) ),
                       se( weight0*(IPO_D1-IPO_D0-mean(weight0*(IPO_D1-IPO_D0))) ),
                       se( weight1*(IPO_D1-IPO_D0-mean(weight1*(IPO_D1-IPO_D0))) - weight0*(IPO_D1-IPO_D0-mean(weight0*(IPO_D1-IPO_D0))) ),
                       se( weight1*(data[,Y]-Y_G1)-weight1*(IPO_D0-psi_01)-EIF_dgg(1,1,1)+EIF_dgg(0,1,1) ),
                       se( weight0*(data[,Y]-Y_G0)-weight0*(IPO_D0-psi_00)-EIF_dgg(1,0,0)+EIF_dgg(0,0,0) ),
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

  if (trim==0) {
    output <- list(results=results, results_specific=results_specific)
  } else {
    output <- list(results=results, results_specific=results_specific, dropped=dropped)
  }

  return(output)
}

