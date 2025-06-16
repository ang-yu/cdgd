


#' Perform conditional decomposition via parametric models
#'
#' @param Y Outcome. The name of a numeric variable (can be binary and take values of 0 and 1).
#' @param D Treatment status. The name of a binary numeric variable taking values of 0 and 1.
#' @param G Advantaged group membership. The name of a binary numeric variable taking values of 0 and 1.
#' @param Q Conditional set. A vector of variable names.
#' @param X Confounders. A vector of variable names.
#' @param data A data frame.
#' @param alpha 1-alpha confidence interval.
#' @param trim1 Threshold for trimming the propensity score. When trim1=a, individuals with propensity scores lower than a or higher than 1-a will be dropped.
#' @param trim2 Threshold for trimming the G given Q predictions. When trim2=a, individuals with G given Q predictions lower than a or higher than 1-a will be dropped.
#' @param weight Sampling weights. The name of a numeric variable. If unspecified, equal weights are used. Technically, the weight should be a deterministic function of X only (note that this is different from the unconditional decomposition).
#'
#' @return A dataframe of estimates.
#'
#' @export
#'
#' @examples
#' data(exp_data)
#'
#' results <- cdgd1_pa(
#' Y="outcome",
#' D="treatment",
#' G="group_a",
#' X="confounder",
#' Q="Q",
#' data=exp_data)
#'
#' results



cdgd1_pa <- function(Y,D,G,X,Q,data,alpha=0.05,trim1=0,trim2=0,weight=NULL) {

  data <- as.data.frame(data)

  if ( sum(is.na(data[,c(Y,D,G,X,Q)]))>0 ) {
    stop(
      "There are missing values in key variables.",
      call. = FALSE
    )
  }

  ### treatment model
  DgivenGXQ.Model <- stats::glm(stats::as.formula(paste(D, paste(G,paste(Q,collapse="+"),paste(X,collapse="+"),sep="+"), sep="~")), data=data, family=stats::binomial(link="logit"))

  # treatment predictions
  DgivenGXQ.Pred <- rep(NA, nrow(data))
  DgivenGXQ.Pred <- stats::predict(DgivenGXQ.Model, newdata = data, type="response")

  # trim the sample based on the propensity score
  dropped <- sum(DgivenGXQ.Pred<trim1 | DgivenGXQ.Pred>1-trim1)  # the number of dropped obs
  data <- data[DgivenGXQ.Pred>=trim1 & DgivenGXQ.Pred<=1-trim1, ]
  DgivenGXQ.Pred <- DgivenGXQ.Pred[DgivenGXQ.Pred>=trim1 & DgivenGXQ.Pred<=1-trim1]

  zero_one <- sum(DgivenGXQ.Pred==0)+sum(DgivenGXQ.Pred==1)
  if ( zero_one>0 ) {
    stop(
      paste("D given X, Q, and G are exact 0 or 1 in", zero_one, "cases.", sep=" "),
      call. = FALSE
    )
  }

  ### Estimate p_g(Q)=Pr(G=g | Q)
  GgivenQ.Model <- stats::glm(stats::as.formula(paste(G, paste(Q,collapse="+"), sep="~")), data=data, family=stats::binomial(link="logit"))

  GgivenQ.Pred <- rep(NA, nrow(data))
  GgivenQ.Pred <- stats::predict(GgivenQ.Model, newdata = data, type="response")

  # trim the sample based on the G given Q predictions
  dropped <- dropped + sum(GgivenQ.Pred<trim2 | GgivenQ.Pred>1-trim2)  # update the number of dropped obs
  data <- data[GgivenQ.Pred>=trim2 & GgivenQ.Pred<=1-trim2, ]
  DgivenGXQ.Pred <- DgivenGXQ.Pred[GgivenQ.Pred>=trim2 & GgivenQ.Pred<=1-trim2]
  GgivenQ.Pred <- GgivenQ.Pred[GgivenQ.Pred>=trim2 & GgivenQ.Pred<=1-trim2]

  zero_one <- sum(GgivenQ.Pred==0)+sum(GgivenQ.Pred==1)
  if ( zero_one>0 ) {
    stop(
      paste("G given Q are exact 0 or 1 in", zero_one, "cases.", sep=" "),
      call. = FALSE
    )
  }

  ### outcome regression model
  YgivenDGXQ.Model <- stats::lm(stats::as.formula(paste(Y, paste(paste(D,c(G,Q,X),sep="*"),collapse="+"), sep="~")), data=data)

  # outcome predictions
  YgivenGXQ.Pred_D0 <- YgivenGXQ.Pred_D1 <- rep(NA, nrow(data))

  pred_data <- data
  pred_data[,D] <- 0
  YgivenGXQ.Pred_D0 <- stats::predict(YgivenDGXQ.Model, newdata = pred_data)

  pred_data <- data
  pred_data[,D] <- 1
  YgivenGXQ.Pred_D1 <- stats::predict(YgivenDGXQ.Model, newdata = pred_data)

  ### Estimate E(Y_d | Q,g)
  YgivenGXQ.Pred_D1 <- YgivenGXQ.Pred_D0 <- DgivenGXQ.Pred <- rep(NA, nrow(data))

  pred_data <- data
  pred_data[,D] <- 1
  YgivenGXQ.Pred_D1 <- stats::predict(YgivenDGXQ.Model, newdata = pred_data)

  pred_data <- data
  pred_data[,D] <- 0
  YgivenGXQ.Pred_D0 <- stats::predict(YgivenDGXQ.Model, newdata = pred_data)

  DgivenGXQ.Pred <- stats::predict(DgivenGXQ.Model, newdata = pred_data, type="response")

  zero_one <- sum(DgivenGXQ.Pred==0)+sum(DgivenGXQ.Pred==1)
  if ( zero_one>0 ) {
    stop(
      paste("D given X, Q, and G are exact 0 or 1 in", zero_one, "cases.", sep=" "),
      call. = FALSE
    )
  }

  ### The "IPO" (individual potential outcome) function
  # For each d and g value, we have IE(d,g)=\frac{\one(D=d)}{\pi(d,X,g)}[Y-\mu(d,X,g)]+\mu(d,X,g)
  # We stabilize the weight by dividing the sample average of estimated weights
  IPO_D0 <- (1-data[,D])/(1-DgivenGXQ.Pred)/mean((1-data[,D])/(1-DgivenGXQ.Pred))*(data[,Y]-YgivenGXQ.Pred_D0) + YgivenGXQ.Pred_D0
  IPO_D1 <- data[,D]/DgivenGXQ.Pred/mean(data[,D]/DgivenGXQ.Pred)*(data[,Y]-YgivenGXQ.Pred_D1) + YgivenGXQ.Pred_D1

  if (is.null(weight)) {
    weight <- rep(1, nrow(data))
    tr.weight <- rep(1, nrow(data))
  } else {
    weight <- data[,weight]
    tr.weight <- weight/stats::predict(stats::lm(stats::as.formula(paste("weight", paste(paste("data[,G]","data[,Q]",sep="*"),collapse="+"), sep="~"))))
  }
  # tr.weight (transformed weight) is the original weight divided by E(weight|G,Q)

  data_temp <- data[,c(G,Q)]
  data_temp$IPO_D0 <- IPO_D0*tr.weight
  data_temp$IPO_D1 <- IPO_D1*tr.weight
  data_temp[,D] <- data[,D]*tr.weight

  Y0givenGQ.Model <- stats::lm(stats::as.formula(paste("IPO_D0", paste(paste(G,Q,sep="*"),collapse="+"), sep="~")), data=data_temp)
  Y1givenGQ.Model <- stats::lm(stats::as.formula(paste("IPO_D1", paste(paste(G,Q,sep="*"),collapse="+"), sep="~")), data=data_temp)

  Y0givenQ.Pred_G0 <- Y0givenQ.Pred_G1 <- Y1givenQ.Pred_G0 <- Y1givenQ.Pred_G1 <- rep(NA, nrow(data))

  pred_data <- data
  pred_data[,G] <- 1
  Y0givenQ.Pred_G1 <- stats::predict(Y0givenGQ.Model, newdata = pred_data)
  Y1givenQ.Pred_G1 <- stats::predict(Y1givenGQ.Model, newdata = pred_data)

  pred_data <- data
  pred_data[,G] <- 0
  Y0givenQ.Pred_G0 <- stats::predict(Y0givenGQ.Model, newdata = pred_data)
  Y1givenQ.Pred_G0 <- stats::predict(Y1givenGQ.Model, newdata = pred_data)

  ### Estimate E(D | Q,g')
  if (all(weight == 1)) {
    DgivenGQ.Model <- stats::glm(stats::as.formula(paste(D, paste(paste(G,Q,sep="*"),collapse="+"), sep="~")), data=data_temp, family=stats::binomial(link="logit"))

    DgivenQ.Pred_G0 <- DgivenQ.Pred_G1 <- rep(NA, nrow(data))

    pred_data <- data
    pred_data[,G] <- 0
    DgivenQ.Pred_G0 <- stats::predict(DgivenGQ.Model, newdata = pred_data, type="response")

    pred_data <- data
    pred_data[,G] <- 1
    DgivenQ.Pred_G1 <- stats::predict(DgivenGQ.Model, newdata = pred_data, type="response")
  } else {
    DgivenGQ.Model <- stats::lm(stats::as.formula(paste(D, paste(paste(G,Q,sep="*"),collapse="+"), sep="~")), data=data_temp)

    DgivenQ.Pred_G0 <- DgivenQ.Pred_G1 <- rep(NA, nrow(data))

    pred_data <- data
    pred_data[,G] <- 0
    DgivenQ.Pred_G0 <- stats::predict(DgivenGQ.Model, newdata = pred_data)

    pred_data <- data
    pred_data[,G] <- 1
    DgivenQ.Pred_G1 <- stats::predict(DgivenGQ.Model, newdata = pred_data)
  }



  ### The one-step estimate of \xi_{dg}
  weight0 <- (1-data[,G])/(1-mean(data[,G]))*weight/mean((1-data[,G])/(1-mean(data[,G]))*weight)
  weight1 <- data[,G]/mean(data[,G])*weight/mean(data[,G]/mean(data[,G])*weight)
  psi_00 <- mean( weight0*IPO_D0 )
  psi_01 <- mean( weight1*IPO_D0 )

  ### The one-step estimate of \xi_{dgg'g''}
  # There are 8 possible dgg'g'' combinations, so we define a function first
  psi_dggg <- function(d,g1,g2,g3) {
    if (d==0 & g1==0) {
      IPO_arg <- IPO_D0
      YdgivenQ.Pred_arg <- Y0givenQ.Pred_G0
      g1givenQ.Pred_arg <- 1-GgivenQ.Pred}
    if (d==1 & g1==0) {
      IPO_arg <- IPO_D1
      YdgivenQ.Pred_arg <- Y1givenQ.Pred_G0
      g1givenQ.Pred_arg <- 1-GgivenQ.Pred}
    if (d==0 & g1==1) {
      IPO_arg <- IPO_D0
      YdgivenQ.Pred_arg <- Y0givenQ.Pred_G1
      g1givenQ.Pred_arg <- GgivenQ.Pred}
    if (d==1 & g1==1) {
      IPO_arg <- IPO_D1
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

    weight_g3 <- weight/mean(as.numeric(data[,G]==g3)/mean(data[,G]==g3)*weight)

    psi_dggg <- mean( weight_g3*as.numeric(data[,G]==g3)/mean(data[,G]==g3)*YdgivenQ.Pred_arg*DgivenQ.Pred_arg +
                        weight_g3*tr.weight*as.numeric(data[,G]==g1)/mean(data[,G]==g3)*g3givenQ.Pred_arg/g1givenQ.Pred_arg/stab1*(IPO_arg-YdgivenQ.Pred_arg)*DgivenQ.Pred_arg +
                        weight_g3*tr.weight*as.numeric(data[,G]==g2)/mean(data[,G]==g3)*g3givenQ.Pred_arg/g2givenQ.Pred_arg/stab2*(data[,D]-DgivenQ.Pred_arg)*YdgivenQ.Pred_arg )

    return(psi_dggg)
  }

  ### point estimates
  Y_G0 <- mean(weight0*data[,Y])       # mean outcome estimate for group 0
  Y_G1 <- mean(weight1*data[,Y])       # mean outcome estimate for group 1
  total <- Y_G1-Y_G0

  baseline <- psi_01-psi_00
  cond_prevalence <- psi_dggg(1,0,1,0)-psi_dggg(0,0,1,0)-psi_dggg(1,0,0,0)+psi_dggg(0,0,0,0)
  cond_effect <- psi_dggg(1,1,1,1)-psi_dggg(0,1,1,1)-psi_dggg(1,0,1,1)+psi_dggg(0,0,1,1)
  Q_dist <- psi_dggg(1,0,1,1)-psi_dggg(0,0,1,1)-psi_dggg(1,0,1,0)+psi_dggg(0,0,1,0)
  cond_selection <- total-baseline-cond_prevalence-cond_effect-Q_dist

  cond_Jackson_reduction <- psi_00+psi_dggg(1,0,1,0)-psi_dggg(0,0,1,0)-Y_G0

  ### standard error estimates
  se <- function(x) {sqrt( mean(x^2)/nrow(data) )}
  total_se <- se( weight1*(data[,Y]-Y_G1) - weight0*(data[,Y]-Y_G0) )
  baseline_se <- se( weight1*(IPO_D0-psi_01) - weight0*(IPO_D0-psi_00) )

  EIF_dggg <- function(d,g1,g2,g3) {
    if (d==0 & g1==0) {
      IPO_arg <- IPO_D0
      YdgivenQ.Pred_arg <- Y0givenQ.Pred_G0
      g1givenQ.Pred_arg <- 1-GgivenQ.Pred}
    if (d==1 & g1==0) {
      IPO_arg <- IPO_D1
      YdgivenQ.Pred_arg <- Y1givenQ.Pred_G0
      g1givenQ.Pred_arg <- 1-GgivenQ.Pred}
    if (d==0 & g1==1) {
      IPO_arg <- IPO_D0
      YdgivenQ.Pred_arg <- Y0givenQ.Pred_G1
      g1givenQ.Pred_arg <- GgivenQ.Pred}
    if (d==1 & g1==1) {
      IPO_arg <- IPO_D1
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

    weight_g3 <- weight/mean(as.numeric(data[,G]==g3)/mean(data[,G]==g3)*weight)

    return(
      weight_g3*as.numeric(data[,G]==g3)/mean(data[,G]==g3)*(YdgivenQ.Pred_arg*DgivenQ.Pred_arg-psi_dggg(d,g1,g2,g3)) +
        weight_g3*tr.weight*as.numeric(data[,G]==g1)/mean(data[,G]==g3)*g3givenQ.Pred_arg/g1givenQ.Pred_arg/stab1*(IPO_arg-YdgivenQ.Pred_arg)*DgivenQ.Pred_arg +
        weight_g3*tr.weight*as.numeric(data[,G]==g2)/mean(data[,G]==g3)*g3givenQ.Pred_arg/g2givenQ.Pred_arg/stab2*(data[,D]-DgivenQ.Pred_arg)*YdgivenQ.Pred_arg
    )
  }

  cond_prevalence_se <- se( EIF_dggg(1,0,1,0)-EIF_dggg(0,0,1,0)-EIF_dggg(1,0,0,0)+EIF_dggg(0,0,0,0) )
  cond_effect_se <- se( EIF_dggg(1,1,1,1)-EIF_dggg(0,1,1,1)-EIF_dggg(1,0,1,1)+EIF_dggg(0,0,1,1) )
  Q_dist_se <- se( EIF_dggg(1,0,1,1)-EIF_dggg(0,0,1,1)-EIF_dggg(1,0,1,0)+EIF_dggg(0,0,1,0) )
  cond_selection_se <- se( weight1*(data[,Y]-Y_G1) - weight0*(data[,Y]-Y_G0) -
                             ( weight1*(IPO_D0-psi_01) - weight0*(IPO_D0-psi_00) ) -
                             ( EIF_dggg(1,0,1,0)-EIF_dggg(0,0,1,0)-EIF_dggg(1,0,0,0)+EIF_dggg(0,0,0,0) ) -
                             ( EIF_dggg(1,1,1,1)-EIF_dggg(0,1,1,1)-EIF_dggg(1,0,1,1)+EIF_dggg(0,0,1,1) ) -
                             ( EIF_dggg(1,0,1,1)-EIF_dggg(0,0,1,1)-EIF_dggg(1,0,1,0)+EIF_dggg(0,0,1,0) ))

  cond_Jackson_reduction_se <- se( weight0*(IPO_D0-psi_00)+EIF_dggg(1,0,1,0)-EIF_dggg(0,0,1,0)-weight0*(data[,Y]-Y_G0) )

  ### output results
  point <- c(total,
             baseline,
             cond_prevalence,
             cond_effect,
             cond_selection,
             Q_dist,
             cond_Jackson_reduction)
  se <- c(total_se,
          baseline_se,
          cond_prevalence_se,
          cond_effect_se,
          cond_selection_se,
          Q_dist_se,
          cond_Jackson_reduction_se)
  p_value <- (1-stats::pnorm(abs(point/se)))*2
  CI_lower <- point - stats::qnorm(1-alpha/2)*se
  CI_upper <- point + stats::qnorm(1-alpha/2)*se
  names <- c("total",
             "baseline",
             "conditional prevalence",
             "conditional effect",
             "conditional selection",
             "Q distribution",
             "conditional Jackson reduction")

  output <- as.data.frame(cbind(point,se,p_value,CI_lower,CI_upper))
  rownames(output) <- names

  if (trim1==0 & trim2==0) {
    output <- output
  } else {
    output <- list(results=output, dropped=dropped)
  }

  return(output)
}
