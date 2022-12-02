
#' Perform conditional decomposition using nuisances predicted beforehand
#'
#' This function gives user full control over the estimation of the nuisance terms.
#' For the unconditional decomposition, six nuisance terms need to be estimated.
#' The nuisance terms should be estimated using cross-fitting if Donsker class is not assumed.
#'
#' @param Y Outcome. The name of a continuous variable.
#' @param D Treatment status. The name of a binary numeric variable taking values of 0 and 1.
#' @param G Advantaged group membership. The name of a binary numeric variable taking values of 0 and 1.
#' @param YgivenXQ.Pred_D0G0 A numeric vector of predicted Y values given X, D=0, and G=0. Vector length=nrow(data).
#' @param YgivenXQ.Pred_D1G0 A numeric vector of predicted Y values given X, D=1, and G=0. Vector length=nrow(data).
#' @param YgivenXQ.Pred_D0G1 A numeric vector of predicted Y values given X, D=0, and G=1. Vector length=nrow(data).
#' @param YgivenXQ.Pred_D1G1 A numeric vector of predicted Y values given X, D=1, and G=1. Vector length=nrow(data).
#' @param DgivenXQ.Pred_G0 A numeric vector of predicted D values given X and G=1. Vector length=nrow(data).
#' @param DgivenXQ.Pred_G1 A numeric vector of predicted D values given X and G=0. Vector length=nrow(data).
#' @param Y0givenQ.Pred_G0 A numeric vector of predicted Y(0) values given Q and G=0. Vector length=nrow(data).
#' @param Y0givenQ.Pred_G1 A numeric vector of predicted Y(0) values given Q and G=0. Vector length=nrow(data).
#' @param Y1givenQ.Pred_G0 A numeric vector of predicted Y(1) values given Q and G=0. Vector length=nrow(data).
#' @param Y1givenQ.Pred_G1 A numeric vector of predicted Y(1) values given Q and G=1. Vector length=nrow(data).
#' @param DgivenQ.Pred_G0 A numeric vector of predicted D values given Q and G=0. Vector length=nrow(data).
#' @param DgivenQ.Pred_G1 A numeric vector of predicted D values given Q and G=0. Vector length=nrow(data).
#' @param GgivenQ.Pred A numeric vector of predicted G values given Q. Vector length=nrow(data).
#' @param data A data frame.
#' @param alpha 1-alpha confidence interval.
#'
#' @return A dataframe of estimates.
#'
#' @export
#'
#' @examples
#' data(exp_data)
#'
#' Y="outcome"
#' D="treatment"
#' G="group_a"
#' X="confounder"
#' Q="Q"
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
#' message <- utils::capture.output( YgivenDGXQ.Model.sample1 <- caret::train(stats::as.formula(paste(Y, paste(D,G,Q,paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample1,], method="ranger",
#'                                                                            trControl=caret::trainControl(method="cv"),
#'                                                                            tuneGrid=expand.grid(mtry=c(2,4),splitrule=c("variance"),min.node.size=c(50,100))) )
#' message <- utils::capture.output( YgivenDGXQ.Model.sample2 <- caret::train(stats::as.formula(paste(Y, paste(D,G,Q,paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample2,], method="ranger",
#'                                                                            trControl=caret::trainControl(method="cv"),
#'                                                                            tuneGrid=expand.grid(mtry=c(2,4),splitrule=c("variance"),min.node.size=c(50,100))) )
#'
#' ### propensity score model
#' data[,D] <- as.factor(data[,D])
#' levels(data[,D]) <- c("D0","D1")  # necessary for caret implementation of ranger
#'
#' message <- utils::capture.output( DgivenGXQ.Model.sample1 <- caret::train(stats::as.formula(paste(D, paste(G,Q,paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample1,], method="ranger",
#'                                                                           trControl=caret::trainControl(method="cv", classProbs=TRUE),
#'                                                                           tuneGrid=expand.grid(mtry=c(1,2),splitrule=c("gini"),min.node.size=c(50,100))) )
#' message <- utils::capture.output( DgivenGXQ.Model.sample2 <- caret::train(stats::as.formula(paste(D, paste(G,Q,paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample2,], method="ranger",
#'                                                                           trControl=caret::trainControl(method="cv", classProbs=TRUE),
#'                                                                           tuneGrid=expand.grid(mtry=c(1,2),splitrule=c("gini"),min.node.size=c(50,100))) )
#' data[,D] <- as.numeric(data[,D])-1
#'
#' ### cross-fitted predictions
#' YgivenXQ.Pred_D0G0 <- YgivenXQ.Pred_D1G0 <- YgivenXQ.Pred_D0G1 <- YgivenXQ.Pred_D1G1 <- DgivenXQ.Pred_G0 <- DgivenXQ.Pred_G1 <- rep(NA, nrow(data))
#'
#' pred_data <- data
#' pred_data[,D] <- 0
#' pred_data[,G] <- 0
#' YgivenXQ.Pred_D0G0[sample1] <- stats::predict(YgivenDGXQ.Model.sample1, newdata = pred_data[sample2,])
#' YgivenXQ.Pred_D0G0[sample2] <- stats::predict(YgivenDGXQ.Model.sample2, newdata = pred_data[sample1,])
#'
#' pred_data <- data
#' pred_data[,D] <- 1
#' pred_data[,G] <- 0
#' YgivenXQ.Pred_D1G0[sample1] <- stats::predict(YgivenDGXQ.Model.sample1, newdata = pred_data[sample2,])
#' YgivenXQ.Pred_D1G0[sample2] <- stats::predict(YgivenDGXQ.Model.sample2, newdata = pred_data[sample1,])
#'
#' pred_data <- data
#' pred_data[,D] <- 0
#' pred_data[,G] <- 1
#' YgivenXQ.Pred_D0G1[sample1] <- stats::predict(YgivenDGXQ.Model.sample1, newdata = pred_data[sample2,])
#' YgivenXQ.Pred_D0G1[sample2] <- stats::predict(YgivenDGXQ.Model.sample2, newdata = pred_data[sample1,])
#'
#' pred_data <- data
#' pred_data[,D] <- 1
#' pred_data[,G] <- 1
#' YgivenXQ.Pred_D1G1[sample1] <- stats::predict(YgivenDGXQ.Model.sample1, newdata = pred_data[sample2,])
#' YgivenXQ.Pred_D1G1[sample2] <- stats::predict(YgivenDGXQ.Model.sample2, newdata = pred_data[sample1,])
#'
#' pred_data <- data
#' pred_data[,G] <- 0
#' DgivenXQ.Pred_G0[sample1] <- stats::predict(DgivenGXQ.Model.sample1, newdata = pred_data[sample2,], type="prob")[,2]
#' DgivenXQ.Pred_G0[sample2] <- stats::predict(DgivenGXQ.Model.sample2, newdata = pred_data[sample1,], type="prob")[,2]
#'
#' pred_data <- data
#' pred_data[,G] <- 1
#' DgivenXQ.Pred_G1[sample1] <- stats::predict(DgivenGXQ.Model.sample1, newdata = pred_data[sample2,], type="prob")[,2]
#' DgivenXQ.Pred_G1[sample2] <- stats::predict(DgivenGXQ.Model.sample2, newdata = pred_data[sample1,], type="prob")[,2]
#'
#' ### Estimate E(Y_d | Q,g)
#' YgivenGXQ.Pred_D1_ncf <- YgivenGXQ.Pred_D0_ncf <- rep(NA, nrow(data)) # ncf stands for non-cross-fitted
#'
#' pred_data <- data
#' pred_data[,D] <- 1
#' YgivenGXQ.Pred_D1_ncf[sample1] <- stats::predict(YgivenDGXQ.Model.sample1, newdata = pred_data[sample1,])
#' YgivenGXQ.Pred_D1_ncf[sample2] <- stats::predict(YgivenDGXQ.Model.sample2, newdata = pred_data[sample2,])
#'
#' pred_data <- data
#' pred_data[,D] <- 0
#' YgivenGXQ.Pred_D0_ncf[sample1] <- stats::predict(YgivenDGXQ.Model.sample1, newdata = pred_data[sample1,])
#' YgivenGXQ.Pred_D0_ncf[sample2] <- stats::predict(YgivenDGXQ.Model.sample2, newdata = pred_data[sample2,])
#'
#' data_temp <- data[,c(G,Q)]
#' data_temp$YgivenGXQ.Pred_D0_ncf <- YgivenGXQ.Pred_D0_ncf
#' data_temp$YgivenGXQ.Pred_D1_ncf <- YgivenGXQ.Pred_D1_ncf
#'
#' message <- utils::capture.output( Y0givenGQ.Model.sample1 <- caret::train(stats::as.formula(paste("YgivenGXQ.Pred_D0_ncf", paste(G,Q,sep="+"), sep="~")), data=data_temp[sample1,], method="ranger",
#'                                                                           trControl=caret::trainControl(method="cv"),
#'                                                                           tuneGrid=expand.grid(mtry=1,splitrule=c("variance"),min.node.size=c(25,50))) )
#' message <- utils::capture.output( Y0givenGQ.Model.sample2 <- caret::train(stats::as.formula(paste("YgivenGXQ.Pred_D0_ncf", paste(G,Q,sep="+"), sep="~")), data=data_temp[sample2,], method="ranger",
#'                                                                           trControl=caret::trainControl(method="cv"),
#'                                                                           tuneGrid=expand.grid(mtry=1,splitrule=c("variance"),min.node.size=c(25,50))) )
#' message <- utils::capture.output( Y1givenGQ.Model.sample1 <- caret::train(stats::as.formula(paste("YgivenGXQ.Pred_D1_ncf", paste(G,Q,sep="+"), sep="~")), data=data_temp[sample1,], method="ranger",
#'                                                                           trControl=caret::trainControl(method="cv"),
#'                                                                           tuneGrid=expand.grid(mtry=1,splitrule=c("variance"),min.node.size=c(25,50))) )
#' message <- utils::capture.output( Y1givenGQ.Model.sample2 <- caret::train(stats::as.formula(paste("YgivenGXQ.Pred_D1_ncf", paste(G,Q,sep="+"), sep="~")), data=data_temp[sample2,], method="ranger",
#'                                                                           trControl=caret::trainControl(method="cv"),
#'                                                                           tuneGrid=expand.grid(mtry=1,splitrule=c("variance"),min.node.size=c(25,50))) )
#'
#' Y0givenQ.Pred_G0 <- Y0givenQ.Pred_G1 <- Y1givenQ.Pred_G0 <- Y1givenQ.Pred_G1 <- rep(NA, nrow(data))
#'
#' pred_data <- data
#' pred_data[,G] <- 1
#' Y0givenQ.Pred_G1[sample1] <- stats::predict(Y0givenGQ.Model.sample1, newdata = pred_data[sample2,])  # cross-fitting is used
#' Y0givenQ.Pred_G1[sample2] <- stats::predict(Y0givenGQ.Model.sample2, newdata = pred_data[sample1,])
#' Y1givenQ.Pred_G1[sample1] <- stats::predict(Y1givenGQ.Model.sample1, newdata = pred_data[sample2,])
#' Y1givenQ.Pred_G1[sample2] <- stats::predict(Y1givenGQ.Model.sample2, newdata = pred_data[sample1,])
#'
#' pred_data <- data
#' pred_data[,G] <- 0
#' Y0givenQ.Pred_G0[sample1] <- stats::predict(Y0givenGQ.Model.sample1, newdata = pred_data[sample2,])  # cross-fitting is used
#' Y0givenQ.Pred_G0[sample2] <- stats::predict(Y0givenGQ.Model.sample2, newdata = pred_data[sample1,])
#' Y1givenQ.Pred_G0[sample1] <- stats::predict(Y1givenGQ.Model.sample1, newdata = pred_data[sample2,])
#' Y1givenQ.Pred_G0[sample2] <- stats::predict(Y1givenGQ.Model.sample2, newdata = pred_data[sample1,])
#'
#' ### Estimate E(D | Q,g')
#' data[,D] <- as.factor(data[,D])
#' levels(data[,D]) <- c("D0","D1")  # necessary for caret implementation of ranger
#'
#' message <- utils::capture.output( DgivenGQ.Model.sample1 <- caret::train(stats::as.formula(paste(D, paste(G,Q,sep="+"), sep="~")), data=data[sample1,], method="ranger",
#'                                                                          trControl=caret::trainControl(method="cv", classProbs=TRUE),
#'                                                                          tuneGrid=expand.grid(mtry=1,splitrule=c("gini"),min.node.size=c(25,50))) )
#' message <- utils::capture.output( DgivenGQ.Model.sample2 <- caret::train(stats::as.formula(paste(D, paste(G,Q,sep="+"), sep="~")), data=data[sample2,], method="ranger",
#'                                                                          trControl=caret::trainControl(method="cv", classProbs=TRUE),
#'                                                                          tuneGrid=expand.grid(mtry=1,splitrule=c("gini"),min.node.size=c(25,50))) )
#' data[,D] <- as.numeric(data[,D])-1
#'
#' DgivenQ.Pred_G0 <- DgivenQ.Pred_G1 <- rep(NA, nrow(data))
#'
#' pred_data <- data
#' pred_data[,G] <- 0
#' DgivenQ.Pred_G0[sample1] <- stats::predict(DgivenGQ.Model.sample1, newdata = pred_data[sample2,], type="prob")[,2]
#' DgivenQ.Pred_G0[sample2] <- stats::predict(DgivenGQ.Model.sample2, newdata = pred_data[sample1,], type="prob")[,2]
#'
#' pred_data <- data
#' pred_data[,G] <- 1
#' DgivenQ.Pred_G1[sample1] <- stats::predict(DgivenGQ.Model.sample1, newdata = pred_data[sample2,], type="prob")[,2]
#' DgivenQ.Pred_G1[sample2] <- stats::predict(DgivenGQ.Model.sample2, newdata = pred_data[sample1,], type="prob")[,2]
#'
#' ### Estimate p_g(Q)=Pr(G=g | Q)
#' data[,G] <- as.factor(data[,G])
#' levels(data[,G]) <- c("G0","G1")  # necessary for caret implementation of ranger
#'
#' message <- utils::capture.output( GgivenQ.Model.sample1 <- caret::train(stats::as.formula(paste(G, paste(Q,sep="+"), sep="~")), data=data[sample1,], method="ranger",
#'                                                                         trControl=caret::trainControl(method="cv", classProbs=TRUE),
#'                                                                         tuneGrid=expand.grid(mtry=1,splitrule=c("gini"),min.node.size=c(25,50))) )
#' message <- utils::capture.output( GgivenQ.Model.sample2 <- caret::train(stats::as.formula(paste(G, paste(Q,sep="+"), sep="~")), data=data[sample2,], method="ranger",
#'                                                                         trControl=caret::trainControl(method="cv", classProbs=TRUE),
#'                                                                         tuneGrid=expand.grid(mtry=1,splitrule=c("gini"),min.node.size=c(25,50))) )
#' data[,G] <- as.numeric(data[,G])-1
#'
#' GgivenQ.Pred <- rep(NA, nrow(data))
#' GgivenQ.Pred[sample1] <- stats::predict(DgivenGQ.Model.sample1, newdata = data[sample2,], type="prob")[,2]
#' GgivenQ.Pred[sample2] <- stats::predict(DgivenGQ.Model.sample2, newdata = data[sample1,], type="prob")[,2]
#'
#'
#' results <- cdgd1_manual(Y=Y,D=D,G=G,
#'                         YgivenXQ.Pred_D0G0=YgivenXQ.Pred_D0G0,
#'                         YgivenXQ.Pred_D1G0=YgivenXQ.Pred_D1G0,
#'                         YgivenXQ.Pred_D0G1=YgivenXQ.Pred_D0G1,
#'                         YgivenXQ.Pred_D1G1=YgivenXQ.Pred_D1G1,
#'                         DgivenXQ.Pred_G0=DgivenXQ.Pred_G0,
#'                         DgivenXQ.Pred_G1=DgivenXQ.Pred_G1,
#'                         Y0givenQ.Pred_G0=Y0givenQ.Pred_G0,
#'                         Y0givenQ.Pred_G1=Y0givenQ.Pred_G1,
#'                         Y1givenQ.Pred_G0=Y1givenQ.Pred_G0,
#'                         Y1givenQ.Pred_G1=Y1givenQ.Pred_G1,
#'                         DgivenQ.Pred_G0=DgivenQ.Pred_G0,
#'                         DgivenQ.Pred_G1=DgivenQ.Pred_G1,
#'                         GgivenQ.Pred=GgivenQ.Pred,
#'                         data,alpha=0.05)
#'
#' results


cdgd1_manual <- function(Y,D,G,
                     YgivenXQ.Pred_D0G0,YgivenXQ.Pred_D1G0,YgivenXQ.Pred_D0G1,YgivenXQ.Pred_D1G1,DgivenXQ.Pred_G0,DgivenXQ.Pred_G1,
                     Y0givenQ.Pred_G0,Y0givenQ.Pred_G1,Y1givenQ.Pred_G0,Y1givenQ.Pred_G1,DgivenQ.Pred_G0,DgivenQ.Pred_G1,GgivenQ.Pred,
                     data,alpha=0.05) {

  data <- as.data.frame(data)

  zero_one <- sum(DgivenXQ.Pred_G0==0)+sum(DgivenXQ.Pred_G1==0)+
    sum(DgivenXQ.Pred_G0==1)+sum(DgivenXQ.Pred_G1==1)
  if ( zero_one>0 ) {
    stop(
      paste("D given X and are exact 0 or 1 in", zero_one, "cases.", sep=" "),
      call. = FALSE
    )
  }

  zero_one <- sum(GgivenQ.Pred==0)+sum(GgivenQ.Pred==1)
  if ( zero_one>0 ) {
    stop(
      paste("G given Q are exact 0 or 1 in", zero_one, "cases.", sep=" "),
      call. = FALSE
    )
  }

  ### The "IPO" (individual potential outcome) function
  # For each d and g value, we have IE(d,g)=\frac{\one(D=d)}{\pi(d,X,g)}[Y-\mu(d,X,g)]+\mu(d,X,g)
  # We stablize the weight by dividing the sample average of estimated weights

  IPO_D0G0 <- (1-data[,D])/(1-DgivenXQ.Pred_G0)/mean((1-data[,D])/(1-DgivenXQ.Pred_G0))*(data[,Y]-YgivenXQ.Pred_D0G0) + YgivenXQ.Pred_D0G0
  IPO_D1G0 <- data[,D]/DgivenXQ.Pred_G0/(mean(data[,D]/DgivenXQ.Pred_G0))*(data[,Y]-YgivenXQ.Pred_D1G0) + YgivenXQ.Pred_D1G0
  IPO_D0G1 <- (1-data[,D])/(1-DgivenXQ.Pred_G0)/(mean((1-data[,D])/(1-DgivenXQ.Pred_G0)))*(data[,Y]-YgivenXQ.Pred_D0G1) + YgivenXQ.Pred_D0G1
  IPO_D1G1 <- data[,D]/DgivenXQ.Pred_G0/mean(data[,D]/DgivenXQ.Pred_G0)*(data[,Y]-YgivenXQ.Pred_D1G1) + YgivenXQ.Pred_D1G1

  ### The one-step estimate of \xi_{dg}
  psi_00 <- mean( (1-data[,G])/(1-mean(data[,G]))*IPO_D0G0 )
  psi_01 <- mean( data[,G]/mean(data[,G])*IPO_D0G1 )
  # Note that this is basically DML2. We could also use DML1:
  #psi_00_S1 <- mean( (1-data[sample1,G])/(1-mean(data[sample1,G]))*IPO_D0G0[sample1] )     # sample 1 estimate
  #psi_00_S2 <- mean( (1-data[sample2,G])/(1-mean(data[sample2,G]))*IPO_D0G0[sample2] )     # sample 2 estimate
  #psi_00 <- (1/2)*(psi_00_S1+psi_00_S2)
  #psi_01_S1 <- mean( data[sample1,G]/mean(data[sample1,G])*IPO_D0G1[sample1] )     # sample 1 estimate
  #psi_01_S2 <- mean( data[sample1,G]/mean(data[sample1,G])*IPO_D0G1[sample2] )     # sample 2 estimate
  #psi_01 <- (1/2)*(psi_01_S1+psi_01_S2)

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
    # Note that this is basically DML2. We could also use DML1:
    #psi_dggg_S1 <- mean( as.numeric(data[sample1,G]==g3)/mean(data[sample1,G]==g3)*YdgivenQ.Pred_arg[sample1]*DgivenQ.Pred_arg[sample1] +
    #                       as.numeric(data[sample1,G]==g1)/mean(data[sample1,G]==g3)*g3givenQ.Pred_arg[sample1]/g1givenQ.Pred_arg[sample1]*(IPO_arg[sample1]-YdgivenQ.Pred_arg[sample1])*DgivenQ.Pred_arg[sample1] +
    #                       as.numeric(data[sample1,G]==g2)/mean(data[sample1,G]==g3)*g3givenQ.Pred_arg[sample1]/g2givenQ.Pred_arg[sample1]*(data[sample1,D]-DgivenQ.Pred_arg[sample1])*YdgivenQ.Pred_arg[sample1] )
    #psi_dggg_S2 <- mean( as.numeric(data[sample2,G]==g3)/mean(data[sample2,G]==g3)*YdgivenQ.Pred_arg[sample2]*DgivenQ.Pred_arg[sample2] +
    #                       as.numeric(data[sample2,G]==g1)/mean(data[sample2,G]==g3)*g3givenQ.Pred_arg[sample2]/g1givenQ.Pred_arg[sample2]*(IPO_arg[sample2]-YdgivenQ.Pred_arg[sample2])*DgivenQ.Pred_arg[sample2] +
    #                       as.numeric(data[sample2,G]==g2)/mean(data[sample2,G]==g3)*g3givenQ.Pred_arg[sample2]/g2givenQ.Pred_arg[sample2]*(data[sample2,D]-DgivenQ.Pred_arg[sample2])*YdgivenQ.Pred_arg[sample2] )
    #psi_dggg <- (1/2)*(psi_dggg_S1+psi_dggg_S2)

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

  cond_Jackson_reduction <- psi_00+psi_dggg(1,0,1,0)-psi_dggg(0,0,1,0)-mean((1-data[,G])/(1-mean(data[,G]))*data[,Y])

  ### standard error estimates
  se <- function(x) {sqrt( mean(x^2)/nrow(data) )}
  total_se <- se( data[,G]/mean(data[,G])*(data[,Y]-Y_G1) - (1-data[,G])/(1-mean(data[,G]))*(data[,Y]-Y_G0) )
  baseline_se <- se( data[,G]/mean(data[,G])*(IPO_D0G1-psi_01) - (1-data[,G])/(1-mean(data[,G]))*(IPO_D0G0-psi_00) )
  # Alternatively, we could use
  # se( c( data[sample1,G]/mean(data[sample1,G])*(IPO_D0G1[sample1]-psi_01) - (1-data[sample1,G])/(1-mean(data[sample1,G]))*(IPO_D0G0[sample1]-psi_00),
  #         data[sample2,G]/mean(data[sample2,G])*(IPO_D0G1[sample2]-psi_01) - (1-data[sample2,G])/(1-mean(data[sample2,G]))*(IPO_D0G0[sample2]-psi_00) ) )
  # But there isn't a theoretically strong reason to prefer one over the other.

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

  cond_Jackson_reduction_se <- se( (1-data[,G])/(1-mean(data[,G]))*(IPO_D0G0-psi_00)+EIF_dggg(1,0,1,0)-EIF_dggg(0,0,1,0)-(1-data[,G])/(1-mean(data[,G]))*(data[,Y]-Y_G0) )

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
  CI_lower <- point - stats::qnorm(1-alpha/2)*se
  CI_upper <- point + stats::qnorm(1-alpha/2)*se
  names <- c("total",
             "baseline",
             "conditional prevalence",
             "conditional effect",
             "conditional selection",
             "Q distribution",
             "conditional Jackson reduction")

  output <- as.data.frame(cbind(names, point,se,CI_lower,CI_upper))

  return(output)
}