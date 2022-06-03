
#' Only produce point estimates of cdgd
#'
#' @param Y Outcome. The name of a continuous variable.
#' @param D Treatment status. The name of a binary numeric variable taking values of 0 and 1.
#' @param G1 Group 1 membership. The name of a binary factor variable taking values of 0 and 1.
#' @param G2 Group 2 membership. The name of a binary factor variable taking values of 0 and 1.
#' @param Q Conditional set. The vector of the names of numeric variables.
#' @param X Confounders. The vector of the names of numeric variables.
#' @param data A data frame.
#' @param weight Survey weights. The name of a numeric variable.
#' @param t Threshold of propensity score censoring. Propensity scores larger than 1-t or smaller than t will be censored.
#' @param algorithm The ML algorithm for modelling. "nnet" for neural network and "ranger" for random forests.
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
#' X=c("confounder","Q"),
#' Q="Q",
#' data=exp_data,
#' t=0.05,
#' algorithm="nnet")
#'
#' results0



cdgd0 <- function(Y,D,G1,G2,Q=NULL,X,data,weight=NULL,t=0.05,algorithm) {

  if (is.null(weight)) {
    data$weight=rep(1, nrow(data))
    weight <- "weight"
  }

  data <- as.data.frame(data)

  # estimate the nuisance functions within each group, so that the final estimates are independent across groups
  G1_index <- data[,G1]==1
  G2_index <- data[,G2]==1

  mainsample_G1 <- sample(nrow(data[G1_index,]), floor(nrow(data[G1_index,])/2), replace=F)
  auxsample_G1 <- setdiff(1:nrow(data[G1_index,]), mainsample_G1)
  mainsample_G2 <- sample(nrow(data[G2_index,]), floor(nrow(data[G2_index,])/2), replace=F)
  auxsample_G2 <- setdiff(1:nrow(data[G2_index,]), mainsample_G2)

  YgivenX.Pred_D0 <- YgivenX.Pred_D1 <- DgivenX.Pred <- rep(NA, nrow(data))

  if (algorithm=="nnet") {
    message <- utils::capture.output( YgivenX.Model.Aux_G1 <- caret::train(stats::as.formula(paste(Y, paste(D,paste(X,collapse="+"),sep="+"), sep="~")), data=data[G1_index,][mainsample_G1,], method="nnet",
                                                                           preProc=c("center","scale"), trControl=caret::trainControl(method="none"), linout=TRUE,
                                                                           tuneGrid=expand.grid(size=2,decay=0.02)) )
    message <- utils::capture.output( YgivenX.Model.Main_G1 <- caret::train(stats::as.formula(paste(Y, paste(D,paste(X,collapse="+"),sep="+"), sep="~")), data=data[G1_index,][auxsample_G1,], method="nnet",
                                                                            preProc=c("center","scale"), trControl=caret::trainControl(method="none"), linout=TRUE,
                                                                            tuneGrid=expand.grid(size=2,decay=0.02)) )
    message <- utils::capture.output( YgivenX.Model.Aux_G2 <- caret::train(stats::as.formula(paste(Y, paste(D,paste(X,collapse="+"),sep="+"), sep="~")), data=data[G2_index,][mainsample_G2,], method="nnet",
                                                                           preProc=c("center","scale"), trControl=caret::trainControl(method="none"), linout=TRUE,
                                                                           tuneGrid=expand.grid(size=2,decay=0.02)) )
    message <- utils::capture.output( YgivenX.Model.Main_G2 <- caret::train(stats::as.formula(paste(Y, paste(D,paste(X,collapse="+"),sep="+"), sep="~")), data=data[G2_index,][auxsample_G2,], method="nnet",
                                                                            preProc=c("center","scale"), trControl=caret::trainControl(method="none"), linout=TRUE,
                                                                            tuneGrid=expand.grid(size=2,decay=0.02)) )
  }
  if (algorithm=="ranger") {
    message <- utils::capture.output( YgivenX.Model.Aux_G1 <- caret::train(stats::as.formula(paste(Y, paste(D,paste(X,collapse="+"),sep="+"), sep="~")), data=data[G1_index,][mainsample_G1,], method="ranger",
                                                                           trControl=caret::trainControl(method="cv"),
                                                                           tuneGrid=expand.grid(mtry=floor(sqrt(length(X))),splitrule="variance",min.node.size=c(5,10,100))) )
    message <- utils::capture.output( YgivenX.Model.Main_G1 <- caret::train(stats::as.formula(paste(Y, paste(D,paste(X,collapse="+"),sep="+"), sep="~")), data=data[G1_index,][auxsample_G1,], method="ranger",
                                                                            trControl=caret::trainControl(method="cv"),
                                                                            tuneGrid=expand.grid(mtry=floor(sqrt(length(X))),splitrule="variance",min.node.size=c(5,10,100))) )
    message <- utils::capture.output( YgivenX.Model.Aux_G2 <- caret::train(stats::as.formula(paste(Y, paste(D,paste(X,collapse="+"),sep="+"), sep="~")), data=data[G2_index,][mainsample_G2,], method="ranger",
                                                                           trControl=caret::trainControl(method="cv"),
                                                                           tuneGrid=expand.grid(mtry=floor(sqrt(length(X))),splitrule="variance",min.node.size=c(5,10,100))) )
    message <- utils::capture.output( YgivenX.Model.Main_G2 <- caret::train(stats::as.formula(paste(Y, paste(D,paste(X,collapse="+"),sep="+"), sep="~")), data=data[G2_index,][auxsample_G2,], method="ranger",
                                                                            trControl=caret::trainControl(method="cv"),
                                                                            tuneGrid=expand.grid(mtry=floor(sqrt(length(X))),splitrule="variance",min.node.size=c(5,10,100))) )
  }


  pred_data <- data
  pred_data[,colnames(pred_data)%in%D] <- 0
  YgivenX.Pred_D0[G1_index][mainsample_G1] <- stats::predict(YgivenX.Model.Aux_G1, newdata = pred_data[G1_index,][mainsample_G1,])
  YgivenX.Pred_D0[G1_index][auxsample_G1] <- stats::predict(YgivenX.Model.Main_G1, newdata = pred_data[G1_index,][auxsample_G1,])
  YgivenX.Pred_D0[G2_index][mainsample_G2] <- stats::predict(YgivenX.Model.Aux_G2, newdata = pred_data[G2_index,][mainsample_G2,])
  YgivenX.Pred_D0[G2_index][auxsample_G2] <- stats::predict(YgivenX.Model.Main_G2, newdata = pred_data[G2_index,][auxsample_G2,])

  pred_data <- data
  pred_data[,colnames(pred_data)%in%D] <- 1
  YgivenX.Pred_D1[G1_index][mainsample_G1] <- stats::predict(YgivenX.Model.Aux_G1, newdata = pred_data[G1_index,][mainsample_G1,])
  YgivenX.Pred_D1[G1_index][auxsample_G1] <- stats::predict(YgivenX.Model.Main_G1, newdata = pred_data[G1_index,][auxsample_G1,])
  YgivenX.Pred_D1[G2_index][mainsample_G2] <- stats::predict(YgivenX.Model.Aux_G2, newdata = pred_data[G2_index,][mainsample_G2,])
  YgivenX.Pred_D1[G2_index][auxsample_G2] <- stats::predict(YgivenX.Model.Main_G2, newdata = pred_data[G2_index,][auxsample_G2,])


  data[,D] <- as.factor(data[,D])
  levels(data[,D]) <- c("D0","D1")  # necessary for caret implementation of ranger

  if (algorithm=="nnet") {
    message <- utils::capture.output( DgivenX.Model.Aux_G1 <- caret::train(stats::as.formula(paste(D, paste(X,collapse="+"), sep="~")), data=data[G1_index,][auxsample_G1,], method="nnet",
                                                                           preProc=c("center","scale"), trControl=caret::trainControl(method="none"), linout=FALSE,
                                                                           tuneGrid=expand.grid(size=2,decay=0.02)) )
    message <- utils::capture.output( DgivenX.Model.Main_G1 <- caret::train(stats::as.formula(paste(D, paste(X,collapse="+"), sep="~")), data=data[G1_index,][mainsample_G1,], method="nnet",
                                                                            preProc=c("center","scale"), trControl=caret::trainControl(method="none"), linout=FALSE,
                                                                            tuneGrid=expand.grid(size=2,decay=0.02)) )
    message <- utils::capture.output( DgivenX.Model.Aux_G2 <- caret::train(stats::as.formula(paste(D, paste(X,collapse="+"), sep="~")), data=data[G2_index,][auxsample_G2,], method="nnet",
                                                                           preProc=c("center","scale"), trControl=caret::trainControl(method="none"), linout=FALSE,
                                                                           tuneGrid=expand.grid(size=2,decay=0.02)) )
    message <- utils::capture.output( DgivenX.Model.Main_G2 <- caret::train(stats::as.formula(paste(D, paste(X,collapse="+"), sep="~")), data=data[G2_index,][mainsample_G2,], method="nnet",
                                                                            preProc=c("center","scale"), trControl=caret::trainControl(method="none"), linout=FALSE,
                                                                            tuneGrid=expand.grid(size=2,decay=0.02)) )
  }
  if (algorithm=="ranger") {
    message <- utils::capture.output( DgivenX.Model.Aux_G1 <- caret::train(stats::as.formula(paste(D, paste(X,collapse="+"), sep="~")), data=data[G1_index,][auxsample_G1,], method="ranger",
                                                                           trControl=caret::trainControl(method="cv", classProbs=TRUE),
                                                                           tuneGrid=expand.grid(mtry=floor(sqrt(length(X))),splitrule="gini",min.node.size=c(1,10,100))) )
    message <- utils::capture.output( DgivenX.Model.Main_G1 <- caret::train(stats::as.formula(paste(D, paste(X,collapse="+"), sep="~")), data=data[G1_index,][mainsample_G1,], method="ranger",
                                                                            trControl=caret::trainControl(method="cv", classProbs=TRUE),
                                                                            tuneGrid=expand.grid(mtry=floor(sqrt(length(X))),splitrule="gini",min.node.size=c(1,10,100))) )
    message <- utils::capture.output( DgivenX.Model.Aux_G2 <- caret::train(stats::as.formula(paste(D, paste(X,collapse="+"), sep="~")), data=data[G2_index,][auxsample_G2,], method="ranger",
                                                                           trControl=caret::trainControl(method="cv", classProbs=TRUE),
                                                                           tuneGrid=expand.grid(mtry=floor(sqrt(length(X))),splitrule="gini",min.node.size=c(1,10,100))) )
    message <- utils::capture.output( DgivenX.Model.Main_G2 <- caret::train(stats::as.formula(paste(D, paste(X,collapse="+"), sep="~")), data=data[G2_index,][mainsample_G2,], method="ranger",
                                                                            trControl=caret::trainControl(method="cv", classProbs=TRUE),
                                                                            tuneGrid=expand.grid(mtry=floor(sqrt(length(X))),splitrule="gini",min.node.size=c(1,10,100))) )
  }

  DgivenX.Pred[G1_index][mainsample_G1] <- stats::predict(DgivenX.Model.Aux_G1, newdata=data[G1_index,][mainsample_G1,], type="prob")[,2]
  DgivenX.Pred[G1_index][auxsample_G1] <- stats::predict(DgivenX.Model.Main_G1, newdata=data[G1_index,][auxsample_G1,], type="prob")[,2]
  DgivenX.Pred[G2_index][mainsample_G2] <- stats::predict(DgivenX.Model.Aux_G2, newdata=data[G2_index,][mainsample_G2,], type="prob")[,2]
  DgivenX.Pred[G2_index][auxsample_G2] <- stats::predict(DgivenX.Model.Main_G2, newdata=data[G2_index,][auxsample_G2,], type="prob")[,2]

  data[,D] <- as.numeric(data[,D])-1

  DgivenX.Pred[DgivenX.Pred<=t] <- t
  DgivenX.Pred[DgivenX.Pred>=1-t] <- 1-t

  Y0_i <- ATT_i <- ATE_i <- wht <- rep(NA, nrow(data))

  wht[G1_index] <- data[,weight][G1_index]/mean(data[,weight][G1_index])
  wht[G2_index] <- data[,weight][G2_index]/mean(data[,weight][G2_index])

  Y0_i[G1_index] <- ( YgivenX.Pred_D0 + (1-data[,D])*(data[,Y]-YgivenX.Pred_D0)/(1-DgivenX.Pred) )[G1_index]*wht[G1_index]
  Y0_i[G2_index] <- ( YgivenX.Pred_D0 + (1-data[,D])*(data[,Y]-YgivenX.Pred_D0)/(1-DgivenX.Pred) )[G2_index]*wht[G2_index]

  ATT_i[G1_index] <- ( (data[,D]-(1-data[,D])*DgivenX.Pred/(1-DgivenX.Pred))*(data[,Y]-YgivenX.Pred_D0) )[G1_index]/mean(data[,D][G1_index]*wht[G1_index])*wht[G1_index]
  ATT_i[G2_index] <- ( (data[,D]-(1-data[,D])*DgivenX.Pred/(1-DgivenX.Pred))*(data[,Y]-YgivenX.Pred_D0) )[G2_index]/mean(data[,D][G2_index]*wht[G2_index])*wht[G2_index]

  ATE_i[G1_index] <- ( YgivenX.Pred_D1 + data[,D]*(data[,Y]-YgivenX.Pred_D1)/DgivenX.Pred - ( YgivenX.Pred_D0 + (1-data[,D])*(data[,Y]-YgivenX.Pred_D0)/(1-DgivenX.Pred) ) )[G1_index]*wht[G1_index]
  ATE_i[G2_index] <- ( YgivenX.Pred_D1 + data[,D]*(data[,Y]-YgivenX.Pred_D1)/DgivenX.Pred - ( YgivenX.Pred_D0 + (1-data[,D])*(data[,Y]-YgivenX.Pred_D0)/(1-DgivenX.Pred) ) )[G2_index]*wht[G2_index]

  total <- mean(data[,Y][G1_index]*wht[G1_index])-mean(data[,Y][G2_index]*wht[G2_index])
  baseline <- mean(Y0_i[G1_index])-mean(Y0_i[G2_index])
  prevalence <- mean(ATE_i[G2_index])*(mean(data[,D][G1_index]*wht[G1_index])-mean(data[,D][G2_index]*wht[G2_index]))
  effect <- mean(data[,D][G1_index]*wht[G1_index])*(mean(ATE_i[G1_index])-mean(ATE_i[G2_index]))
  selection <- (mean(ATT_i[G1_index])-mean(ATE_i[G1_index]))*mean(data[,D][G1_index]*wht[G1_index])-
    (mean(ATT_i[G2_index])-mean(ATE_i[G2_index]))*mean(data[,D][G2_index]*wht[G2_index])

  ### conditional prevalence ###
  if (!is.null(Q)) {
    data_cond <- cbind(ATE_i, data[,D], data[,Q])
    data_cond <- as.data.frame(data_cond)
    colnames(data_cond)[1:2] <- c("tau","D")
    Q_names <- colnames(data_cond)[3:ncol(data_cond)]
    data_cond$D <- as.factor(data_cond$D)
    levels(data_cond$D) <- c("D0","D1")  # necessary for caret implementation of ranger

    TaugivenQ.Pred_G1_G1 <- TaugivenQ.Pred_G2_G1 <- DgivenQ.Pred_G1_G1 <- DgivenQ.Pred_G2_G1 <- rep(NA, sum(G1_index))
    TaugivenQ.Pred_G1_G2 <- TaugivenQ.Pred_G2_G2 <- DgivenQ.Pred_G1_G2 <- DgivenQ.Pred_G2_G2 <- rep(NA, sum(G2_index))

    if (algorithm=="nnet") {
      message <- utils::capture.output( TaugivenQ.Model_G1 <- caret::train(stats::as.formula(paste("tau", paste(Q_names,sep="+"), sep="~")), data=data_cond[G1_index,], method="nnet",
                                                                           preProc=c("center","scale"), trControl=caret::trainControl(method="none"), linout=TRUE,
                                                                           tuneGrid=expand.grid(size=2,decay=0.02)) )

      message <- utils::capture.output( TaugivenQ.Model_G2 <- caret::train(stats::as.formula(paste("tau", paste(Q_names,sep="+"), sep="~")), data=data_cond[G2_index,], method="nnet",
                                                                           preProc=c("center","scale"), trControl=caret::trainControl(method="none"), linout=TRUE,
                                                                           tuneGrid=expand.grid(size=2,decay=0.02)) )
    }
    if (algorithm=="ranger") {
      message <- utils::capture.output( TaugivenQ.Model_G1 <- caret::train(stats::as.formula(paste("tau", paste(Q_names,sep="+"), sep="~")), data=data_cond[G1_index,], method="ranger",
                                                                           trControl=caret::trainControl(method="cv"),
                                                                           tuneGrid=expand.grid(mtry=floor(sqrt(length(Q_names))),splitrule="variance",min.node.size=c(5,10,100))) )

      message <- utils::capture.output( TaugivenQ.Model_G2 <- caret::train(stats::as.formula(paste("tau", paste(Q_names,sep="+"), sep="~")), data=data_cond[G2_index,], method="ranger",
                                                                           trControl=caret::trainControl(method="cv"),
                                                                           tuneGrid=expand.grid(mtry=floor(sqrt(length(Q_names))),splitrule="variance",min.node.size=c(5,10,100))) )
    }


    TaugivenQ.Pred_G1_G1 <- stats::predict(TaugivenQ.Model_G1, newdata = data_cond[G1_index,])
    TaugivenQ.Pred_G2_G1 <- stats::predict(TaugivenQ.Model_G2, newdata = data_cond[G1_index,])
    TaugivenQ.Pred_G1_G2 <- stats::predict(TaugivenQ.Model_G1, newdata = data_cond[G2_index,])
    TaugivenQ.Pred_G2_G2 <- stats::predict(TaugivenQ.Model_G2, newdata = data_cond[G2_index,])

    #plotLowess(data_cond[G2_index,]$tau ~ data_cond[G2_index,]$V3)
    #plot(data_cond[G2_index,]$V3, TaugivenQ.Pred_G2_G2)

    if (algorithm=="nnet") {
      message <- utils::capture.output( DaugivenQ.Model_G1 <- caret::train(stats::as.formula(paste("D", paste(Q_names,sep="+"), sep="~")), data=data_cond[G1_index,], method="nnet",
                                                                           preProc=c("center","scale"), trControl=caret::trainControl(method="none"), linout=FALSE,
                                                                           tuneGrid=expand.grid(size=2,decay=0.02)), weights=wht[G1_index] )

      message <- utils::capture.output( DaugivenQ.Model_G2 <- caret::train(stats::as.formula(paste("D", paste(Q_names,sep="+"), sep="~")), data=data_cond[G2_index,], method="nnet",
                                                                           preProc=c("center","scale"), trControl=caret::trainControl(method="none"), linout=FALSE,
                                                                           tuneGrid=expand.grid(size=2,decay=0.02)), weights=wht[G2_index] )
    }
    if (algorithm=="ranger") {
      message <- utils::capture.output( DaugivenQ.Model_G1 <- caret::train(stats::as.formula(paste("D", paste(Q_names,sep="+"), sep="~")), data=data_cond[G1_index,], method="ranger",
                                                                           trControl=caret::trainControl(method="cv", classProbs=TRUE),
                                                                           tuneGrid=expand.grid(mtry=floor(sqrt(length(Q_names))),splitrule="gini",min.node.size=c(1,10,100))), weights=wht[G1_index] )

      message <- utils::capture.output( DaugivenQ.Model_G2 <- caret::train(stats::as.formula(paste("D", paste(Q_names,sep="+"), sep="~")), data=data_cond[G2_index,], method="ranger",
                                                                           trControl=caret::trainControl(method="cv", classProbs=TRUE),
                                                                           tuneGrid=expand.grid(mtry=floor(sqrt(length(Q_names))),splitrule="gini",min.node.size=c(1,10,100))), weights=wht[G2_index] )
    }


    DgivenQ.Pred_G1_G1 <- stats::predict(DaugivenQ.Model_G1, newdata = data_cond[G1_index,], type="prob")[,2]
    DgivenQ.Pred_G2_G1 <- stats::predict(DaugivenQ.Model_G2, newdata = data_cond[G1_index,], type="prob")[,2]
    DgivenQ.Pred_G1_G2 <- stats::predict(DaugivenQ.Model_G1, newdata = data_cond[G2_index,], type="prob")[,2]
    DgivenQ.Pred_G2_G2 <- stats::predict(DaugivenQ.Model_G2, newdata = data_cond[G2_index,], type="prob")[,2]

    cond_prevalence <- mean((DgivenQ.Pred_G1_G2-DgivenQ.Pred_G2_G2)*TaugivenQ.Pred_G2_G2*wht[G2_index])
    cond_effect <- mean((TaugivenQ.Pred_G1_G1-TaugivenQ.Pred_G2_G1)*DgivenQ.Pred_G1_G1*wht[G1_index])
    Q_dist <- mean(DgivenQ.Pred_G1_G1*TaugivenQ.Pred_G2_G1*wht[G1_index]) - mean(DgivenQ.Pred_G1_G2*TaugivenQ.Pred_G2_G2*wht[G2_index])
    cond_selection <- total-baseline-cond_prevalence-cond_effect-Q_dist
  }

  ###  ###

  if (is.null(Q)) {
    point <- c(
      mean(data[,Y][G1_index]*wht[G1_index]),
      mean(data[,Y][G2_index]*wht[G2_index]),
      mean(Y0_i[G1_index]),
      mean(Y0_i[G2_index]),
      mean(data[,D][G1_index]*wht[G1_index]),
      mean(data[,D][G2_index]*wht[G2_index]),
      mean(ATE_i[G1_index]),
      mean(ATE_i[G2_index]),
      mean(ATT_i[G1_index]),
      mean(ATT_i[G2_index]),
      total,
      baseline,
      mean(data[,D][G1_index]*wht[G1_index])-mean(data[,D][G2_index]*wht[G2_index]),
      mean(ATE_i[G1_index])-mean(ATE_i[G2_index]),
      mean(ATT_i[G1_index])-mean(ATT_i[G2_index]),
      total,
      baseline,
      prevalence,
      effect,
      selection
    )

    item <- c(
      "mean_Y_G1",
      "mean_Y_G2",
      "mean_Y0_G1",
      "mean_Y0_G2",
      "mean_D_G1",
      "mean_D_G2",
      "ATE_G1",
      "ATE_G2",
      "ATT_G1",
      "ATT_G2",
      "diff_mean_Y",
      "diff_mean_Y0",
      "diff_mean_D",
      "diff_ATE",
      "dff_ATT",
      "total",
      "baseline",
      "prevalence",
      "effect",
      "selection")
  } else {
    point <- c(
      mean(data[,Y][G1_index]*wht[G1_index]),
      mean(data[,Y][G2_index]*wht[G2_index]),
      mean(Y0_i[G1_index]),
      mean(Y0_i[G2_index]),
      mean(data[,D][G1_index]*wht[G1_index]),
      mean(data[,D][G2_index]*wht[G2_index]),
      mean(ATE_i[G1_index]),
      mean(ATE_i[G2_index]),
      mean(ATT_i[G1_index]),
      mean(ATT_i[G2_index]),
      total,
      baseline,
      mean(data[,D][G1_index]*wht[G1_index])-mean(data[,D][G2_index]*wht[G2_index]),
      mean(ATE_i[G1_index])-mean(ATE_i[G2_index]),
      mean(ATT_i[G1_index])-mean(ATT_i[G2_index]),
      total,
      baseline,
      prevalence,
      effect,
      selection,
      cond_prevalence,
      cond_effect,
      cond_selection,
      Q_dist
    )

    item <- c(
      "mean_Y_G1",
      "mean_Y_G2",
      "mean_Y0_G1",
      "mean_Y0_G2",
      "mean_D_G1",
      "mean_D_G2",
      "ATE_G1",
      "ATE_G2",
      "ATT_G1",
      "ATT_G2",
      "diff_mean_Y",
      "diff_mean_Y0",
      "diff_mean_D",
      "diff_ATE",
      "dff_ATT",
      "total",
      "baseline",
      "prevalence",
      "effect",
      "selection",
      "cond_prevalence",
      "cond_effect",
      "cond_selection",
      "Q_dist")
  }


  output <- as.data.frame(cbind(item, point))
  output$point <- as.numeric(output$point)

  return(output)
}
