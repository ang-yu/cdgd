


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
#' @param t Threshold of propensity score censoring. Propensity scores larger than (1-t)th quantile or smaller than tth quantile are censored.
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



cdgd0 <- function(Y,D,G,Q=NULL,X,data,t=0.025,algorithm) {

  if (!requireNamespace("caret", quietly=TRUE)) {
    stop(
      "Package \"caret\" must be installed to use this function.",
      call. = FALSE
    )
  }

  data <- as.data.frame(data)

### estimate the nuisance functions within each group, so that the final estimates are independent across groups
  sample1 <- sample(nrow(data), floor(nrow(data)/2), replace=F)
  sample2 <- setdiff(1:nrow(data), sample1)

### outcome regression model
  if (algorithm=="nnet") {
    if (!requireNamespace("nnet", quietly=TRUE)) {
      stop(
        "Package \"nnet\" must be installed to use this function.",
        call. = FALSE
      )
    }
    message <- utils::capture.output( YgivenDGX.Model.sample1 <- caret::train(stats::as.formula(paste(Y, paste(D,G,paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample1,], method="nnet",
                                                                             preProc=c("center","scale"), trControl=caret::trainControl(method="cv"), linout=TRUE,
                                                                             tuneGrid=expand.grid(size=c(1,2),decay=c(0,0.1,0.2,0.3)) ))
    message <- utils::capture.output( YgivenDGX.Model.sample2 <- caret::train(stats::as.formula(paste(Y, paste(D,G,paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample2,], method="nnet",
                                                                             preProc=c("center","scale"), trControl=caret::trainControl(method="cv"), linout=TRUE,
                                                                             tuneGrid=expand.grid(size=c(1,2),decay=c(0,0.1,0.2,0.3)) ))
  }
  if (algorithm=="ranger") {
    if (!requireNamespace("ranger", quietly=TRUE)) {
      stop(
        "Package \"ranger\" must be installed to use this function.",
        call. = FALSE
      )
    }
    message <- utils::capture.output( YgivenDGX.Model.sample1 <- caret::train(stats::as.formula(paste(Y, paste(D,G,paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample1,], method="ranger",
                                                                             trControl=caret::trainControl(method="cv"),
                                                                             tuneGrid=expand.grid(mtry=c(5,10,15,20),splitrule=c("variance"),min.node.size=c(50,100,200))) )
    message <- utils::capture.output( YgivenDGX.Model.sample2 <- caret::train(stats::as.formula(paste(Y, paste(D,G,paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample2,], method="ranger",
                                                                             trControl=caret::trainControl(method="cv"),
                                                                             tuneGrid=expand.grid(mtry=c(5,10,15,20),splitrule=c("variance"),min.node.size=c(50,100,200))) )
  }
  if (algorithm=="gbm") {
    if (!requireNamespace("gbm", quietly=TRUE)) {
      stop(
        "Package \"gbm\" must be installed to use this function.",
        call. = FALSE
      )
    }
    message <- utils::capture.output( YgivenDGX.Model.sample1 <- caret::train(stats::as.formula(paste(Y, paste(D,G,paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample1,], method="gbm",
                                                                             trControl=caret::trainControl(method="cv"),
                                                                             tuneGrid=expand.grid(n.trees=c(25,75,125), interaction.depth=c(1,2,4), shrinkage=0.1, n.minobsinnode=c(5,10,15))) )
    message <- utils::capture.output( YgivenDGX.Model.sample2 <- caret::train(stats::as.formula(paste(Y, paste(D,G,paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample2,], method="gbm",
                                                                             trControl=caret::trainControl(method="cv"),
                                                                             tuneGrid=expand.grid(n.trees=c(25,75,125), interaction.depth=c(1,2,4), shrinkage=0.1, n.minobsinnode=c(5,10,15))) )
  }

### propensity score model
  data[,D] <- as.factor(data[,D])
  levels(data[,D]) <- c("D0","D1")  # necessary for caret implementation of ranger

  if (algorithm=="nnet") {
    message <- utils::capture.output( DgivenGX.Model.sample1 <- caret::train(stats::as.formula(paste(D, paste(G,paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample1,], method="nnet",
                                                                            preProc=c("center","scale"), trControl=caret::trainControl(method="cv"), linout=FALSE,
                                                                            tuneGrid=expand.grid(size=c(1,2),decay=c(0,0.1,0.2,0.3)) ))
    message <- utils::capture.output( DgivenGX.Model.sample2 <- caret::train(stats::as.formula(paste(D, paste(G,paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample2,], method="nnet",
                                                                            preProc=c("center","scale"), trControl=caret::trainControl(method="cv"), linout=FALSE,
                                                                            tuneGrid=expand.grid(size=c(1,2),decay=c(0,0.1,0.2,0.3)) ))
  }
  if (algorithm=="ranger") {
    message <- utils::capture.output( DgivenGX.Model.sample1 <- caret::train(stats::as.formula(paste(D, paste(G,paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample1,], method="ranger",
                                                                            trControl=caret::trainControl(method="cv", classProbs=TRUE),
                                                                            tuneGrid=expand.grid(mtry=c(5,10,15,20),splitrule=c("gini"),min.node.size=c(50,100,200))) )
    message <- utils::capture.output( DgivenGX.Model.sample2 <- caret::train(stats::as.formula(paste(D, paste(G,paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample2,], method="ranger",
                                                                            trControl=caret::trainControl(method="cv", classProbs=TRUE),
                                                                            tuneGrid=expand.grid(mtry=c(5,10,15,20),splitrule=c("gini"),min.node.size=c(50,100,200))) )
  }
  if (algorithm=="gbm") {
    message <- utils::capture.output( DgivenGX.Model.sample1 <- caret::train(stats::as.formula(paste(D, paste(G,paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample1,], method="gbm",
                                                                            trControl=caret::trainControl(method="cv"),
                                                                            tuneGrid=expand.grid(n.trees=c(25,75,125), interaction.depth=c(1,2,4), shrinkage=0.1, n.minobsinnode=c(5,10,15))) )
    message <- utils::capture.output( DgivenGX.Model.sample2 <- caret::train(stats::as.formula(paste(D, paste(G,paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample2,], method="gbm",
                                                                            trControl=caret::trainControl(method="cv"),
                                                                            tuneGrid=expand.grid(n.trees=c(25,75,125), interaction.depth=c(1,2,4), shrinkage=0.1, n.minobsinnode=c(5,10,15))) )
  }

  data[,D] <- as.numeric(data[,D])-1

### cross-fitted predictions
  YgivenX.Pred_D0G0 <- YgivenX.Pred_D1G0 <- YgivenX.Pred_D0G1 <- YgivenX.Pred_D1G1 <- DgivenX.Pred_G0 <- DgivenX.Pred_G1 <- rep(NA, nrow(data))

  pred_data <- data
  pred_data[,D] <- 0
  pred_data[,G] <- 0
  YgivenX.Pred_D0G0[sample1] <- stats::predict(YgivenDGX.Model.sample1, newdata = pred_data[sample2,])
  YgivenX.Pred_D0G0[sample2] <- stats::predict(YgivenDGX.Model.sample2, newdata = pred_data[sample1,])

  pred_data <- data
  pred_data[,D] <- 1
  pred_data[,G] <- 0
  YgivenX.Pred_D1G0[sample1] <- stats::predict(YgivenDGX.Model.sample1, newdata = pred_data[sample2,])
  YgivenX.Pred_D1G0[sample2] <- stats::predict(YgivenDGX.Model.sample2, newdata = pred_data[sample1,])

  pred_data <- data
  pred_data[,D] <- 0
  pred_data[,G] <- 1
  YgivenX.Pred_D0G1[sample1] <- stats::predict(YgivenDGX.Model.sample1, newdata = pred_data[sample2,])
  YgivenX.Pred_D0G1[sample2] <- stats::predict(YgivenDGX.Model.sample2, newdata = pred_data[sample1,])

  pred_data <- data
  pred_data[,D] <- 1
  pred_data[,G] <- 1
  YgivenX.Pred_D1G1[sample1] <- stats::predict(YgivenDGX.Model.sample1, newdata = pred_data[sample2,])
  YgivenX.Pred_D1G1[sample2] <- stats::predict(YgivenDGX.Model.sample2, newdata = pred_data[sample1,])

  pred_data <- data
  pred_data[,G] <- 0
  DgivenX.Pred_G0[sample1] <- stats::predict(DgivenGX.Model.sample1, newdata = pred_data[sample2,], type="prob")[,2]
  DgivenX.Pred_G0[sample2] <- stats::predict(DgivenGX.Model.sample2, newdata = pred_data[sample1,], type="prob")[,2]

  pred_data <- data
  pred_data[,G] <- 1
  DgivenX.Pred_G1[sample1] <- stats::predict(DgivenGX.Model.sample1, newdata = pred_data[sample2,], type="prob")[,2]
  DgivenX.Pred_G1[sample2] <- stats::predict(DgivenGX.Model.sample2, newdata = pred_data[sample1,], type="prob")[,2]

### The "IPO" (individual potential outcome) function
  # For each d and g value, we have IE(d,g)=\frac{\one(D=d)}{\pi(d,X,g)}[Y-\mu(d,X,g)]+\mu(d,X,g)
  # We stablize the weight by dividing the sample average of estimated weights

  IPO_D0G0 <- (1-data[,D])/(1-DgivenX.Pred_G0)/mean((1-data[,D])/(1-DgivenX.Pred_G0))*(data[,Y]-YgivenX.Pred_D0G0) + YgivenX.Pred_D0G0
  IPO_D1G0 <- data[,D]/DgivenX.Pred_G0/(mean(data[,D]/DgivenX.Pred_G0))*(data[,Y]-YgivenX.Pred_D1G0) + YgivenX.Pred_D1G0
  IPO_D0G1 <- (1-data[,D])/(1-DgivenX.Pred_G0)/(mean((1-data[,D])/(1-DgivenX.Pred_G0)))*(data[,Y]-YgivenX.Pred_D0G1) + YgivenX.Pred_D0G1
  IPO_D1G1 <- data[,D]/DgivenX.Pred_G0/mean(data[,D]/DgivenX.Pred_G0)*(data[,Y]-YgivenX.Pred_D1G1) + YgivenX.Pred_D1G1

### The cross-fitted substitution estimates of \xi_{dg} and \xi_{dgg'}
  xi_D0G0 <- xi_D0G1 <- rep(NA, nrow(data))
  xi_D0G0[sample1] <- mean(YgivenX.Pred_D0G0[intersect(which(data[,G]==0),sample2)])
  xi_D0G0[sample2] <- mean(YgivenX.Pred_D0G0[intersect(which(data[,G]==0),sample1)])
  xi_D0G1[sample1] <- mean(YgivenX.Pred_D0G0[intersect(which(data[,G]==1),sample2)])
  xi_D0G1[sample2] <- mean(YgivenX.Pred_D0G0[intersect(which(data[,G]==1),sample1)])

### The one-step estimate of \xi_{dg} and \xi_{dgg'}


  Y0_i <- ATT_i <- ATE_i <- rep(NA, nrow(data))


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
