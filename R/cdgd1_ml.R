


#' Perform conditional decomposition via machine learning
#'
#' @param Y Outcome. The name of a continuous variable.
#' @param D Treatment status. The name of a binary numeric variable taking values of 0 and 1.
#' @param G Advantaged group membership. The name of a binary numeric variable taking values of 0 and 1.
#' @param Q Conditional set. The vector of the names of numeric variables.
#' @param X Confounders. The vector of the names of numeric variables.
#' @param data A data frame.
#' @param algorithm The ML algorithm for modelling. "nnet" for neural network, "ranger" for random forests, "gbm" for generalized boosted models.
#' @param alpha 1-alpha confidence interval.
#' @param trim1 Threshold for trimming the propensity score. When trim1=a, individuals with propensity scores lower than a or higher than 1-a will be dropped.
#' @param trim2 Threshold for trimming the G given Q predictions. When trim2=a, individuals with G given Q predictions lower than a or higher than 1-a will be dropped.
#'
#' @return A dataframe of estimates.
#'
#' @export
#'
#' @examples
#' # This example will take a minute to run.
#' \donttest{
#' data(exp_data)
#'
#' set.seed(1)
#'
#' results <- cdgd1_ml(
#' Y="outcome",
#' D="treatment",
#' G="group_a",
#' X="confounder",
#' Q="Q",
#' data=exp_data,
#' algorithm="gbm")
#'
#' results}



cdgd1_ml <- function(Y,D,G,X,Q,data,algorithm,alpha=0.05,trim1=0,trim2=0) {

  if (!requireNamespace("caret", quietly=TRUE)) {
    stop(
      "Package \"caret\" must be installed to use this function.",
      call. = FALSE
    )
  }

  data <- as.data.frame(data)

### estimate the nuisance functions with cross-fitting
  sample1 <- sample(nrow(data), floor(nrow(data)/2), replace=FALSE)
  sample2 <- setdiff(1:nrow(data), sample1)

### propensity score model
  data[,D] <- as.factor(data[,D])
  levels(data[,D]) <- c("D0","D1")  # necessary for caret implementation of ranger

  if (algorithm=="nnet") {
    message <- utils::capture.output( DgivenGXQ.Model.sample1 <- caret::train(stats::as.formula(paste(D, paste(G,paste(Q,collapse="+"),paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample1,], method="nnet",
                                                                              preProc=c("center","scale"), trControl=caret::trainControl(method="cv"), linout=FALSE ))
    message <- utils::capture.output( DgivenGXQ.Model.sample2 <- caret::train(stats::as.formula(paste(D, paste(G,paste(Q,collapse="+"),paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample2,], method="nnet",
                                                                              preProc=c("center","scale"), trControl=caret::trainControl(method="cv"), linout=FALSE ))
  }
  if (algorithm=="ranger") {
    message <- utils::capture.output( DgivenGXQ.Model.sample1 <- caret::train(stats::as.formula(paste(D, paste(G,paste(Q,collapse="+"),paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample1,], method="ranger",
                                                                              trControl=caret::trainControl(method="cv", classProbs=TRUE)) )
    message <- utils::capture.output( DgivenGXQ.Model.sample2 <- caret::train(stats::as.formula(paste(D, paste(G,paste(Q,collapse="+"),paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample2,], method="ranger",
                                                                              trControl=caret::trainControl(method="cv", classProbs=TRUE)) )
  }
  if (algorithm=="gbm") {
    message <- utils::capture.output( DgivenGXQ.Model.sample1 <- caret::train(stats::as.formula(paste(D, paste(G,paste(Q,collapse="+"),paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample1,], method="gbm",
                                                                              trControl=caret::trainControl(method="cv")) )
    message <- utils::capture.output( DgivenGXQ.Model.sample2 <- caret::train(stats::as.formula(paste(D, paste(G,paste(Q,collapse="+"),paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample2,], method="gbm",
                                                                              trControl=caret::trainControl(method="cv")) )
  }

  data[,D] <- as.numeric(data[,D])-1

  # propensity score predictions
  DgivenGXQ.Pred <- rep(NA, nrow(data))

  DgivenGXQ.Pred[sample2] <- stats::predict(DgivenGXQ.Model.sample1, newdata = data[sample2,], type="prob")[,2]
  DgivenGXQ.Pred[sample1] <- stats::predict(DgivenGXQ.Model.sample2, newdata = data[sample1,], type="prob")[,2]

  DgivenGXQ.Pred_ncf <- rep(NA, nrow(data)) # ncf stands for non-cross-fitted

  DgivenGXQ.Pred_ncf[sample1] <- stats::predict(DgivenGXQ.Model.sample1, newdata = data[sample1,], type="prob")[,2]
  DgivenGXQ.Pred_ncf[sample2] <- stats::predict(DgivenGXQ.Model.sample2, newdata = data[sample2,], type="prob")[,2]

  # trim the sample based on the propensity score
  dropped <- sum(DgivenGXQ.Pred<trim1 | DgivenGXQ.Pred>1-trim1 | DgivenGXQ.Pred_ncf<trim1 | DgivenGXQ.Pred_ncf>1-trim1)  # the number of dropped obs

  data$sample1 <- 1:nrow(data) %in% sample1    # sample 1 indicator for the new data
  keeping_ind <- DgivenGXQ.Pred>=trim1 & DgivenGXQ.Pred<=1-trim1 & DgivenGXQ.Pred_ncf>=trim1 & DgivenGXQ.Pred_ncf<=1-trim1 # must have this for cdgd1_ml
  data <- data[keeping_ind, ]
  sample1 <- which(data$sample1)
  sample2 <- setdiff(1:nrow(data), sample1)

  DgivenGXQ.Pred <- DgivenGXQ.Pred[keeping_ind]
  DgivenGXQ.Pred_ncf <- DgivenGXQ.Pred_ncf[keeping_ind]

  zero_one <- sum(DgivenGXQ.Pred==0 | DgivenGXQ.Pred==1 | DgivenGXQ.Pred_ncf==0 | DgivenGXQ.Pred_ncf==1)
  if ( zero_one>0 ) {
    stop(
      paste("D given X, Q, and G are exact 0 or 1 in", zero_one, "cases.", sep=" "),
      call. = FALSE
    )
  }

  ### outcome regression model
  if (algorithm=="nnet") {
    if (!requireNamespace("nnet", quietly=TRUE)) {
      stop(
        "Package \"nnet\" must be installed to use this function.",
        call. = FALSE
      )
    }
    message <- utils::capture.output( YgivenDGXQ.Model.sample1 <- caret::train(stats::as.formula(paste(Y, paste(D,G,paste(Q,collapse="+"),paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample1,], method="nnet",
                                                                               preProc=c("center","scale"), trControl=caret::trainControl(method="cv"), linout=TRUE ))
    message <- utils::capture.output( YgivenDGXQ.Model.sample2 <- caret::train(stats::as.formula(paste(Y, paste(D,G,paste(Q,collapse="+"),paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample2,], method="nnet",
                                                                               preProc=c("center","scale"), trControl=caret::trainControl(method="cv"), linout=TRUE ))
  }
  if (algorithm=="ranger") {
    if (!requireNamespace("ranger", quietly=TRUE)) {
      stop(
        "Package \"ranger\" must be installed to use this function.",
        call. = FALSE
      )
    }
    message <- utils::capture.output( YgivenDGXQ.Model.sample1 <- caret::train(stats::as.formula(paste(Y, paste(D,G,paste(Q,collapse="+"),paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample1,], method="ranger",
                                                                               trControl=caret::trainControl(method="cv")) )
    message <- utils::capture.output( YgivenDGXQ.Model.sample2 <- caret::train(stats::as.formula(paste(Y, paste(D,G,paste(Q,collapse="+"),paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample2,], method="ranger",
                                                                               trControl=caret::trainControl(method="cv")) )
  }
  if (algorithm=="gbm") {
    if (!requireNamespace("gbm", quietly=TRUE)) {
      stop(
        "Package \"gbm\" must be installed to use this function.",
        call. = FALSE
      )
    }
    message <- utils::capture.output( YgivenDGXQ.Model.sample1 <- caret::train(stats::as.formula(paste(Y, paste(D,G,paste(Q,collapse="+"),paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample1,], method="gbm",
                                                                               trControl=caret::trainControl(method="cv")) )
    message <- utils::capture.output( YgivenDGXQ.Model.sample2 <- caret::train(stats::as.formula(paste(Y, paste(D,G,paste(Q,collapse="+"),paste(X,collapse="+"),sep="+"), sep="~")), data=data[sample2,], method="gbm",
                                                                               trControl=caret::trainControl(method="cv")) )
  }

### outcome predictions
  YgivenGXQ.Pred_D0 <- YgivenGXQ.Pred_D1 <- rep(NA, nrow(data))

  pred_data <- data
  pred_data[,D] <- 0
  YgivenGXQ.Pred_D0[sample2] <- stats::predict(YgivenDGXQ.Model.sample1, newdata = pred_data[sample2,])
  YgivenGXQ.Pred_D0[sample1] <- stats::predict(YgivenDGXQ.Model.sample2, newdata = pred_data[sample1,])

  pred_data <- data
  pred_data[,D] <- 1
  YgivenGXQ.Pred_D1[sample2] <- stats::predict(YgivenDGXQ.Model.sample1, newdata = pred_data[sample2,])
  YgivenGXQ.Pred_D1[sample1] <- stats::predict(YgivenDGXQ.Model.sample2, newdata = pred_data[sample1,])

### Estimate p_g(Q)=Pr(G=g | Q)
  data[,G] <- as.factor(data[,G])
  levels(data[,G]) <- c("G0","G1")  # necessary for caret implementation of ranger

  if (algorithm=="nnet") {
    message <- utils::capture.output( GgivenQ.Model.sample1 <- caret::train(stats::as.formula(paste(G, paste(Q,collapse="+"), sep="~")), data=data[sample1,], method="nnet",
                                                                            preProc=c("center","scale"), trControl=caret::trainControl(method="cv"), linout=FALSE ))
    message <- utils::capture.output( GgivenQ.Model.sample2 <- caret::train(stats::as.formula(paste(G, paste(Q,collapse="+"), sep="~")), data=data[sample2,], method="nnet",
                                                                            preProc=c("center","scale"), trControl=caret::trainControl(method="cv"), linout=FALSE ))
  }
  if (algorithm=="ranger") {
    message <- utils::capture.output( GgivenQ.Model.sample1 <- caret::train(stats::as.formula(paste(G, paste(Q,collapse="+"), sep="~")), data=data[sample1,], method="ranger",
                                                                            trControl=caret::trainControl(method="cv", classProbs=TRUE)) )
    message <- utils::capture.output( GgivenQ.Model.sample2 <- caret::train(stats::as.formula(paste(G, paste(Q,collapse="+"), sep="~")), data=data[sample2,], method="ranger",
                                                                            trControl=caret::trainControl(method="cv", classProbs=TRUE)) )
  }
  if (algorithm=="gbm") {
    message <- utils::capture.output( GgivenQ.Model.sample1 <- caret::train(stats::as.formula(paste(G, paste(Q,collapse="+"), sep="~")), data=data[sample1,], method="gbm",
                                                                            trControl=caret::trainControl(method="cv")) )
    message <- utils::capture.output( GgivenQ.Model.sample2 <- caret::train(stats::as.formula(paste(G, paste(Q,collapse="+"), sep="~")), data=data[sample2,], method="gbm",
                                                                            trControl=caret::trainControl(method="cv")) )
  }

  data[,G] <- as.numeric(data[,G])-1

  GgivenQ.Pred <- rep(NA, nrow(data))
  GgivenQ.Pred[sample2] <- stats::predict(GgivenQ.Model.sample1, newdata = data[sample2,], type="prob")[,2]
  GgivenQ.Pred[sample1] <- stats::predict(GgivenQ.Model.sample2, newdata = data[sample1,], type="prob")[,2]

  # trim the sample based on the G given Q predictions
  dropped <- dropped + sum(GgivenQ.Pred<trim2 | GgivenQ.Pred>1-trim2)  # update the number of dropped obs

  data$sample1 <- 1:nrow(data) %in% sample1    # sample 1 indicator for the new data
  data <- data[GgivenQ.Pred>=trim2 & GgivenQ.Pred<=1-trim2, ]
  sample1 <- which(data$sample1)
  sample2 <- setdiff(1:nrow(data), sample1)

  YgivenGXQ.Pred_D0 <- YgivenGXQ.Pred_D0[GgivenQ.Pred>=trim2 & GgivenQ.Pred<=1-trim2]
  YgivenGXQ.Pred_D1 <- YgivenGXQ.Pred_D1[GgivenQ.Pred>=trim2 & GgivenQ.Pred<=1-trim2]
  DgivenGXQ.Pred <- DgivenGXQ.Pred[GgivenQ.Pred>=trim2 & GgivenQ.Pred<=1-trim2]
  DgivenGXQ.Pred_ncf <- DgivenGXQ.Pred_ncf[GgivenQ.Pred>=trim2 & GgivenQ.Pred<=1-trim2]
  GgivenQ.Pred <- GgivenQ.Pred[GgivenQ.Pred>=trim2 & GgivenQ.Pred<=1-trim2]

  zero_one <- sum(GgivenQ.Pred==0)+sum(GgivenQ.Pred==1)
  if ( zero_one>0 ) {
    stop(
      paste("G given Q are exact 0 or 1 in", zero_one, "cases.", sep=" "),
      call. = FALSE
    )
  }

### Estimate E(Y_d | Q,g)
  YgivenGXQ.Pred_D1_ncf <- YgivenGXQ.Pred_D0_ncf <- rep(NA, nrow(data)) # ncf stands for non-cross-fitted

  pred_data <- data
  pred_data[,D] <- 1
  YgivenGXQ.Pred_D1_ncf[sample1] <- stats::predict(YgivenDGXQ.Model.sample1, newdata = pred_data[sample1,])
  YgivenGXQ.Pred_D1_ncf[sample2] <- stats::predict(YgivenDGXQ.Model.sample2, newdata = pred_data[sample2,])

  pred_data <- data
  pred_data[,D] <- 0
  YgivenGXQ.Pred_D0_ncf[sample1] <- stats::predict(YgivenDGXQ.Model.sample1, newdata = pred_data[sample1,])
  YgivenGXQ.Pred_D0_ncf[sample2] <- stats::predict(YgivenDGXQ.Model.sample2, newdata = pred_data[sample2,])

  # IPOs for modelling E(Y_d | Q,g)
  IPO_D0_ncf <- (1-data[,D])/(1-DgivenGXQ.Pred_ncf)/mean((1-data[,D])/(1-DgivenGXQ.Pred_ncf))*(data[,Y]-YgivenGXQ.Pred_D0_ncf) + YgivenGXQ.Pred_D0_ncf
  IPO_D1_ncf <- data[,D]/DgivenGXQ.Pred_ncf/mean(data[,D]/DgivenGXQ.Pred_ncf)*(data[,Y]-YgivenGXQ.Pred_D1_ncf) + YgivenGXQ.Pred_D1_ncf

  data_temp <- data[,c(G,Q)]
  data_temp$IPO_D0_ncf <- IPO_D0_ncf
  data_temp$IPO_D1_ncf <- IPO_D1_ncf

  if (algorithm=="nnet") {
    message <- utils::capture.output( Y0givenGQ.Model.sample1 <- caret::train(stats::as.formula(paste("IPO_D0_ncf", paste(G,paste(Q,collapse="+"),sep="+"), sep="~")), data=data_temp[sample1,], method="nnet",
                                                                              preProc=c("center","scale"), trControl=caret::trainControl(method="cv"), linout=TRUE ))
    message <- utils::capture.output( Y0givenGQ.Model.sample2 <- caret::train(stats::as.formula(paste("IPO_D0_ncf", paste(G,paste(Q,collapse="+"),sep="+"), sep="~")), data=data_temp[sample2,], method="nnet",
                                                                              preProc=c("center","scale"), trControl=caret::trainControl(method="cv"), linout=TRUE ))
    message <- utils::capture.output( Y1givenGQ.Model.sample1 <- caret::train(stats::as.formula(paste("IPO_D1_ncf", paste(G,paste(Q,collapse="+"),sep="+"), sep="~")), data=data_temp[sample1,], method="nnet",
                                                                              preProc=c("center","scale"), trControl=caret::trainControl(method="cv"), linout=TRUE ))
    message <- utils::capture.output( Y1givenGQ.Model.sample2 <- caret::train(stats::as.formula(paste("IPO_D1_ncf", paste(G,paste(Q,collapse="+"),sep="+"), sep="~")), data=data_temp[sample2,], method="nnet",
                                                                              preProc=c("center","scale"), trControl=caret::trainControl(method="cv"), linout=TRUE ))
  }
  if (algorithm=="ranger") {
    message <- utils::capture.output( Y0givenGQ.Model.sample1 <- caret::train(stats::as.formula(paste("IPO_D0_ncf", paste(G,paste(Q,collapse="+"),sep="+"), sep="~")), data=data_temp[sample1,], method="ranger",
                                                                              trControl=caret::trainControl(method="cv")) )
    message <- utils::capture.output( Y0givenGQ.Model.sample2 <- caret::train(stats::as.formula(paste("IPO_D0_ncf", paste(G,paste(Q,collapse="+"),sep="+"), sep="~")), data=data_temp[sample2,], method="ranger",
                                                                              trControl=caret::trainControl(method="cv")) )
    message <- utils::capture.output( Y1givenGQ.Model.sample1 <- caret::train(stats::as.formula(paste("IPO_D1_ncf", paste(G,paste(Q,collapse="+"),sep="+"), sep="~")), data=data_temp[sample1,], method="ranger",
                                                                              trControl=caret::trainControl(method="cv")) )
    message <- utils::capture.output( Y1givenGQ.Model.sample2 <- caret::train(stats::as.formula(paste("IPO_D1_ncf", paste(G,paste(Q,collapse="+"),sep="+"), sep="~")), data=data_temp[sample2,], method="ranger",
                                                                              trControl=caret::trainControl(method="cv")) )
  }
  if (algorithm=="gbm") {
    message <- utils::capture.output( Y0givenGQ.Model.sample1 <- caret::train(stats::as.formula(paste("IPO_D0_ncf", paste(G,paste(Q,collapse="+"),sep="+"), sep="~")), data=data_temp[sample1,], method="gbm",
                                                                              trControl=caret::trainControl(method="cv")) )
    message <- utils::capture.output( Y0givenGQ.Model.sample2 <- caret::train(stats::as.formula(paste("IPO_D0_ncf", paste(G,paste(Q,collapse="+"),sep="+"), sep="~")), data=data_temp[sample2,], method="gbm",
                                                                              trControl=caret::trainControl(method="cv")) )
    message <- utils::capture.output( Y1givenGQ.Model.sample1 <- caret::train(stats::as.formula(paste("IPO_D1_ncf", paste(G,paste(Q,collapse="+"),sep="+"), sep="~")), data=data_temp[sample1,], method="gbm",
                                                                              trControl=caret::trainControl(method="cv")) )
    message <- utils::capture.output( Y1givenGQ.Model.sample2 <- caret::train(stats::as.formula(paste("IPO_D1_ncf", paste(G,paste(Q,collapse="+"),sep="+"), sep="~")), data=data_temp[sample2,], method="gbm",
                                                                              trControl=caret::trainControl(method="cv")) )
  }

  Y0givenQ.Pred_G0 <- Y0givenQ.Pred_G1 <- Y1givenQ.Pred_G0 <- Y1givenQ.Pred_G1 <- rep(NA, nrow(data))

  pred_data <- data
  pred_data[,G] <- 1
  Y0givenQ.Pred_G1[sample2] <- stats::predict(Y0givenGQ.Model.sample1, newdata = pred_data[sample2,])  # cross-fitting is used
  Y0givenQ.Pred_G1[sample1] <- stats::predict(Y0givenGQ.Model.sample2, newdata = pred_data[sample1,])
  Y1givenQ.Pred_G1[sample2] <- stats::predict(Y1givenGQ.Model.sample1, newdata = pred_data[sample2,])
  Y1givenQ.Pred_G1[sample1] <- stats::predict(Y1givenGQ.Model.sample2, newdata = pred_data[sample1,])

  pred_data <- data
  pred_data[,G] <- 0
  Y0givenQ.Pred_G0[sample2] <- stats::predict(Y0givenGQ.Model.sample1, newdata = pred_data[sample2,])  # cross-fitting is used
  Y0givenQ.Pred_G0[sample1] <- stats::predict(Y0givenGQ.Model.sample2, newdata = pred_data[sample1,])
  Y1givenQ.Pred_G0[sample2] <- stats::predict(Y1givenGQ.Model.sample1, newdata = pred_data[sample2,])
  Y1givenQ.Pred_G0[sample1] <- stats::predict(Y1givenGQ.Model.sample2, newdata = pred_data[sample1,])

### The "IPO" (individual potential outcome) function
  # For each d and g value, we have IE(d,g)=\frac{\one(D=d)}{\pi(d,X,g)}[Y-\mu(d,X,g)]+\mu(d,X,g)
  # We stabilize the weight by dividing the sample average of estimated weights

  IPO_D0 <- (1-data[,D])/(1-DgivenGXQ.Pred)/mean((1-data[,D])/(1-DgivenGXQ.Pred))*(data[,Y]-YgivenGXQ.Pred_D0) + YgivenGXQ.Pred_D0
  IPO_D1 <- data[,D]/DgivenGXQ.Pred/mean(data[,D]/DgivenGXQ.Pred)*(data[,Y]-YgivenGXQ.Pred_D1) + YgivenGXQ.Pred_D1

### Estimate E(D | Q,g')
  data[,D] <- as.factor(data[,D])
  levels(data[,D]) <- c("D0","D1")  # necessary for caret implementation of ranger

  if (algorithm=="nnet") {
    message <- utils::capture.output( DgivenGQ.Model.sample1 <- caret::train(stats::as.formula(paste(D, paste(G,paste(Q,collapse="+"),sep="+"), sep="~")), data=data[sample1,], method="nnet",
                                                                             preProc=c("center","scale"), trControl=caret::trainControl(method="cv"), linout=FALSE ))
    message <- utils::capture.output( DgivenGQ.Model.sample2 <- caret::train(stats::as.formula(paste(D, paste(G,paste(Q,collapse="+"),sep="+"), sep="~")), data=data[sample2,], method="nnet",
                                                                             preProc=c("center","scale"), trControl=caret::trainControl(method="cv"), linout=FALSE ))
  }
  if (algorithm=="ranger") {
    message <- utils::capture.output( DgivenGQ.Model.sample1 <- caret::train(stats::as.formula(paste(D, paste(G,paste(Q,collapse="+"),sep="+"), sep="~")), data=data[sample1,], method="ranger",
                                                                              trControl=caret::trainControl(method="cv", classProbs=TRUE)) )
    message <- utils::capture.output( DgivenGQ.Model.sample2 <- caret::train(stats::as.formula(paste(D, paste(G,paste(Q,collapse="+"),sep="+"), sep="~")), data=data[sample2,], method="ranger",
                                                                              trControl=caret::trainControl(method="cv", classProbs=TRUE)) )
  }
  if (algorithm=="gbm") {
    message <- utils::capture.output( DgivenGQ.Model.sample1 <- caret::train(stats::as.formula(paste(D, paste(G,paste(Q,collapse="+"),sep="+"), sep="~")), data=data[sample1,], method="gbm",
                                                                              trControl=caret::trainControl(method="cv")) )
    message <- utils::capture.output( DgivenGQ.Model.sample2 <- caret::train(stats::as.formula(paste(D, paste(G,paste(Q,collapse="+"),sep="+"), sep="~")), data=data[sample2,], method="gbm",
                                                                              trControl=caret::trainControl(method="cv")) )
  }

  data[,D] <- as.numeric(data[,D])-1

  DgivenQ.Pred_G0 <- DgivenQ.Pred_G1 <- rep(NA, nrow(data))

  pred_data <- data
  pred_data[,G] <- 0
  DgivenQ.Pred_G0[sample2] <- stats::predict(DgivenGQ.Model.sample1, newdata = pred_data[sample2,], type="prob")[,2]
  DgivenQ.Pred_G0[sample1] <- stats::predict(DgivenGQ.Model.sample2, newdata = pred_data[sample1,], type="prob")[,2]

  pred_data <- data
  pred_data[,G] <- 1
  DgivenQ.Pred_G1[sample2] <- stats::predict(DgivenGQ.Model.sample1, newdata = pred_data[sample2,], type="prob")[,2]
  DgivenQ.Pred_G1[sample1] <- stats::predict(DgivenGQ.Model.sample2, newdata = pred_data[sample1,], type="prob")[,2]

### The one-step estimate of \xi_{dg}
  psi_00 <- mean( (1-data[,G])/(1-mean(data[,G]))*IPO_D0 )
  psi_01 <- mean( data[,G]/mean(data[,G])*IPO_D0 )
  # Note that this is basically DML2. We could also use DML1:
  #psi_00_S1 <- mean( (1-data[sample1,G])/(1-mean(data[sample1,G]))*IPO_D0[sample1] )     # sample 1 estimate
  #psi_00_S2 <- mean( (1-data[sample2,G])/(1-mean(data[sample2,G]))*IPO_D0[sample2] )     # sample 2 estimate
  #psi_00 <- (1/2)*(psi_00_S1+psi_00_S2)
  #psi_01_S1 <- mean( data[sample1,G]/mean(data[sample1,G])*IPO_D0[sample1] )     # sample 1 estimate
  #psi_01_S2 <- mean( data[sample1,G]/mean(data[sample1,G])*IPO_D0[sample2] )     # sample 2 estimate
  #psi_01 <- (1/2)*(psi_01_S1+psi_01_S2)

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
  baseline_se <- se( data[,G]/mean(data[,G])*(IPO_D0-psi_01) - (1-data[,G])/(1-mean(data[,G]))*(IPO_D0-psi_00) )
  # Alternatively, we could use
  # se( c( data[sample1,G]/mean(data[sample1,G])*(IPO_D0[sample1]-psi_01) - (1-data[sample1,G])/(1-mean(data[sample1,G]))*(IPO_D0[sample1]-psi_00),
  #         data[sample2,G]/mean(data[sample2,G])*(IPO_D0[sample2]-psi_01) - (1-data[sample2,G])/(1-mean(data[sample2,G]))*(IPO_D0[sample2]-psi_00) ) )
  # But there isn't a theoretically strong reason to prefer one over the other.

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
                         ( data[,G]/mean(data[,G])*(IPO_D0-psi_01) - (1-data[,G])/(1-mean(data[,G]))*(IPO_D0-psi_00) ) -
                         ( EIF_dggg(1,0,1,0)-EIF_dggg(0,0,1,0)-EIF_dggg(1,0,0,0)+EIF_dggg(0,0,0,0) ) -
                         ( EIF_dggg(1,1,1,1)-EIF_dggg(0,1,1,1)-EIF_dggg(1,0,1,1)+EIF_dggg(0,0,1,1) ) -
                         ( EIF_dggg(1,0,1,1)-EIF_dggg(0,0,1,1)-EIF_dggg(1,0,1,0)+EIF_dggg(0,0,1,0) ))

  cond_Jackson_reduction_se <- se( (1-data[,G])/(1-mean(data[,G]))*(IPO_D0-psi_00)+EIF_dggg(1,0,1,0)-EIF_dggg(0,0,1,0)-(1-data[,G])/(1-mean(data[,G]))*(data[,Y]-Y_G0) )

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
