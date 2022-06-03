
#' Decomposition based on conditional ignorability
#'
#' @param Y Outcome. The name of a continuous variable.
#' @param D Treatment status. The name of a binary numeric variable taking values of 0 and 1.
#' @param G1 Group 1 membership. The name of a binary factor variable taking values of 0 and 1.
#' @param G2 Group 2 membership. The name of a binary factor variable taking values of 0 and 1.
#' @param Q Conditional set. The vector of the names of numeric variables.
#' @param X Confounders. The vector of the names of numeric variables.
#' @param data A data frame.
#' @param weight Survey weights. The name of a numeric variable.
#' @param alpha (1-alpha) confidence intervals.
#' @param k Number of monte carlo simulation.
#' @param t Threshold of propensity score censoring. Propensity scores larger than 1-t or smaller than t will be censored.
#' @param algorithm The ML algorithm for modelling. "nnet" for neural network and "ranger" for random forests.
#'
#' @return A data frame of point estimates and confidence intervals
#' @export
#'
#' @examples
#' data(exp_data)
#'
#' set.seed(1)
#'
#' results <- cdgd(
#' Y="outcome",
#' D="treatment",
#' G1="group_a",
#' G2="group_b",
#' X=c("confounder","Q"),
#' Q="Q",
#' data=exp_data,
#' t=0.05,
#' algorithm="nnet",
#' alpha=0.05,
#' k=20)
#'
#' results


cdgd <-  function(Y,D,G1,G2,Q=NULL,X,data,weight=NULL,alpha=0.05,k=500,t=0.05,algorithm) {

  item_point <- cdgd0(Y=Y,D=D,G1=G1,G2=G2,Q=Q,X=X,data=data,weight=weight,t=t,algorithm=algorithm)

  data_boot <- replicate(k, data[sample(1:nrow(data), nrow(data), replace=TRUE),], simplify=FALSE)
  output_boot <- sapply(data_boot, function (x) {
    cdgd0(Y=Y,D=D,G1=G1,G2=G2,Q=Q,X=X,data=x,weight=weight,t=t,algorithm=algorithm)[,2]
  }
    )

  se <- sqrt( (1/k)*rowSums(  (output_boot-rowMeans(output_boot))^2  ) )

  upper <- sapply(1:nrow(output_boot), function(x) output_boot[x,][order(as.numeric(output_boot[x,]))][floor(k*(1-alpha/2))])
  lower <- sapply(1:nrow(output_boot), function(x) output_boot[x,][order(as.numeric(output_boot[x,]))][ceiling(k*(alpha/2))])

  output_df <- as.data.frame(cbind(item_point, se, lower, upper))

  return(output_df)

}

