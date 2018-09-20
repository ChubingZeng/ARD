#' Automatic Relevance Determination Prediction
#'
#' \code{predict.ARD} returns predicted value
#' @param ARD.result ARD object
#' @param data.to.predict prediction data
#' @return predicted value
#' @examples
#' predict.ARD(ARD.result,x_test)

predict.ARD<-function(ARD.result,data.to.predict){
        return(data.to.predict%*%ARD.result$coefficients)
}
