#' @title Transfer learning in group testing scenario
#'
#' @description
#' Print a summary of GT_tran.
#'
#' @param x An object of class `GT_tran`.
#' @param ... Additional arguments passed to the print method.
#'
#' @import stats
#' @examples
#' tran = GT_tran(500,5,0.95,0.95,41,
#' c(0.009627190,0.009347203,0.008905579,0.009290185,0.010037513),
#' c(9.083533e-07,2.552419e-07,6.768645e-07,5.780736e-07,8.428058e-07))
#' print(tran)
#' @export

print.GT_tran = function(x, ...){
  if(!inherits(x,'GT_tran')){
    stop("Object must be of class 'GT_tran'")
  }
  cat('* Model: Transfer learning in group testing * \n')
  if(is.null(x$data)){
    cat('No information will be transferred! \n')
    cat('The estimation of p is :', x$p_hat)
  }else{
    cat('the transferred setting is :', sort(as.numeric(rownames(x$data))), '\n')
    cat('The estimation of p is :', x$p_hat)
  }
}

