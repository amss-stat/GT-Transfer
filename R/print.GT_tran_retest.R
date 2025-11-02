#' @title Transfer learning in group testing with retesting information
#'
#' @description
#' Print a summary of GT_tran_retest.
#'
#' @param x An object of class `GT_tran_retest`.
#' @param ... Additional arguments passed to the print method.
#'
#' @import stats
#' @examples
#' tran_retest = GT_tran_retest(500,5,0.95,0.95,0.999,0.999,41,c(23,18,0,0,0,0),
#' c(0.009627190,0.009347203,0.008905579,0.009290185,0.010037513),
#' c(9.083533e-07,2.552419e-07,6.768645e-07,5.780736e-07,8.428058e-07))
#' print(tran_retest)
#' @export

print.GT_tran_retest = function(x,...){
  if(!inherits(x,'GT_tran_retest')){
    stop("Object must be of class 'GT_tran_retest'")
  }
  cat('* Model: Transfer learning in group testing with retesting information * \n')
  if(is.null(x$data)){
    cat('No information will be transferred! \n')
    cat('The estimation of p is :', x$p_grave)
  }else{
    cat('The transferred setting is :', sort(as.numeric(rownames(x$data))), '\n')
    cat('The estimation of p is :', x$p_breve)
  }
}
