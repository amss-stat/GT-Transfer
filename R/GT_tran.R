#' @title Transfer learning in group testing scenario
#'
#' @description
#' GT_tran is used to estimate the disease prevalence in group testing with auxiliary summary statistics.
#'
#' @param m Numeric, the group number
#' @param k Numeric, the group size
#' @param pi_0 Numeric, the specificity in group level
#' @param pi_1 Numeric, the sensitivity in group level
#' @param m_1 Integer, Number of groups which has been test positive
#' @param p_tilde_j Vector, the estimation of p from J source datasets
#' @param v_tilde_j Vector, the asymptotic variance of p_tilde_j
#' @param screen A logical value (T or F) indicating whether to screen the source data: If screen = TRUE, the source data will be screened using the method described in the article, based on subsequent parameters w, b1, and b2; If screen = F, no screening will be performed, and all source data will be used directly for transfer.
#' @param w The weight which measure the importance of type 1 error and type 2 error in the first strategy of screening method. Default is 1.
#' @param b1 Numeric, a small number to evaluate the type 1 error (Default is 0.01)
#' @param b2 Numeric, a small number to evaluate the type 2 error (Default is 0.01)
#' @param eta Numeric, the tolerent parameter of Condition 1 (Default is `\eqn{0.5 + log_m J - log_m \max\limits_{j}\{ \tilde{v}_j/\tilde{v} \} }`)
#' @param lambda Numeric, the parameter to evaluate the importance between target information and source information (Default is `\eqn{ \sqrt{\frac{s \tilde{v}}{\max\limits_{j}\{\tilde{v}_j\}}} }`)
#' @return A list contains all the statistics and the estimation after transfer learning.
#' \item{m}{The group number}
#' \item{k}{The group size}
#' \item{pi_0}{The specificity in group level}
#' \item{pi_1}{The sensitivity in group level}
#' \item{m_1}{Number of groups which has been test positive}
#' \item{g_tilde}{The estimation of probability that a group reported to be positive without transfer leaning}
#' \item{p_tilde}{The estimation of disease prevalence without transfer leaning}
#' \item{v_tilde}{The variance of p_tilde}
#' \item{lambda}{The parameter to evaluate the importance between target information and source information}
#' \item{g_hat}{The estimation of probability that a group reported to be positive after transfer leaning}
#' \item{p_hat}{The estimation of disease prevalence after transfer leaning}
#' \item{v_hat}{The variance of p_hat}
#' @import stats
#' @examples
#' tran = GT_tran(500,5,0.95,0.95,41,
#' c(0.009627190,0.009347203,0.008905579,0.009290185,0.010037513),
#' c(9.083533e-07,2.552419e-07,6.768645e-07,5.780736e-07,8.428058e-07))
#' tran
#' @export

GT_tran = function(m,k,pi_0,pi_1,m_1,p_tilde_j,v_tilde_j,screen = T,w = 1,b1 = 0.01,b2 = 0.01,eta = NULL,lambda = NULL){
  if(length(p_tilde_j) != length(v_tilde_j)){
    stop("Error: The length of p_tilde_j and v_tilde_j must be same")
  }
  result = list()
  class(result) = c("GT_tran", "list")
  result$m = m
  result$k = k
  result$pi_0 = pi_0
  result$pi_1 = pi_1
  result$m_1 = m_1

  g_tilde=m_1/m
  p_tilde=1-((pi_1-min(pi_1,max(1-pi_0,g_tilde)))/(pi_0+pi_1-1))^(1/k)

  result$g_tilde = g_tilde
  result$p_tilde = p_tilde
  if(p_tilde==1){
    v_tilde=0
  }else{
    v_tilde=g_tilde*(1-g_tilde)*(1-p_tilde)^(2-2*k)/((pi_1+pi_0-1)^2*k^2*m)
  }
  result$v_tilde = v_tilde

  g_tilde_j = c()
  for (j in 1:length(p_tilde_j)){
    g_tilde_j[j]=pi_1-(pi_0+pi_1-1)*(1-p_tilde_j[j])^k
  }
  data=data.frame(p_tilde_j,v_tilde_j,g_tilde_j)
  data=data[order(data$v_tilde_j),]

  if(v_tilde>0){
    if(screen ==T){
      if(is.null(eta)){
        taus=c()
        for(s in 1:length(p_tilde_j)){
          taus[s]=data$v_tilde_j[s]/(v_tilde*s)
        }
        eta = 0.5 + (log(dim(data)[1]) - min(taus))/log(m)
      }

      for(j in dim(data)[1]:1){
        m_eta = m^(-eta)
        alpha1 = function(z) 2 - 2*pnorm(z/sqrt(v_tilde+data$v_tilde_j[j]))
        alpha2 = function(z) 2*(pnorm((z-m_eta)/sqrt(v_tilde+data$v_tilde_j[j])) - pnorm((-m_eta)/sqrt(v_tilde+data$v_tilde_j[j])))
        L = function(z) w*alpha1(z) + (1-w)*alpha2(z)
        z_min = max(c(0,m_eta+sqrt(v_tilde+data$v_tilde_j[j])*qnorm((b2/2)+dnorm(-m_eta/sqrt(v_tilde+data$v_tilde_j[j])))))
        z_max = qnorm(1-b1/2)*sqrt(v_tilde+data$v_tilde_j[j])
        if(z_min < z_max){
          Z = optimize(L,interval=c(z_min,z_max),tol = (z_max-z_min)/10000)
          Z = Z$minimum
        }else{
          Z = z_max
        }

        if(abs(data$p_tilde_j[j] - p_tilde)>Z){
          data = data[-j,]
        }
      }
    }


    if(dim(data)[1]!=0){
      taus=c()
      for(s in 1:dim(data)[1]){
        taus[s]=data$v_tilde_j[s]/(v_tilde*s)
      }
      data$taus = taus

      if(screen ==T){
        if(min(taus)<1){
          data=data[1:which.min(taus),]
        }else{
          data = NULL
        }
      }

      result$data = data

      if(is.null(data)){
        result$lambda = 0
        result$g_hat = g_tilde
        result$p_hat = p_tilde
        result$v_hat = v_tilde
        return(result)
      }else{
        if(is.null(lambda)){
          lambda=1/sqrt(data$taus[dim(data)[1]])
        }
        loss1=function(g,lambda) - (log(g))*m_1/m - (log(1-g))*(m-m_1)/m + lambda*(sum(g*log(g) +(1-g)*(log(1-g)) - g*log(data$g_tilde_j) - (1-g)*log(1-data$g_tilde_j)))/dim(data)[1]

        g_hat=optimize(loss1,interval=c(0,1),lambda=lambda)
        p_hat=1-((pi_1-min(pi_1,max(1-pi_0,g_hat$minimum)))/(pi_0+pi_1-1))^(1/k)
        result$lambda = lambda
        result$g_hat = g_hat$minimum
        result$p_hat = p_hat
        result$v_hat = ( (1-result$p_hat)^(2-2*result$k) * result$g_hat^2 * (1-result$g_hat)^2  * ( 1/(result$m*result$g_hat*(1-result$g_hat)) + (result$lambda^2) * (result$k^2) * ((result$pi_0+result$pi_1-1)^2) * sum(result$data$v_tilde_j * ((1-result$data$p_tilde_j)^(2*result$k-2)) / ((result$data$g_tilde_j^2)*(1-result$data$g_tilde_j)^2)) / (dim(result$data)[1]^2))) / (result$k*(1+result$lambda)*(result$pi_0+result$pi_1-1))^2
        return(result)
      }

    }else{
      p_hat = p_tilde
      result$lambda = 0
      result$g_hat = g_tilde
      result$p_hat = p_tilde
      result$v_hat = v_tilde
      return(result)
    }
  }else{
    loss2=function(g) sum(g*log(g) +(1-g)*(log(1-g)) - g*log(data$g_tilde_j) - (1-g)*log(1-data$g_tilde_j))

    g_hat=optimize(loss2,interval=c(0,1))
    p_hat=1-((pi_1-min(pi_1,max(1-pi_0,g_hat$minimum)))/(pi_0+pi_1-1))^(1/k)
    result$data = data
    result$lambda = Inf
    result$g_hat = g_hat$minimum
    result$p_hat = p_hat
    result$v_hat = ( (1-result$p_hat)^(2-2*result$k) * result$g_hat^2 * (1-result$g_hat)^2  * ((result$k^2) * ((result$pi_0+result$pi_1-1)^2) * sum(result$data$v_tilde_j * ((1-result$data$p_tilde_j)^(2*result$k-2)) / ((result$data$g_tilde_j^2)*(1-result$data$g_tilde_j)^2)) / (dim(result$data)[1]^2))) / (result$k*(result$pi_0+result$pi_1-1))^2
    return(result)
  }
}
