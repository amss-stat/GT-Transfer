#' @title Transfer learning in group testing with retesting information
#'
#' @description
#' GT_tran_retest is used to estimate the disease prevalence in group testing with auxiliary summary statistics and retesting information.
#'
#' @param m Numeric, group number
#' @param k Numeric, group size
#' @param pi_0 Numeric, the specificity in group level
#' @param pi_1 Numeric, the sensitivity in group level
#' @param S0 Numeric, the specificity in individual level
#' @param S1 Numeric, the sensitivity  in individual level
#' @param m_1 Numeric, number of groups which has been test positive
#' @param m_1i Vector, number of positive groups which has i individuals are tested to be positive. Note that the sum of m_1i must be equal to m1, and the length of m_1i must be k+1
#' @param p_tilde_j Vector, the estimation of p from source data
#' @param v_tilde_j Vector, the asymptotic variance of p_tilde_j
#' @param w Numeric, the weight which measure the importance of type 1 error and type 2 error in the first strategy of screening method. Default is 1.
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
#' \item{m_1i}{Number of positive groups which has i individuals are tested to be positive}
#' \item{g_tilde}{The estimation of probability that a group reported to be positive without transfer leaning and retesting}
#' \item{p_tilde}{The estimation of disease prevalence without transfer leaning and retesting}
#' \item{v_tilde}{The variance of p_tilde}
#' \item{g_grave}{The estimation of probability that a group reported to be positive with retesting information}
#' \item{p_grave}{The estimation of disease prevalence with retesting information}
#' \item{v_grave}{The variance of p_grave}
#' \item{lambda}{The parameter to evaluate the importance between target information and source information}
#' \item{g_breve}{The estimation of probability that a group reported to be positive after transfer leaning and introducing retesting information}
#' \item{p_breve}{The estimation of disease prevalence after transfer leaning and introducing retesting information}
#' \item{v_breve}{The variance of p_hat}
#' @import stats
#' @examples
#' tran_retest = GT_tran_retest(500,5,0.95,0.95,0.999,0.999,41,c(23,18,0,0,0,0),
#' c(0.009627190,0.009347203,0.008905579,0.009290185,0.010037513),
#' c(9.083533e-07,2.552419e-07,6.768645e-07,5.780736e-07,8.428058e-07))
#' tran_retest
#' @export

GT_tran_retest = function(m,k,pi_0,pi_1,S0,S1,m_1,m_1i,p_tilde_j,v_tilde_j,w = 1,b1 = 0.01,b2 = 0.01,eta = NULL,lambda = NULL){
  P1_i=function(g){
    if(length(g)==1){
      P1_i=rep(0,k+1)
      for(I in 0:k){
        s=0
        for(a in 0:I){
          for(b in 0:(k-I)){
            s= s + pi_1 * factorial(I) * factorial(k-I) * ((1-((pi_1-g)/(pi_0+pi_1-1))^(1/k))^(a+b)) * ((pi_1-g)/(pi_0+pi_1-1))^((k-a-b)/k) *(S1^a) * ((1-S1)^b) * (S0^(k-I-b)) *((1-S0)^(I-a)) / ( factorial(a) * factorial(I-a) * factorial(b) * factorial(k-I-b) )
          }
        }
        P1_i[I+1]=( (pi_1-g) * (S0^(k-I)) * ((1-S0)^I) * (1-pi_0) /(pi_0+pi_1-1) + s - pi_1 * ((pi_1-g)/(pi_0+pi_1-1)) * ((S0)^(k-I)) * ((1-S0)^I)) *factorial(k)/(factorial(I)*factorial(k-I))
      }
      return(P1_i)
    }else{
      result=c()
      for (g_1 in g){
        P1_i=rep(0,k+1)
        for(I in 0:k){
          s=0
          for(a in 0:I){
            for(b in 0:(k-I)){
              s= s + pi_1 * factorial(I) * factorial(k-I) * ((1-((pi_1-g_1)/(pi_0+pi_1-1))^(1/k))^(a+b)) * ((pi_1-g_1)/(pi_0+pi_1-1))^((k-a-b)/k) *(S1^a) * ((1-S1)^b) * (S0^(k-I-b)) *((1-S0)^(I-a)) / ( factorial(a) * factorial(I-a) * factorial(b) * factorial(k-I-b) )
            }
          }
          P1_i[I+1]=( (pi_1-g_1) * (S0^(k-I)) * ((1-S0)^I) * (1-pi_0) /(pi_0+pi_1-1) + s - pi_1 * ((pi_1-g_1)/(pi_0+pi_1-1)) * ((S0)^(k-I)) * ((1-S0)^I)) *factorial(k)/(factorial(I)*factorial(k-I))
        }
        result=c(result,P1_i)
      }
      return(result)
    }
  }
  d_P1_i = function(g){
    if(length(g)==1){
      d_P1_i=rep(0,k+1)
      for(I in 0:k){
        s=0
        for(a in 0:I){
          if((k-I)>=max(c(1-a,0))){
            for(b in max(c(1-a,0)):(k-I)){
              s = s + pi_1 * factorial(I) * factorial(k-I) * (S1^a) * ((1-S1)^b) * (S0^(k-I-b)) * ((1-S0)^(I-a)) * ((1-((pi_1-g)/(pi_0+pi_1-1))^(1/k))^(a+b-1)) * (pi_1-g)^((-a-b)/k) * (1/(pi_0+pi_1-1))^(1-(a+b)/k) * (k*((pi_1-g)/(pi_0+pi_1-1))^(1/k) -k+a+b)  / ( factorial(a) * factorial(I-a) * factorial(b) * factorial(k-I-b) * k )
            }
          }
        }
        d_P1_i[I+1]=(-(S0^(k-I)) * ((1-S0)^I) * (1-pi_0) /(pi_0+pi_1-1) + s) *factorial(k)/(factorial(I)*factorial(k-I))
      }
      return(d_P1_i)
    }else{
      out=c()
      for(g_1 in g){
        d_P1_i=rep(0,k+1)
        for(I in 0:k){
          s=0
          for(a in 0:I){
            if((k-I)>=max(c(1-a,0))){
              for(b in max(c(1-a,0)):(k-I)){
                s = s + pi_1 * factorial(I) * factorial(k-I) * (S1^a) * ((1-S1)^b) * (S0^(k-I-b)) * ((1-S0)^(I-a)) * ((1-((pi_1-g_1)/(pi_0+pi_1-1))^(1/k))^(a+b-1)) * (pi_1-g_1)^((-a-b)/k) * (1/(pi_0+pi_1-1))^(1-(a+b)/k) * (k*((pi_1-g_1)/(pi_0+pi_1-1))^(1/k) -k+a+b)  / ( factorial(a) * factorial(I-a) * factorial(b) * factorial(k-I-b) * k )
              }
            }
          }
          d_P1_i[I+1]=(-(S0^(k-I)) * ((1-S0)^I) * (1-pi_0) /(pi_0+pi_1-1) + s) *factorial(k)/(factorial(I)*factorial(k-I))
        }
        out=c(out,d_P1_i)
      }
      return(out)
    }
  }
  sF=function(D,k){
    result=c()
    for(kk in 1:(length(D)/k)){
      result=c(result,sum(D[((kk-1)*k+1):(kk*k)]))
    }
    return(result)
  }
  if(length(p_tilde_j) != length(v_tilde_j)){
    stop("Error: The length of p_tilde_j and v_tilde_j must be same")
  }
  if(sum(m_1i)!=m_1){
    stop("Error: The retest sample size isn't equal to m_1")
  }
  if(length(m_1i)!=k+1){
    stop("Error: The retest sample size isn't equal to m_1")
  }
  result = list()
  class(result) = c("GT_tran_retest", "list")
  result$m = m
  result$k = k
  result$pi_0 = pi_0
  result$pi_1 = pi_1
  result$m_1 = m_1
  result$m_1i = m_1i
  result$S0 = S0
  result$S1 = S1

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

  loss1=function(g) -sum(log(P1_i(g))*m_1i) - (m - m_1)*(log(1-g))
  g_grave=optimize(loss1,interval=c(1-pi_0,pi_1))
  p_grave=1-((pi_1-min(pi_1,max(1-pi_0,g_grave$minimum)))/(pi_0+pi_1-1))^(1/k)
  result$g_grave = g_grave$minimum
  result$p_grave = p_grave
  result$v_grave = ((1-result$p_breve)^(2-2*result$k)) * ( (result$g_breve/(result$m*(1-result$g_breve)) +sum( (d_P1_i(result$g_breve)^2) / (result$m*P1_i(result$g_breve)) ) - 3/result$m)) / (1/(1-result$g_breve) + sum( (d_P1_i(result$g_breve)^2) / (P1_i(result$g_breve)) ) * result$k*(result$pi_0+result$pi_1-1))^2

  g_tilde_j = c()
  for (j in 1:length(p_tilde_j)){
    g_tilde_j[j]=pi_1-(pi_0+pi_1-1)*(1-p_tilde_j[j])^k
  }
  data=data.frame(p_tilde_j,v_tilde_j,g_tilde_j)
  data=data[order(data$v_tilde_j),]

  if(v_tilde>0){
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
    if(dim(data)[1]!=0){
      taus=c()
      for(s in 1:dim(data)[1]){
        taus[s]=data$v_tilde_j[s]/(v_tilde*s)
      }
      data$taus = taus
      if(min(taus)<1){
        data=data[1:which.min(taus),]
      }else{
        data = NULL
      }
      result$data = data

      if(is.null(data)){
        loss1=function(g)  -sum(log(P1_i(g))*m_1i) - (m - m_1)*(log(1-g))
        g_breve=optimize(loss1,interval=c(1-pi_0,pi_1))
        p_breve=1-((pi_1-min(pi_1,max(1-pi_0,g_breve$minimum)))/(pi_0+pi_1-1))^(1/k)
        result$lambda = 0
        result$g_breve = g_breve$minimum
        result$p_breve = p_breve
        result$v_breve = ((1-result$p_breve)^(2-2*result$k)) * ( (result$g_breve/(result$m*(1-result$g_breve)) +sum( (d_P1_i(result$g_breve)^2) / (result$m*P1_i(result$g_breve)) ) - 3/result$m)) / (1/(1-result$g_breve) + sum( (d_P1_i(result$g_breve)^2) / (P1_i(result$g_breve)) ) * result$k*(result$pi_0+result$pi_1-1))^2
        return(result)
      }else{
        if(is.null(lambda)){
          lambda=1/sqrt(min(taus))
        }
        loss2=function(g,lambda) - sum(log(P1_i(g))*m_1i/m) - (log(1-g))*(m-m_1)/m + lambda*sum(log(P1_i(g))*P1_i(g)) +lambda*(1-g)*(log(1-g)) - lambda*sum(P1_i(g)*log(P1_i(data$g_tilde_j)))/dim(data)[1] - lambda*sum((1-g)*log(1-data$g_tilde_j))/dim(data)[1]
        g_breve=optimize(loss2,interval=c(1-pi_0,pi_1),lambda=lambda)
        p_breve=1-((pi_1-min(pi_1,max(1-pi_0,g_breve$minimum)))/(pi_0+pi_1-1))^(1/k)
        result$lambda = lambda
        result$g_breve = g_breve$minimum
        result$p_breve = p_breve
        result$v_breve = ((1-result$p_breve)^(2-2*result$k)) * ( result$g_breve/(result$m*(1-result$g_breve)) - 3/result$m +sum( (d_P1_i(result$g_breve)^2) / (result$m*P1_i(result$g_breve)) ) + (result$lambda^2)* (sum( (( 1/(1-result$data$g_tilde_j) + sF(d_P1_i(result$data$g_tilde_j) * d_P1_i(result$g_breve)/P1_i(result$data$g_tilde_j) ,result$k+1) )^2) * (result$k^2) * ((result$pi_0+result$pi_1-1)^2) * ((1-result$data$p_tilde_j)^(2*result$k-2)) * result$data$v_tilde_j ))/dim(result$data)[1]^2 ) / (( (1+result$lambda)/(1-result$g_breve) + sum( (result$lambda+1)*(d_P1_i(result$g_breve)^2)/P1_i(result$g_breve)) )*k*(result$pi_0+result$pi_1-1))^2
        return(result)
      }
    }else{
      loss1=function(g)  -sum(log(P1_i(g))*m_1i) - (m - m_1)*(log(1-g))
      g_breve=optimize(loss1,interval=c(1-pi_0,pi_1))
      p_breve=1-((pi_1-min(pi_1,max(1-pi_0,g_breve$minimum)))/(pi_0+pi_1-1))^(1/k)
      result$lambda = 0
      result$g_breve = g_breve$minimum
      result$p_breve = p_breve
      result$v_breve = ((1-result$p_breve)^(2-2*result$k)) * ( (result$g_breve/(result$m*(1-result$g_breve)) +sum( (d_P1_i(result$g_breve)^2) / (result$m*P1_i(result$g_breve)) ) - 3/result$m)) / (1/(1-result$g_breve) + sum( (d_P1_i(result$g_breve)^2) / (P1_i(result$g_breve)) ) * result$k*(result$pi_0+result$pi_1-1))^2
      return(result)
    }
  }else{
    loss3=function(g)  sum((1-g)*(log(1-g) - log(1-data$g_tilde_j)) + P1_i(g)*(log(P1_i(g)) - log(P1_i(data$g_tilde_j))))
    g_breve=optimize(loss3,interval=c(1-pi_0,pi_1))
    p_breve=1-((pi_1-min(pi_1,max(1-pi_0,g_breve$minimum)))/(pi_0+pi_1-1))^(1/k)
    result$data = data
    result$lambda = Inf
    result$g_breve = g_breve$minimum
    result$p_breve = p_breve
    result$v_breve = ((1-result$p_breve)^(2-2*result$k)) * ( (sum( (( 1/(1-result$data$g_tilde_j) + sF(d_P1_i(result$data$g_tilde_j) * d_P1_i(result$g_breve)/P1_i(result$data$g_tilde_j) ,result$k+1) )^2) * (result$k^2) * ((result$pi_0+result$pi_1-1)^2) * ((1-result$data$p_tilde_j)^(2*result$k-2)) * result$data$v_tilde_j ))/dim(result$data)[1]^2 ) / (( 1/(1-result$g_breve) + sum( (d_P1_i(result$g_breve)^2)/P1_i(result$g_breve)) )*k*(result$pi_0+result$pi_1-1))^2
    return(result)
  }
}

