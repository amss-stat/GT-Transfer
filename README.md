# GT-Transfer
We present an integrated multiple knowledge framework to improve the prevalence estimation for group testing data.
Unlike the existing methods, the proposed method built on the public source summary statistics because of the privacy and ethnicity restrictions.
Moreover, we extend the method to retesting design which can further imporve the estimation ability.

## Installation
 ```
library(devtools)
install_github("amss-stat/GT-Transfer")
 ```
## Quick start
 ```
library(GTT)
set.seed(4321)

p=0.01
k=5
m=500
pi_1=0.95
pi_0=0.95
S1=0.999
S0=0.999

d=matrix(rbinom(k*m,size=1,prob=p), nrow=k)
gresult=c()
for (h in 1:m){
  if(sum(d[,h])>0){
    gresult=c(gresult, rbinom(1,1,pi_1))
  }else{
    gresult=c(gresult, rbinom(1,1,1-pi_0))
  }
}
m_1=sum(gresult)
  
m_1i=rep(0,k+1)
retest=which(gresult==1)
R=0
for (r in retest){
  R=0
  for(d_1 in d[,r]){
    if(d_1>0){
      R=R+rbinom(1,1,S1)
    }else{
      R=R+rbinom(1,1,1-S0)
    }
  }
  m_1i[R+1]=m_1i[R+1]+1
}

g_tilde_j=c()
p_tilde_j=c()
v_tilde_j=c()
j=1
while(j<6){
  k_j=sample(c(5,10),size=1)
  m_j=sample(2000:5000,size=1,replace = T)
  pi_1_j=runif(1, min = 0.9, max = 0.999)
  pi_0_j=runif(1, min = 0.9, max = 0.999)
  p_j=p
  g_j=pi_1_j-(pi_0_j+pi_1_j-1)*(1-p_j)^k_j
  x=rbinom(m_j,1,g_j)
  if(sum(x)/m_j > 1-pi_0_j & sum(x)/m_j < pi_1_j){
    p_tilde_j[j]=1-((pi_1_j-sum(x)/m_j)/(pi_0_j+pi_1_j-1))^(1/k_j)
    v_tilde_j[j]=(sum(x)/m_j)*(1-sum(x)/m_j)*(1-p_tilde_j[j])^(2-2*k_j)/((pi_1_j+pi_0_j-1)^2*k_j^2*m_j)
    j=j+1
  }
}

tran = GT_tran(m,k,pi_0,pi_1,m_1,p_tilde_j,v_tilde_j,w = 1,eta = NULL)
print(tran)

tran_retest = GT_tran_retest(m,k,pi_0,pi_1,S0,S1,m_1,m_1i,p_tilde_j,v_tilde_j,w = 1,eta = NULL)
print(tran_retest)
 ```
