#'Gives estimated values of SSI
#'@param x1 is sample from Exp Population 1.
#'@param x2 is sample from Exp Population 2.
#'@param x3 is sample from Exp Population 3.
#'@param n is number of sample
#'@export
CLEXP3EST<-function(x1,x2,x3,n){
  k=3


  #library(devtools)
  #library(testthat)

  P_MLE<-rep(1,k)  ########################
  P_U<- rep(1,k)    #########################
  P_improved_U<- rep(1,k)  ########################
  P_UMVUE<- rep(1,k)  ########################

  P_BSE<- rep(1,k)  ########################
  P_improved_BSE<- rep(1,k)  ########################
  P_Bayes<- rep(1,k)  ########################

  #library(tibble)
  #library(gsl)
  #library(cubature) # load the package "cubature"
  t<-rep(1,k)
  t[1]<-mean(x1)-min(x1,x2,x3)
  t[2]<-mean(x2)-min(x1,x2,x3)
  t[3]<-mean(x3)-min(x1,x2,x3)


  formu<-function(a){
    formu_alpha<-rep(1,k)
    formu_original<-rep(1,k)
    for (i in 1:k) {
      for (j in 1:k) {
        formu_alpha[j]<-a[i]/(a[i]+a[j])

      }
      formu_original[i]<-1-((1/k)*(sum(formu_alpha)))

    }
    b<-(formu_original)
    b
  }
  sum_tulika<-sum(1/t)

  #################################  MLE #############################################
  mu_MLE<-min(x1,x2,x3)  ########################
  alpha_MLE<- rep(1,k)  ########################
  t_MLE<-rep(1,k)  ########################
  for (i in 1:k) {  ########################
    for (j in 1:k) {  ########################
      alpha_MLE[j]<-t[i]/(t[i]+t[j])  ########################
      ########################
    }  ########################
    t_MLE[i]<-1-((1/k)*(sum(alpha_MLE)))  ########################

  }  ########################


  #                               theta                              #  ########################
  theta_U<-rep(1,k)  ########################
  theta_improved_U<-rep(1,k)
  alpha_U<-rep(1,k)
  t_U<-rep(1,k)
  alpha_improved_U<-rep(1,k)
  t_improved_U<-rep(1,k)
  theta_improved_BSE<-rep(1,k)

  for (i in 1:k) {
    theta_U[i]<-t[i]+(1/((n-1)*sum_tulika))
    theta_improved_U[i]<-((n/(n+1))*t[i])+(1/((n-1)*sum_tulika))
    theta_improved_BSE[i]<-((n-1)/n)*theta_U[i]
  }
  ########################

  ########################





  #############################################################################################
  for (i in 1:k) {
    for (j in 1:k) {
      alpha_U[j]<-theta_U[i]/(theta_U[i]+theta_U[j])
      alpha_improved_U[j]<-theta_improved_U[i]/(theta_improved_U[i]+theta_improved_U[j])

    }
    t_U[i]<-1-((1/k)*(sum(alpha_U)))
    t_improved_U[i]<-1-((1/k)*(sum(alpha_improved_U)))

  }
  ######################################### BSE ####################################

  theta_BSE<-rep(1,k)
  theta_BSE[1]<-mean(x1)-min(x1)
  theta_BSE[2]<-mean(x2)-min(x2)
  theta_BSE[3]<-mean(x3)-min(x3)
  for (i in 1:k) {
    P_MLE[i]<-formu(t)[i]
    P_U[i]<- formu(theta_U)[i]
    P_improved_U[i]<-formu(theta_improved_U)[i]
    P_BSE[i]<-formu(theta_BSE)[i]
    P_improved_BSE[i]<-formu(theta_improved_BSE)[i]
  }














  ########################################### UMVUE TOTAL ##########################################
  mult<-function(a,b){
    ((a*b)^(n-1))/(n-1) # 2,1 for p_21
  }
  prodx<-function(a,b,c){
    (n-1)^2/((a^(n-1))*(b^(n-1))*c*sum_tulika)
  } # a=x1, b=x2,c=x3 for p_21,

  prod2<-function(a,b){
    ((a*b)^(n-1))/(n-1)^2
  }



  firstterm<-function(b){
    1/(b*sum_tulika)  # 2 for 21
  }





  d_matrix<-matrix(1,k,k)
  for(i in 1:k){
    for(j in 1:k){
      d_matrix[i,j]<-t[i]-t[j]
    }
  }



  int_21<-function(x){
    (x^(n-2))*(x+t[2]-t[1])^(n-1)
  }
  int_31<-function(x){
    (x^(n-2))*(x+t[3]-t[1])^(n-1)
  }
  int_12<-function(x){
    (x^(n-2))*(x+t[1]-t[2])^(n-1)
  }
  int_32<-function(x){
    (x^(n-2))*(x+t[3]-t[2])^(n-1)
  }
  int_13<-function(x){
    (x^(n-2))*(x+t[1]-t[3])^(n-1)
  }
  int_23<-function(x){
    (x^(n-2))*(x+t[2]-t[3])^(n-1)
  }

  dint_21<-function(x){
    (x^(n-2))*(x-d_matrix[1,2])^(n-1)
  }
  dint_31<-function(x){
    (x^(n-2))*(x-d_matrix[1,3])^(n-1)
  }
  dint_12<-function(x){
    (x^(n-2))*(x-d_matrix[2,1])^(n-1)
  }
  dint_32<-function(x){
    (x^(n-2))*(x-d_matrix[2,3])^(n-1)
  }
  dint_13<-function(x){
    (x^(n-2))*(x-d_matrix[3,1])^(n-1)
  }
  dint_23<-function(x){
    (x^(n-2))*(x-d_matrix[3,2])^(n-1)
  }


  if(d_matrix[1,2]>0){
    s21<- ((prodx(t[1],t[2],t[3]))*(  prod2(t[1],t[2])  - ((integrate(dint_21,d_matrix[1,2],t[1])$value)/(n-1))  ))+firstterm(t[2])
  } else
    s21<- ( (prodx(t[1],t[2],t[3]))*(   prod2(t[1],t[2]) - ( (integrate(int_21,0,t[1])$value)/(n-1)   )      ))    +firstterm(t[2])



  if(d_matrix[1,3]>0){
    s31<- ((prodx(t[1],t[3],t[2]))*(  prod2(t[1],t[3])  - ((integrate(dint_31,d_matrix[1,3],t[1])$value)/(n-1))  ))+firstterm(t[3])
  } else
    s31<- ( (prodx(t[1],t[3],t[2]))*(   prod2(t[1],t[3]) - ( (integrate(int_31,0,t[1])$value)/(n-1)   )      ))    +firstterm(t[3])




  if(d_matrix[2,1]>0){
    s12<- ((prodx(t[2],t[1],t[3]))*(  prod2(t[2],t[1])  - ((integrate(dint_12,d_matrix[2,1],t[2])$value)/(n-1))  ))+firstterm(t[1])
  } else
    s12<- ( (prodx(t[2],t[1],t[3]))*(   prod2(t[2],t[1]) - ( (integrate(int_12,0,t[2])$value)/(n-1)   )      ))    +firstterm(t[1])

  if(d_matrix[2,3]>0){
    s32<- ((prodx(t[2],t[3],t[1]))*(  prod2(t[2],t[3])  - ((integrate(dint_32,d_matrix[2,3],t[2])$value)/(n-1))  ))+firstterm(t[3])
  } else
    s32<- ( (prodx(t[2],t[3],t[1]))*(   prod2(t[2],t[3]) - ( (integrate(int_32,0,t[2])$value)/(n-1)   )      ))    +firstterm(t[3])


  if(d_matrix[3,1]>0){
    s13<- ((prodx(t[3],t[1],t[2]))*(  prod2(t[3],t[1])  - ((integrate(dint_13,d_matrix[3,1],t[3])$value)/(n-1))  ))+firstterm(t[1])
  } else
    s13<- ( (prodx(t[3],t[1],t[2]))*(   prod2(t[3],t[1]) - ( (integrate(int_13,0,t[3])$value)/(n-1)   )      ))    +firstterm(t[1])
if(d_matrix[3,2]>0){
  s23<- ((prodx(t[3],t[2],t[1]))*(  prod2(t[3],t[2])  - ((integrate(dint_23,d_matrix[3,2],t[3])$value)/(n-1))  ))+firstterm(t[2])
} else
  s23<- ( (prodx(t[3],t[2],t[1]))*(   prod2(t[3],t[2]) - ( (integrate(int_23,0,t[3])$value)/(n-1)   )      ))    +firstterm(t[2])




u1 <-(5/6)-(1/3)*(s21+s31)
u2 <-(5/6)-(1/3)*(s12+s32)
u3 <-(5/6)-(1/3)*(s13+s23)


##################################### Bayes ###############################################
s<-n*t
constant<-function(a,b,c){
  (gamma(n+3)*(gamma(n+2))^2)/((a^(n+3))*((b*c)^(n+2)))
}
K_Tulika<-1/(   constant(s[1],s[2],s[3])+constant(s[2],s[1],s[3])+constant(s[3],s[2],s[1]))
K_Tulika

#####################i=2 ,l=1#################################

# x[1] corresponds to l
# x[2] corresponds to i


first_term_bayes1<-function(x){
  ((x[1]+x[2])^(-1))*(x[1]^(n+1))*(x[2]^(n+2))*(exp(-s[1]*x[1]))*(exp(-s[2]*x[2]))
}

term11<-(adaptIntegrate(first_term_bayes1, lowerLimit = c(0, 0), upperLimit = c(Inf, Inf))$integral)*gamma(n+3)*K_Tulika/(s[3]^(n+3))
term21<-((gamma(n+3))*(gamma(n+2))^2)*K_Tulika/((s[2]^(n+3))*((s[1]*s[3])^(n+2)))
l1i2 <-term11+term21

############################ i=3, l=1
second_term_bayes1<-function(x){
  ((x[1]+x[2])^(-1))*(x[1]^(n+1))*(x[2]^(n+2))*(exp(-s[1]*x[1]))*(exp(-s[3]*x[2]))
}

term31<-(adaptIntegrate(second_term_bayes1, lowerLimit = c(0, 0), upperLimit = c(Inf, Inf))$integral)*gamma(n+3)*K_Tulika/(s[2]^(n+3))
term41<-((gamma(n+3))*(gamma(n+2))^2)*K_Tulika/((s[3]^(n+3))*((s[1]*s[2])^(n+2)))
l1i3 <-term31+term41
bayes1 <- (5/6)-(l1i2 +l1i3 )/3




#####################i=1 ,l=2#################################
# x[1] corresponds to l
# x[2] corresponds to i



first_term_bayes2<-function(x){
  ((x[1]+x[2])^(-1))*(x[1]^(n+1))*(x[2]^(n+2))*(exp(-s[2]*x[1]))*(exp(-s[1]*x[2]))
}

term12<-(adaptIntegrate(first_term_bayes2, lowerLimit = c(0, 0), upperLimit = c(Inf, Inf))$integral)*gamma(n+3)*K_Tulika/(s[3]^(n+3))
term22<-((gamma(n+3))*(gamma(n+2))^2)*K_Tulika/((s[1]^(n+3))*((s[2]*s[3])^(n+2)))
term22
term12
l2i1 <-term12+term22

############################ i=3, l=2

second_term_bayes2<-function(x){
  ((x[1]+x[2])^(-1))*(x[1]^(n+1))*(x[2]^(n+2))*(exp(-s[2]*x[1]))*(exp(-s[3]*x[2]))
}

term32<-(adaptIntegrate(second_term_bayes2, lowerLimit = c(0, 0), upperLimit = c(Inf, Inf))$integral)*gamma(n+3)*K_Tulika/(s[1]^(n+3))
term42<-((gamma(n+3))*(gamma(n+2))^2)*K_Tulika/((s[3]^(n+3))*((s[1]*s[2])^(n+2)))
l2i3 <-term32+term42
bayes2 <- (5/6)-(l2i3 +l2i1 )/3


#####################i=1 ,l=3#################################
# x[1] corresponds to l
# x[2] corresponds to i



first_term_bayes3<-function(x){
  ((x[1]+x[2])^(-1))*(x[1]^(n+1))*(x[2]^(n+2))*(exp(-s[3]*x[1]))*(exp(-s[1]*x[2]))
}

term13<-(adaptIntegrate(first_term_bayes3, lowerLimit = c(0, 0), upperLimit = c(Inf, Inf))$integral)*gamma(n+3)*K_Tulika/(s[2]^(n+3))
term23<-((gamma(n+3))*(gamma(n+2))^2)*K_Tulika/((s[1]^(n+3))*((s[2]*s[3])^(n+2)))
l3i1 <-term13+term23
l3i1
############################ i=2, l=3

second_term_bayes3<-function(x){
  ((x[1]+x[2])^(-1))*(x[1]^(n+1))*(x[2]^(n+2))*(exp(-s[3]*x[1]))*(exp(-s[2]*x[2]))
}

term33<-(adaptIntegrate(second_term_bayes3, lowerLimit = c(0, 0), upperLimit = c(Inf, Inf))$integral)*gamma(n+3)*K_Tulika/(s[1]^(n+3))
term43<-((gamma(n+3))*(gamma(n+2))^2)*K_Tulika/((s[2]^(n+3))*((s[1]*s[3])^(n+2)))
l3i2 <-term33+term43
bayes3 <- (5/6)-(l3i1 +l3i2 )/3



######################################################################




P_UMVUE[1]<-u1   ########################
P_UMVUE[2]<-u2   ########################
P_UMVUE[3]<-u3   ########################

P_Bayes[1]<-bayes1
P_Bayes[2]<-bayes2
P_Bayes[3]<-bayes3

comp_data <- data.frame(


  P_i_for_i = c (1:3),


  MLE = c(P_MLE[1], P_MLE[2], P_MLE[3]),
  UMVUE = c(P_UMVUE[1], P_UMVUE[2], P_UMVUE[3]),
  Gen_Bayes=c(P_Bayes[1],P_Bayes[2],P_Bayes[3]),
  Plug_in_UMVUE=c(P_U[1],P_U[2],P_U[3]),
  Plug_in_Improved_UMVUE=c(P_improved_U[1],P_improved_U[2],P_improved_U[3])
  #AE= c(P_BSE[1],P_BSE[2],P_BSE[3])
)

comp_data











}
