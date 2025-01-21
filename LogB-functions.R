# log-Bilal Distribuição


# PROBABILITY DENSITY FUNCTION

# theta <- (mu / (mu + 24))^(-1/2) - 5

dLB <- function(y, mu = 0.5) {
  
  fx1 <- (12 / ((mu / (mu + 24))^(-1/2) - 5)) * y^(4 / ((mu / (mu + 24))^(-1/2) - 5) - 1) * (1 - y^(2 / ((mu / (mu + 24))^(-1/2) - 5)))
  
  return(fx1)
}


integrate(dLB,0,1) 


# CUMULATIVE DISTRIBUTION FUNCTION 

pLB <- function(y, mu = 0.5) {
  
  cdf <- 3 * y^(2 / ((mu / (mu + 24))^(-1/2) - 5)) - 2 * y^(3 / ((mu / (mu + 24))^(-1/2) - 5))
  
  return(cdf)
}

pLB(0.25)


integrate(dLB, 0, 0.5)



# QUANTILE FUNCTION

# qLB<-function(u=0.5, mu=0.5)
# {
#   theta <- (mu / (mu + 24))^(-1/2) - 5
#   
#   erro <- 2 * sqrt(u^2 - u) - 2 * u + 1
#   
#   q <- (2 / theta * (1 + exp^(1/3) + 1 / exp^(1/3)))^theta
#   
#   return(q)
#   
# }
# u=pLB(.82)
# qLB(u,mu=.5) #In sqrt(u^2 - u) : NaNs produzidos


#ctrl+shift C para comentar e descomentar tudo marcado


#Como existem apenas dois numeros no Reais (0,1), busca outro metodo para gerar numeros aleatorios
  

# METHOD FOR RANDOM GENERATION

# numbers<-function(n, phi, theta){
#   inverse <- function(f, lower, upper){
#     function(y){
#       uniroot(function(x){f(x) - y}, lower = lower, upper = upper, tol=1e-3)$root
#     }
#   }
#   qLB <- inverse(pLB, lower = 0, upper = 1)
#   uniform_random <- runif(n, min = 0, max = 1)
#   y<-sapply(uniform_random, qLB)
#   return(y)
# }
# 
# numbers(100)


inverse <- function(f, lower, upper){
  function(y){
    uniroot(function(x){f(x) - y}, lower = lower, upper = upper, tol=1e-3)$root
  }
}

qLB <- function(u,mu){
  f<- function(y) {3 * y^(2 / ((mu / (mu + 24))^(-1/2) - 5)) - 2 * y^(3 / ((mu / (mu + 24))^(-1/2) - 5))}
  q<-inverse(f, lower = 0, upper = 1)
  return(q(u))
}

# y=.6
# mu=.15
# u=pLB(y,mu)
# 
# # qLB(u,rep(mu,n))
# # 
# # rUGO <- function(n, mu=.5) {
# #   u <- runif(n)
# #   y <- sapply(u,mu, qLB)
# #   return(y)
# # }
# 
# u <- runif(1)
# qLB(u,.8)


