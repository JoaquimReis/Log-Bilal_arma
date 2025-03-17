# log-Bilal Distribuição


# PROBABILITY DENSITY FUNCTION

dLB <- function(y, mu = 0.75) {
  
  fx1 <- (12 / ((mu / (mu + 24))^(-1/2) - 5)) * y^(4 / ((mu / (mu + 24))^(-1/2) - 5) - 1) * (1 - y^(2 / ((mu / (mu + 24))^(-1/2) - 5)))
  
  return(fx1)
}


#integrate(dLB,0,1) 


# CUMULATIVE DISTRIBUTION FUNCTION 

pLB <- function(y, mu = 0.75) {
  
  cdf <- 3 * y^(4 / ((mu / (mu + 24))^(-1/2) - 5)) - 2 * y^(6 / ((mu / (mu + 24))^(-1/2) - 5))
  
  return(cdf)
}

# pLB(.77)
# 
# integrate(dLB, 0, .77) 


# METHOD FOR RANDOM GENERATION

inverse <- function(f, lower, upper){
  function(y){
    uniroot(function(x){f(x) - y}, lower = lower, upper = upper, tol=1e-10)$root
  }
}

qLB <- function(u,mu){
  f<- function(y) {3 * y^(4 / ((mu / (mu + 24))^(-1/2) - 5)) - 2 * y^(6 / ((mu / (mu + 24))^(-1/2) - 5))
  }
  q<-inverse(f, lower = 0, upper = 1)
  return(q(u))
}

# y=.6
# mu=.15
# u=pLB(y,mu)
