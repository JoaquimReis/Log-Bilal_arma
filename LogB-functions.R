# log-Bilal Distribuição


# PROBABILITY DENSITY FUNCTION

dLB <- function(y, theta= .15)
{
  
  fx1 <- 6 / theta * y^(2 / theta-1) * (1 - y^(1/ theta))
  
  return(fx1)
  
}  

integrate(dLB,0,1) 


# CUMULATIVE DISTRIBUTION FUNCTION 

pLB <-  function(y, theta = 0.15) 
{
  cdf <- 3*y^(2/theta)-2*y^(3/theta)
  
  return(cdf)
}

pLB(.9999,2)
integrate(dLB, 0, .75,theta=2)


# QUANTILE FUNCTION

qLB<-function(u,theta=0.5)
{
  q <- 
    return(q)
}

u=pLB(.82)
qLB(u,theta = 0.5)


# INVERSION METHOD FOR RANDOM GENERATION

rLB <- function(n, alpha = 0.5) {
  u <- runif(n)
  y <- 1-(sqrt((1-u)^(1/alpha)))
  return(y)
}

rLB(100)