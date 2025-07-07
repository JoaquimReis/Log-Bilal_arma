# log-Bilal Distribuição

# ---------------------------------------------------------------------------
# FUNÇÃO DE DENSIDADE DE PROBABILIDADE (PDF)
# ---------------------------------------------------------------------------

dLB <- function(y, mu = 0.75) {
  
  fx1 <- (12 / ((mu / (mu + 24))^(-1/2) - 5)) * y^(4 / ((mu / (mu + 24))^(-1/2) - 5) - 1) * (1 - y^(2 / ((mu / (mu + 24))^(-1/2) - 5)))
  
  return(fx1)
}


#integrate(dLB,0,1) 


# ---------------------------------------------------------------------------
# FUNÇÃO DE DISTRIBUIÇÃO ACUMULADA (CDF)
# ---------------------------------------------------------------------------
pLB <- function(y, mu = 0.75) {
  
  cdf <- 3 * y^(4 / ((mu / (mu + 24))^(-1/2) - 5)) - 2 * y^(6 / ((mu / (mu + 24))^(-1/2) - 5))
  
  return(cdf)
}

# pLB(.77)
# 
# integrate(dLB, 0, .77) 


# Função inversa (quantílica) Log-Bilal via método da inversão numérica
#------------------------------------------------------------------------------------------
inverse <- function(f, lower, upper){
  function(y){
    uniroot(function(x){f(x) - y}, lower = lower, upper = upper, tol = 1e-10)$root
  }
}

qLB <- function(u, mu){
  if (any(u <= 0 | u >= 1)) stop(paste("u must be between 0 and 1", "\n", ""))
  
  f<- function(y) {3 * y^(4 / ((mu / (mu + 24))^(-1/2) - 5)) - 2 * y^(6 / ((mu / (mu + 24))^(-1/2) - 5))
  }
  q<-inverse(f, lower = 0, upper = 1)
  return(q(u))
}

# y=.7
# mu=.15
# u=pLB(y,mu)
# qLB(u,mu)

#------------------------------------------------------------------------------------------
# Função de geração de números aleatórios - Distribuição Log-Bilal
#------------------------------------------------------------------------------------------
rLB <- function(n, mu) {
  u <- runif(n)
  y <- mapply(function(ui, mui) qLB(ui, mu = mui), u, mu)
  return(y)
}

# Exemplo de uso:
# y <- rLB(1000, mu = 0.6)
# hist(y, breaks = 30, main = "Amostras da distribuição Log-Bilal")
