# SIMULACAO DE MONTE CARLO - ARMA(1,1) SEM ESTIMATIVA DE MU

rm(list = ls())

source("simu.LogBarma.R")
source("LogBarmafit.R")

alpha = 1.5
phi = 0.2  # AR
theta = 0.4 # MA
true_values = c(alpha, phi, theta)  # alpha, phi, theta
vn = c(70, 150, 300)#, 500, 1000)
R = 50
z = 1.96

ar1 = 1
ma1 = 1


start_time <- Sys.time()

system.time({
  for (n in vn) {
    set.seed(5)
    # Matriz de resultados (3 colunas: alpha, phi, theta)
    estim <- ICi <- ICs <- err <- matrix(NA, nrow = R, ncol = length(true_values))
    
    # Contadores
    calpha <- cphi <- ctheta <- 0 
    bug <- 0  # Inicializa o contador de erros
    
    pb <- txtProgressBar(min = 0, max = R, style = 3)
    
    for (i in 1:R) {
      y <- simu.LogBarma(n, phi = phi, theta = theta, alpha = alpha, freq = 12, link = "logit")
      fit1 <- try(LogBarma.fit(y, ma=ma1, ar=ar1), silent = TRUE)
      
      if (!inherits(fit1, "try-error") && fit1$conv == 0 && nrow(fit1$model) == 3) {
        estim[i, ] <- fit1$model[, 1] 
        err[i, ] <- fit1$model[, 2]
        
        if (!any(is.na(estim[i, ])) && !any(is.na(err[i, ]))) {
          ICi[i, ] <- estim[i, ] - (z * err[i, ])
          ICs[i, ] <- estim[i, ] + (z * err[i, ])
          
          if (ICi[i, 1] <= alpha && ICs[i, 1] >= alpha) calpha <- calpha + 1
          if (ICi[i, 2] <= phi && ICs[i, 2] >= phi) cphi <- cphi + 1
          if (ICi[i, 3] <= theta && ICs[i, 3] >= theta) ctheta <- ctheta + 1
        }
      } else {
        bug <- bug + 1  # Incrementa contador de falhas
      }
      setTxtProgressBar(pb, i)
    }

          
          # Verifica se os parâmetros estão dentro de intervalos válidos
          if (ICi[i, 1] <= alpha && ICs[i, 1] >= alpha && 
              ICi[i, 2] >= -1 && ICs[i, 2] <= 1 && 
              ICi[i, 3] >= -1 && ICs[i, 3] <= 1) {
            if (ICi[i, 1] <= alpha && ICs[i, 1] >= alpha) calpha <- calpha + 1
            if (ICi[i, 2] <= phi && ICs[i, 2] >= phi) cphi <- cphi + 1
            if (ICi[i, 3] <= theta && ICs[i, 3] >= theta) ctheta <- ctheta + 1
          }
        }
      } else {
        bug <- bug + 1  # Incrementa contador de falhas
      }
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # Estatísticas
    m <- apply(estim, 2, mean, na.rm = TRUE)
    bias <- true_values - m
    biasP <- (bias / true_values) * 100
    erro <- apply(estim, 2, sd, na.rm = TRUE)
    MSE <- apply(estim, 2, var, na.rm = TRUE) + bias^2
    TC <- c(calpha, cphi, ctheta) / R
    
    results <- rbind(m, bias, biasP, erro, MSE, TC)
    rownames(results) <- c("Mean", "Bias", "RB%", "SE", "MSE", "TC")
    colnames(results) <- c("alpha", "phi", "theta")
    
    print(c("Tamanho da Amostra:", n))
    print(round(results, 4))
    print(warnings())
  }
})  

end_time <- Sys.time()
execution_time <- end_time - start_time
print(paste("Tempo total de execução:", round(as.numeric(execution_time, units = "secs")), "segundos"))
