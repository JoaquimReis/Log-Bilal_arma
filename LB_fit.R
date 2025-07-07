# ---------------------------------------------------------------------------
# FUNÇÃO PARA AJUSTE DO MODELO LOG-BILAL ARMA
# ---------------------------------------------------------------------------
LogBarma.fit<-function (y, ar = NA, ma = NA, link = "logit",
                        h=1, diag=0,X = NA,X_hat=NA)
{
  source("LogB-functions.r")
  if (min(y) <= 0 || max(y) >= 1)
    stop("OUT OF RANGE (0,1)!")
  
  if(is.ts(y)==T)  freq<-frequency(y) else stop("data can be a time-series object")
  
  z<-c()
  maxit1<-50
  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  m <- max(p,q,na.rm=T)
  y1 <- y[(m+1):n]
  p1 <- length(ar)  
  q1 <- length(ma)
  error <- rep(0,n) 
  eta <- rep(NA,n)
  
  y_prev <- c(rep(NA,(n+h))) 
  
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("logit", "probit", "cloglog"))){
    stats <- make.link(linktemp)
  }  else {
    stop(paste(linktemp, "link not available, available links are \"logit\", ","\"probit\" and \"cloglog\""))
  }
  
  link = linktemp 
  linkfun = stats$linkfun
  linkinv = stats$linkinv 
  mu.eta = stats$mu.eta 
  diflink = function(t) 1/(stats$mu.eta(stats$linkfun(t)))
  ynew = linkfun(y) 
  ynew_ar <- suppressWarnings(matrix(ynew,(n-1),max(p,1,na.rm=T)))
  
  ###########################################################
  if(any(is.na(ar)) == F) {
    names_phi <- c(paste("phi", ar, sep = ""))
    Z <- suppressWarnings(matrix(ynew, (n-1), p1)[m:(n-1),])} else {
      ar = p1<-0; Z <- NA  
    } 
  
  if(any(is.na(ma)) == F) {
    names_theta <- c(paste("theta", ma, sep = ""))
  } else ma = q1 <- 0 
  
  if(any(is.na(X)) == F){
    names_beta<-c(paste("beta", 1 : ncol(as.matrix(X)), sep = ""))
    Xm <- X[(m+1):n, ]     
    k = ncol(X)
  } else {
    k = 0 
    X <- matrix(rep(0,n), nrow = n)
    Xm <- NA
  }
  
  # Recorrences
  q_1 <- max(q1, 1)
  
  k_i <- q1/q_1                                 
  
  Xstart <- (cbind(rep(1, (n-m)), Xm, Z))         
  Xstart <- matrix(apply(Xstart, 1, na.omit),nrow = (n-m),byrow = T)
  ols <- lm.fit(Xstart, ynew[(m+1) : n])$coef
  initial <- c(rep(0, k+p1+q1+1))
  initial[1 : (k+p1+1)] <- ols
  
  loglik <- function(z) 
  {
    alpha <- z[1]
    if(k==0)  beta = as.matrix(0) else beta = as.matrix(z[2:(k+1)])
    if(p1==0) {phi = as.matrix(0);ar=1} else phi = as.matrix(z[(k+2):(k+p1+1)]) 
    if(q1==0) theta = as.matrix(0) else  theta = as.matrix(z[(k+p1+2):(k+p1+q1+1)])
    
    Xbeta <- X%*%beta
    Xbeta_ar <- suppressWarnings(matrix(Xbeta, (n-1), max(p, 1, na.rm = T)))
    
    for(i in (m+1):n)
    {
      eta[i] <- alpha + Xbeta[i] + (ynew_ar[(i-1), ar] - Xbeta_ar[(i-1), ar])%*%phi + t(theta)%*%error[i-ma]
      error[i] <- ynew[i] - eta[i] 
    }
    mu <- linkinv(eta[(m+1):n])
    
    ll <-log(12 / ((mu / (mu + 24))^(-1/2) - 5)) + 
      (4 / ((mu / (mu + 24))^(-1/2) - 5) - 1) * log(y1) + 
      log(1 - y1^(2 / ((mu / (mu + 24))^(-1/2) - 5)))
    sum(ll)
  } 
  
  #############################################
  
  opt<-optim(initial, loglik, 
             method = "BFGS", hessian = TRUE,
             control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
  
  
  if (opt$conv != 0)
  {
    warning("FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
    opt<-optim(initial, loglik, 
               method = "BFGS", hessian = TRUE,
               control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
    if (opt$conv != 0)
    {
      warning("FUNCTION DID NOT CONVERGE NEITHER WITH NUMERICAL GRADIENT!")
    }else{
      warning("IT WORKS WITH NUMERICAL GRADIENT!")
    }
  }
  
  z$conv <- opt$conv
  coef <- (opt$par)[1:(p1+q1+k+1)]
  alpha <- coef[1]
  if(k==0) beta=names_beta=NULL else z$beta <- coef[2:(k+1)]
  if(p1==0) phi=names_phi=NULL else z$phi <- coef[(k+2):(k+p1+1)]
  if(q1==0) theta=names_theta=NULL else z$theta <- coef[(k+p1+2):(k+p1+q1+1)]
  
  names_par <- c("alpha",names_beta,names_phi,names_theta)
  names(coef)<-names_par
  z$coeff <- coef
  J_inv <- solve(-(opt$hessian))
  z$stderror<-sqrt(diag(J_inv))
  z$zstat <- z$coef/z$stderror
  z$pvalues <- 2*(1 - pnorm(abs(z$zstat)) )
  z$loglik <- opt$value
  z$counts <- as.numeric(opt$counts[1])
  
  if(any(is.na(X)==F))
  { 
    z$k<- (p1+q1+k+2)
    z$aic <- -2*(z$loglik)+2*(z$k)
    z$bic <- -2*(z$loglik)+log(n)*(z$k)
    z$hq <- -2*(z$loglik)+log(log(n))*(z$k)
  }else{
    z$k<- (p1+q1+2)
    z$aic <- -2*(z$loglik)+2*(z$k)
    z$bic <- -2*(z$loglik)+log(n)*(z$k)
    z$hq <- -2*(z$loglik)+log(log(n))*(z$k)
  }
  
  model_presentation <- cbind(round(z$coef,4),round(z$stderror,4),round(z$zstat,4),round(z$pvalues,4))
  colnames(model_presentation)<-c("Estimate","Std. Error","z value","Pr(>|z|)")
  z$model <- model_presentation
  
  return(z)
}


#------------------------------------------------------------------------------
# Exemplo de uso
#------------------------------------------------------------------------------
# Ajuste o caminho para suas funções:
# source("simu.LogBarma.R")

# Exemplo com (p1=1, q1=1) e nenhuma covariável => total de 3 parâmetros:
#set.seed(2)
#y <- simu.LogBarma(100,phi=0.2, theta=0.3, alpha=2)
#fit <- LogBarma.fit(y, ma=1, ar=1,)
# fit$model
