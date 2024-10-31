geees <- function(Time, Status, Ibeta,Igamma,X, Z, id, model, corstr, stdz, esmax, eps) {
  Kn <- length(id)
  K <- length(unique(id))
  n <- as.vector(table(id))
  Z1 <- Z
  X1 <- X
  if (stdz) {
    for (i in 2:ncol(Z1)) {
      Z[, i] <- (Z[, i] - mean(Z[, i]))/sd(Z[, i])
    }
    
    
    for (i in 1:ncol(X1)) {
      X[, i] <- (X[, i] - mean(X[, i]))/sd(X[, i])
    }
    
  }
  ppmt2 <- c(Igamma , Ibeta)
  gamma1 <- Igamma
  beta2 <- Ibeta
  w <- Status
  
  Survival <- basesurv(Time, Status, X, beta = beta2, w)
  Lambda <- Survival$bcumhaz
  survivaL <- Survival$uncuresurv
  KK1 <- 0
  repeat {
  
    w <- Status + ((1 - Status) * exp(gamma1 %*% t(Z)) * survivaL)/(1 + exp(gamma1 %*% t(Z)) * survivaL)
    w <- as.vector(w)

    Gamma1 <- eval(parse(text = paste("geese", "(", "w~Z[,-1,drop=FALSE]", ",id=", "id", ",family = binomial", ",corstr='", 
                                      corstr, "'", ")", sep = "") ))
    gamma1<- Gamma1$beta
    
    Beta1 <- geebt(Status, Lambda, X, beta = beta2 , w, id, corstr,eps,esmax)
    beta1 <- Beta1$beta
    
    Survival <- basesurv(Time, Status, X, beta = beta1, w)
    Lambda <- Survival$bcumhaz
    survivaL <- Survival$uncuresurv
    ppmt1 <- c(gamma1, beta1)   

    if (any(abs(ppmt1 - ppmt2) > eps) && (KK1 < esmax)) {
      ppmt2 <- ppmt1
      KK1 <- KK1 + 1
    } else break
  }
  w <- Status + ((1 - Status) * exp(gamma1 %*% t(Z)) * survivaL)/(1 + exp(gamma1 %*% t(Z)) * survivaL)
  w <- as.vector(w)
  if (stdz) {
    gamma1 <- c(gamma1[1] - sum((gamma1[-1] * apply(Z1[, -1, drop = FALSE], 2, mean)/apply(Z1[, -1, drop = FALSE], 2, sd))), gamma1[-1]/apply(Z1[, 
                                                                                                                                                 -1, drop = FALSE], 2, sd))
    beta1 <- beta1/apply(X1, 2, sd)
  }
  
    alpha <- Gamma1$alpha 
    phi  <- Gamma1$gamma
    rho <- Beta1$rho
    pphi <- Beta1$pphi
   
  
  convergence <-  (Beta1$convergence) & (KK1 < esmax)
  es <- list(gamma = gamma1, beta = beta1, Lambda = Lambda, gcor = alpha, gphi = phi, bcor = rho, 
             bphi = pphi, w = w, convergence = convergence)
}