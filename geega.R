geega <- function(w, Z, gamma, id, corstr,eps,itermax) {
  gamma1 <- gamma
  id <- id
  K <- length(unique(id))
  n <- as.vector(table(id))
  P <- exp(gamma1 %*% t(Z))/(1 + exp(gamma1 %*% t(Z)))
  PP <- P * (1 - P)
  SK1 <- 1
  repeat {
    if (corstr == "independence") {
      phi <- 1
      alpha <- 0
    }else{ 
    rr <- as.vector((w - P)/sqrt(PP))
    rr1 <- 0
    phi <- sum(rr^2)/(sum(n) - dim(Z)[2])
    rrm <- matrix(0, ncol = K, nrow = max(n))
    for (i in 1:K) {
      rrm[(1:n[i]), i] <- rr[id == i]
    }
    rr <- rrm
    rr <- t(rr)
     if(corstr == "exchangeable") {
      for (i in 1:K) {
        if (n[i] == 1) {
          rr1 <- rr1 + rr[i, 1]
        } else {
          for (j in 1:(n[i] - 1)) rr1 <- rr1 + rr[i, j] * sum(rr[i, (j + 1):n[i]])
        }
      }
      alpha <- (phi^(-1)) * rr1/(sum(n * (n - 1))/2 - dim(Z)[2])
      
    }else if(corstr == "AR1"){
      
      for (i in 1:K) {
        if (n[i] == 1)
          rr1 <- rr1 + rr[i, 1] else {
            for (j in 1:(n[i] - 1)) rr1 <- rr1 + rr[i, j] * rr[i, j+1]
          }
      }
      alpha <- (phi^(-1)) * rr1/(sum((n - 1)) - dim(Z)[2]  )
    }
   }

    SK <- 1
    repeat {
      D <- diag(c(PP[id == 1], 0))[1:n[1], 1:n[1]] %*% Z[id == 1, ]
      for (i in 2:K) {
        D <- rbind(D, diag(c(PP[id == i], 0))[1:n[i], 1:n[i]] %*% Z[id == i, ])
      }
      S <- w - P

      if(corstr=="AR1"){
        exponent <- abs(matrix(1:n[1] - 1, nrow = n[1], ncol = n[1], byrow = TRUE) - (1:n[1] - 1))
        R1 <- alpha^exponent
      }else  if(corstr=="exchangeable" | corstr=="independence"){
        R1 <- matrix(alpha, n[1], n[1])
        diag(R1) <- 1
      }
      V <- sqrt(diag(c(PP[id == 1], 0))[1:n[1], 1:n[1]]) %*% R1 %*% sqrt(diag(c(PP[id == 1], 0))[1:n[1], 1:n[1]]) * phi
      for (i in 2:K) {
        if(corstr=="AR1"){
          exponent <- abs(matrix(1:n[i] - 1, nrow = n[i], ncol = n[i], byrow = TRUE) - (1:n[i] - 1))
          R1 <- alpha^exponent
        }else  if(corstr=="exchangeable" | corstr=="independence"){
          R1 <- matrix(alpha, n[i], n[i])
          diag(R1) <- 1
        }
        V <- bdiag(V, sqrt(diag(c(PP[id == i], 0))[1:n[i], 1:n[i]]) %*% R1 %*% sqrt(diag(c(PP[id == i], 0))[1:n[i], 1:n[i]]) * phi)
      }
      V <- as.matrix(V)
      ZU <- gamma %*% t(D) + S
      newgamma <- t(ginv(t(D) %*% ginv(V) %*% D) %*% t(D) %*% ginv(V) %*% t(ZU))
      if (any(abs(newgamma - gamma) > eps) & SK <= itermax) {
        gamma <- newgamma
        SK <- SK + 1
        P <- exp(gamma %*% t(Z))/(1 + exp(gamma %*% t(Z)))
        PP <- P * (1 - P)
      } else break
    }

    if (any(abs(gamma - gamma1) > eps) & SK1 <= itermax) {
      gamma1 <- gamma
      P <- exp(gamma1 %*% t(Z))/(1 + exp(gamma1 %*% t(Z)))
      PP <- P * (1 - P)
      SK1 <- SK1 + 1
    } else break
  }
  convergence <- (SK1 <= itermax) & (SK <= itermax)
  list(gamma = gamma1, alpha = alpha, phi = phi, convergence =convergence )
}