geevarest <- function(Time, Status, X, Z, id, gamma, beta, gphi, gcor, bphi, bcor, Lambda, w) {
    K <- length(unique(id))
    n <- as.vector(table(id))

        t2 <- Time
        c1 <- Status
        t11 <- sort(Time)
        c11 <- Status[order(Time)]
        tt1 <- unique(t11[c11 == 1])
        kk <- length(table(t11[c11 == 1]))
        gamascale <- gphi
        betascale <- bphi
        gamacorr <- gcor
        betacorr <- bcor
        newppmt2c <- as.vector(gamma)
        newppmt2s <- as.vector(beta)
        gg1 <- w
        z1 <- Z
        p1 <- exp(z1 %*% newppmt2c)/(1 + exp(z1 %*% newppmt2c))
        g1 <- gg1
        ABC <- rep(0, K)
        VA <- matrix(0, dim(z1)[2], dim(z1)[2])
        for (v in 1:(dim(z1)[2])) {
            for (w1 in 1:(dim(z1)[2])) {
                for (i in 1:K) {
                  R1 <- matrix(gamacorr, n[i], n[i])
                  diag(R1) <- 1
                  IR1 <- solve(R1)
                  B1 <- matrix(0, n[i], n[i])
                  z11 <- matrix(z1[id == i, ], nrow = n[i], )
                  A1 <- t(z11[, v])
                  p11 <- p1[id == i]
                  g11 <- g1[id == i]
                  pp11 <- p11 * (1 - p11)
                  BB <- (pp11^(1/2)) %*% ((t(pp11))^(-1/2)) * IR1
                  for (s in 1:n[i]) {
                    for (l in 1:n[i]) {
                      B1[s, l] <- (1/2) * (z11[s, w1] * (1 - 2 * p11[s]) - z11[l, w1] * (1 - 2 * p11[l])) * BB[s, l]
                    }
                  }
                  C1 <- g11 - p11
                  D1 <- BB
                  E1 <- z11[, w1] * pp11
                  ABC[i] <- A1 %*% (B1 %*% C1 - D1 %*% E1)
                }
                VA[v, w1] <- sum(ABC) * (gamascale^(-1))
                ABC <- rep(0, K)
            }
        }
        sdm_gamma <- VA
        be <- as.vector(newppmt2s)
        z2 <- X
        c2 <- Status
        mu2 <- exp(z2 %*% be)
        Lambda0 <- Lambda
        ABC1 <- rep(0, K)
        VA1 <- matrix(0, dim(z2)[2], dim(z2)[2])
        for (v in 1:dim(z2)[2]) {
            for (w1 in 1:dim(z2)[2]) {
                for (i in 1:K) {
                  Q1 <- matrix(betacorr, n[i], n[i])
                  diag(Q1) <- 1
                  IQ1 <- solve(Q1)
                  B2 <- matrix(0, n[i], n[i])
                  z22 <- matrix(z2[id == i, ], nrow = n[i], )
                  A2 <- t(z22[, v])
                  c22 <- c2[id == i]
                  g22 <- g1[id == i]
                  mu22 <- mu2[id == i]
                  Lambda22 <- Lambda0[id == i]
                  BB1 <- (mu22^(1/2)) %*% ((t(mu22))^(-1/2)) * IQ1
                  for (s in 1:n[i]) {
                    for (l in 1:n[i]) {
                      B2[s, l] <- (1/2) * (z22[s, w1] - z22[l, w1]) * BB1[s, l]
                    }
                  }
                  C2 <- c22/Lambda22 - mu22
                  D2 <- BB1
                  E2 <- z22[, w1] * mu22
                  G2 <- diag(g22 * Lambda22)
                  ABC1[i] <- A2 %*% (B2 %*% G2 %*% C2 - D2 %*% G2 %*% E2)
                }
                VA1[v, w1] <- sum(ABC1) * (betascale^(-1))
                ABC1 <- rep(0, K)
            }
        }
        gSS <- sort(unique(Lambda))
        if(gSS[1] == 1e-08)
        { gSS <- gSS[-1] }
        gS <- c(gSS[1], gSS[2:kk] - gSS[1:(kk - 1)])
        BBC <- matrix(0, kk, dim(z2)[2])
        xxxx <- z2
        gg1 <- g1
        for (j in 1:dim(z2)[2]) {
            for (s in 1:(kk)) {
                BCm <- gS[s] * exp(xxxx[(c1 == 1) & (t2 == tt1[s]), ,drop = FALSE] %*% be)
                BBC[s, j] <- sum(exp(xxxx[(c1 == 1) & (t2 == tt1[s]), ,drop = FALSE] %*% be) * (exp(-BCm) + BCm * exp(-BCm) - 1)/((1 - exp(-BCm))^2) * xxxx[(c1 == 1) & (t2 == tt1[s]), j]) + sum(gg1[t2 >= 
                  tt1[s]] * exp(xxxx[t2 >= tt1[s], ,drop = FALSE] %*% be) * xxxx[t2 >= tt1[s], j])
            }
        }
        CCC <- rep(0, (kk))
        for (s in 1:(kk)) {
            CCm <- gS[s] * exp(xxxx[(c1 == 1) & (t2 == tt1[s]), ,drop = FALSE] %*% be)
            CCC[s] <- -sum(exp(2 * (xxxx[(c1 == 1) & (t2 == tt1[s]), ,drop = FALSE] %*% be) - CCm)/(1 - exp(-CCm))^2)
        }
        BC <- matrix(0, dim(z2)[2], kk)
        for (r in 1:dim(z2)[2]) {
            for (s in 1:(kk)) {
                elem <- 0
                for (i in 1:K) {
                  mu22 <- mu2[id == i]
                  xxx1 <- xxxx[id == i, r]
                  t21 <- t2[id == i]
                  g22 <- g1[id == i]
                  Q1 <- matrix(betacorr, n[i], n[i])
                  diag(Q1) <- 1
                  IQ1 <- solve(Q1)
                  for (j in 1:n[i]) {
                    if (t21[j] >= tt1[s]) 
                    elem <- elem + sum(xxx1 * ((mu22)^(1/2)) * ((mu22[j])^(-1/2)) * IQ1[, j]) * g22[j] * mu22[j] * (betascale^(-1))
                  }
                }
                BC[r, s] <- -elem
            }
        }
        sdm_betalpha <- rbind(cbind(VA1, BC), cbind(-BBC, diag(CCC)))
        sdm <- bdiag(sdm_gamma, sdm_betalpha)
        sdm <- as.matrix(sdm)
        PP <- p1 * (1 - p1)
        zzzz <- t(z1)
        ept <- g1 * (1 - g1)
        R1 <- matrix(gamacorr, n[1], n[1])
        diag(R1) <- 1
        R1m <- R1
        for (i in 2:K) {
            R1a <- matrix(gamacorr, n[i], n[i])
            diag(R1a) <- 1
            R1m <- bdiag(R1m, R1a)
        }
        R1m <- as.matrix(R1m)
        eptm <- (sqrt(ept) %*% t(sqrt(ept))) * R1m
        D <- diag(PP[id == 1]) %*% diag(rep(1, n[1])) %*% t(zzzz[, id == 1])
        for (i in 2:K) {
            D <- rbind(D, diag(PP[id == i], n[i], n[i]) %*% diag(1, n[i], n[i]) %*% t(zzzz[, id == i]))
        }
        V <- sqrt(diag(PP[id == 1])) %*% R1 %*% sqrt(diag(PP[id == 1])) * gamascale
        for (i in 2:K) {
            R1 <- matrix(gamacorr, n[i], n[i])
            diag(R1) <- 1
            V <- bdiag(V, sqrt(diag(PP[id == i], n[i], n[i])) %*% R1 %*% sqrt(diag(PP[id == i], n[i], n[i])) * gamascale)
        }
        V <- as.matrix(V)
        LV11 <- t(D) %*% ginv(V) %*% eptm %*% t(ginv(V)) %*% D
        mu <- mu2
        xxxx <- z2
        Y1 <- Status
        D1 <- diag(mu[id == 1]) %*% diag(1, n[1], n[1]) %*% t(t(xxxx[id == 1, ]))
        for (i in 2:K) {
            D1 <- rbind(D1, diag(mu[id == i], n[i], n[i]) %*% diag(1, n[i], n[i]) %*% matrix(xxxx[id == i, ], nrow = n[i]))
        }
        S1 <- Y1 - Lambda0 * mu
        eptm1 <- (S1 %*% t(S1)) * eptm
        Q1 <- matrix(betacorr, n[1], n[1])
        diag(Q1) <- 1
        V1 <- sqrt(diag(mu[id == 1])) %*% Q1 %*% sqrt(diag(mu[id == 1])) * betascale
        for (i in 2:K) {
            Q1 <- matrix(betacorr, n[i], n[i])
            diag(Q1) <- 1
            V1 <- bdiag(V1, sqrt(diag(mu[id == i], n[i], n[i])) %*% Q1 %*% sqrt(diag(mu[id == i], n[i], n[i])) * betascale)
        }
        V1 <- as.matrix(V1)
        LV22 <- t(D1) %*% ginv(V1) %*% eptm1 %*% t(ginv(V1)) %*% D1
        cmu22 <- c2[id == 1] - (mu2[id == 1]) * (Lambda0[id == 1])
        R3m <- t(matrix(rep(cmu22, n[1]), n[1], ))
        for (i in 2:K) {
            cmu22 <- c2[id == i] - (mu2[id == i]) * (Lambda0[id == i])
            R3m <- bdiag(R3m, t(matrix(rep(cmu22, n[i]), n[i], )))
        }
        R3m <- as.matrix(R3m)
        eptm2 <- eptm * R3m
        LV12 <- t(D) %*% ginv(V) %*% eptm2 %*% t(ginv(V1)) %*% D1
        cmu22 <- c2[id == 1] - (mu2[id == 1]) * (Lambda0[id == 1])
        R3mr <- (matrix(rep(cmu22, n[1]), n[1], ))
        for (i in 2:K) {
            cmu22 <- c2[id == i] - (mu2[id == i]) * (Lambda0[id == i])
            R3mr <- bdiag(R3mr, (matrix(rep(cmu22, n[i]), n[i], )))
        }
        R3mr <- as.matrix(R3mr)
        eptm3 <- eptm * (R3mr)
        LV21 <- t(D1) %*% ginv(V1) %*% eptm3 %*% t(ginv(V)) %*% D
        LV33 <- matrix(0, kk, kk)
        for (s in 1:kk) {
            for (l in 1:kk) {
                DDD <- 0
                for (i in 1:K) {
                  t21 <- t2[id == i]
                  A33 <- (1:n[i])[t21 >= tt1[s]]
                  B33 <- (1:n[i])[t21 >= tt1[l]]
                  if (length(A33) == 0 || length(B33) == 0) {
                    DDD <- DDD
                  } else {
                    R1 <- matrix(gamacorr, n[i], n[i])
                    diag(R1) <- 1
                    ept1 <- ept[id == i]
                    ept1m <- sqrt(ept1) %*% t(sqrt(ept1)) * R1
                    mu22 <- mu2[id == i]
                    mu22m <- mu22 %*% t(mu22)
                    DDD <- DDD + sum((ept1m * mu22m)[A33, B33])
                  }
                }
                LV33[s, l] <- DDD
            }
        }
        LV13 <- matrix(0, dim(z1)[2], kk)
        for (r in 1:dim(z1)[2]) {
            for (s in 1:kk) {
                EEE <- 0
                for (i in 1:K) {
                  t21 <- t2[id == i]
                  A13 <- (1:n[i])[t21 >= tt1[s]]
                  if (length(A13) == 0) {
                    EEE <- EEE
                  } else {
                    R1 <- matrix(gamacorr, n[i], n[i])
                    diag(R1) <- 1
                    IR1 <- solve(R1)
                    zz13 <- z1[id == i, r]
                    mpp <- PP[id == i]
                    mm1 <- 0
                    for (j in 1:n[i]) {
                      mm1[j] <- sum(zz13 * sqrt(mpp) * IR1[, j])/sqrt(mpp[j])
                    }
                    mm2 <- mu[id == i]
                    ept1 <- ept[id == i]
                    ept1m <- sqrt(ept1) %*% t(sqrt(ept1)) * R1
                    Newi <- mm1 %*% t(-(mm2)) * ept1m
                    EEE <- EEE + sum(Newi[, A13])
                  }
                }
                LV13[r, s] <- EEE
            }
        }
        LV23 <- matrix(0, dim(z2)[2], kk)
        for (r in 1:dim(z2)[2]) {
            for (s in 1:kk) {
                FFF <- 0
                for (i in 1:K) {
                  t21 <- t2[id == i]
                  A23 <- (1:n[i])[t21 >= tt1[s]]
                  if (length(A23) == 0) {
                    FFF <- FFF
                  } else {
                    R1 <- matrix(gamacorr, n[i], n[i])
                    diag(R1) <- 1
                    R2 <- matrix(betacorr, n[i], n[i])
                    diag(R2) <- 1
                    IR2 <- solve(R2)
                    c21 <- c2[id == i]
                    Lambda01 <- Lambda0[id == i]
                    xx23 <- z2[id == i, r]
                    mmu <- mu[id == i]
                    nn1 <- 0
                    for (j in 1:n[i]) {
                      nn1[j] <- sum(xx23 * sqrt(mmu) * IR2[, j])/sqrt(mmu[j]) * (c21[j] - Lambda01[j] * mmu[j])
                    }
                    ept1 <- ept[id == i]
                    ept1m <- sqrt(ept1) %*% t(sqrt(ept1)) * R1
                    New1i <- nn1 %*% t(-mmu) * ept1m
                    FFF <- FFF + sum(New1i[, A23])
                  }
                }
                LV23[r, s] <- FFF
            }
        } 
        LV <- matrix(0, dim(z1)[2] + dim(z2)[2] + kk, dim(z1)[2] + dim(z2)[2] + kk)
        LV[1:dim(z1)[2], 1:dim(z1)[2]] <- LV11
        LV[(dim(z1)[2] + 1):(dim(z1)[2] + dim(z2)[2]), (dim(z1)[2] + 1):(dim(z1)[2] + dim(z2)[2])] <- LV22
        LV[1:dim(z1)[2], (dim(z1)[2] + 1):(dim(z1)[2] + dim(z2)[2])] <- LV12
        LV[(dim(z1)[2] + 1):(dim(z1)[2] + dim(z2)[2]), 1:dim(z1)[2]] <- LV21
        LV[1:dim(z1)[2], (dim(z1)[2] + dim(z2)[2] + 1):(dim(z1)[2] + dim(z2)[2] + kk)] <- LV13
        LV[(dim(z1)[2] + 1):(dim(z1)[2] + dim(z2)[2]), (dim(z1)[2] + dim(z2)[2] + 1):(dim(z1)[2] + dim(z2)[2] + kk)] <- LV23
        LV[(dim(z1)[2] + dim(z2)[2] + 1):(dim(z1)[2] + dim(z2)[2] + kk), 1:dim(z1)[2]] <- t(LV13)
        LV[(dim(z1)[2] + dim(z2)[2] + 1):(dim(z1)[2] + dim(z2)[2] + kk), (dim(z1)[2] + 1):(dim(z1)[2] + dim(z2)[2])] <- t(LV23)
        LV[(dim(z1)[2] + dim(z2)[2] + 1):(dim(z1)[2] + dim(z2)[2] + kk), (dim(z1)[2] + dim(z2)[2] + 1):(dim(z1)[2] + dim(z2)[2] + kk)] <- LV33 
        lsdm <- -sdm - LV
        fdv <- rep(0, dim(z1)[2] + dim(z2)[2] + kk)
        fdm <- matrix(0, dim(z1)[2] + dim(z2)[2] + kk, dim(z1)[2] + dim(z2)[2] + kk)
        for (i in 1:K) {
            z11 <- z1[id == i, ]
            p11 <- p1[id == i]
            g11 <- g1[id == i]
            pp11 <- p11 * (1 - p11)
            pp11m <- diag(pp11)
            C1 <- g11 - p11
            R1 <- matrix(gamacorr, n[i], n[i])
            diag(R1) <- 1
            fdv[1:dim(z1)[2]] <- t(pp11m %*% z11) %*% ginv(sqrt(pp11m) %*% R1 %*% sqrt(pp11m) * (gamascale)) %*% C1
            z22 <- z2[id == i, ]
            mu22 <- mu2[id == i]
            mu22m <- diag(mu22, n[i], n[i])
            Lambda01 <- Lambda0[id == i]
            G2 <- diag(g11 * Lambda01)
            c22 <- c2[id == i]
            C2 <- (c22/Lambda01) - mu22
            Q1 <- matrix(betacorr, n[i], n[i])
            diag(Q1) <- 1
            fdv[(dim(z1)[2] + 1):(dim(z1)[2] + dim(z2)[2])] <- t(mu22m %*% z22) %*% ginv(sqrt(mu22m) %*% Q1 %*% sqrt(mu22m) * (betascale)) %*% G2 %*% C2
            t22 <- t2[id == i]
            eqalpha <- rep(0, kk)
            for (j in 1:kk) {
                Aalpha <- (1:n[i])[t22 == tt1[j]]
                Balpha <- (1:n[i])[t22 >= tt1[j]]
                if (length(Balpha) == 0) {
                  eqalpha[j] <- 0
                }
                if (length(Balpha) != 0 & length(Aalpha) == 0) {
                  eqalpha[j] <- -sum(g11[Balpha] * mu22[Balpha])
                } else eqalpha[j] <- sum(mu22[Aalpha]/(1 - exp(-gS[j] * mu22[Aalpha]))) - sum(g11[Balpha] * mu22[Balpha])
            }
            fdv[(dim(z1)[2] + dim(z2)[2] + 1):(dim(z1)[2] + dim(z2)[2] + kk)] <- eqalpha
            fdm1 <- fdv %*% t(fdv)
            fdm <- fdm + fdm1
        }
        vcm <- ginv(lsdm) %*% fdm %*% t(ginv(lsdm))
        var_gamma <- diag(vcm)[1:dim(z1)[2]]
        var_beta <- diag(vcm)[(dim(z1)[2] + 1):(dim(z1)[2] + dim(z2)[2])]
   
        sd_gamma <- sqrt(var_gamma)
        sd_beta <- sqrt(var_beta)

    return(list(varga = var_gamma, varbe = var_beta))
}