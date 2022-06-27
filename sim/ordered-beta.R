library(tidyverse)

rorderedbeta <- function(n, mu=0.5, phi=1, k0=-1, k1=1) {
    pi.0 <- 1 - plogis(qlogis(mu) - k0)
    pi.1 <- plogis(qlogis(mu) - k1)
    pi.01 <- plogis(qlogis(mu) - k0) - plogis(qlogis(mu) - k1)

    I <- sample.int(3, n, prob=c(pi.0, pi.1, pi.01), replace=TRUE)
    R <- cbind(0, 1, rbeta(n, mu*phi, (1-mu)*phi))
    return(R[cbind(1:n, I)])
}

porderedbeta <- function(x, mu=0.5, phi=1, k0=-1, k1=1) {
    pi.0 <- 1 - plogis(qlogis(mu) - k0)
    pi.1 <- plogis(qlogis(mu) - k1)
    pi.01 <- plogis(qlogis(mu) - k0) - plogis(qlogis(mu) - k1)

    if (x == 0)
        return(pi.0)
    else if (x == 1)
        return(pi.1)
    else
        return(pi.01 * dbeta(x, mu*phi, (1-mu)*phi))
}


MU <- .75
PHI <- 10
K0 <- -5
K1 <- 10.5
PI.0 <- 1 - plogis(qlogis(MU) - K0)
PI.1 <- plogis(qlogis(MU) - K1)
PI.01 <- plogis(qlogis(MU) - K0) - plogis(qlogis(MU) - K1)
X <- rorderedbeta(10000000, mu=MU, phi=PHI, k0=K0, k1=K1)
hist(X, breaks=100)

## Estimates of mean and variance
mean(X) - (PI.1 + PI.01*MU)
var(X) - (PI.1*(1-PI.1) - 2*PI.1*PI.01*MU + PI.01*MU*(1 + MU*PHI - PI.01*MU*PHI - PI.01*MU)/(PHI + 1))
var(X) - (PI.1*(1-PI.1) - 2*PI.1*PI.01*MU + PI.01*MU*(1-MU)/(PHI+1) + PI.01*(1-PI.01)*MU**2)

PI.1*(1-PI.1) - 2*PI.1*PI.01*MU + PI.01*MU*(1 + MU*PHI - PI.01*MU*PHI - PI.01*MU)/(PHI + 1)

PI.1*(1-PI.1) - 2*PI.1*PI.01*MU + PI.01*MU*(1-MU)/(PHI+1) + PI.01*(1-PI.01)*MU*(MU*PHI + 1)/(PHI+1)
