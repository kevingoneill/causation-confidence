library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(patchwork)
library(infotheo)


## sample n 0s/1s with probability prob
rbernoulli <- function(n, prob) rbinom(n, 1, prob)

## Create a sequence, but shifting 0/1 to close values
## to avoid NAs in calculations
prob_seq <- function(..., delta=1e-6) {
    c <- seq(...)
    c[c == 0] <- delta
    c[c == 1] <- 1 - delta
    return(c)
}

## Normalize a relative to b
softmax <- function(a, b) exp(a) / (exp(a) + exp(b))
softmaxv <- function(a, b) c(softmax(a,b), softmax(b,a))

softmax.gradient <- function(a, b) c(softmax(a,b)*(1-softmax(b,a)),
                                     -softmax(a,b)*softmax(b,a))

## probability weighting
weight.prob <- function(p, gamma=0.6) {
    p^gamma / (p^gamma + (1-p)^gamma)^(1/gamma)
}


normalize <- function(x, lower=0, upper=1) {
    maximum <- max(x, na.rm=TRUE)
    minimum <- min(x, na.rm=TRUE)

    if (maximum == minimum) {
        x
    } else {
        (x - minimum) * ((upper-lower) / (maximum-minimum)) + lower
    }
}


## helper function to plot predicted/simulated values
pred_plot <- function(df, var='K', ylab='Causal Strength', trans='identity',
                      ylim=c(NA, NA), facet=TRUE) {
    p <- ggplot(df) +
        aes(x=p_C, group=p_A, color=p_A) + aes_string(y=var) +
        geom_line(size=1) +
        scale_color_viridis(name='P(A)') +
        scale_x_continuous('P(C)', limits=0:1,
                           breaks=c(0, .5, 1), labels=c('0', '.5', '1')) +
        scale_y_continuous(ylab, trans=trans, limits=ylim) +
        theme_bw()

    if (facet) p <- p + facet_grid(method ~ structure)

    return(p)
}
