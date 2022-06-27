library(dplyr)
library(tidyr)
library(ggplot2)
library(cmdstanr)
library(modelr)
library(bayestestR)
library(tidybayes)
library(bayesplot)
library(rstantools)
library(viridis)
library(ggthemes)
library(patchwork)
library(ggquiver)
library(latex2exp)
options(mc.cores=parallel::detectCores())

####################################################################################################
##                                     Model Predictions
####################################################################################################
## Create a sequence, but shifting 0/1 to close values
## to avoid NAs in calculations
prob_seq <- function(..., delta=1e-6) {
    c <- seq(...)
    c[c == 0] <- delta
    c[c == 1] <- 1 - delta
    return(c)
}

## linearly normalize x to the range [lower, upper]
normalize <- function(x, lower=0, upper=1) {
    maximum <- max(x, na.rm=TRUE)
    minimum <- min(x, na.rm=TRUE)

    if (maximum == minimum) {
        x
    } else {
        (x - minimum) * ((upper-lower) / (maximum-minimum)) + lower
    }
}

## Derived causal judgments from each model
K_deltaP <- function(p_C, p_A, structure)
    ifelse(structure=='Conjunctive', p_A, 1-p_A)
K_PPC <- function(p_C, p_A, structure)
    ifelse(structure=='Conjunctive', p_A, 1)
K_SP <- function(p_C, p_A, structure)
    ifelse(structure=='Conjunctive', (1-p_C)*p_A, (1-p_C) * (1-p_A))
K_Icard <- function(p_C, p_A, structure)
    ifelse(structure=='Conjunctive', p_C*p_A-p_C+1, p_C)
K_Quillien <- function(p_C, p_A, structure)
    ifelse(structure=='Conjunctive',
           sqrt((1-p_C)*p_A / (1-p_C*p_A)),
           sqrt((1-p_A)*p_C / (p_C + p_A - p_C*p_A)))

df.pred <- expand_grid(structure=c('Conjunctive', 'Disjunctive'),
                       p_C=prob_seq(.1, 1, .1), p_A=prob_seq(.1, 1, .1)) %>%
    mutate(p_E=ifelse(structure=='Conjunctive', p_C*p_A, p_C+p_A-p_C*p_A),
           deltaP=    K_deltaP(p_C, p_A, structure),
           PPC=       K_PPC(p_C, p_A, structure),
           SP=        K_SP(p_C, p_A, structure),
           Icard=     K_Icard(p_C, p_A, structure),
           Quillien=  K_Quillien(p_C, p_A, structure)) %>%
    pivot_longer(deltaP:Quillien, names_to='model', values_to='K') %>%
    group_by(model) %>%
    mutate(model=factor(model, levels=c('deltaP', 'PPC', 'SP', 'Icard', 'Quillien')),
           norm=ifelse(model=='Quillien', sqrt(p_C*(1-p_C)/(p_E*(1-p_E))), 1),
           K.raw=K/norm,
           K.var=norm^2*K.raw*(1-K.raw),
           K.sd=sqrt(K.var),
           K.cv=K.sd/K,
           K.entropy=-(K.raw*log(K.raw) + (1-K.raw)*log(1-K.raw + 1e-12)),
           p_C=ifelse(1 - p_C < 0.09, 1, p_C),
           p_A=ifelse(1 - p_A < 0.09, 1, p_A))

df.pred.wide <- df.pred %>% select(-norm) %>%
    pivot_wider(names_from='model',
                values_from=c('K', 'K.raw', 'K.var', 'K.sd', 'K.entropy', 'K.cv')) %>%
    mutate(structure=tolower(substring(structure, 1, 3)),
           p_C=factor(p_C, ordered=TRUE),
           p_A=factor(p_A, ordered=TRUE))



####################################################################################################
##                                     Load Data
####################################################################################################
judgments <- read.csv('data/data-processed.csv') %>% tibble %>%
    mutate(valence=factor(ifelse(vignette %in% c('casino', 'ducks', 'fruit'), '+', '-'),
                          levels=c('+', '-')),
           p_C=factor(p_C, ordered=TRUE),
           p_A=factor(p_A, ordered=TRUE)) %>%
    left_join(df.pred.wide) %>%
    mutate(structure=factor(structure, labels=c('Conjunctive', 'Disjunctive')),
           vignette=factor(vignette),
           p_C=as.numeric(as.character(p_C)),
           p_A=as.numeric(as.character(p_A)))

N <- length(unique(judgments$id))
writeLines(sprintf('Age: %.2f (%.2f)', mean(judgments$age, na.rm=TRUE), sd(judgments$age, na.rm=TRUE)))
judgments %>% select(sex) %>% table

judgments <- filter(judgments, attn_check=='Yes.') %>% select(-attn_check)
writeLines(sprintf('Age: %.2f (%.2f)', mean(judgments$age, na.rm=TRUE), sd(judgments$age, na.rm=TRUE)))
judgments %>% select(sex) %>% table
N2 <- length(unique(judgments$id))
writeLines(sprintf('N = %d (pre), %d (post), %.2f percent dropout', N, N2, 100*(1 - N2/N)))


## set up data for Stan
grid <- judgments %>% data_grid(p_C, p_A) %>%
    mutate(grid_index=row_number())

judgments <- judgments %>% left_join(grid, by=c('p_C', 'p_A'))

data.stan <- list(prior_only=TRUE,
                  N=nrow(judgments), G=100, Dx=2, Dy=2,
                  C=length(levels(judgments$structure)),
                  V=length(levels(judgments$vignette)),
                  g=judgments$grid_index,
                  c=as.integer(judgments$structure),
                  v=as.integer(judgments$vignette),
                  x=select(grid, -grid_index),
                  y=select(judgments, cause, conf))
str(data.stan)

df.pred.wide <- df.pred.wide %>%
    mutate(structure=factor(structure, labels=c('Conjunctive', 'Disjunctive')),
           p_C=as.numeric(as.character(p_C)),
           p_A=as.numeric(as.character(p_A)))


####################################################################################################
##                                          GP prior
####################################################################################################
model <- cmdstan_model('gp.stan')
if (file.exists('gp-prior.rds')) {
    gp.prior <- readRDS('gp-prior.rds')
} else {
    data.stan$prior_only <- TRUE
    gp.prior <- model$sample(data=data.stan)
    gp.prior$save_object('gp-prior.rds')
}

gp.prior$summary(c('rho', 'rho_tilde', 'alpha', 'alpha_tilde', 'cutpoints', 'Omega', 'Omega_tilde')) %>% print(n=132)
mcmc_dens(gp.prior$draws(c('rho', 'rho_tilde', 'alpha', 'alpha_tilde', 'cutpoints')))

draws.prior <- gp.prior %>%
    spread_draws(e_pred[structure, grid_index, resp], var_pred[structure, grid_index, resp]) %>%
    full_join(gp.prior %>%
              spread_draws(e_pred_grad[structure, dx, grid_index, resp],
                           var_pred_grad[structure, dx, grid_index, resp]) %>%
              mutate(dx=factor(dx, levels=1:2, labels=c('pC', 'pA'))) %>%
              pivot_wider(names_from=dx, values_from=c(e_pred_grad, var_pred_grad))) %>%
    left_join(grid, by='grid_index') %>%
    mutate(structure=factor(levels(judgments$structure)[structure]),
           resp=factor(resp, levels=1:2, labels=c('Causal Judgment', 'Confidence'))) %>%
    group_by(structure, p_C, p_A, resp) %>%
    select(-grid_index) %>%
    rename_with(~ paste0('prior_', .), starts_with(c('e_pred', 'var_pred')))

draws.prior.vignette <- gp.prior %>%
    spread_draws(e_pred_tilde[structure, vignette, grid_index, resp], var_pred_tilde[structure, vignette, grid_index, resp]) %>%
    rename(e_pred=e_pred_tilde, var_pred=var_pred_tilde) %>%
    full_join(gp.prior %>%
              spread_draws(e_pred_tilde_grad[structure, vignette, dx, grid_index, resp],
                           var_pred_tilde_grad[structure, vignette, dx, grid_index, resp]) %>%
              mutate(dx=factor(dx, levels=1:2, labels=c('pC', 'pA'))) %>%
              rename(e_pred_grad=e_pred_tilde_grad, var_pred_grad=var_pred_tilde_grad) %>%
              pivot_wider(names_from=dx, values_from=c(e_pred_grad, var_pred_grad))) %>%
    left_join(grid, by='grid_index') %>%
    mutate(structure=factor(levels(judgments$structure)[structure]),
           vignette=factor(levels(judgments$vignette)[vignette]),
           resp=factor(resp, levels=1:2, labels=c('Causal Judgment', 'Confidence'))) %>%
    group_by(structure, vignette, p_C, p_A, resp) %>%
    select(-grid_index) %>%
    rename_with(~ paste0('prior_', .), starts_with(c('e_pred', 'var_pred')))


####################################################################################################
##                                          GP posterior
####################################################################################################
if (file.exists('gp.rds')) {
    gp.fit <- readRDS('gp.rds')
} else {
    data.stan$prior_only <- FALSE
    gp.fit <- model$sample(data=data.stan)
    gp.fit$save_object('gp.rds')
}

gp.fit$summary(c('rho', 'rho_tilde', 'alpha', 'alpha_tilde', 'cutpoints', 'Omega', 'Omega_tilde')) %>% print(n=132)
mcmc_dens(gp.fit$draws(c('rho', 'rho_tilde', 'alpha', 'alpha_tilde', 'cutpoints')))

draws <- gp.fit %>%
    spread_draws(e_pred[structure, grid_index, resp], var_pred[structure, grid_index, resp]) %>%
    full_join(gp.fit %>%
              spread_draws(e_pred_grad[structure, dx, grid_index, resp],
                           var_pred_grad[structure, dx, grid_index, resp]) %>%
              mutate(dx=factor(dx, levels=1:2, labels=c('pC', 'pA'))) %>%
              pivot_wider(names_from=dx, values_from=c(e_pred_grad, var_pred_grad))) %>%
    left_join(grid, by='grid_index') %>%
    mutate(structure=factor(levels(judgments$structure)[structure]),
           resp=factor(resp, levels=1:2, labels=c('Causal Judgment', 'Confidence'))) %>%
    group_by(structure, p_C, p_A, resp) %>%
    select(-grid_index) %>%
    full_join(draws.prior) %>%
    group_by(structure, p_C, p_A, resp) %>%
    mutate(P_e_pred_grad_pC=as.numeric(pd_to_p(pd(e_pred_grad_pC))),
           P_e_pred_grad_pA=as.numeric(pd_to_p(pd(e_pred_grad_pA))),
           P_var_pred_grad_pC=as.numeric(pd_to_p(pd(var_pred_grad_pC))),
           P_var_pred_grad_pA=as.numeric(pd_to_p(pd(var_pred_grad_pA))),
           log10_BF_e_pred_grad_pC=bf_pointnull(e_pred_grad_pC, prior_e_pred_grad_pC)$log_BF / log(10),
           log10_BF_e_pred_grad_pA=bf_pointnull(e_pred_grad_pA, prior_e_pred_grad_pA)$log_BF / log(10),
           log10_BF_var_pred_grad_pC=bf_pointnull(var_pred_grad_pC, prior_var_pred_grad_pC)$log_BF / log(10),
           log10_BF_var_pred_grad_pA=bf_pointnull(var_pred_grad_pA, prior_var_pred_grad_pA)$log_BF / log(10))

draws.vignette <- gp.fit %>%
    spread_draws(e_pred_tilde[structure, vignette, grid_index, resp], var_pred_tilde[structure, vignette, grid_index, resp]) %>%
    rename(e_pred=e_pred_tilde, var_pred=var_pred_tilde) %>%
    full_join(gp.fit %>%
              spread_draws(e_pred_tilde_grad[structure, vignette, dx, grid_index, resp],
                           var_pred_tilde_grad[structure, vignette, dx, grid_index, resp]) %>%
              mutate(dx=factor(dx, levels=1:2, labels=c('pC', 'pA'))) %>%
              rename(e_pred_grad=e_pred_tilde_grad, var_pred_grad=var_pred_tilde_grad) %>%
              pivot_wider(names_from=dx, values_from=c(e_pred_grad, var_pred_grad))) %>%
    left_join(grid, by='grid_index') %>%
    mutate(structure=factor(levels(judgments$structure)[structure]),
           vignette=factor(levels(judgments$vignette)[vignette]),
           resp=factor(resp, levels=1:2, labels=c('Causal Judgment', 'Confidence'))) %>%
    group_by(structure, p_C, p_A, resp, vignette) %>%
    select(-grid_index) %>%
    full_join(draws.prior.vignette) %>%
    group_by(structure, p_C, p_A, resp, vignette) %>%
    mutate(P_e_pred_grad_pC=as.numeric(pd_to_p(pd(e_pred_grad_pC))),
           P_e_pred_grad_pA=as.numeric(pd_to_p(pd(e_pred_grad_pA))),
           P_var_pred_grad_pC=as.numeric(pd_to_p(pd(var_pred_grad_pC))),
           P_var_pred_grad_pA=as.numeric(pd_to_p(pd(var_pred_grad_pA))),
           log10_BF_e_pred_grad_pC=bf_pointnull(e_pred_grad_pC, prior_e_pred_grad_pC)$log_BF / log(10),
           log10_BF_e_pred_grad_pA=bf_pointnull(e_pred_grad_pA, prior_e_pred_grad_pA)$log_BF / log(10),
           log10_BF_var_pred_grad_pC=bf_pointnull(var_pred_grad_pC, prior_var_pred_grad_pC)$log_BF / log(10),
           log10_BF_var_pred_grad_pA=bf_pointnull(var_pred_grad_pA, prior_var_pred_grad_pA)$log_BF / log(10))



####################################################################################################
##                                 GP Fit Checks
####################################################################################################

## Posterior predictive distribution (per data point)
post_pred <- gp.fit$draws('y_hat', format='draws_df') %>%
    pivot_longer(starts_with('y_hat'), names_pattern='y_hat\\[(.*),(.*)\\]', names_to=c('.row', 'resp')) %>%
    mutate(.row=as.integer(.row),
           resp=factor(resp, levels=c('1','2'), labels=c('cause_hat', 'conf_hat'))) %>%
    pivot_wider(names_from=resp) %>%
    left_join(judgments %>% select(vignette:conf) %>% mutate(.row=row_number()), by='.row') %>%
    mutate(cause_resid=cause-cause_hat,
           conf_resid=conf-conf_hat)
## Expectation of posterior predictive distribution (per data point)
draws.e_pred <- judgments %>% select(structure, p_C, p_A, vignette, cause, conf) %>%
    mutate(.row=row_number()) %>%
    expand_grid(.draw=unique(draws$.draw)) %>%
    left_join(draws.vignette %>%
              select(structure, p_C, p_A, vignette, resp, e_pred, prior_e_pred, .draw) %>%
              pivot_wider(names_from=resp, values_from=e_pred:prior_e_pred) %>%
              rename(cause_hat=`e_pred_Causal Judgment`,
                     prior_cause_hat=`prior_e_pred_Causal Judgment`,
                     conf_hat=e_pred_Confidence,
                     prior_conf_hat=prior_e_pred_Confidence))

draws.e_pred.summary <- draws.e_pred %>%
    group_by(.row, vignette, structure, p_C, p_A) %>%
    median_hdci(cause, cause_hat, prior_cause_hat, conf, conf_hat, prior_conf_hat) %>%
    select(-cause.lower, -cause.upper, -conf.lower, -conf.upper)


## Fit plot
s <- sample.int(max(post_pred$.draw), 100)
p.fit.cause <- ggplot(draws.e_pred.summary, aes(x=cause, y=cause_hat, ymin=cause_hat.lower, ymax=cause_hat.upper)) +
    geom_abline(intercept=0, slope=1) +
    geom_pointinterval(point_size=.5, point_alpha=.33, interval_size=.01, interval_alpha=.05, stroke=0) +
    stat_smooth(method='lm') +
    theme_bw(base_size=14) +
    coord_fixed(xlim=c(0, 1), ylim=c(0, 1)) +
    xlab('Causal Judgment') + ylab('Expectation of\nPredicted Causal Judgment')
p.fit.conf <- ggplot(draws.e_pred.summary, aes(x=conf, y=conf_hat, ymin=conf_hat.lower, ymax=conf_hat.upper)) +
    geom_abline(intercept=0, slope=1) +
    geom_pointinterval(point_size=.5, point_alpha=.33, interval_size=.01, interval_alpha=.05, stroke=0) +
    stat_smooth(method='lm') +
    theme_bw(base_size=14) +
    coord_fixed(xlim=c(0, 1), ylim=c(0, 1)) +
    xlab('Confidence Rating') + ylab('Expectation of\nPredicted Confidence Rating')

p.fit.cause | p.fit.conf
ggsave('plots/fit.png', width=8, height=4, dpi=1000)

## Posterior Predictive Check
p.ppc.cause <- post_pred %>% select(.row, .draw, cause_hat) %>%
    pivot_wider(names_from=.draw, values_from=cause_hat) %>%
    select(-.row) %>% as.matrix() %>% t %>%
    ppc_ecdf_overlay(y=judgments$cause, yrep=.[s,]) +
    xlab('Causal Judgment') + ylab('Cumulative Probability')
p.ppc.conf <- post_pred %>% select(.row, .draw, conf_hat) %>%
    pivot_wider(names_from=.draw, values_from=conf_hat) %>%
    select(-.row) %>% as.matrix() %>% t %>%
    ppc_ecdf_overlay(y=judgments$conf, yrep=.[s,]) +
    xlab('Confidence Rating') + ylab('Cumulative Probability')

(p.ppc.cause | p.ppc.conf) + plot_layout(guides='collect')
ggsave('plots/ppc.png', width=8, height=4)


## Get R^2 values
cause_hat <- draws.e_pred %>% select(.row, .draw, cause_hat) %>%
    pivot_wider(names_from=.draw, values_from=cause_hat) %>%
    select(-.row) %>% as.matrix() %>% t
conf_hat <- draws.e_pred %>% select(.row, .draw, conf_hat) %>%
    pivot_wider(names_from=.draw, values_from=conf_hat) %>%
    select(-.row) %>% as.matrix() %>% t
draws.R2 <- tibble(.draw=as.integer(1:4000),
                   R2_cause=bayes_R2(cause_hat, judgments$cause),
                   R2_conf=bayes_R2(conf_hat, judgments$conf)) %>%
    pivot_longer(R2_cause:R2_conf, names_to='resp', values_to='R Squared') %>%
    mutate(resp=factor(resp, levels=c('R2_conf', 'R2_cause'), labels=c('Confidence', 'Causal Judgment')))

draws.R2 %>% group_by(resp) %>% median_hdi(`R Squared`)

draws.R2 %>%
    ggplot(aes(y=resp, x=`R Squared`)) +
    stat_halfeye() +
    theme_bw(base_size=14) +
    theme(axis.title.y=element_blank())
ggsave('plots/r_squared.png', width=8, height=4)




####################################################################################################
##                                   Plots
####################################################################################################
## Fixed effects
p.cause.mean <- draws %>%
    filter(resp=='Causal Judgment') %>%
    median_hdci(e_pred, e_pred_grad_pC, e_pred_grad_pA, P_e_pred_grad_pC, P_e_pred_grad_pA, log10_BF_e_pred_grad_pC, log10_BF_e_pred_grad_pA) %>%
    mutate(e_pred_grad_pC=ifelse(log10_BF_e_pred_grad_pC > 1 & P_e_pred_grad_pC < .05, e_pred_grad_pC, 0),
           e_pred_grad_pA=ifelse(log10_BF_e_pred_grad_pA > 1 & P_e_pred_grad_pA < .05, e_pred_grad_pA, 0)) %>%
    ggplot(aes(x=p_C, y=p_A)) +
    geom_raster(aes(fill=e_pred)) +
    geom_quiver(aes(u=e_pred_grad_pC, v=e_pred_grad_pA), color='white', size=.25, center=TRUE) +
    scale_fill_viridis(name='Mean\nCausal\nJudgment', option='magma') +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of Focal Cause') + ylab('Probability of Alternate Cause') +
    facet_grid(structure ~ .) +
    theme_classic() + coord_fixed()
p.cause.var <- draws %>%
    filter(resp=='Causal Judgment') %>%
    median_hdci(var_pred, var_pred_grad_pC, var_pred_grad_pA, P_var_pred_grad_pC, P_var_pred_grad_pA, log10_BF_var_pred_grad_pC, log10_BF_var_pred_grad_pA) %>%
    mutate(var_pred_grad_pC=ifelse(log10_BF_var_pred_grad_pC > 1 & P_var_pred_grad_pC < .05, var_pred_grad_pC, 0),
           var_pred_grad_pA=ifelse(log10_BF_var_pred_grad_pA > 1 & P_var_pred_grad_pA < .05, var_pred_grad_pA, 0)) %>%
    ggplot(aes(x=p_C, y=p_A)) +
    geom_raster(aes(fill=var_pred)) +
    geom_quiver(aes(u=var_pred_grad_pC, v=var_pred_grad_pA), color='white', size=.25, center=TRUE) +
    scale_fill_viridis(name='Variance\nCausal\nJudgment', option='magma', limits=c(0, NA)) +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of Focal Cause') + ylab('Probability of Alternate Cause') +
    facet_grid(structure ~ .) +
    theme_classic() + coord_fixed()
p.conf.mean <- draws %>%
    filter(resp=='Confidence') %>%
    median_hdci(e_pred, e_pred_grad_pC, e_pred_grad_pA, P_e_pred_grad_pC, P_e_pred_grad_pA, log10_BF_e_pred_grad_pC, log10_BF_e_pred_grad_pA) %>%
    mutate(e_pred_grad_pC=ifelse(log10_BF_e_pred_grad_pC > 1 & P_e_pred_grad_pC < .05, e_pred_grad_pC, 0),
           e_pred_grad_pA=ifelse(log10_BF_e_pred_grad_pA > 1 & P_e_pred_grad_pA < .05, e_pred_grad_pA, 0)) %>%
    ggplot(aes(x=p_C, y=p_A)) +
    geom_raster(aes(fill=e_pred)) +
    geom_quiver(aes(u=e_pred_grad_pC, v=e_pred_grad_pA), color='white', size=.25, center=TRUE) +
    scale_fill_viridis(name='Mean\nConfidence', option='magma') +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of Focal Cause') + ylab('Probability of Alternate Cause') +
    facet_grid(structure ~ .) +
    theme_classic() + coord_fixed()
p.conf.var <- draws %>%
    filter(resp=='Confidence') %>%
    median_hdci(var_pred, var_pred_grad_pC, var_pred_grad_pA, P_var_pred_grad_pC, P_var_pred_grad_pA, log10_BF_var_pred_grad_pC, log10_BF_var_pred_grad_pA) %>%
    mutate(var_pred_grad_pC=ifelse(log10_BF_var_pred_grad_pC > 1 & P_var_pred_grad_pC < .05, var_pred_grad_pC, 0),
           var_pred_grad_pA=ifelse(log10_BF_var_pred_grad_pA > 1 & P_var_pred_grad_pA < .05, var_pred_grad_pA, 0)) %>%
    ggplot(aes(x=p_C, y=p_A)) +
    geom_raster(aes(fill=var_pred)) +
    geom_quiver(aes(u=var_pred_grad_pC, v=var_pred_grad_pA), color='white', size=.25, center=TRUE) +
    scale_fill_viridis(name='Variance\nConfidence', option='magma', limits=c(0, NA)) +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of Focal Cause') + ylab('Probability of Alternate Cause') +
    facet_grid(structure ~ .) +
    theme_classic() + coord_fixed()

p.cause.mean + p.cause.var + p.conf.mean + p.conf.var
ggsave('plots/gp-mean-var.png', width=8, height=8)

## Vignette-level effects
p.cause.mean.vignette <- draws.vignette %>%
    filter(resp=='Causal Judgment') %>%
    median_hdci(e_pred, e_pred_grad_pC, e_pred_grad_pA, P_e_pred_grad_pC, P_e_pred_grad_pA, log10_BF_e_pred_grad_pC, log10_BF_e_pred_grad_pA) %>%
    mutate(e_pred_grad_pC=ifelse(log10_BF_e_pred_grad_pC > 1 & P_e_pred_grad_pC < .05, e_pred_grad_pC, 0),
           e_pred_grad_pA=ifelse(log10_BF_e_pred_grad_pA > 1 & P_e_pred_grad_pA < .05, e_pred_grad_pA, 0)) %>%
    ggplot(aes(x=p_C, y=p_A)) +
    geom_raster(aes(fill=e_pred)) +
    geom_quiver(aes(u=e_pred_grad_pC, v=e_pred_grad_pA), color='white', size=.25, center=TRUE) +
    scale_fill_viridis(name='Mean\nCausal\nJudgment', option='magma') +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of Focal Cause') + ylab('Probability of Alternate Cause') +
    facet_grid(structure ~ vignette) +
    theme_classic() + coord_fixed()
p.cause.var.vignette <- draws.vignette %>%
    filter(resp=='Causal Judgment') %>%
    median_hdci(var_pred, var_pred_grad_pC, var_pred_grad_pA, P_var_pred_grad_pC, P_var_pred_grad_pA, log10_BF_var_pred_grad_pC, log10_BF_var_pred_grad_pA) %>%
    mutate(var_pred_grad_pC=ifelse(log10_BF_var_pred_grad_pC > 1 & P_var_pred_grad_pC < .05, var_pred_grad_pC, 0),
           var_pred_grad_pA=ifelse(log10_BF_var_pred_grad_pA > 1 & P_var_pred_grad_pA < .05, var_pred_grad_pA, 0)) %>%
    ggplot(aes(x=p_C, y=p_A)) +
    geom_raster(aes(fill=var_pred)) +
    geom_quiver(aes(u=var_pred_grad_pC, v=var_pred_grad_pA), color='white', size=.25, center=TRUE) +
    scale_fill_viridis(name='Variance\nCausal\nJudgment', option='magma', limits=c(0, NA)) +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of Focal Cause') + ylab('Probability of Alternate Cause') +
    facet_grid(structure ~ vignette) +
    theme_classic() + coord_fixed()
p.conf.mean.vignette <- draws.vignette %>%
    filter(resp=='Confidence') %>%
    median_hdci(e_pred, e_pred_grad_pC, e_pred_grad_pA, P_e_pred_grad_pC, P_e_pred_grad_pA, log10_BF_e_pred_grad_pC, log10_BF_e_pred_grad_pA) %>%
    mutate(e_pred_grad_pC=ifelse(log10_BF_e_pred_grad_pC > 1 & P_e_pred_grad_pC < .05, e_pred_grad_pC, 0),
           e_pred_grad_pA=ifelse(log10_BF_e_pred_grad_pA > 1 & P_e_pred_grad_pA < .05, e_pred_grad_pA, 0)) %>%
    ggplot(aes(x=p_C, y=p_A)) +
    geom_raster(aes(fill=e_pred)) +
    geom_quiver(aes(u=e_pred_grad_pC, v=e_pred_grad_pA), color='white', size=.25, center=TRUE) +
    scale_fill_viridis(name='Mean\nConfidence', option='magma') +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of Focal Cause') + ylab('Probability of Alternate Cause') +
    facet_grid(structure ~ vignette) +
    theme_classic() + coord_fixed()
p.conf.var.vignette <- draws.vignette %>%
    filter(resp=='Confidence') %>%
    median_hdci(var_pred, var_pred_grad_pC, var_pred_grad_pA, P_var_pred_grad_pC, P_var_pred_grad_pA, log10_BF_var_pred_grad_pC, log10_BF_var_pred_grad_pA) %>%
    mutate(var_pred_grad_pC=ifelse(log10_BF_var_pred_grad_pC > 1 & P_var_pred_grad_pC < .05, var_pred_grad_pC, 0),
           var_pred_grad_pA=ifelse(log10_BF_var_pred_grad_pA > 1 & P_var_pred_grad_pA < .05, var_pred_grad_pA, 0)) %>%
    ggplot(aes(x=p_C, y=p_A)) +
    geom_raster(aes(fill=var_pred)) +
    geom_quiver(aes(u=var_pred_grad_pC, v=var_pred_grad_pA), color='white', size=.25, center=TRUE) +
    scale_fill_viridis(name='Variance\nConfidence', option='magma', limits=c(0, NA)) +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of Focal Cause') + ylab('Probability of Alternate Cause') +
    facet_grid(structure ~ vignette) +
    theme_classic() + coord_fixed()

p.cause.mean.vignette / p.cause.var.vignette / p.conf.mean.vignette / p.conf.var.vignette
ggsave('plots/gp-mean-var-vignette.png', width=8, height=12)

## One-dimensional plots for uncertainty
draws %>%
    filter(p_A %in% c(.1, .5, 1)) %>%
    median_hdi(e_pred) %>%    
    ggplot(aes(x=p_C, y=e_pred, fill=as.factor(p_A))) +
    geom_line(aes(color=as.factor(p_A)), size=2) +
    geom_ribbon(aes(ymin=.lower, ymax=.upper), alpha=0.1) +
    scale_color_viridis(discrete=TRUE, name='Probability of\nAlternate Cause') +
    scale_fill_viridis(discrete=TRUE, name='Probability of\nAlternate Cause') +
    xlab('Probability of Focal Cause') + ylab('Estimated Marginal Mean') +
    facet_grid(resp ~ structure, scales='free_y') +
    theme_classic()
ggsave('plots/gp-mean.png', width=6, height=4)

draws %>%
    filter(p_A %in% c(.1, .5, 1)) %>%
    median_hdi(var_pred) %>%    
    ggplot(aes(x=p_C, y=var_pred, fill=as.factor(p_A))) +
    geom_line(aes(color=as.factor(p_A)), size=2) +
    geom_ribbon(aes(ymin=.lower, ymax=.upper), alpha=0.1) +
    scale_color_viridis(discrete=TRUE, name='Probability of\nAlternate Cause') +
    scale_fill_viridis(discrete=TRUE, name='Probability of\nAlternate Cause') +
    xlab('Probability of Focal Cause') + ylab('Estimated Marginal Variance') +
    facet_grid(resp ~ structure, scales='free_y') +
    theme_classic()
ggsave('plots/gp-variance.png', width=6, height=4)




####################################################################################################
##                                   Compare Models to GP
####################################################################################################
df.pred <- df.pred %>%
    mutate(model=factor(model, levels=c('deltaP', 'PPC', 'SP', 'Icard', 'Quillien'),
                        labels=c('Delta P', 'Power PC', 'Crediting Causality',
                                 'Necessity-Sufficiency', 'CES'))) %>%
    group_by(model)

p.cause.pred <- ggplot(df.pred, aes(x=p_C, y=p_A, fill=K)) +
    geom_raster() +
    facet_grid(structure ~ model) +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of Focal Cause') + ylab('Probability of Alternate Cause') +
    scale_fill_viridis(name=TeX('$\\kappa_{C\\rightarrow E}$'), option='magma', limits=c(0,1)) +
    theme_classic() + coord_fixed()
(p.cause.mean | p.cause.pred) + plot_layout(widths=c(1/6, 5/6)) + plot_annotation(tag_levels='A')
ggsave('plots/gp-model-cause.pdf', width=12, height=4)

p.conf.pred <- df.pred %>%
    group_by(model) %>%
    mutate(K.var=1-normalize(K.var)) %>%
    ggplot(aes(x=p_C, y=p_A, fill=K.var)) +
    geom_raster() +
    facet_grid(structure ~ model) +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of Focal Cause') + ylab('Probability of Alternate Cause') +
    scale_fill_viridis(name=TeX('$Var(\\kappa_{C\\rightarrow E}$)'), option='magma') +
    theme_classic() + coord_fixed()
(p.conf.mean | p.conf.pred) + plot_layout(widths=c(1/6, 5/6)) + plot_annotation(tag_levels='A')
ggsave('plots/gp-model-var.pdf', width=12, height=4)

p.conf.pred <- df.pred %>%
    group_by(model) %>%
    mutate(K.sd=1-normalize(K.sd)) %>%
    ggplot(aes(x=p_C, y=p_A, fill=K.sd)) +
    geom_raster() +
    facet_grid(structure ~ model) +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of Focal Cause') + ylab('Probability of Alternate Cause') +
    scale_fill_viridis(name=TeX('$\\sigma_{\\kappa_{C\\rightarrow E}$}'), option='magma') +
    theme_classic() + coord_fixed()
(p.conf.mean | p.conf.pred) + plot_layout(widths=c(1/6, 5/6)) + plot_annotation(tag_levels='A')
ggsave('plots/gp-model-sd.pdf', width=12, height=4)

p.conf.pred <- df.pred %>%
    group_by(model) %>%
    mutate(K.cv=1-normalize(K.cv)) %>%
    ggplot(aes(x=p_C, y=p_A, fill=K.cv)) +
    geom_raster() +
    facet_grid(structure ~ model) +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of Focal Cause') + ylab('Probability of Alternate Cause') +
    scale_fill_viridis(name=TeX('$CV(\\kappa_{C\\rightarrow E}$)'), option='magma') +
    theme_classic() + coord_fixed()
(p.conf.mean | p.conf.pred) + plot_layout(widths=c(1/6, 5/6)) + plot_annotation(tag_levels='A')
ggsave('plots/gp-model-cv.pdf', width=12, height=4)

p.conf.pred <- df.pred %>%
    group_by(model) %>%
    mutate(K.entropy=1-normalize(K.entropy)) %>%
    ggplot(aes(x=p_C, y=p_A, fill=K.entropy)) +
    geom_raster() +
    facet_grid(structure ~ model) +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of Focal Cause') + ylab('Probability of Alternate Cause') +
    scale_fill_viridis(name=TeX('$H(\\kappa_{C\\rightarrow E}$)'), option='magma') +
    theme_classic() + coord_fixed()
(p.conf.mean | p.conf.pred) + plot_layout(widths=c(1/6, 5/6)) + plot_annotation(tag_levels='A')
ggsave('plots/gp-model-entropy.pdf', width=12, height=4)




## Obtain correlations between models & data
draws.model <- draws %>% ungroup %>%
    select(structure, p_C, p_A, resp, .chain, .iteration, .draw, e_pred, prior_e_pred) %>%
    mutate(p_C=as.character(p_C), p_A=as.character(p_A)) %>%
    pivot_wider(names_from=resp, values_from=c(e_pred, prior_e_pred)) %>%
    rename(cause=`e_pred_Causal Judgment`,
           confidence=`e_pred_Confidence`,
           prior_cause=`prior_e_pred_Causal Judgment`,
           prior_confidence=`prior_e_pred_Confidence`) %>%
    left_join(df.pred.wide %>% mutate(p_C=as.character(p_C), p_A=as.character(p_A))) %>%
    group_by(.draw, structure) %>%
    summarize(across(K_deltaP:K_Quillien, ~cor(cause, .), .names='posterior_{.col}'),
              across(K.var_deltaP:K.var_Quillien, ~cor(confidence, .), .names='posterior_{.col}'),
              across(K.sd_deltaP:K.sd_Quillien, ~cor(confidence, .), .names='posterior_{.col}'),
              across(K.cv_deltaP:K.cv_Quillien, ~cor(confidence, .), .names='posterior_{.col}'),
              across(K.entropy_deltaP:K.entropy_Quillien, ~cor(confidence, .), .names='posterior_{.col}'),
              across(K_deltaP:K_Quillien, ~cor(prior_cause, .), .names='prior_{.col}'),              
              across(K.var_deltaP:K.var_Quillien, ~cor(prior_confidence, .), .names='prior_{.col}'),
              across(K.sd_deltaP:K.sd_Quillien, ~cor(prior_confidence, .), .names='prior_{.col}'),
              across(K.cv_deltaP:K.cv_Quillien, ~cor(prior_confidence, .), .names='prior_{.col}'),
              across(K.entropy_deltaP:K.entropy_Quillien, ~cor(prior_confidence, .), .names='prior_{.col}')) %>%
    pivot_longer(posterior_K_deltaP:prior_K.entropy_Quillien, names_to=c('distribution', 'variable', 'model'), names_sep='_') %>%
    ## reverse correlations for confidence
    mutate(value=ifelse(variable=='K', value, -value),
           model=factor(model, levels=c('deltaP', 'PPC', 'SP', 'Icard', 'Quillien'),
                        labels=c('Delta P', 'Power PC', 'Crediting\nCausality', 'Necessity\nSufficiency', 'Counterfactual\nEffect Size'))) %>%
    replace_na(list(value=0)) %>%
    pivot_wider(names_from=distribution) %>%
    group_by(variable, structure, model) %>%
    mutate(P=as.numeric(pd_to_p(pd(posterior))),
           log10_BF=ifelse(all(posterior==0), 0, bf_pointnull(posterior, prior)$log_BF / log(10)))

## Print model summaries
draws.model %>%
    filter(variable=='K' | variable=='K.sd') %>%
    median_hdi(posterior, P, log10_BF) %>%
    arrange(variable, structure, desc(posterior)) %>%
    select(-P.lower, -P.upper, -log10_BF.lower, -log10_BF.upper) %>%
    print(width=100)

draws.model %>%
    filter(variable=='K' | variable=='K.sd') %>%
    ggplot(aes(x=posterior, y=model, group=structure, fill=structure)) +
    geom_vline(xintercept=0, linetype='dashed') +
    stat_halfeye(aes(side=ifelse(structure=='Conjunctive', 'top', 'bottom')),
                 position=position_dodge(width=-.5),
                 point_interval='median_hdi',
                 normalize='panels') +
    scale_fill_manual(name='Causal Structure', values=colorblind_pal()(3)[-1]) +
    facet_wrap(~variable, labeller=labeller(variable=c('K'='Causal Judgment', 'K.sd'='Confidence'))) +
    ylab('Model') + xlab('Correlation') +
    theme_classic() +
    theme(legend.position='bottom')
ggsave('plots/model_comparisons.pdf', width=10, height=5)


draws.model %>%
    filter(variable!='K' & variable!='K.sd') %>%
    mutate(variable=factor(variable, levels=c('K.var', 'K.entropy', 'K.cv'))) %>%
    ggplot(aes(x=posterior, y=model, group=structure, fill=structure)) +
    geom_vline(xintercept=0, linetype='dashed') +
    stat_halfeye(aes(side=ifelse(structure=='Conjunctive', 'top', 'bottom')),
                 position=position_dodge(width=-.5),
                 point_interval='median_hdi',
                 normalize='panels') +
    scale_fill_manual(name='Causal Structure', values=colorblind_pal()(3)[-1]) +
    facet_wrap(~variable, labeller=labeller(variable=c('K.var'='Variance', 'K.cv'='Coefficient of Variation',
                                                       'K.entropy'='Entropy'))) +
    ylab('Model') + xlab('Correlation') +
    theme_classic() +
    theme(legend.position='bottom')
ggsave('plots/model_comparisons_confidence.pdf', width=10, height=5)
