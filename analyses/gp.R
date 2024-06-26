library(tidyverse)
library(cmdstanr)
library(bayestestR)
library(tidybayes)
library(bayesplot)
library(ggdist)
library(viridis)
library(ggthemes)
library(patchwork)
library(ggquiver)
library(latex2exp)
library(multidplyr)
library(loo)
library(posterior)

####################################################################################################
##                                     Load Data
####################################################################################################
judgments <- read_csv('data/data-processed.csv') %>%
    mutate(valence=factor(ifelse(vignette %in% c('casino', 'ducks', 'fruit'), '+', '-'),
                          levels=c('+', '-')),
           structure=factor(structure, labels=c('Conjunctive', 'Disjunctive')),
           vignette=factor(vignette))

N <- length(unique(judgments$id))
writeLines(sprintf('Age: %.2f (%.2f)', mean(judgments$age, na.rm=TRUE), sd(judgments$age, na.rm=TRUE)))
judgments %>% select(sex) %>% table

judgments <- filter(judgments, attn_check=='Yes.') %>% select(-attn_check)
writeLines(sprintf('Age: %.2f (%.2f)', mean(judgments$age, na.rm=TRUE), sd(judgments$age, na.rm=TRUE)))
judgments %>% select(sex) %>% table
N2 <- length(unique(judgments$id))
writeLines(sprintf('N = %d (pre), %d (post), %.2f percent dropout', N, N2, 100*(1 - N2/N)))


## set up data for Stan
grid <- judgments %>% distinct(p_C, p_A) %>%
    arrange(p_C, p_A) %>%
    mutate(grid_index=row_number())
judgments <- judgments %>% left_join(grid, by=c('p_C', 'p_A'))

data.stan <- list(N=nrow(judgments), G=100, Dx=2, Dy=2,
                  C=length(levels(judgments$structure)),
                  V=length(levels(judgments$vignette)),
                  g=judgments$grid_index,
                  c=as.integer(judgments$structure),
                  v=as.integer(judgments$vignette),
                  x=select(grid, -grid_index),
                  y=select(judgments, cause, conf))
str(data.stan)


####################################################################################################
##                                          GP prior
####################################################################################################
model <- cmdstan_model('gp.stan')
if (file.exists('gp-prior.rds')) {
    gp.prior <- readRDS('gp-prior.rds')
} else {
    gp.prior <- model$sample(data=c(data.stan, prior_only=TRUE), parallel_chains=4)
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
    gp.fit <- model$sample(data=c(data.stan, prior_only=FALSE), parallel_chains=4)
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
    mutate(P_e_pred_grad_pC=as.numeric(pd_to_p(as.numeric(pd(e_pred_grad_pC)))),
           P_e_pred_grad_pA=as.numeric(pd_to_p(as.numeric(pd(e_pred_grad_pA)))),
           P_var_pred_grad_pC=as.numeric(pd_to_p(as.numeric(pd(var_pred_grad_pC)))),
           P_var_pred_grad_pA=as.numeric(pd_to_p(as.numeric(pd(var_pred_grad_pA)))),
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


## reduce memory overhead
rm(draws.prior, draws.prior.vignette)
gc()

## Print out significance tests for gradients
draws %>%
    select(.draw, contains('pred_grad')) %>%
    select(-starts_with('prior')) %>%
    pivot_longer(contains('pred_grad'), names_to=c('statistic', 'variable', 'gradient'), names_pattern='(P|log10_BF)*_*(e|var)_pred_grad_(p.)') %>%
    mutate(statistic=ifelse(statistic=='', 'value', statistic)) %>%
    pivot_wider(names_from=statistic) %>%
    mutate(BF = 10^log10_BF) %>%
    select(-log10_BF) %>%
    ##filter(P < .05, BF > 10) %>%   ## comment out this line to get non-significant trends
    group_by(structure, resp, p_C, p_A, variable, gradient, P, BF) %>%
    median_hdi(value) %>%
    group_by(structure, resp, variable, gradient, direction=sign(value)) %>%
    filter(abs(value) == max(abs(value))) %>%
    arrange(resp, structure, variable, gradient, direction) %>%
    knitr::kable('simple')




####################################################################################################
##                                 GP Fit Checks
####################################################################################################

## Calculate variance explained with & without vignette-level effects
draws.R2 <- judgments %>%
    mutate(.row=row_number()) %>%
    select(structure, p_C, p_A, cause, conf) %>%
    expand_grid(.draw=unique(draws$.draw)) %>%
    left_join(draws %>%
              select(structure, p_C, p_A, resp, e_pred, var_pred, .draw) %>%
              pivot_wider(names_from=resp, values_from=e_pred:var_pred) %>%
              rename(e_pred_cause=`e_pred_Causal Judgment`,
                     var_pred_cause=`var_pred_Causal Judgment`,
                     e_pred_confidence=`e_pred_Confidence`,
                     var_pred_confidence=`var_pred_Confidence`)) %>%
    group_by(.draw) %>%
    summarize(R2_cause=var(e_pred_cause)/(var(e_pred_cause) + mean(var_pred_cause)),
              R2_confidence=var(e_pred_confidence)/(var(e_pred_confidence) + mean(var_pred_confidence))) %>%
    pivot_longer(R2_cause:R2_confidence, names_to='.variable', names_prefix='R2_', values_to='R2')
draws.R2.vignette <- judgments %>%
    mutate(.row=row_number()) %>%
    select(structure, p_C, p_A, vignette, cause, conf) %>%
    expand_grid(.draw=unique(draws.vignette$.draw)) %>%
    left_join(draws.vignette %>%
              select(structure, p_C, p_A, vignette, resp, e_pred, var_pred, .draw) %>%
              pivot_wider(names_from=resp, values_from=e_pred:var_pred) %>%
              rename(e_pred_cause=`e_pred_Causal Judgment`,
                     var_pred_cause=`var_pred_Causal Judgment`,
                     e_pred_confidence=`e_pred_Confidence`,
                     var_pred_confidence=`var_pred_Confidence`)) %>%
    group_by(.draw) %>%
    summarize(R2_cause=var(e_pred_cause)/(var(e_pred_cause) + mean(var_pred_cause)),
              R2_confidence=var(e_pred_confidence)/(var(e_pred_confidence) + mean(var_pred_confidence))) %>%
    pivot_longer(R2_cause:R2_confidence, names_to='.variable', names_prefix='R2_', values_to='R2')


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
draws.e_pred.summary

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
draws.R2 %>% group_by(.variable) %>% median_hdi(R2)

draws.R2 %>%
    ggplot(aes(y=.variable, x=R2)) +
    stat_halfeye() +
    scale_x_continuous(name='R Squared', lim=c(0, NA)) +
    scale_y_discrete(labels=c('Causal\nJudgment', 'Confidence')) +
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
    geom_quiver(aes(u=e_pred_grad_pC, v=e_pred_grad_pA), color='white', linewidth=.25, center=TRUE) +
    scale_fill_viridis(name='Mean\nCausal\nJudgment', option='magma') +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of Candidate Cause') + ylab('Probability of Alternate Cause') +
    facet_grid(structure ~ .) +
    theme_classic() + coord_fixed(expand=FALSE)
p.cause.var <- draws %>%
    filter(resp=='Causal Judgment') %>%
    median_hdci(var_pred, var_pred_grad_pC, var_pred_grad_pA, P_var_pred_grad_pC, P_var_pred_grad_pA, log10_BF_var_pred_grad_pC, log10_BF_var_pred_grad_pA) %>%
    mutate(var_pred_grad_pC=ifelse(log10_BF_var_pred_grad_pC > 1 & P_var_pred_grad_pC < .05, var_pred_grad_pC, 0),
           var_pred_grad_pA=ifelse(log10_BF_var_pred_grad_pA > 1 & P_var_pred_grad_pA < .05, var_pred_grad_pA, 0)) %>%
    ggplot(aes(x=p_C, y=p_A)) +
    geom_raster(aes(fill=var_pred)) +
    geom_quiver(aes(u=var_pred_grad_pC, v=var_pred_grad_pA), color='white', linewidth=.25, center=TRUE) +
    scale_fill_viridis(name='Variance\nCausal\nJudgment', option='magma', limits=c(0, NA)) +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of Candidate Cause') + ylab('Probability of Alternate Cause') +
    facet_grid(structure ~ .) +
    theme_classic() + coord_fixed(expand=FALSE)
p.conf.mean <- draws %>%
    filter(resp=='Confidence') %>%
    median_hdci(e_pred, e_pred_grad_pC, e_pred_grad_pA, P_e_pred_grad_pC, P_e_pred_grad_pA, log10_BF_e_pred_grad_pC, log10_BF_e_pred_grad_pA) %>%
    mutate(e_pred_grad_pC=ifelse(log10_BF_e_pred_grad_pC > 1 & P_e_pred_grad_pC < .05, e_pred_grad_pC, 0),
           e_pred_grad_pA=ifelse(log10_BF_e_pred_grad_pA > 1 & P_e_pred_grad_pA < .05, e_pred_grad_pA, 0)) %>%
    ggplot(aes(x=p_C, y=p_A)) +
    geom_raster(aes(fill=e_pred)) +
    geom_quiver(aes(u=e_pred_grad_pC, v=e_pred_grad_pA), color='white', linewidth=.25, center=TRUE) +
    scale_fill_viridis(name='Mean\nConfidence', option='magma') +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of Candidate Cause') + ylab('Probability of Alternate Cause') +
    facet_grid(structure ~ .) +
    theme_classic() + coord_fixed(expand=FALSE)
p.conf.var <- draws %>%
    filter(resp=='Confidence') %>%
    median_hdci(var_pred, var_pred_grad_pC, var_pred_grad_pA, P_var_pred_grad_pC, P_var_pred_grad_pA, log10_BF_var_pred_grad_pC, log10_BF_var_pred_grad_pA) %>%
    mutate(var_pred_grad_pC=ifelse(log10_BF_var_pred_grad_pC > 1 & P_var_pred_grad_pC < .05, var_pred_grad_pC, 0),
           var_pred_grad_pA=ifelse(log10_BF_var_pred_grad_pA > 1 & P_var_pred_grad_pA < .05, var_pred_grad_pA, 0)) %>%
    ggplot(aes(x=p_C, y=p_A)) +
    geom_raster(aes(fill=var_pred)) +
    geom_quiver(aes(u=var_pred_grad_pC, v=var_pred_grad_pA), color='white', linewidth=.25, center=TRUE) +
    scale_fill_viridis(name='Variance\nConfidence', option='magma', limits=c(0, NA)) +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of Candidate Cause') + ylab('Probability of Alternate Cause') +
    facet_grid(structure ~ .) +
    theme_classic() + coord_fixed(expand=FALSE)

p.cause.mean + p.cause.var + p.conf.mean + p.conf.var
ggsave('plots/gp-mean-var.png', width=8, height=8)
ggsave('plots/gp-mean-var.pdf', width=8, height=8)

(p.cause.var | p.conf.var) + plot_annotation(tag_levels='A')
ggsave('plots/gp-var.pdf', width=8, height=4)

## Vignette-level effects
p.cause.mean.vignette <- draws.vignette %>%
    filter(resp=='Causal Judgment') %>%
    median_hdci(e_pred, e_pred_grad_pC, e_pred_grad_pA, P_e_pred_grad_pC, P_e_pred_grad_pA, log10_BF_e_pred_grad_pC, log10_BF_e_pred_grad_pA) %>%
    mutate(e_pred_grad_pC=ifelse(log10_BF_e_pred_grad_pC > 1 & P_e_pred_grad_pC < .05, e_pred_grad_pC, 0),
           e_pred_grad_pA=ifelse(log10_BF_e_pred_grad_pA > 1 & P_e_pred_grad_pA < .05, e_pred_grad_pA, 0)) %>%
    ggplot(aes(x=p_C, y=p_A)) +
    geom_raster(aes(fill=e_pred)) +
    geom_quiver(aes(u=e_pred_grad_pC, v=e_pred_grad_pA), color='white', linewidth=.25, center=TRUE) +
    scale_fill_viridis(name='Mean\nCausal\nJudgment', option='magma') +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of Candidate Cause') + ylab('Probability of Alternate Cause') +
    facet_grid(structure ~ vignette) +
    theme_classic() + coord_fixed(expand=FALSE)
p.cause.var.vignette <- draws.vignette %>%
    filter(resp=='Causal Judgment') %>%
    median_hdci(var_pred, var_pred_grad_pC, var_pred_grad_pA, P_var_pred_grad_pC, P_var_pred_grad_pA, log10_BF_var_pred_grad_pC, log10_BF_var_pred_grad_pA) %>%
    mutate(var_pred_grad_pC=ifelse(log10_BF_var_pred_grad_pC > 1 & P_var_pred_grad_pC < .05, var_pred_grad_pC, 0),
           var_pred_grad_pA=ifelse(log10_BF_var_pred_grad_pA > 1 & P_var_pred_grad_pA < .05, var_pred_grad_pA, 0)) %>%
    ggplot(aes(x=p_C, y=p_A)) +
    geom_raster(aes(fill=var_pred)) +
    geom_quiver(aes(u=var_pred_grad_pC, v=var_pred_grad_pA), color='white', linewidth=.25, center=TRUE) +
    scale_fill_viridis(name='Variance\nCausal\nJudgment', option='magma', limits=c(0, NA)) +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of Candidate Cause') + ylab('Probability of Alternate Cause') +
    facet_grid(structure ~ vignette) +
    theme_classic() + coord_fixed(expand=FALSE)
p.conf.mean.vignette <- draws.vignette %>%
    filter(resp=='Confidence') %>%
    median_hdci(e_pred, e_pred_grad_pC, e_pred_grad_pA, P_e_pred_grad_pC, P_e_pred_grad_pA, log10_BF_e_pred_grad_pC, log10_BF_e_pred_grad_pA) %>%
    mutate(e_pred_grad_pC=ifelse(log10_BF_e_pred_grad_pC > 1 & P_e_pred_grad_pC < .05, e_pred_grad_pC, 0),
           e_pred_grad_pA=ifelse(log10_BF_e_pred_grad_pA > 1 & P_e_pred_grad_pA < .05, e_pred_grad_pA, 0)) %>%
    ggplot(aes(x=p_C, y=p_A)) +
    geom_raster(aes(fill=e_pred)) +
    geom_quiver(aes(u=e_pred_grad_pC, v=e_pred_grad_pA), color='white', linewidth=.25, center=TRUE) +
    scale_fill_viridis(name='Mean\nConfidence', option='magma') +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of Candidate Cause') + ylab('Probability of Alternate Cause') +
    facet_grid(structure ~ vignette) +
    theme_classic() + coord_fixed(expand=FALSE)
p.conf.var.vignette <- draws.vignette %>%
    filter(resp=='Confidence') %>%
    median_hdci(var_pred, var_pred_grad_pC, var_pred_grad_pA, P_var_pred_grad_pC, P_var_pred_grad_pA, log10_BF_var_pred_grad_pC, log10_BF_var_pred_grad_pA) %>%
    mutate(var_pred_grad_pC=ifelse(log10_BF_var_pred_grad_pC > 1 & P_var_pred_grad_pC < .05, var_pred_grad_pC, 0),
           var_pred_grad_pA=ifelse(log10_BF_var_pred_grad_pA > 1 & P_var_pred_grad_pA < .05, var_pred_grad_pA, 0)) %>%
    ggplot(aes(x=p_C, y=p_A)) +
    geom_raster(aes(fill=var_pred)) +
    geom_quiver(aes(u=var_pred_grad_pC, v=var_pred_grad_pA), color='white', linewidth=.25, center=TRUE) +
    scale_fill_viridis(name='Variance\nConfidence', option='magma', limits=c(0, NA)) +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of Candidate Cause') + ylab('Probability of Alternate Cause') +
    facet_grid(structure ~ vignette) +
    theme_classic() + coord_fixed(expand=FALSE)

p.cause.mean.vignette / p.cause.var.vignette / p.conf.mean.vignette / p.conf.var.vignette
ggsave('plots/gp-mean-var-vignette.pdf', width=8, height=12)

## One-dimensional plots for uncertainty
draws %>%
    filter(p_A %in% c(.1, .5, 1)) %>%
    median_hdi(e_pred) %>%    
    ggplot(aes(x=p_C, y=e_pred, fill=as.factor(p_A))) +
    geom_line(aes(color=as.factor(p_A)), linewidth=2) +
    geom_ribbon(aes(ymin=.lower, ymax=.upper), alpha=0.1) +
    scale_color_viridis(discrete=TRUE, name='Probability of\nAlternate Cause') +
    scale_fill_viridis(discrete=TRUE, name='Probability of\nAlternate Cause') +
    xlab('Probability of Candidate Cause') + ylab('Estimated Marginal Mean') +
    facet_grid(resp ~ structure, scales='free_y') +
    theme_classic()
ggsave('plots/gp-mean.png', width=6, height=4)

draws %>%
    filter(p_A %in% c(.1, .5, 1)) %>%
    median_hdi(var_pred) %>%    
    ggplot(aes(x=p_C, y=var_pred, fill=as.factor(p_A))) +
    geom_line(aes(color=as.factor(p_A)), linewidth=2) +
    geom_ribbon(aes(ymin=.lower, ymax=.upper), alpha=0.1) +
    scale_color_viridis(discrete=TRUE, name='Probability of\nAlternate Cause') +
    scale_fill_viridis(discrete=TRUE, name='Probability of\nAlternate Cause') +
    xlab('Probability of Candidate Cause') + ylab('Estimated Marginal Variance') +
    facet_grid(resp ~ structure, scales='free_y') +
    theme_classic()
ggsave('plots/gp-variance.png', width=6, height=4)



draws.cor <- draws %>% select(.draw, ends_with('pred')) %>%
    pivot_wider(names_from=resp, values_from=ends_with('pred')) %>%
    rename(epred_cause=`e_pred_Causal Judgment`,
           epred_conf=`e_pred_Confidence`,
           varpred_cause=`var_pred_Causal Judgment`,
           varpred_conf=`var_pred_Confidence`,
           prior_epred_cause=`prior_e_pred_Causal Judgment`,
           prior_epred_conf=`prior_e_pred_Confidence`,
           prior_varpred_cause=`prior_var_pred_Causal Judgment`,
           prior_varpred_conf=`prior_var_pred_Confidence`) %>%
    group_by(.draw) %>%
    summarize(posterior_r_epred_cause_varpred_cause=cor(epred_cause, varpred_cause),
              posterior_r_epred_conf_varpred_conf=cor(epred_conf, varpred_conf),
              posterior_r_epred_cause_epred_conf=cor(epred_cause, epred_conf),
              posterior_r_varpred_cause_epred_conf=cor(varpred_cause, epred_conf),
              posterior_r_epred_cause_varpred_conf=cor(epred_cause, varpred_conf),
              posterior_r_varpred_cause_varpred_conf=cor(varpred_cause, varpred_conf),
              prior_r_epred_cause_varpred_cause=cor(prior_epred_cause, prior_varpred_cause),
              prior_r_epred_conf_varpred_conf=cor(prior_epred_conf, prior_varpred_conf),
              prior_r_epred_cause_epred_conf=cor(prior_epred_cause, prior_epred_conf),
              prior_r_varpred_cause_epred_conf=cor(prior_varpred_cause, prior_epred_conf),
              prior_r_epred_cause_varpred_conf=cor(prior_epred_cause, prior_varpred_conf),
              prior_r_varpred_cause_varpred_conf=cor(prior_varpred_cause, prior_varpred_conf)) %>%
    mutate(across(posterior_r_epred_cause_varpred_cause:prior_r_varpred_cause_varpred_conf, ~ replace_na(.x, 0)),  ## replace NA correlations with 0
           P_r_epred_cause_varpred_cause=as.numeric(pd_to_p(pd(posterior_r_epred_cause_varpred_cause))),
           P_r_epred_conf_varpred_conf=as.numeric(pd_to_p(pd(posterior_r_epred_conf_varpred_conf))),
           P_r_epred_cause_epred_conf=as.numeric(pd_to_p(pd(posterior_r_epred_cause_epred_conf))),
           P_r_varpred_cause_epred_conf=as.numeric(pd_to_p(pd(posterior_r_varpred_cause_epred_conf))),
           P_r_epred_cause_varpred_conf=as.numeric(pd_to_p(pd(posterior_r_epred_cause_varpred_conf))),
           P_r_varpred_cause_varpred_conf=as.numeric(pd_to_p(pd(posterior_r_varpred_cause_varpred_conf))),
           log10_BF_r_epred_cause_varpred_cause=bf_pointnull(posterior_r_epred_cause_varpred_cause, prior_r_epred_cause_varpred_cause)$log_BF / log(10),
           log10_BF_r_epred_conf_varpred_conf=bf_pointnull(posterior_r_epred_conf_varpred_conf, prior_r_epred_conf_varpred_conf)$log_BF / log(10),
           log10_BF_r_epred_cause_epred_conf=bf_pointnull(posterior_r_epred_cause_epred_conf, prior_r_epred_cause_epred_conf)$log_BF / log(10),
           log10_BF_r_varpred_cause_epred_conf=bf_pointnull(posterior_r_varpred_cause_epred_conf, prior_r_varpred_cause_epred_conf)$log_BF / log(10),
           log10_BF_r_epred_cause_varpred_conf=bf_pointnull(posterior_r_epred_cause_varpred_conf, prior_r_epred_cause_varpred_conf)$log_BF / log(10),
           log10_BF_r_varpred_cause_varpred_conf=bf_pointnull(posterior_r_varpred_cause_varpred_conf, prior_r_varpred_cause_varpred_conf)$log_BF / log(10)) %>%
    pivot_longer(-.draw, names_to=c('dist', 'var1', 'var2'), names_pattern='(.*)_r_(.*_.*)_(.*_.*)') %>%
    pivot_wider(names_from='dist')


draws.cor %>%
    group_by(var1, var2, P, BF=10^log10_BF) %>%
    median_hdci %>%
    arrange(desc(BF)) %>%
    print(width=200)


draws.cor %>%
    mutate(var1=str_remove(var1, '(epred|pred)_'),
           var2=str_remove(var2, '(epred|pred)_'),
           var=factor(paste0('r(', var1, ', ', var2, ')'),
                      levels=rev(c('r(cause, varcause)', 'r(conf, varconf)', 'r(cause, conf)',
                                   'r(varcause, conf)', 'r(varcause, varconf)', 'r(cause, varconf)')))) %>%
    ggplot(aes(y=var, x=posterior)) + #, fill=structure)) +
    geom_vline(xintercept=0, linetype='dashed') +
    stat_slab(aes(x=prior), ##, side=ifelse(structure=='Conjunctive', 'top', 'bottom')),
              position=position_dodge(width=.5), scale=.5, fill='black', alpha=.05,
              normalize='groups') +    
    stat_halfeye(#aes(side=ifelse(structure=='Conjunctive', 'top', 'bottom')),
        fill=colorblind_pal()(2)[-1],
        position=position_dodge(width=.5), scale=.5,
                 point_interval='median_hdi', normalize='groups') +
    #scale_fill_manual(name='Causal Structure', values=colorblind_pal()(2)[-1]) +
    theme_classic() +
    theme(axis.title=element_blank())
ggsave('plots/correlations.png', width=6, height=6)



####################################################################################################
##                              Fit Counterfactual Sampling Models
####################################################################################################
CAUSE_MEASURES <- c('Delta P', 'Power PC', 'Crediting Causality', 'Necessity-Sufficiency', 'CES')
CONFIDENCE_MEASURES <- c('Variance', 'Standard Deviation', 'Coefficient of Variation', 'Entropy')

## helper function to abbreviate the names of the different models
abbrv <- function(s) {
    str_to_lower(str_remove_all(s, '[a-z[:space:]\\-]'))
}

data.models.stan <- list(N=nrow(judgments), P=10,
                         S=length(levels(judgments$structure)),
                         pc=as.integer(as.factor(judgments$p_C)),
                         pa=as.integer(as.factor(judgments$p_A)),
                         s=as.integer(judgments$structure),
                         cause=judgments$cause,
                         confidence=judgments$conf)
str(data.models.stan)

cf_model <- cmdstan_model('counterfactual_sampling_model.stan')


if (length(list.files('models/')) > 0) {
    cf_models <- expand_grid(strength=factor(CAUSE_MEASURES, levels=CAUSE_MEASURES),
                             precision=factor(CONFIDENCE_MEASURES, levels=CONFIDENCE_MEASURES)) %>%
        group_by(strength, precision) %>%
        mutate(fname=paste0(abbrv(strength), '_', abbrv(precision), '.rds'),
               prior=map(fname, function(s) readRDS(paste0('models/prior_', s))),
               posterior=map(fname, function(s) readRDS(paste0('models/', s))))
} else {
    ## Set up a cluster to fit models in parallel
    cluster <- new_cluster(parallel::detectCores())
    cluster_library(cluster, c('tidyverse', 'cmdstanr'))
    cluster_copy(cluster, c('cf_model', 'data.models.stan', 'abbrv'))
    
    cf_models <- expand_grid(strength=factor(CAUSE_MEASURES, levels=CAUSE_MEASURES),
                             precision=factor(CONFIDENCE_MEASURES, levels=CONFIDENCE_MEASURES)) %>%
        group_by(strength, precision) %>%
        partition(cluster) %>%
        mutate(fname=paste0(abbrv(strength), '_', abbrv(precision), '.rds'),
               prior=map2(strength, precision,
                          function(s, p) cf_model$sample(c(data.models.stan, prior_only=TRUE, cause_measure=as.integer(s),
                                                           confidence_measure=as.integer(p)), adapt_delta=.9)),
               posterior=map2(strength, precision,
                              function(s, p) cf_model$sample(c(data.models.stan, prior_only=FALSE, cause_measure=as.integer(s),
                                                               confidence_measure=as.integer(p)), adapt_delta=.9))) %>%
        collect()
    
    ## save model objects
    dir.create('models')
    pwalk(list(cf_models$fname, cf_models$prior), function(s, m) m$save_object(paste0('models/prior_', s)))
    pwalk(list(cf_models$fname, cf_models$posterior), function(s, m) m$save_object(paste0('models/', s)))
}

## check for convergence issues
lapply(cf_models$posterior, function(m) print(m$cmdstan_diagnose()))

## extract model draws
cf_models <- cf_models %>%
    mutate(probs.posterior=map(posterior, ~ spread_draws(., pC[structure, p], pA[structure, p]) %>%
                                              mutate(structure=factor(levels(judgments$structure)[structure]),
                                                     p=p/10)),
           draws.posterior=map(posterior, ~ spread_draws(., epred_cause[structure, p_C, p_A],
                                                         epred_confidence[structure, p_C, p_A],
                                                         cause_hat[structure, p_C, p_A],
                                                         confidence_hat[structure, p_C, p_A]) %>%
                                              mutate(structure=factor(levels(judgments$structure)[structure]),
                                                     p_C=p_C/10,
                                                     p_A=p_A/10)),
           loo_cause=map(posterior, ~ loo(.$draws('log_lik_cause'), r_eff=relative_eff(exp(.$draws('log_lik_cause'))))),
           loo_confidence=map(posterior, ~ loo(.$draws('log_lik_confidence'), r_eff=relative_eff(exp(.$draws('log_lik_confidence'))))),
           loo=map(posterior, ~ loo(.$draws('log_lik_cause') + .$draws('log_lik_confidence'),
                                    r_eff=relative_eff(exp(.$draws('log_lik_cause') + .$draws('log_lik_confidence'))))),
           R2=map(posterior, ~ spread_draws(., R2_cause, R2_confidence)))


cf_models %>%
    unnest(R2) %>%
    pivot_longer(R2_cause:R2_confidence, names_to='.variable', names_prefix='R2_', values_to='R2') %>%
    left_join(draws.R2 %>% rename(descriptive_R2=R2), by=c('.draw', '.variable')) %>%
    mutate(R2_ratio = R2/descriptive_R2,
           .variable=ifelse(.variable=='cause', 'Causal Judgment', 'Confidence')) %>%
    group_by(strength, precision, .variable) %>%
    median_hdci(R2_ratio) %>%
    print(n=40)

cf_models %>%
    unnest(R2) %>%
    pivot_longer(R2_cause:R2_confidence, names_to='.variable', names_prefix='R2_', values_to='R2') %>%
    left_join(draws.R2 %>% rename(descriptive_R2=R2), by=c('.draw', '.variable')) %>%
    mutate(R2_ratio = R2/descriptive_R2,
           .variable=ifelse(.variable=='cause', 'Causal Judgment', 'Confidence')) %>%
    ggplot(aes(y=strength, x=R2_ratio, group=precision, color=precision)) +
    geom_vline(xintercept=0:1, linetype='dashed') +
    stat_pointinterval(point_interval=median_hdi, .width=.95, position=position_dodge(.5)) +
    scale_y_discrete(name='Measure of Difference-Making',
                     labels=c('Delta P', 'Power PC', 'Crediting\nCausality', 'Necessity\nSufficiency', 'Counterfactual\nEffect Size')) +
    scale_color_discrete(name='Measure of\nPrecision') +
    scale_x_continuous(name='R Squared Ratio\n(Relative to Descriptive Model Excluding Vignette-Level Effects)') +
    ##coord_cartesian(xlim=0:1) +
    facet_wrap(.variable ~ ., scales='free_x') +
    theme_classic(base_size=18) +
    theme(panel.grid.major=element_line(linewidth=.5))
ggsave('plots/csm_r2.pdf', width=12, height=6)


## compute ELPD_LOO for descriptive model
## (since pareto k values are too high, these estimates are unreliable)
ll_cause <- loo(gp.fit$draws('log_lik')[,,1:nrow(judgments)],
                r_eff=relative_eff(exp(gp.fit$draws('log_lik')[,,1:nrow(judgments)])))
ll_confidence <- loo(gp.fit$draws('log_lik')[,,(nrow(judgments)+1):(2*nrow(judgments))],
                     r_eff=relative_eff(exp(gp.fit$draws('log_lik')[,,(nrow(judgments)+1):(2*nrow(judgments))])))

rbind(do.call(loo_compare, cf_models$loo_cause),
      do.call(loo_compare, cf_models$loo_confidence)) %>%
    as_tibble(rownames=NA) %>%
    rownames_to_column() %>%
    mutate(.variable=rep(c('Causal Judgment', 'Confidence'), each=nrow(cf_models)),
           strength=cf_models$strength[as.numeric(substring(rowname, 6))],
           precision=cf_models$precision[as.numeric(substring(rowname, 6))]) %>%
    relocate(strength, precision, .variable) %>%
    select(-rowname) %>%
    write_csv('counterfactual_sampling_model_ll_diff.csv')

read_csv('counterfactual_sampling_model_ll_diff.csv') %>%
    mutate(elpd=sprintf('%.2f (%.2f)', elpd_loo, se_elpd_loo),
           diff=sprintf('%.2f (%.2f)', elpd_diff, se_diff)) %>%
    select(.variable, strength, precision, elpd, diff) %>%
    write_delim('counterfactual_sampling_model_ll_diff.txt', delim='&', eol=' \\\\\n')

## Model comparisons using LOO-IC
do.call(loo_compare, cf_models$loo) %>%
    as_tibble(rownames=NA) %>%
    rownames_to_column() %>%
    mutate(strength=cf_models$strength[as.numeric(substring(rowname, 6))],
           precision=cf_models$precision[as.numeric(substring(rowname, 6))]) %>%
    relocate(strength, precision) %>%
    select(-rowname) %>%
    ggplot(aes(y=strength, x=elpd_diff, group=precision, color=precision)) +
    geom_vline(xintercept=0, linetype='dashed') +
    geom_pointrange(aes(xmin=elpd_diff - 2*se_diff, xmax=elpd_diff + 2*se_diff), position=position_dodge(.5)) +
    scale_y_discrete(name='Measure of Difference-Making',
                     labels=c('Delta P', 'Power PC', 'Crediting\nCausality', 'Necessity\nSufficiency', 'Counterfactual\nEffect Size')) +
    scale_color_discrete(name='Measure of\nPrecision') + 
    xlab('Expected Log Pointwise Predictive Density Difference\n(relative to best performing model)') +
    theme_classic(base_size=18)
ggsave('plots/csm_loo.pdf', width=10, height=6)

rbind(do.call(loo_compare, cf_models$loo_cause),
      do.call(loo_compare, cf_models$loo_confidence)) %>%
    as_tibble(rownames=NA) %>%
    rownames_to_column() %>%
    mutate(.variable=rep(c('Causal Judgment', 'Confidence'), each=nrow(cf_models)),
           strength=cf_models$strength[as.numeric(substring(rowname, 6))],
           precision=cf_models$precision[as.numeric(substring(rowname, 6))]) %>%
    relocate(strength, precision, .variable) %>%
    select(-rowname) %>%
    ggplot(aes(y=strength, x=elpd_diff, group=precision, color=precision)) +
    geom_vline(xintercept=0, linetype='dashed') +
    geom_pointrange(aes(xmin=elpd_diff - 2*se_diff, xmax=elpd_diff + 2*se_diff), position=position_dodge(.5)) +
    scale_y_discrete(name='Measure of Difference-Making',
                     labels=c('Delta P', 'Power PC', 'Crediting\nCausality', 'Necessity\nSufficiency', 'Counterfactual\nEffect Size')) +
    scale_color_discrete(name='Measure of\nPrecision') +
    xlab('Expected Log Pointwise Predictive Density Difference\n(relative to best performing model)') +
    facet_wrap(.variable ~ ., scales='free_x') +
    theme_classic(base_size=18)
ggsave('plots/csm_loo2.pdf', width=12, height=6)


## comparisons separately by precision measure
cf_models %>%
    expand_grid(.variable=c('Causal Judgment', 'Confidence')) %>%
    group_by(precision, .variable) %>%
    summarize(loo=ifelse(first(.variable) =='Causal Judgment',
                         list(rownames_to_column(as_tibble(loo_compare(loo_cause), rownames=NA))),
                         list(rownames_to_column(as_tibble(loo_compare(loo_confidence), rownames=NA))))) %>%
    unnest(c(loo)) %>%
    mutate(strength=factor(CAUSE_MEASURES[as.numeric(substring(rowname, 6))], levels=CAUSE_MEASURES)) %>%
    ggplot(aes(y=strength, x=elpd_diff, group=precision, color=precision)) +
    geom_vline(xintercept=0, linetype='dashed') +
    geom_pointrange(aes(xmin=elpd_diff - 2*se_diff, xmax=elpd_diff + 2*se_diff), position=position_dodge(.5)) +
    #scale_y_discrete(name='Measure of Difference-Making',
    #                 labels=c('Delta P', 'Power PC', 'Crediting\nCausality', 'Necessity\nSufficiency', 'Counterfactual\nEffect Size')) +
    scale_color_discrete(name='Measure of\nPrecision') +
    xlab('Expected Log Pointwise Predictive Density Difference\n(relative to best performing model)') +
    facet_wrap(.variable ~ ., scales='free_x') +
    theme_classic(base_size=18)


## plot sampling probability calibration
for (p in CONFIDENCE_MEASURES) {
    print(p)
    p.plot <- cf_models %>%
        filter(precision == p) %>%
        unnest(probs.posterior) %>%
        pivot_longer(c(pC, pA), names_to='.variable', values_to='.value') %>%
        group_by(strength, precision, structure, p, .variable) %>%
        median_hdci(.value) %>%
        ggplot(aes(x=p, y=.value, group=.variable)) +
        geom_abline(slope=1, intercept=0, linetype='dashed') +
        geom_ribbon(aes(ymin=.lower, ymax=.upper, fill=.variable), alpha=.15) +
        geom_line(aes(color=.variable), linewidth=1.5) +
        facet_grid(structure ~ strength) +
        scale_x_continuous(name='Objective Probability',
                           breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
        scale_y_continuous(name='Estimated Counterfactual Sampling Probability',
                           breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
        scale_color_discrete(name='', labels=c('Candidate\nCause', 'Alternate\nCause')) +
        scale_fill_discrete(name='', labels=c('Candidate\nCause', 'Alternate\nCause')) +
        coord_fixed(xlim=c(0, 1), ylim=c(0, 1), expand=FALSE) +
        theme_classic() +
        theme(panel.grid.major=element_line(color='black', linewidth=.05),
              panel.background = element_rect(color='black', linewidth=1))
    ggsave(paste0('plots/csm_probabilities_', abbrv(p), '.pdf'), plot=p.plot, width=12, height=4)
}


## calculate axis limits for plotting
range.cause <- draws %>% filter(resp=='Causal Judgment') %>% median_hdi(e_pred) %>% pull(e_pred) %>% range
range.confidence <- draws %>% filter(resp=='Confidence') %>% median_hdi(e_pred) %>% pull(e_pred) %>% range

## plot predicted causal judgments and confidence
for (p in CONFIDENCE_MEASURES) {
    print(p)
    p.draws <- cf_models %>%
        filter(precision == p) %>%
        unnest(draws.posterior) %>%        
        group_by(strength, precision, structure, p_C, p_A) %>%
        median_hdci(epred_cause, epred_confidence)
    
    csm.cause.plot <- ggplot(p.draws, aes(x=p_C, y=p_A, fill=epred_cause)) +
        geom_raster() +
        scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
        scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
        xlab('Probability of Candidate Cause') + ylab('Probability of Alternate Cause') +
        facet_grid(structure ~ strength) +
        theme_classic() + coord_fixed(expand=FALSE)
    csm.confidence.plot <- ggplot(p.draws, aes(x=p_C, y=p_A, fill=epred_confidence)) +
        geom_raster() +
        scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
        scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
        xlab('Probability of Candidate Cause') + ylab('Probability of Alternate Cause') +
        facet_grid(structure ~ strength) +
        theme_classic() + coord_fixed(expand=FALSE)
    
    (((p.cause.mean | csm.cause.plot) &
      scale_fill_viridis(name='Mean\nCausal\nJudgment', option='magma',
                         limits=c(min(range.cause[1], p.draws$epred_cause),
                                  max(range.cause[2], p.draws$epred_cause)))) +
     plot_layout(widths=c(1/6, 5/6), guides='collect')) /
        (((p.conf.mean | csm.confidence.plot) &
          scale_fill_viridis(name='Mean\nConfidence', option='magma',
                             limits=c(min(range.confidence[1], p.draws$epred_confidence),
                                      max(range.confidence[2], p.draws$epred_confidence)))) +
         plot_layout(widths=c(1/6, 5/6), guides='collect')) +
        plot_annotation(tag_levels='A')
    
    ggsave(paste0('plots/csm_fit_', abbrv(p), '.pdf'), width=12, height=8)
    ggsave(paste0('plots/csm_fit_', abbrv(p), '.png'), width=12, height=8)
}


## plot predicted causal judgments
for (p in CONFIDENCE_MEASURES) {
    print(p)
    p.draws <- cf_models %>%
        filter(precision == p) %>%
        unnest(draws.posterior) %>%        
        group_by(strength, precision, structure, p_C, p_A) %>%
        median_hdci(epred_cause)
    p.plot <- ggplot(p.draws, aes(x=p_C, y=p_A, fill=epred_cause)) +
        geom_raster() +
        scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
        scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
        xlab('Probability of Candidate Cause') + ylab('Probability of Alternate Cause') +
        facet_grid(structure ~ strength) +
        theme_classic() + coord_fixed(expand=FALSE)
    ggsave(paste0('plots/csm_cause_', abbrv(p), '.pdf'), width=12, height=4,
           plot=((p.cause.mean | p.plot) &
                 scale_fill_viridis(name='Mean\nCausal\nJudgment', option='magma',
                                    limits=c(min(range.cause[1], p.draws$epred_cause),
                                             max(range.cause[2], p.draws$epred_cause)))) +
               plot_layout(widths=c(1/6, 5/6), guides='collect') + plot_annotation(tag_levels='A'))
}

## plot predicted confidence
for (p in CONFIDENCE_MEASURES) {
    print(p)
    p.draws <- cf_models %>%
        filter(precision == p) %>%
        unnest(draws.posterior) %>%        
        group_by(strength, precision, structure, p_C, p_A) %>%
        median_hdci(epred_confidence)
    p.plot <- ggplot(p.draws, aes(x=p_C, y=p_A, fill=epred_confidence)) +
        geom_raster() +
        scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
        scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
        xlab('Probability of Candidate Cause') + ylab('Probability of Alternate Cause') +
        facet_grid(structure ~ strength) +
        theme_classic() + coord_fixed(expand=FALSE)
    ggsave(paste0('plots/csm_confidence_', abbrv(p), '.pdf'), width=12, height=4,
           plot=((p.conf.mean | p.plot) &
                 scale_fill_viridis(name='Mean\nConfidence', option='magma',
                                    limits=c(min(range.confidence[1], p.draws$epred_confidence),
                                             max(range.confidence[2], p.draws$epred_confidence)))) +
               plot_layout(widths=c(1/6, 5/6), guides='collect') + plot_annotation(tag_levels='A'))
}


p.data <- ggplot(judgments, aes(x=cause, y=conf)) +
    geom_point(shape=16, alpha=.025, size=.25, data=judgments) +
    scale_x_continuous(name='Causal Judgment', breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(name='Confidence', breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    coord_fixed(xlim=c(0, 1), ylim=c(0, 1)) +
    facet_grid(structure ~ .) +
    theme_classic()

p.data <- draws.quad %>%
    group_by(cause) %>%
    median_hdci(e_pred) %>%
    ggplot(aes(x=cause, y=e_pred)) +
    geom_point(aes(y=conf), shape=16, alpha=.025, size=.25, data=judgments) +
    geom_ribbon(aes(ymin=.lower, ymax=.upper), alpha=.25, fill='red') +
    geom_line(linewidth=1, color='red') +
    scale_x_continuous(name='Causal Judgment', breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(name='Confidence', breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    coord_fixed(xlim=c(0, 1), ylim=c(0, 1)) +
    facet_grid(structure ~ .) +
    theme_classic()

for (p in CONFIDENCE_MEASURES) {
    print(p)
    p.plot <- cf_models %>%
        filter(precision == p) %>%
        unnest(draws.posterior) %>%
        ggplot(aes(x=cause_hat, y=confidence_hat)) +
        geom_point(shape=16, size=.25, alpha=.05, data=function (d) filter(d, .draw %in% sample(1:4000, 50))) +
        ##stat_smooth(color='red', fill='red') +
        ##geom_line(aes(x=epred_cause, y=epred_confidence), linewidth=.75, color='red',
        ##          data=function (d) median_qi(group_by(d, strength, precision, structure, p_C, p_A), epred_cause, epred_confidence)) +
        scale_x_continuous(name='Simulated Causal Judgment', breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
        scale_y_continuous(name='Simulated Confidence', breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
        facet_grid(structure ~ strength) +
        coord_fixed(xlim=0:1, ylim=0:1) +
        theme_classic()
    ggsave(paste0('plots/csm_ppred_', abbrv(p), '.pdf'), width=12, height=4,
           plot=(p.data | p.plot) + plot_layout(widths=c(1/6, 5/6), guides='collect') + plot_annotation(tag_levels='A'))
}


####################################################################################################
##               Test for Quadratic Relationship b/t Causal Judgments & Confidence
####################################################################################################
data.quad <- list(N=nrow(judgments), K=3, X=model.matrix(conf ~ poly(cause, 2, raw=TRUE), judgments),
                  V=length(unique(judgments$vignette)), v=judgments$vignette, K_V=3,
                  X_V=model.matrix(conf ~ poly(cause, 2, raw=TRUE), judgments),
                  y=judgments$conf, N_pred=length(seq(0, 1, .01)),
                  X_pred=model.matrix(~ poly(cause, 2, raw=TRUE), tibble(cause=seq(0, 1, .01))),
                  X_V_pred=model.matrix(~ poly(cause, 2, raw=TRUE), tibble(cause=seq(0, 1, .01))))

m.quad <- cmdstan_model('ordbeta_vignette.stan')
prior.quad <- m.quad$sample(c(data.quad, prior_only=TRUE), parallel_chains=4)
if (file.exists('quadratic.rds')) {
    fit.quad <- readRDS('quadratic.rds')
} else {
    fit.quad <- m.quad$sample(c(data.quad, prior_only=FALSE), parallel_chains=4)
    fit.quad$save_object('quadratic.rds')
}


## extract model coefficients
draws.coef <- left_join(spread_draws(fit.quad, b[n]),
          spread_draws(prior.quad, b[n]) %>% rename(prior_b=b)) %>%
    group_by(n) %>%
    mutate(log10_BF=bf_pointnull(b, prior_b)$log_BF / log(10),
           P=pd_to_p(as.numeric(pd(b))))
draws.coef %>%
    group_by(BF=10^log10_BF, P, .add=TRUE) %>%
    median_hdi(b)


draws.quad <- fit.quad %>%
  spread_draws(e_pred[n]) %>%
  mutate(cause=(n-1)/100)

s <- sample(draws.quad$.draw, 100)
p.quad <- draws.quad %>%
  group_by(cause) %>%
  median_hdci(e_pred) %>%
  ggplot(aes(x=cause, y=e_pred)) +
  geom_point(aes(y=conf), shape=16, alpha=.025, size=.5, data=judgments) +
  ##geom_line(aes(group=.draw), data=filter(draws.quad, .draw %in% s), alpha=.5) +
  geom_ribbon(aes(ymin=.lower, ymax=.upper), alpha=.25, fill='red') +
  geom_line(linewidth=1, color='red') +
  scale_x_continuous(name='Causal Judgment', breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
  scale_y_continuous(name='Confidence', breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
  coord_fixed(xlim=c(0, 1), ylim=c(0, 1)) +
  theme_classic(base_size=18)
p.quad
ggsave('plots/quadratic.png', width=7, height=6)

p.cause.var + p.conf.var + p.quad +
    plot_annotation(tag_levels='A') +
    plot_layout(widths=c(.25, .25, .5))
ggsave('plots/fig8.pdf', width=10, height=3.5)


predictions.quad <- tibble(p=seq(0, 1, .0001)) %>%
    mutate(var=p*(1-p),
           sd=sqrt(var),
           entropy=-p*log(p) - (1-p)*log(1-p),
           cv=sd/p) %>%
    pivot_longer(var:cv, names_to='measure') %>%
    group_by(measure) %>%
    mutate(value=1 - (value - min(value, na.rm=TRUE)) / max(value, na.rm=TRUE))


p.quad.pred <- predictions.quad %>%
    ggplot(aes(x=p, y=value, color=measure)) +
    geom_line(linewidth=1.5) +
    scale_x_continuous(name='Causal Judgment',
                       labels=c('0', '.25', '.5', '.75', '1'), expand=c(0, 0)) +
    ylab('Predicted Confidence\n(Normalized)') +
    scale_color_discrete(name='Measure', labels=c('CV', 'Entropy', 'SD', 'Variance')) +
    theme_classic(base_size=18) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    coord_fixed()


(p.quad.pred | p.quad) + plot_annotation(tag_levels='A') 
ggsave('plots/cause_conf.pdf', width=11, height=5)










############################## plots for defense #####################################
p.cause.con <- draws %>%
    filter(resp=='Causal Judgment', structure=='Conjunctive') %>%
    median_hdci(e_pred, e_pred_grad_pC, e_pred_grad_pA, P_e_pred_grad_pC, P_e_pred_grad_pA, log10_BF_e_pred_grad_pC, log10_BF_e_pred_grad_pA) %>%
    mutate(e_pred_grad_pC=ifelse(log10_BF_e_pred_grad_pC > 1 & P_e_pred_grad_pC < .05, e_pred_grad_pC, 0),
           e_pred_grad_pA=ifelse(log10_BF_e_pred_grad_pA > 1 & P_e_pred_grad_pA < .05, e_pred_grad_pA, 0)) %>%
    ggplot(aes(x=p_C, y=p_A)) +
    geom_raster(aes(fill=e_pred)) +
    scale_fill_viridis(name='Mean\nCausal\nJudgment', option='magma') +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of   ') + ylab('Probability of   ') +
    theme_classic() + coord_fixed(expand=FALSE)
p.cause.dis <- draws %>%
    filter(resp=='Causal Judgment', structure=='Disjunctive') %>%
    median_hdci(e_pred, e_pred_grad_pC, e_pred_grad_pA, P_e_pred_grad_pC, P_e_pred_grad_pA, log10_BF_e_pred_grad_pC, log10_BF_e_pred_grad_pA) %>%
    mutate(e_pred_grad_pC=ifelse(log10_BF_e_pred_grad_pC > 1 & P_e_pred_grad_pC < .05, e_pred_grad_pC, 0),
           e_pred_grad_pA=ifelse(log10_BF_e_pred_grad_pA > 1 & P_e_pred_grad_pA < .05, e_pred_grad_pA, 0)) %>%
    ggplot(aes(x=p_C, y=p_A)) +
    geom_raster(aes(fill=e_pred)) +
    scale_fill_viridis(name='Mean\nCausal\nJudgment', option='magma') +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of   ') + ylab('Probability of   ') +
    theme_classic() + coord_fixed(expand=FALSE)
p.confidence.con <- draws %>%
    filter(resp=='Confidence', structure=='Conjunctive') %>%
    median_hdci(e_pred, e_pred_grad_pC, e_pred_grad_pA, P_e_pred_grad_pC, P_e_pred_grad_pA, log10_BF_e_pred_grad_pC, log10_BF_e_pred_grad_pA) %>%
    mutate(e_pred_grad_pC=ifelse(log10_BF_e_pred_grad_pC > 1 & P_e_pred_grad_pC < .05, e_pred_grad_pC, 0),
           e_pred_grad_pA=ifelse(log10_BF_e_pred_grad_pA > 1 & P_e_pred_grad_pA < .05, e_pred_grad_pA, 0)) %>%
    ggplot(aes(x=p_C, y=p_A)) +
    geom_raster(aes(fill=e_pred)) +
    scale_fill_viridis(name='Mean\nConfidence', option='magma') +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of   ') + ylab('Probability of   ') +
    theme_classic() + coord_fixed(expand=FALSE)
p.confidence.dis <- draws %>%
    filter(resp=='Confidence', structure=='Disjunctive') %>%
    median_hdci(e_pred, e_pred_grad_pC, e_pred_grad_pA, P_e_pred_grad_pC, P_e_pred_grad_pA, log10_BF_e_pred_grad_pC, log10_BF_e_pred_grad_pA) %>%
    mutate(e_pred_grad_pC=ifelse(log10_BF_e_pred_grad_pC > 1 & P_e_pred_grad_pC < .05, e_pred_grad_pC, 0),
           e_pred_grad_pA=ifelse(log10_BF_e_pred_grad_pA > 1 & P_e_pred_grad_pA < .05, e_pred_grad_pA, 0)) %>%
    ggplot(aes(x=p_C, y=p_A)) +
    geom_raster(aes(fill=e_pred)) +
    scale_fill_viridis(name='Mean\nConfidence', option='magma') +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of   ') + ylab('Probability of   ') +
    theme_classic() + coord_fixed(expand=FALSE)


ns.draws <- cf_models %>%
    filter(strength=='Necessity-Sufficiency', precision == 'Standard Deviation') %>%
    unnest(draws.posterior) %>%
    group_by(strength, precision, structure, p_C, p_A) %>%
    median_hdci(epred_cause, epred_confidence)

ns.cause.con <- ns.draws %>%
    filter(structure == 'Conjunctive') %>%
    ggplot(aes(x=p_C, y=p_A, fill=epred_cause)) +
    geom_raster() +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of   ') + ylab('Probability of   ') +
    theme_classic() + coord_fixed(expand=FALSE)
ns.cause.dis <- ns.draws %>%
    filter(structure == 'Disjunctive') %>%
    ggplot(aes(x=p_C, y=p_A, fill=epred_cause)) +
    geom_raster() +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of   ') + ylab('Probability of   ') +
    theme_classic() + coord_fixed(expand=FALSE)
ns.confidence.con <- ns.draws %>%
    filter(structure == 'Conjunctive') %>%
    ggplot(aes(x=p_C, y=p_A, fill=epred_confidence)) +
    geom_raster() +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of   ') + ylab('Probability of   ') +
    theme_classic() + coord_fixed(expand=FALSE)
ns.confidence.dis <- ns.draws %>%
    filter(structure == 'Disjunctive') %>%
    ggplot(aes(x=p_C, y=p_A, fill=epred_confidence)) +
    geom_raster() +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of   ') + ylab('Probability of   ') +
    theme_classic() + coord_fixed(expand=FALSE)



((p.cause.con | ns.cause.con) &
 scale_fill_viridis(name='Mean\nCausal\nJudgment', option='magma',
                    limits=c(min(range.cause[1], ns.draws$epred_cause),
                             max(range.cause[2], ns.draws$epred_cause))) &
 theme_classic(base_size=18) &
 theme(axis.text=element_text(color='white', size=48))) +
    plot_layout(guides='collect')
ggsave('plots/ns_cause_con.png', width=12, height=5, dpi=200)

((p.cause.dis | ns.cause.dis) &
 scale_fill_viridis(name='Mean\nCausal\nJudgment', option='magma',
                    limits=c(min(range.cause[1], ns.draws$epred_cause),
                             max(range.cause[2], ns.draws$epred_cause))) &
 theme_classic(base_size=18) &
 theme(axis.text=element_text(color='white', size=48))) +
    plot_layout(guides='collect')
ggsave('plots/ns_cause_dis.png', width=12, height=5, dpi=200)

((p.confidence.con | ns.confidence.con) &
 scale_fill_viridis(name='Mean\nConfidence', option='magma',
                    limits=c(min(range.confidence[1], ns.draws$epred_confidence),
                             max(range.confidence[2], ns.draws$epred_confidence))) &
 theme_classic(base_size=18) &
 theme(axis.text=element_text(color='white', size=48))) +
    plot_layout(guides='collect')
ggsave('plots/ns_confidence_con.png', width=12, height=5, dpi=200)

((p.confidence.dis | ns.confidence.dis) &
 scale_fill_viridis(name='Mean\nConfidence', option='magma',
                    limits=c(min(range.confidence[1], ns.draws$epred_confidence),
                             max(range.confidence[2], ns.draws$epred_confidence))) &
 theme_classic(base_size=18) &
 theme(axis.text=element_text(color='white', size=48))) +
    plot_layout(guides='collect')
ggsave('plots/ns_confidence_dis.png', width=12, height=5, dpi=200)
