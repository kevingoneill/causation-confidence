library(dplyr)
library(tidyr)
library(ggplot2)
library(rstan)
library(modelr)
library(bayestestR)
library(tidybayes)
library(viridis)
library(scico)
library(patchwork)
library(ggquiver)
library(latex2exp)

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



## Define simulation values
p_C <- prob_seq(.1, 1, .1)    # normality of candidate cause
p_A <- prob_seq(.1, 1, .1)   # normality of alternate cause
structure <- c('Conjunctive', 'Disjunctive')

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

df.pred <- expand_grid(structure, p_C, p_A) %>%
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
    mutate(grid_index=row_number()) %>%
    expand_grid(d=0:2) %>%
    mutate(grid_index = grid_index + d*max(grid_index),
           d=factor(d, labels=c('Function', 'dP(C)', 'dP(A)'))) %>%
    relocate(grid_index)

judgments <- judgments %>%
    left_join(grid %>% filter(d=='Function') %>% select(-d),
              by=c('p_C', 'p_A'))


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
if (file.exists('gp-prior.rds')) {
    gp.prior <- readRDS('gp-prior.rds')
} else {
    data.stan$prior_only <- TRUE
    gp.prior <- stan(file='gp-gradient.stan', data=data.stan, cores=4)
    saveRDS(gp.prior, 'gp-prior.rds')
}


prior_draws <- gp.prior %>%
    spread_draws(rho[structure, dx], alpha[structure, resp], phi[resp], cutpoints[resp, idx],
                 f[structure, grid_index, resp]) %>%
    mutate(resp=factor(resp, levels=1:2, labels=c('Causal Judgment', 'Confidence')),
           structure=levels(judgments$structure)[structure],
           dx=factor(dx, levels=1:2, labels=c('p(C)', 'p(A)')),
           idx=idx-1) %>%
    left_join(grid, by='grid_index') %>%
    group_by(structure, p_C, p_A, resp) %>%
    select(-grid_index) %>%
    pivot_wider(names_from=idx, values_from=cutpoints, names_prefix='cutpoint.') %>%
    pivot_wider(names_from=dx, values_from=rho, names_prefix='rho.') %>%
    mutate(p.0=1 - plogis(f-cutpoint.0),
           p.cont=plogis(f-cutpoint.0) - plogis(f-cutpoint.1),
           p.1=plogis(f-cutpoint.1),
           mu=ifelse(d=='Function', p.0*0 + p.cont*plogis(f) + p.1*1, f))



####################################################################################################
##                                          GP posterior
####################################################################################################
if (file.exists('gp-gradient.rds')) {
    gp.fit <- readRDS('gp-gradient.rds')
} else {
    data.stan$prior_only <- FALSE
    gp.fit <- stan(file='gp-gradient.stan', data=data.stan, cores=4)
    saveRDS(gp.fit, 'gp-gradient.rds')
}

print(gp.fit, pars=c('rho', 'alpha', 'phi', 'cutpoints', 'Omega'), prob=c(.025, .5, .975))
plot(gp.fit, pars=c('rho', 'alpha', 'phi', 'cutpoints', 'Omega'))


draws <- gp.fit %>%
    spread_draws(rho[structure, dx], alpha[structure, resp], phi[resp], cutpoints[resp, idx],
                 f[structure, grid_index, resp]) %>%
    mutate(resp=factor(resp, levels=1:2, labels=c('Causal Judgment', 'Confidence')),
           structure=levels(judgments$structure)[structure],
           dx=factor(dx, levels=1:2, labels=c('p(C)', 'p(A)')),
           idx=idx-1) %>%
    left_join(grid, by='grid_index') %>%
    group_by(structure, p_C, p_A, resp) %>%
    select(-grid_index) %>%
    pivot_wider(names_from=idx, values_from=cutpoints, names_prefix='cutpoint.') %>%
    pivot_wider(names_from=dx, values_from=rho, names_prefix='rho.') %>%
    mutate(p.0=1 - plogis(f-cutpoint.0),
           p.cont=plogis(f-cutpoint.0) - plogis(f-cutpoint.1),
           p.1=plogis(f-cutpoint.1),
           mu=ifelse(d=='Function', p.0*0 + p.cont*plogis(f) + p.1*1, f))

## join prior/posterior into one DF
draws <- prior_draws %>%
    rename(prior_alpha=alpha,
           prior_phi=phi,
           prior_f=f,
           prior_cutpoint.0=cutpoint.0,
           prior_cutpoint.1=cutpoint.1,
           `prior_rho.p(C)`=`rho.p(C)`,
           `prior_rho.p(A)`=`rho.p(A)`,
           prior_p.0=p.0,
           prior_p.cont=p.cont,
           prior_p.1=p.1,
           prior_mu=mu) %>%
    full_join(draws) %>%
    group_by(structure, p_C, p_A, resp, d)

## Function draws
draws %>%
    filter(d=='Function', p_A %in% c(.1, .5, 1)) %>%
    median_hdi(mu) %>%    
    ggplot(aes(x=p_C, y=mu, fill=as.factor(p_A))) +
    geom_line(aes(color=as.factor(p_A)), size=2) +
    geom_ribbon(aes(ymin=.lower, ymax=.upper), alpha=0.1) +
    scale_color_viridis(discrete=TRUE, name='Probability of\nAlternate Cause') +
    scale_fill_viridis(discrete=TRUE, name='Probability of\nAlternate Cause') +
    xlab('Probability of Focal Cause') + ylab('Estimated Marginal Mean') +
    facet_grid(resp ~ structure, scales='free_y') +
    theme_classic()
ggsave('plots/gp-beta.png', width=6, height=4)






####################################################################################################
##                                   Compare Models to GP
####################################################################################################

df.pred <- df.pred %>%
    mutate(model=factor(model, levels=c('deltaP', 'PPC', 'SP', 'Icard', 'Quillien'),
                        labels=c('Delta P', 'Power PC', 'Crediting Causality',
                                 'Necessity-Sufficiency', 'CES'))) %>%
    group_by(model)


p.cause.data <- draws %>%
    filter(resp=='Causal Judgment') %>%
    median_hdi(mu) %>%
    rename(.value=mu) %>%
    pivot_wider(names_from=d, values_from=c(.value, .lower, .upper), names_glue='{d}{.value}') %>%
    mutate(`dP(C).value`=ifelse(`dP(C).lower` > 0 | `dP(C).upper` < 0, `dP(C).value`, 0),
           `dP(A).value`=ifelse(`dP(A).lower` > 0 | `dP(A).upper` < 0, `dP(A).value`, 0)) %>%
    ggplot(aes(x=p_C, y=p_A, u=`dP(C).value`, v=`dP(A).value`)) +
    geom_raster(aes(fill=Function.value)) +
    geom_quiver(color='white', size=.25, center=TRUE) +
    scale_fill_viridis(name='Mean\nCausal\nJudgment', option='magma') +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of Focal Cause') + ylab('Probability of Alternate Cause') +
    facet_grid(structure ~ resp, labeller=labeller(resp=c('Causal Judgment'='Data'))) +
    theme_classic() + coord_fixed()

p.cause.pred <- ggplot(df.pred, aes(x=p_C, y=p_A, fill=K)) +
    geom_raster() +
    facet_grid(structure ~ model) +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of Focal Cause') + ylab('Probability of Alternate Cause') +
    scale_fill_viridis(name=TeX('$\\kappa_{C\\rightarrow E}$'), option='magma', limits=c(0,1)) +
    theme_classic() + coord_fixed()

(p.cause.data | p.cause.pred) + plot_layout(widths=c(1/6, 5/6))
ggsave('plots/gp-model-cause.png', width=12, height=4)



p.conf.data <- draws %>%
    filter(resp=='Confidence') %>%
    median_hdci(mu) %>%
    rename(.value=mu) %>%
    pivot_wider(names_from=d, values_from=c(.value, .lower, .upper), names_glue='{d}{.value}') %>%
    mutate(`dP(C).value`=ifelse(`dP(C).lower` > 0 | `dP(C).upper` < 0, `dP(C).value`, 0),
           `dP(A).value`=ifelse(`dP(A).lower` > 0 | `dP(A).upper` < 0, `dP(A).value`, 0)) %>%
    ggplot(aes(x=p_C, y=p_A, u=`dP(C).value`, v=`dP(A).value`)) +
    geom_raster(aes(fill=Function.value)) +
    geom_quiver(color='white', size=.25, center=TRUE) +
    scale_fill_viridis(name='Mean\nConfidence', option='magma') +
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1), labels=c('0', '.25', '.5', '.75', '1')) +
    xlab('Probability of Focal Cause') + ylab('Probability of Alternate Cause') +
    facet_grid(structure ~ resp, labeller=labeller(resp=c('Confidence'='Data'))) +
    theme_classic() + coord_fixed()

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
(p.conf.data | p.conf.pred) + plot_layout(widths=c(1/6, 5/6))
ggsave('plots/gp-model-var.png', width=12, height=4)


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
(p.conf.data | p.conf.pred) + plot_layout(widths=c(1/6, 5/6))
ggsave('plots/gp-model-sd.png', width=12, height=4)


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
(p.conf.data | p.conf.pred) + plot_layout(widths=c(1/6, 5/6))
ggsave('plots/gp-model-cv.png', width=12, height=4)


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
(p.conf.data | p.conf.pred) + plot_layout(widths=c(1/6, 5/6))
ggsave('plots/gp-model-entropy.png', width=12, height=4)





draws.mu <- draws %>% ungroup %>%
    filter(d == 'Function') %>%
    select(structure, p_C, p_A, resp, .chain, .iteration, .draw, mu) %>%
    mutate(p_C=as.character(p_C), p_A=as.character(p_A)) %>%
    pivot_wider(names_from=resp, values_from=mu) %>%
    left_join(df.pred.wide %>% mutate(p_C=as.character(p_C), p_A=as.character(p_A)))

draws.model <- draws.mu %>%
    group_by(.draw, structure) %>%
    summarize(across(K_deltaP:K_Quillien, ~cor(`Causal Judgment`, .)),
              across(K.var_deltaP:K.var_Quillien, ~cor(Confidence, .)),
              across(K.sd_deltaP:K.sd_Quillien, ~cor(Confidence, .)),
              across(K.cv_deltaP:K.cv_Quillien, ~cor(Confidence, .)),
              across(K.entropy_deltaP:K.entropy_Quillien, ~cor(Confidence, .))) %>%
    pivot_longer(K_deltaP:K.entropy_Quillien, names_to=c('variable', 'model'), names_sep='_') %>%
    ## reverse correlations for confidence
    mutate(value=ifelse(variable=='K', value, -value),
           model=factor(model, levels=c('deltaP', 'PPC', 'SP', 'Icard', 'Quillien'),
                        labels=c('Delta P', 'Power PC', 'Crediting Causality', 'Necessity-Sufficiency', 'Counterfactual Effect Size'))) %>%
    replace_na(list(value=0))

draws.model %>%
    filter(variable=='K' | variable=='K.sd') %>%
    group_by(variable, structure, model) %>%
    median_hdi(value) %>%
    arrange(variable, structure, desc(value))

draws.model %>%
    filter(variable=='K' | variable=='K.sd') %>%
    ggplot(aes(x=value, y=model, group=structure, fill=structure)) +
    geom_vline(xintercept=0, linetype='dashed') +
    stat_halfeye(aes(side=ifelse(structure=='Conjunctive', 'top', 'bottom')),
                 position=position_dodge(width=-.5),
                 point_interval='median_hdi',
                 normalize='panels') +
    scale_fill_discrete(name='Causal Structure') +
    facet_wrap(~variable, labeller=labeller(variable=c('K'='Causal Judgment', 'K.sd'='Confidence'))) +
    ylab('Model') + xlab('Correlation') +
    theme_classic() +
    theme(legend.position='bottom')
ggsave('plots/model_comparisons.png', width=10, height=5)


draws.model %>%
    filter(variable!='K' & variable!='K.sd') %>%
    mutate(variable=factor(variable, levels=c('K.var', 'K.entropy', 'K.cv'))) %>%
    ggplot(aes(x=value, y=model, group=structure, fill=structure)) +
    geom_vline(xintercept=0, linetype='dashed') +
    stat_halfeye(aes(side=ifelse(structure=='Conjunctive', 'top', 'bottom')),
                 position=position_dodge(width=-.5),
                 point_interval='median_hdi',
                 normalize='panels') +
    scale_fill_discrete(name='Causal Structure') +
    facet_wrap(~variable, labeller=labeller(variable=c('K.var'='Variance', 'K.cv'='Coefficient of Variation',
                                                       'K.entropy'='Entropy'))) +
    ylab('Model') + xlab('Correlation') +
    theme_classic() +
    theme(legend.position='bottom')
ggsave('plots/model_comparisons_confidence.png', width=10, height=5)


