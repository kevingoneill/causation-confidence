source('utils.R')

## Define simulation values
N <- 10000
n <- 1:N         # number of samples
p_C <- prob_seq(0, 1, 0.01, delta=10/N)   # normality of candidate cause
p_A <- prob_seq(0, 1, 0.01, , delta=10/N)   # normality of alternate cause
structure <- c('Conjunctive', 'Disjunctive')

DIR <- 'plots/icard/icard-'

## Conjunctive
##  p(E) = p(C)*p(A)
##  p(E|C) = p(A)
##  p(-E|-C,A) = 1
##  K = p(C)p(E|C) + p(-C)p(-E|-C) = p(C)p(A) + p(-C)
##
## Disjunctive
##  p(E) = p(C) + p(A) - p(C)*p(A)
##  p(E|C) = 1
##  p(-E|-C,A) = 0
##  K = p(C)p(E|C) + p(-C)p(-E|-C) = p(C)

##  Make analytic predictions based on Icard et al. (2017)
df.pred <- expand_grid(p_C, p_A, structure, method='derived') %>%
    mutate(N=max(n),                  ## number of CF samples
           K=ifelse(structure=='Conjunctive',
                    p_C*p_A-p_C+1, p_C), ## causal strength estimate
           K.var=K * (1-K),
           K.sd=sqrt(K.var),
           K.se=sqrt(K.var / N),  ## SE of causal strength estimates
           K.cv=K.sd/K,
           K.entropy=-(K*log(K) + (1-K)*log(1-K)))

## Estimate those values using sampling
df.sim <- expand_grid(p_C, p_A, structure, n, method='simulated') %>%
    group_by(p_C, p_A, structure, method) %>%
    mutate(Cn=rbernoulli(n(), weight.prob(p_C)),      ## sample which events occur in sample n
           An=rbernoulli(n(), weight.prob(p_A)),
           En=ifelse(structure=='Conjunctive', pmin(Cn, An), pmax(Cn, An)),
           SSn=ifelse(structure=='Conjunctive', An, 1),
           NSn=ifelse(structure=='Conjunctive', 1, 0),
           Kn=Cn*SSn + (1-Cn)*NSn) %>% ## determine whether C caused E in sample n
    summarize(K=mean(Kn), ## estimate causal strength & certainty
              N=max(n),
              K.sd=sd(Kn),
              K.var=var(Kn),
              K.cv=sd(Kn)/mean(Kn),
              K.entropy=entropy(Kn))

## Simulated values match analytically derived values
## note- this correlation should increase with N
cor.test(df.pred$K, df.sim$K)
cor.test(df.pred$K.sd, df.sim$K.sd)
cor.test(df.pred$K.var, df.sim$K.var)
cor.test(df.pred$K.cv, df.sim$K.cv)
cor.test(df.pred$K.entropy, df.sim$K.entropy)



## Mean causal judgment
pred_plot(bind_rows(df.pred, df.sim))
ggsave(paste0(DIR, 'mean.png'))
df.pred %>% filter(structure=='Conjunctive') %>% pred_plot(facet=FALSE)
ggsave(paste0(DIR, 'mean-conjunctive.png'))
df.pred %>% filter(structure=='Disjunctive') %>% pred_plot(facet=FALSE)
ggsave(paste0(DIR, 'mean-disjunctive.png'))

## SD causal judgment
pred_plot(bind_rows(df.pred, df.sim), var='K.sd', ylab='Standard Deviation')
ggsave(paste0(DIR, 'sd.png'))
df.pred %>% filter(structure=='Conjunctive') %>%
    pred_plot(var='K.sd', ylab='Standard Deviation')
ggsave(paste0(DIR, 'sd-conjunctive.png'))
df.pred %>% filter(structure=='Disjunctive') %>%
    pred_plot(var='K.sd', ylab='Standard Deviation')
ggsave(paste0(DIR, 'sd-disjunctive.png'))

## Var causal judgment
pred_plot(bind_rows(df.pred, df.sim), var='K.var', ylab='Variance')
ggsave(paste0(DIR, 'var.png'))
df.pred %>% filter(structure=='Conjunctive') %>%
    pred_plot(var='K.var', ylab='Variance')
ggsave(paste0(DIR, 'var-conjunctive.png'))
df.pred %>% filter(structure=='Disjunctive') %>%
    pred_plot(var='K.var', ylab='Variance')
ggsave(paste0(DIR, 'var-disjunctive.png'))

## CV causal judgment
pred_plot(bind_rows(df.pred, df.sim), var='K.cv',
          ylab='Coefficient of Variation', ylim=c(0, 10))
ggsave(paste0(DIR, 'cv.png'))
df.pred %>% filter(structure=='Conjunctive') %>%
    pred_plot(var='K.cv', ylab='Coefficient of Variation', ylim=c(0, 10))
ggsave(paste0(DIR, 'cv-conjunctive.png'))
df.pred %>% filter(structure=='Disjunctive') %>%
    pred_plot(var='K.cv', ylab='Coefficient of Variation', ylim=c(0, 10))
ggsave(paste0(DIR, 'cv-disjunctive.png'))

## Entropy causal judgment
pred_plot(bind_rows(df.pred, df.sim), var='K.entropy', ylab='Entropy')
ggsave(paste0(DIR, 'entropy.png'))
df.pred %>% filter(structure=='Conjunctive') %>%
    pred_plot(var='K.entropy', ylab='Entropy')
ggsave(paste0(DIR, 'entropy-conjunctive.png'))
df.pred %>% filter(structure=='Disjunctive') %>%
    pred_plot(var='K.entropy', ylab='Entropy')
ggsave(paste0(DIR, 'entropy-disjunctive.png'))


bind_rows(df.pred, df.sim) %>%
    mutate(p_C=round(p_C*100)/100,
           p_A=round(p_A*100)/100) %>%
    ggplot() +
    aes(x=p_C, y=p_A, z=K, fill=K) +
    geom_tile() + coord_fixed() +
    geom_contour(color='grey50', alpha=0.5, bins=10) +
    scale_fill_viridis(name='Causal Strength', option='magma') +
    scale_x_continuous('P(C)', limits=0:1,
                       breaks=c(0, .5, 1), labels=c('0', '.5', '1')) +
    scale_y_continuous('P(A)', limits=0:1,
                       breaks=c(0, .5, 1), labels=c('0', '.5', '1')) +
    theme_classic() +
    facet_grid(method ~ structure)


bind_rows(df.pred, df.sim) %>%
    mutate(p_C=round(p_C*100)/100,
           p_A=round(p_A*100)/100) %>%
    ggplot() +
    aes(x=p_C, y=p_A, z=K.var, fill=K.var) +
    geom_tile() + coord_fixed() +
    geom_contour(color='grey50', alpha=0.5, bins=10) +
    scale_fill_viridis(name='Var(Causal Strength)', option='magma') +
    scale_x_continuous('P(C)', limits=0:1,
                       breaks=c(0, .5, 1), labels=c('0', '.5', '1')) +
    scale_y_continuous('P(A)', limits=0:1,
                       breaks=c(0, .5, 1), labels=c('0', '.5', '1')) +
    theme_classic() +
    facet_grid(method ~ structure)
