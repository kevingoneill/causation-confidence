source('utils.R')

## Define simulation values
N <- 10000
n <- 1:N         # number of samples
p_C <- prob_seq(0, 1, 0.01, delta=10/N)   # normality of candidate cause
p_A <- prob_seq(0, 1, 0.1, , delta=10/N)   # normality of alternate cause
structure <- c('Conjunctive', 'Disjunctive')

DIR <- 'plots/power-pc/power-pc-'

## Conjunctive
##  p(E) = p(C)*p(A)
##  p(E|C) = p(A)
##  p(E|-C) = 0
##  p(-E|-C) = 1
##  PPC = (p(E|C) - p(E|-C)) / p(-E|-C) = (p(A) - 0) / 1 = p(A)
##
## Disjunctive
##  p(E) = p(C) + p(A) - p(C)*p(A)
##  p(E|C) = 1
##  p(E|-C) = p(A)
##  p(-E|-C) = 1 - p(A)
##  PPC = (p(E|C) - p(E|-C)) / p(-E|-C) = (1 - p(A)) / (1 - p(A)) = 1


##  Make analytic predictions based on Icard et al. (2017)
df.pred <- expand_grid(p_C, p_A, structure, method='derived') %>%
    mutate(N=max(n),                  ## number of CF samples
           K=ifelse(structure=='Conjunctive', p_A, 1), ## causal strength estimate
           K.var=ifelse(structure=='Conjunctive', p_A*(1-p_A), 0),
           K.sd=sqrt(K.var),
           K.cv=K.sd/K,
           K.entropy=-(K*log(K) + (1-K)*log(1-K)))

## Estimate those values using sampling
df.sim <- expand_grid(p_C, p_A, structure, n, method='simulated') %>%
    group_by(p_C, p_A, structure, method) %>%
    mutate(Cn=rbernoulli(n(), p_C),      ## sample which events occur in sample n
           An=rbernoulli(n(), p_A),
           En=ifelse(structure=='Conjunctive', pmin(Cn, An), pmax(Cn, An)),
           En.given.Cn=ifelse(structure=='Conjunctive', An, 1),
           En.given.notCn=ifelse(structure=='Conjunctive', 0, An),
           notEn.given.notCn=ifelse(structure=='Conjunctive', 1, 1-An),
           Kn=ifelse(notEn.given.notCn == 0, 1,
           (En.given.Cn - En.given.notCn)/notEn.given.notCn)) %>%
    summarize(K=mean(Kn, na.rm=TRUE), ## estimate causal strength & certainty
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
