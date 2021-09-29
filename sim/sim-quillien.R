source('utils.R')

## Define simulation values
e <- exp(1)
N <- 100000
n <- 1:N         # number of samples
p_C <- prob_seq(0, 1, 0.01, delta=0.005)   # normality of candidate cause
p_A <- prob_seq(0, 1, 0.1, , delta=0.005)   # normality of alternate cause
structure <- c('Conjunctive', 'Disjunctive')

DIR <- 'plots/quillien/quillien-'

##  Make analytic predictions based on Quillien (2020)
df.pred <- expand_grid(p_C, p_A, structure, method='derived') %>%
    mutate(N=max(n),                  ## number of CF samples
           p_E=ifelse(structure=='Conjunctive', p_C*p_A, p_C+p_A-p_C*p_A),
           varC=p_C * (1-p_C),
           varA=p_A * (1-p_A),
           varE=p_E * (1-p_E),
           K=ifelse(structure=='Conjunctive',
                    sqrt((1-p_C)*p_A / (1-p_C*p_A + 1e-6)),
                    sqrt((1-p_A)*p_C / (p_C + p_A - p_C*p_A))),
           Ka=ifelse(structure=='Conjunctive',
                     sqrt((1-p_A)*p_C / (1-p_A*p_C)),
                     sqrt((1-p_C)*p_A / (p_A + p_C + p_A*p_C))),
           normC=sqrt(varC/varE),
           K.var=ifelse(structure=='Conjunctive',
           (1-p_C)*(1-p_A)/(1-p_C*p_A),
           p_C*p_A / (p_C+p_A-p_C*p_A)),
           K.se=sqrt(K.var / N),
           K.cv=sqrt(K.var)/K,
           K.entropy=-p_A*log(p_A) - (1-p_A)*log(1-p_A))
           
## Estimate those values using sampling
df.sim <- expand_grid(p_C, p_A, structure, n, method='simulated') %>%
    group_by(p_C, p_A, structure) %>%
    mutate(Cn=rbernoulli(n(), weight.prob(p_C, gamma=.5)),      ## sample whether C occurs in sample n
           An=rbernoulli(n(), weight.prob(p_A, gamma=.5)),      ## sample whether A occurs in sample n
           Cn.twin=1-Cn,                 ## intervene on Cn
           An.twin=1-An,                 ## intervene on An
           ## determine the effect in all 3 worlds
           En=ifelse(structure=='Conjunctive',
                     pmin(Cn, An), pmax(Cn, An)),
           En.twinC=ifelse(structure=='Conjunctive',
                           pmin(Cn.twin, An), pmax(Cn.twin, An)),
           En.twinA=ifelse(structure=='Conjunctive',
                           pmin(Cn, An.twin), pmax(Cn, An.twin)),
           ## determine the specific causal effects
           normC=sd(Cn) / sd(En),
           normA=sd(An) / sd(En),
           Knc=(En.twinC-En)/(Cn.twin-Cn),
           Kna=(En.twinA-En)/(An.twin-An),
           Knc.norm=softmax(Knc*normC, Kna*normA),
           Kna.norm=softmax(Kna*normA, Knc*normC)) %>% 
    ## estimate causal strength & certainty
    group_by(p_C, p_A, structure, method) %>%
    summarize(N=max(n),
              normC=mean(normC), normA=mean(normA),
              K=mean(Knc*normC), Ka.mean=mean(Kna*normA),
              K.sd=sd(Knc*normC), K.var=var(Knc*normC), K.se=sd(Knc*normC)/sqrt(N))


pred_plot(df.sim, 'K', ylim=0:1)
pred_plot(df.sim, 'K.sd')

pred_plot(df.sim, 'Knorm.mean', ylim=0:1)
pred_plot(df.sim, 'Knorm.sd')


pred_plot(df.sim, 'Knorm2.mean', ylim=0:1)
pred_plot(df.sim, 'Knorm2.sd')



## Simulated values match analytically derived values
## note- this correlation should increase with N
cor.test(df.pred$K, df.sim$K)
cor.test(df.pred$K.se, df.sim$K.se)

cor.test(df.pred$Knorm, df.sim$Knorm)
cor.test(df.pred$Knorm.se, df.sim$Knorm.se)


## Join analytic and simulated predictions into one df
pred_plot(bind_rows(df.pred, df.sim), ylim=0:1)
ggsave('Quillien_mean_predicitons.png')
pred_plot(bind_rows(df.pred, df.sim), var='K.var', ylab='Confidence', trans='reverse')
ggsave('Quillien_confidence_predicitons.png')


pred_plot(bind_rows(df.pred, df.sim), var='Knorm', ylim=0:1)
ggsave('Quillien_norm_mean_predicitons.png')


pred_plot(bind_rows(df.pred, df.sim), var='Knorm.var', ylab='Confidence', trans='reverse')

ggsave('Quillien_norm_confidence_predicitons.png')



pred_plot(bind_rows(df.pred, df.sim), var='K.cv', trans='reverse',
          ylab='Coefficient of Variation', ylim=c(3.25, 0))
ggsave('Quillien_CV_predicitons.png')

pred_plot(bind_rows(df.pred, df.sim), var='K.entropy', trans='reverse',
          ylab='Entropy')
ggsave('Quillien_entropy_predicitons.png')
