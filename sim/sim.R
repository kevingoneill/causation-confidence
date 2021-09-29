source('utils.R')

## Define simulation values
p_C <- prob_seq(0, 1, 0.1, delta=1e-12)   # normality of candidate cause
p_A <- prob_seq(0, 1, 0.1, , delta=1e-12)   # normality of alternate cause
structure <- c('Conjunctive', 'Disjunctive')

## normalize to [0, 1]
normalize <- function(x) if(max(x, na.rm=TRUE) == 0) {x} else {x / max(x, na.rm=TRUE)}

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
           K.var=normalize(norm^2 * K.raw*(1-K.raw)),
           K.sd=normalize(sqrt(K.var)),
           K.cv=K.sd/K,
           K.entropy=-(K.raw*log(K.raw) + (1-K.raw)*log(1-K.raw + 1e-12)))


ggplot(df.pred) +
    aes(x=p_C, y=p_A, fill=K) +
    geom_tile() + coord_fixed() +
    ##geom_quiver(aes(u=dK.dp_C, v=dK.dp_A), center=TRUE) +
    geom_contour(aes(z=K), color='grey50', alpha=0.5, bins=10) +
    scale_fill_viridis(name='Causal Strength', option='magma') +
    scale_x_continuous('P(C)', limits=0:1,
                       breaks=c(0, .5, 1), labels=c('0', '.5', '1')) +
    scale_y_continuous('P(A)', limits=0:1,
                       breaks=c(0, .5, 1), labels=c('0', '.5', '1')) +
    theme_classic() +
    facet_grid(model ~ structure)
ggsave('plots/predictions-causal-strength-2d.png', width=5, height=7.5)

ggplot(df.pred) +
    aes(x=p_C, y=p_A, fill=K.sd) +
    geom_tile() + coord_fixed() +
    geom_contour(aes(z=K.sd), color='grey50', alpha=0.5, bins=5) +
    scale_fill_viridis(name='SD(Causal Strength)', option='magma') +
    scale_x_continuous('P(C)', limits=0:1,
                       breaks=c(0, .5, 1), labels=c('0', '.5', '1')) +
    scale_y_continuous('P(A)', limits=0:1,
                       breaks=c(0, .5, 1), labels=c('0', '.5', '1')) +
    theme_classic() +
    facet_grid(model ~ structure)
ggsave('plots/predictions-sd-2d.png', width=5, height=7.5)

ggplot(df.pred) +
    aes(x=p_C, y=p_A, fill=K.var) +
    geom_tile() + coord_fixed() +
    geom_contour(aes(z=K.var), color='grey50', alpha=0.5, bins=5) +
    scale_fill_viridis(name='Var(Causal Strength)', option='magma') +
    scale_x_continuous('P(C)', limits=0:1,
                       breaks=c(0, .5, 1), labels=c('0', '.5', '1')) +
    scale_y_continuous('P(A)', limits=0:1,
                       breaks=c(0, .5, 1), labels=c('0', '.5', '1')) +
    theme_classic() +
    facet_grid(model ~ structure)
ggsave('plots/predictions-var-2d.png', width=5, height=7.5)





p <- ggplot(df.pred) +
    aes(x=p_C, y=K, group=p_A, color=p_A) +
    geom_line(size=1) + ylab('Causal Strength') +
    scale_color_viridis(name='P(A)') +
    scale_x_continuous('P(C)', limits=0:1,
                       breaks=c(0, .5, 1), labels=c('0', '.5', '1')) +
    theme_bw() +
    facet_grid(model ~ structure) +
    transition_reveal(p_A)

anim_save('plots/prediction-causal-strength.mp4',
          animate(p, renderer=av_renderer(), fps=20))

p <- ggplot(df.pred) +
    aes(x=p_C, y=K.sd, group=p_A, color=p_A) +
    geom_line(size=1) + ylab('Causal Strength') +
    scale_color_viridis(name='P(A)') +
    scale_x_continuous('P(C)', limits=0:1,
                       breaks=c(0, .5, 1), labels=c('0', '.5', '1')) +
    theme_bw() +
    facet_grid(model ~ structure) +
    transition_time(p_A)

anim_save('plots/prediction-sd2.mp4',
          animate(p, renderer=av_renderer(), fps=25))

ggsave('plots/predictions-causal-strength.png', width=5, height=7.5)

ggplot(df.pred) +
    aes(x=p_C, y=K.sd, group=p_A, color=p_A) +
    geom_line(size=1) + ylab('SD(Causal Strength)') +
    scale_color_viridis(name='P(A)') +
    scale_x_continuous('P(C)', limits=0:1,
                       breaks=c(0, .5, 1), labels=c('0', '.5', '1')) +
    theme_bw() +
    facet_grid(model ~ structure, scales='free_y')
ggsave('plots/predictions-sd.png', width=5, height=7.5)

ggplot(df.pred) +
    aes(x=p_C, y=K.var, group=p_A, color=p_A) +
    geom_line(size=1) + ylab('Var(Causal Strength)') +
    scale_color_viridis(name='P(A)') +
    scale_x_continuous('P(C)', limits=0:1,
                       breaks=c(0, .5, 1), labels=c('0', '.5', '1')) +
    theme_bw() +
    facet_grid(model ~ structure, scales='free_y')
ggsave('plots/predictions-var.png', width=5, height=7.5)

ggplot(df.pred) +
    aes(x=p_C, y=K.cv, group=p_A, color=p_A) +
    geom_line(size=1) + ylab('CV(Causal Strength)') +
    coord_cartesian(ylim=c(0, 10)) +
    scale_color_viridis(name='P(A)') +
    scale_x_continuous('P(C)', limits=0:1,
                       breaks=c(0, .5, 1), labels=c('0', '.5', '1')) +
    theme_bw() +
    facet_grid(model ~ structure)

ggsave('plots/predictions-cv.png', width=5, height=7.5)

ggplot(df.pred) +
    aes(x=p_C, y=K.entropy, group=p_A, color=p_A) +
    geom_line(size=1) + ylab('Entropy(Causal Strength)') +
    scale_color_viridis(name='P(A)') +
    scale_x_continuous('P(C)', limits=0:1,
                       breaks=c(0, .5, 1), labels=c('0', '.5', '1')) +
    theme_bw() +
    facet_grid(model ~ structure)
ggsave('plots/predictions-entropy.png', width=5, height=7.5)
