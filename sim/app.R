library(shiny)
source('utils.R')

## Derived causal judgments from each model
deltaP <- function(p_C, p_A, structure)
    ifelse(structure=='Conjunctive', p_A, 1-p_A)
PPC <- function(p_C, p_A, structure)
    ifelse(structure=='Conjunctive', p_A, 1)
SP <- function(p_C, p_A, structure)
    ifelse(structure=='Conjunctive', (1-p_C)*p_A, (1-p_C) * (1-p_A))
Icard <- function(p_C, p_A, structure)
    ifelse(structure=='Conjunctive', p_C*p_A-p_C+1, p_C)
Quillien <- function(p_C, p_A, structure)
    ifelse(structure=='Conjunctive',
           sqrt((1-p_C)*p_A / (1-p_C*p_A)),
           sqrt((1-p_A)*p_C / (p_C + p_A - p_C*p_A)))


model.predictions <- function(kappa=1.0, structure='Conjunctive') {
    expand_grid(pC=prob_seq(0, 1, 0.1, delta=1e-10),
                pA=prob_seq(0, 1, 0.1, delta=1e-10),
                structure, kappa) %>%
        mutate(pC.w=weight.prob(pC, kappa),
               pA.w=weight.prob(pA, kappa),
               pE=ifelse(structure=='Conjunctive', pC.w*pA.w, pC.w+pA.w-pC.w*pA.w),
               deltaP.C=deltaP(pC.w, pA.w, structure),
               deltaP.A=deltaP(pA.w, pC.w, structure),
               PPC.C=PPC(pC.w, pA.w, structure),
               PPC.A=PPC(pA.w, pC.w, structure),
               SP.C=SP(pC.w, pA.w, structure),
               SP.A=SP(pA.w, pC.w, structure),
               Icard.C=Icard(pC.w, pA.w, structure),
               Icard.A=Icard(pA.w, pC.w, structure),
               Quillien.C=Quillien(pC.w, pA.w, structure),
               Quillien.A=Quillien(pA.w, pC.w, structure),
               pC=ifelse(pmin(pC, 1-pC) < .09, round(pC), pC),
               pA=ifelse(pmin(pA, 1-pA) < .09, round(pA), pA)) %>%
        pivot_longer(deltaP.C:Quillien.A, names_to=c('model', 'variable'), names_sep='\\.', values_to='K') %>%
        pivot_wider(names_from='variable', values_from='K', names_prefix='K.') %>%
        group_by(model) %>%
        mutate(model=factor(model, levels=c('deltaP', 'PPC', 'SP', 'Icard', 'Quillien')),
               norm=ifelse(model=='Quillien', sqrt(pC.w*(1-pC.w)/(pE*(1-pE))), 1),
               raw=K.C/norm,
               K.norm=normalize(exp(K.C)/(exp(K.C)+exp(K.A))),
               var=1-normalize(norm^2 * raw*(1-raw)),
               sd=1-normalize(sqrt(norm^2 * raw*(1-raw))),
               cv=sd/K.C,
               entropy=1-normalize(-(raw*log(raw) + (1-raw)*log(1-raw + 1e-12))),
               extremity=normalize(abs(0.5 - K.C)),
               distance=normalize(sqrt(pC^2 + pA^2)),
               log_ratio=normalize(log(K.C/K.A)))
}


## 1 - (1-pC)*pA
## 1 - (pA - pCpA)
## 1 - pA + pCpA
## pCpA + 1 - pA

## pC*pA-pC+1
## pCpA + (1 - pC)

## 1-(1-pC)*(1-pA)
## 1 - (1 - pC - pA + pCpA)
## 1 - 1 + pC + pA - pCpA
## pC + pA - pCpA

## sqrt(pC^2 + pA^2)

ui <- fluidPage(titlePanel('Causal Judgment Model Predictions'),
    sidebarLayout(
        sidebarPanel(radioButtons('causalStructure', 'Causal Structure', choices=c('Conjunctive', 'Disjunctive'), inline=TRUE),
                     radioButtons('plotType', 'Plot Type', choices=c('1D', '2D'), inline=TRUE),
                     plotOutput('weightPlot'),
                     sliderInput('kappa', 'Probability Weight',
                                 -1.5, 1.5, 0.0, step=.1)),
        mainPanel(tabsetPanel(tabPanel('Causal Strength', plotOutput('plotCause')),
                              tabPanel('Causal Strength (normalized)', plotOutput('plotCauseNorm')),
                              tabPanel('Standard Deviation', plotOutput('plotSD')),
                              tabPanel('Variance', plotOutput('plotVar')),
                              tabPanel('Entropy', plotOutput('plotEntropy')),
                              tabPanel('Extremity', plotOutput('plotExtremity')),
                              tabPanel('Distance', plotOutput('plotDistance')),
                              tabPanel('Log Ratio', plotOutput('plotLogRatio'))))))

plot1D <- function(df, dv='K.C', label=ifelse(dv[1]=='K', 'Causal Strength', 'Confidence')) {
    ggplot(df) +
        aes_string(x='pC', y=dv, group='pA', color='pA') +
        geom_line(size=1) +
        scale_x_continuous('P(C)', breaks=c(0, .5, 1), labels=c('0', '.5', '1')) +
        scale_y_continuous(label, breaks=c(0, .5, 1), labels=c('0', '.5', '1')) +
        coord_fixed(expand=FALSE) +
        theme_classic(base_size=18) +
        scale_color_viridis(name='P(A)', option='magma') +
        facet_wrap(~ model)
}

plot2D <- function(df, dv='K.C', label=ifelse(dv[1]=='K', 'Causal\nStrength', 'Confidence')) {
    ggplot(df, aes_string(x='pC', y='pA', z=dv, fill=dv)) +
        geom_tile() + coord_fixed(expand=FALSE) +
        geom_contour(color='grey50', alpha=0.5, bins=5) +
        scale_fill_viridis(name=label, option='magma') +
        scale_x_continuous('P(C)', breaks=c(0, .5, 1), labels=c('0', '.5', '1')) +
        scale_y_continuous('P(A)', breaks=c(0, .5, 1), labels=c('0', '.5', '1')) +
        theme_classic(base_size=18) +
        facet_wrap(~ model)    
}

plot <- function(df, dv='K.C', label=ifelse(dv[1]=='K', 'Causal\nStrength', 'Confidence'), plotType='1D') {
    if (plotType == '1D') {
        plot1D(df, dv, label=label)
    } else {
        plot2D(df, dv, label=label)
    }
}

server <- function(input, output) {
    pred <- reactive({
        model.predictions(exp(input$kappa), input$causalStructure)
    })
    
    output$plotCause <- renderPlot({
        plot(pred(), 'K.C', plotType=input$plotType)
    })

    output$plotCauseNorm <- renderPlot({
        plot(pred(), 'K.norm', plotType=input$plotType)
    })

    output$plotSD <- renderPlot({
        plot(pred(), 'sd', plotType=input$plotType)
    })

    output$plotVar <- renderPlot({
        plot(pred(), 'var', plotType=input$plotType)
    })

    output$plotEntropy <- renderPlot({
        plot(pred(), 'entropy', plotType=input$plotType)
    })

    output$plotExtremity <- renderPlot({
        plot(pred(), 'extremity', plotType=input$plotType)
    })

    output$plotDistance <- renderPlot({
        plot(pred(), 'distance', plotType=input$plotType)
    })

    output$plotLogRatio <- renderPlot({
        plot(pred(), 'log_ratio', plotType=input$plotType)
    })
    
    output$weightPlot <- renderPlot({
        data.frame(p=seq(0, 1, .01))%>%
            mutate(w=weight.prob(p, exp(input$kappa))) %>%
            ggplot(aes(x=p, y=w)) +
            geom_polygon(fill='red', alpha=0.25) +
            geom_abline(intercept=0, slope=1) +
            geom_line(color='red', size=1) +
            scale_x_continuous('Probability', breaks=c(0, .5, 1), labels=c('0', '.5', '1')) +
            scale_y_continuous('Perceived Probability', breaks=c(0, .5, 1), labels=c('0', '.5', '1')) +
            coord_fixed(expand=FALSE) +
            theme_bw(base_size=18)
    })
}

shinyApp(ui = ui, server = server)
