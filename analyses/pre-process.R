library(dplyr)
library(tidyr)
library(stringr)

DIR <- 'data/'
INPUT_FILE <- 'data-raw.csv'
OUTPUT_FILE <- 'data-processed.csv'

read.csv(paste0(DIR, INPUT_FILE)) %>% tibble %>%
    slice(-1, -2) %>%
    select(!(StartDate:Progress) & !(Finished:RecordedDate) & !(RecipientLastName:Instruct8_1) & !(PROLIFIC_PID:SESSION_ID)) %>%
    rename(duration=Duration..in.seconds.,
           id=ResponseId,
           age=Age, sex=Sex, attn_check=AttnCheck) %>%
    mutate_at(vars(ends_with('_1') | ends_with('probC') | ends_with('probA')), as.numeric) %>%
    mutate_at(vars(ends_with('probC') | ends_with('probA')), function(x) x/10) %>%
    mutate(duration=as.numeric(duration)/60,
           age=as.numeric(age)) %>%
    rename_all(function(s) str_remove_all(s, '_1')) %>%
    pivot_longer(cards_con_cause:wheels_dis_conf,
                 names_pattern='(.*)_(.*)_(.*)', names_to=c('vignette', 'structure', 'rating')) %>%
    filter(!is.na(value)) %>%
    pivot_wider(names_from='rating') %>%
    pivot_longer(cards_probC:wheels_probA, names_pattern='(.*)_(.*)', names_to=c('vignette_prob', 'measure')) %>%
    pivot_wider(names_from='measure') %>%
    filter(vignette == vignette_prob) %>%
    select(!vignette_prob) %>%
    rename(p_C=probC, p_A=probA) %>%
    relocate(id, vignette, structure, p_C, p_A, cause, conf, attn_check, duration, age, sex) %>%
    write.csv(paste0(DIR, OUTPUT_FILE), row.names=FALSE)
