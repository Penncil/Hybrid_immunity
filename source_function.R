## source functions
## use BA45 as an example

get_eligibal_cohort <- function(cohort, study_start_date, trial_start_date){
  ############################# eligible criteria #############################
  ## no (re)infection or vaccination from cohort start date to trial start date
  positive_prior = cohort %>% filter((ba45_reinfection_date > as.Date(study_start_date))& (ba45_reinfection_date < as.Date(trial_start_date))) 
  immu_first_prior = cohort %>% filter((ba45_imm_date > as.Date(study_start_date))& (ba45_imm_date < as.Date(trial_start_date)))
  
  cohort1 <- cohort %>% 
    anti_join(positive_prior, by='person_id') %>%
    anti_join(immu_first_prior, by='person_id')
  
  ## age constraint and categorize
  cohort2 <- cohort1 %>%
    mutate(age_trial_start = (as.Date(trial_start_date)-as.Date(birth_date))/365.25) %>%
    filter(age_trial_start >= min_age, age_trial_start < max_age)
  
  return(cohort2)
}

get_vaccinated_elig <- function(cohort, trial_start_date, trial_end_date){
  ## vaccinated during trial period
  cohort_vac <- cohort %>% 
    filter(ba45_imm_date >= as.Date(trial_start_date), ba45_imm_date < as.Date(trial_end_date))
  
  ## exclude patients with (re)infection before vaccination
  positive_prior_vac = cohort_vac %>% filter(ba45_reinfection_date <= as.Date(ba45_imm_date))
  cohort_vac_elig <- cohort_vac %>% anti_join(positive_prior_vac, by='person_id') %>%
    mutate(cohort_entry_date = as.Date(ba45_imm_date)) %>%
    mutate(trt = "vaccinated")
  
  return(cohort_vac_elig)
}


get_unvaccinated_elig <- function(cohort, trial_start_date, trial_end_date){
  ## vaccinated before the trial end
  cohort_vac <- cohort %>% 
    filter(ba45_imm_date < as.Date(trial_end_date))
  
  ## not vaccinated during trial period
  cohort_unvac <- cohort %>% anti_join(cohort_vac, by='person_id')
  
  ## add additional situation that reinfection prior to date of vaccination
  cohort_unvac_c <- cohort %>%
    filter(ba45_imm_date < as.Date(trial_end_date), ba45_imm_date > as.Date(trial_start_date))%>%
    filter(ba45_reinfection_date < as.Date(trial_end_date), ba45_reinfection_date > as.Date(trial_start_date))%>%
    filter(ba45_imm_date > ba45_reinfection_date)
  
  #combine, index date= trial start date
  cohort_unvac <- cohort_unvac %>% full_join(cohort_unvac_c) %>%
    mutate(cohort_entry_date = as.Date(trial_start_date)) %>%
    mutate(trt="unvaccinated")
  
  return(cohort_unvac)
}

