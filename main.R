## main analysis

library(lubridate)
library(dplyr)
library(rlang)
library(MatchIt)
library(survival)
library(sandwich)

source("source_function.R")

# study period
study_start_date <- "2022-06-21"
study_end_date <- "2022-11-30"

min_age = 5
max_age = 12

## main cohort table
cohort <- read.csv("demo_data.csv")

# define the large table
df <- data.frame()

# define length of each trial
n = 7

## For each trial
trial_start_date <- study_start_date
i=1
# test using i=2
# i=2
# trial_start_date <- as.Date(study_start_date) + days(n)
while (as.Date(trial_start_date) < as.Date(study_end_date)){
  cat(paste0("The ", i, "-th trial with start date ", trial_start_date))
  trial_end_date <- min(as.Date(trial_start_date) + days(n), as.Date(study_end_date))
  
  ## eligible criteria
  cohort_elig <- get_eligibal_cohort(cohort, study_start_date, trial_start_date)
  
  ## vaccinated group
  data_vac <- get_vaccinated_elig(cohort_elig, trial_start_date, trial_end_date)
  
  ## unvaccinated group
  data_unvac <- get_unvaccinated_elig(cohort_elig, trial_start_date, trial_end_date)
  
  ############################# exact+ps matching ################################
  mydata <- union(data_vac, data_unvac)
  mydata <- as.data.frame(mydata)
  
  mydata$trt[mydata$trt=='vaccinated'] <- 1
  mydata$trt[mydata$trt=='unvaccinated'] <- 0
  mydata$trt <- as.numeric(mydata$trt)
  
  ##add covariates:
  # the days from last vaccination to trial start date:
  mydata$imm_day_diff <- as.Date(trial_start_date) - as.Date(mydata$pre_ba45_imm_date)
  mydata$imm_day_diff_cat <-ifelse(is.na(mydata$imm_day_diff), "no vaccine", ifelse(mydata$imm_day_diff < days(121), "< 4months", ">=4months"))
  # the days from last infection to trial start date:
  mydata$infection_day_diff <- as.Date(trial_start_date) - as.Date(mydata$last_infection_date)
  
  ## Exact matching:
  ## use exact gender, race, site
  ## vaccine doses prior to index date (0, 1/2, 2+) prior to study start date
  ## indication of pre-delta infection
  ## time since last vaccine date (no vaccine, <4 months, >= 4 months)
  ex.xvars = c("sex_cat", "eth_cat", "site", 
               "imm_day_diff_cat", "last_variant", "imm_prior_count")
  
  ## PS matching:
  ## age, obese
  ## PMCA body systems,
  ## number of outpatient/inaptient/ED visits and unique medications 24 months - 7 days prior to cohort start date (0,1,2,2+),
  ## days from delta infection to trial start date
  ## indicaiton of previous severe covid
  prop.xvars = c("age_trial_start", "obese", colnames(mydata%>%select(starts_with("pmca_"))),
                 "n_ed", "n_inpatient", "n_outpatient", "n_drug", 
                 "infection_day_diff", "person_severity_index")
  formula.prop = as.formula(paste0("trt~", paste0(prop.xvars, collapse = "+")))
  
  set.seed(1)
  m.out = matchit(formula.prop, 
                  data=mydata,
                  caliper =0.2,
                  method="nearest",
                  distance = "glm",
                  ratio=5,
                  verbose=FALSE,
                  estimand="ATT",
                  exact = ex.xvars,
                  distance.options = list(s="lambda.min"))
  m.data = match.data(m.out)
  
  ################### add follow-up time and censoring time for all individuals ######################
  # Split dataframe
  m.treated <- m.data[m.data$trt == 1, ]
  m.untreated <- m.data[m.data$trt == 0, ]
  
  ### add for unvaccinated group
  m.untreated$followup_end_date <- with(m.untreated, pmin(as.Date(ba45_imm_date), as.Date(study_end_date), na.rm = TRUE))
  m.untreated$followup_length <- as.numeric(m.untreated$followup_end_date - m.untreated$cohort_entry_date)
  
  m.untreated <- m.untreated %>%
    mutate(status = case_when((!is.na(ba45_reinfection_date) & is.na(ba45_imm_date)) | (ba45_reinfection_date <= ba45_imm_date)  ~ 1L,
                              TRUE ~ 0L))
  m.untreated$event_date <- with(m.untreated, pmin(as.Date(ba45_imm_date), as.Date(study_end_date), as.Date(ba45_reinfection_date), na.rm = TRUE))
  m.untreated$time_to_event <- as.numeric(m.untreated$event_date - m.untreated$cohort_entry_date)
  
  ### add for vaccinated group
  # first add column of the vaccination date of matched unvaccinated controls
  m.treated <- merge(m.treated, m.untreated %>% 
                       select(ba45_imm_date, subclass) %>% 
                       rename(ba45_imm_date_matched_pair = ba45_imm_date) %>%
                       group_by(subclass) %>%
                       slice_min(order_by = ba45_imm_date_matched_pair, with_ties = F, na_rm = F)%>%
                       ungroup(), by = "subclass") 
  m.treated$followup_end_date <- with(m.treated, pmin(as.Date(ba45_imm_date_matched_pair), as.Date(study_end_date), na.rm = TRUE))
  m.treated$followup_length <- as.numeric(as.Date(m.treated$followup_end_date) - as.Date(m.treated$cohort_entry_date))
  
  m.treated <- m.treated %>%
    mutate(status = case_when(!is.na(ba45_reinfection_date) & is.na(ba45_imm_date_matched_pair)  ~ 1L,
                              TRUE ~ 0L))
  m.treated$event_date <- with(m.treated, pmin(as.Date(ba45_imm_date_matched_pair), as.Date(study_end_date), as.Date(ba45_reinfection_date), na.rm = TRUE))
  m.treated$time_to_event <- as.numeric(as.Date(m.treated$event_date) - as.Date(m.treated$cohort_entry_date))
  
  # Merge back
  m.untreated$ba45_imm_date_matched_pair <- NA
  m.data.updated <- union(m.treated, m.untreated)
  m.data.updated$trial_id <- i
  m.data.updated$trial_start_date <- trial_start_date
  m.data.updated$trial_end_date <- trial_end_date
  m.data.updated$documented_infection <- ifelse(m.data.updated$ba45_reinfection_date <= m.data.updated$followup_end_date, 1, 0) 
  m.data.updated$documented_infection[is.na(m.data.updated$documented_infection)] <- 0
  
  if (i==1){
    df <- rbind(df, m.data.updated)
  }else{
    df <- df %>% union(m.data.updated)
  }
  
  i=i+1
  trial_start_date <- format(as.Date(trial_end_date, origin="1970-01-01"))
}

df$infection_day_diff <- as.numeric(df$infection_day_diff)
df$imm_day_diff <- as.numeric(df$imm_day_diff)

################################# outcome model ####################################
cox <- coxph(Surv(time_to_event, status) ~ trt,
             data = df)

print(1-exp(coef(cox)))
print(1-exp(confint(cox)))


