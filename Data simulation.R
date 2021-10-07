library(lme4)
library(tidyverse)
library(faux)
library(lmerTest)

#DATA SIMULATION 

#Reading in, organizing data
datagraph = read.csv('/Users/maidamuhtar/Documents/Lakarprogrammet/Kurs 8/OCT_Aurora_Tegaderm_Data File ver.3.1_Graphs.csv', sep=';')

#Changing the variables 'Participant', 'Site', 'Tegaderm' into factors in datagraph
datagraph$Participant = factor(datagraph$Participant, levels=c(1,2,3,4,5,6), labels=c('P001','P002','P003','P004','P005','P006'))
datagraph$Site = factor(datagraph$Site, levels=c(1,2,3), labels=c('Forearm','Hypothenar/thenar','Finger Pad')) 
datagraph$Tegaderm = factor(datagraph$Tegaderm, levels=c(0,1), labels=c('Bare','Film'))

datagraphFA = (filter(datagraph, Site == 'Forearm')) 
datagraphFA = arrange(datagraphFA, 'Participant','Site','Tegaderm','Observation_number')

datFA <- subset(datagraphFA, select = -c(X10pc_P, X50pc_P, X50pc_D, X10pc_D, FNo_Maxind, FNo_SD, FNo_Tip, FNo_ED, FNo_Tot, Mx_Ind_10pc, Mx_Ind_50pc, Observation_number, Data_number))

#Creating a model on existing data 
mod2 <- lmer(Mx_Ind ~ Tegaderm + (Tegaderm | Participant), datFA)
summary(mod2)
mod2.sum <- summary(mod2)

summary(datFA)

#Parameters
sub_n <- 6 # number of subjects in this simulation
sub_sd <- 0.20610 # SD for the subjects' random intercept
stim_n <- 60  # number of stimuli in this simulation
grand_i <- 0.7798 # overall mean Mx_Ind
eff_teg <- 0.21890  # mean difference between versions: Bare - Film
error_sd <- 0.06389 # residual (error) SD
sub_slope_sd <- 0.07218 # slope SD
sub_i_slope_cor <- -0.9361 #correlation between intercept and slope ? 

#Subjects 
sub <- faux::rnorm_multi(
  n = sub_n, 
  vars = 2, 
  r = sub_i_slope_cor,
  mu = 0, # means of random intercepts and slopes are always 0
  sd = c(sub_sd, sub_slope_sd),
  varnames = c("sub_i", "sub_teg_slope")
) %>%
  mutate(
    sub_id = 1:sub_n
  )

?rnorm_multi()

#Checking
ggplot(sub, aes(sub_i, sub_teg_slope)) +
  geom_point() +
  geom_smooth(method = lm)

#Stimuli
mod3 <- lmer(Mx_Ind ~ 1 + (1 | Participant) + (1 | Tegaderm), data = datFA) #new model to get the stimuli intercept ? 

stim <- ranef(mod3)$Tegaderm %>%
  as_tibble(rownames = 'Tegaderm') %>%
  rename(stim_i = `(Intercept)`)

view(stim)

summary(mod3)

ggplot(stim, aes(stim_i)) +
  geom_density()

# Putting together subject and sitmuli  
trials <- crossing(
  sub_id = sub$sub_id, # get subject IDs from the sub data table
  Tegaderm = stim$Tegaderm, # get stimulus IDs from the stim data table
  obs_n = c(1, 2, 3, 4, 5),
) %>%
  left_join(sub, by = "sub_id") %>% # includes the intercept and condition for each subject
  left_join(stim, by = "Tegaderm")   # includes the intercept for each stimulus

# Calculating dependent variable (Mx_Ind) and creating a table 
datFA.s <- trials %>%
  mutate(
    # effect-code subject condition and stimulus version
    Tegaderm.e = recode(Tegaderm, "Film" = -0.5, "Bare" = +0.5),
    # calculate trial-specific effects by adding overall effects and slopes
    #version_eff = stim_version_eff + stim_version_slope + sub_version_slope,
    film_eff = eff_teg + sub_teg_slope,
    # calculate error term (normally distributed residual with SD set above)
    err = rnorm(nrow(.), 0, error_sd),
    # calculate DV from intercepts, effects, and error
    Mx_Ind.s = grand_i + sub_i + stim_i + err +
      (Tegaderm.e * film_eff)
  )

#Visualizing and checking data, comparing it to observed data 

ggplot(datFA.s, aes(Tegaderm, Mx_Ind.s, color = Tegaderm)) +
  geom_hline(yintercept = grand_i) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9))

mod.s <- lmer(Mx_Ind.s ~ Tegaderm + (Tegaderm| sub_id), datFA.s)
summary(mod.s)
mod.s_sum <- summary(mod.s)

mod2 <- lmer(Mx_Ind ~ Tegaderm + (Tegaderm | Participant), datFA)
summary(mod2)
mod2.sum <- summary(mod2)

mod.s_sum$ngrps
mod.s_sum$varcor
mod.s_sum$coefficients

mod2.sum$ngrps
mod2.sum$varcor
mod2.sum$coefficients

ggplot(datFA, aes(Tegaderm, Mx_Ind, color = Tegaderm)) +
  geom_hline(yintercept = grand_i) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9))

#Creating a function 
sim_lmer <- function(sub_n = 6, # number of subjects in this simulation
                     sub_sd = 0.20610, # SD for the subjects' random intercept
                     stim_n = 60,  # number of stimuli in this simulation
                     grand_i = 0.7798, # overall mean Mx_Ind
                     eff_teg = 0.21890,  # mean difference between versions: Bare - Film
                     error_sd = 0.06389, # residual (error) SD
                     sub_slope_sd = 0.07218, # slope SD
                     sub_i_slope_cor = -0.9361 #correlation between intercept and slope ? 
){ 
  sub <- faux::rnorm_multi(
    n = sub_n, 
    vars = 2, 
    r = sub_i_slope_cor,
    mu = 0, # means of random intercepts and slopes are always 0
    sd = c(sub_sd, sub_slope_sd),
    varnames = c("sub_i", "sub_teg_slope")
  ) %>%
    mutate(
      sub_id = 1:sub_n
    )
  
  mod3 <- lmer(Mx_Ind ~ 1 + (1 | Participant) + (1 | Tegaderm), data = datFA) #new model to get the stimuli intercept ? 
  
  stim <- ranef(mod3)$Tegaderm %>%
    as_tibble(rownames = 'Tegaderm') %>%
    rename(stim_i = `(Intercept)`)
  
  trials <- crossing(
    sub_id = sub$sub_id, # get subject IDs from the sub data table
    Tegaderm = stim$Tegaderm, # get stimulus IDs from the stim data table
    obs_n = c(1, 2, 3, 4, 5),
  ) %>%
    left_join(sub, by = "sub_id") %>% # includes the intercept and condition for each subject
    left_join(stim, by = "Tegaderm")   # includes the intercept for each stimulus
  
  datFA.s <- trials %>%
    mutate(
      # effect-code subject condition and stimulus version
      Tegaderm.e = recode(Tegaderm, "Film" = -0.5, "Bare" = +0.5),
      # calculate trial-specific effects by adding overall effects and slopes
      #version_eff = stim_version_eff + stim_version_slope + sub_version_slope,
      film_eff = eff_teg + sub_teg_slope,
      # calculate error term (normally distributed residual with SD set above)
      err = rnorm(nrow(.), 0, error_sd),
      # calculate DV from intercepts, effects, and error
      Mx_Ind.s = grand_i + sub_i + stim_i + err +
        (Tegaderm.e * film_eff)
    )
  
  mod.s <- lmer(Mx_Ind.s ~ Tegaderm + (Tegaderm| sub_id), datFA.s)
  mod.s_sum <- summary(mod.s)
  return(mod.s_sum)
}

sim_lmer()

#BIN
#Changing the variables 'Participant', 'Site', 'Tegaderm' into factors in datagraph
datFA.s$sub_id = factor(datFA.s$sub_id, levels=c(1,2,3,4,5,6), labels=c('P001','P002','P003','P004','P005','P006'))
