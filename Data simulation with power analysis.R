library(lme4)
library(tidyverse)
library(faux)
library(lmerTest)

# Primary reference: Barr, Dale J. (2021). Learning statistical models through simulation in R: An interactive textbook. 
# Version 1.0.0. Retrieved from https://psyteachr.github.io/stat-models-v1.

# DATA SIMULATION 

# Reading in, organizing data
datagraph = read.csv('/Users/maidamuhtar/Documents/Lakarprogrammet/Kurs 8/OCT_Aurora_Tegaderm_Data File ver.3.1_Graphs.csv', sep=';')
# Access data through: <iframe src="https://widgets.figshare.com/articles/17041157/embed?show_title=1" width="568" height="351" allowfullscreen frameborder="0"></iframe>

# Changing the variables 'Participant', 'Site', 'Tegaderm' into factors in datagraph
datagraph$Participant = factor(datagraph$Participant, levels=c(1,2,3,4,5,6), labels=c('P001','P002','P003','P004','P005','P006'))
datagraph$Site = factor(datagraph$Site, levels=c(1,2,3), labels=c('Forearm','Hypothenar/thenar','Finger Pad')) 
datagraph$Tegaderm = factor(datagraph$Tegaderm, levels=c(0,1), labels=c('Bare','Film'))

datagraphFA = (filter(datagraph, Site == 'Forearm')) 
datagraphFA = arrange(datagraphFA, 'Participant','Site','Tegaderm','Observation_number')

datFA <- subset(datagraphFA, select = -c(X10pc_P, X50pc_P, X50pc_D, X10pc_D, FNo_Maxind, FNo_SD, FNo_Tip, FNo_ED, FNo_Tot, Mx_Ind_10pc, Mx_Ind_50pc, Observation_number, Data_number))

summary(datFA)

# Creating a mixed effects model on existing data 
mod_datFA <- lmer(Mx_Ind ~ Tegaderm + (Tegaderm | Participant), datFA)
summary(mod_datFA)
mod_datFA_sum <- summary(mod_datFA)

# Parameters
beta_0 <- 0.8924 # mean max indentation for bare skin 
beta_1 <- -0.2189 # mean difference between Film - Bare, i.e. the effect of Film on max indentation
tau_0 <- 0.2061 # by-subject random intercept sd 
tau_1 <- 0.07218 # by-subject random slope sd 
rho <- -0.9361 # correlation between intercept and slope
sigma <- 0.06289 # residual (error) sd 

# Set number of subjects and items
n_subj <- 6 # number of subjects
n_bare <- 30 # number of bare stimuli
n_film <- 30 # number of film stimuli
obs_n <- 5 # observation number, number of repetitions

# Simulate a sample of items
# Total number of items = n_bare + n_film
items <- data.frame(obs_n = seq_len(obs_n),
                    Tegaderm = rep(c("bare","film"), c(n_bare, n_film)))

# Effect-code category
items$X_i <- recode(items$Tegaderm,"bare" = 0, "film" = 1)

# Simulate a sample of subjects
# Calculate random intercept / random slope covariance
covar <- rho * tau_0 * tau_1

# Put values into variance-covariancematrix
cov_mx <- matrix(c(tau_0^2, covar,
                   covar, tau_1^2),
                 nrow = 2, byrow = TRUE)

# Generate the by-subject random effects
subject_rfx <- MASS::mvrnorm(n = n_subj, mu = c(T_0s = 0, T_1s = 0), Sigma = cov_mx) 

# Combine with subject IDs
subjects <- data.frame(subj_id = seq_len(n_subj), subject_rfx)

# Cross subject and item IDs; add an error term
# nrow(.) is the number of rows in the table 
trials.sim <- crossing(subjects, items) %>%
  mutate(e_si = rnorm(nrow(.), mean = 0, sd = sigma)) %>%
  select(subj_id, obs_n, Tegaderm, X_i,
         everything())

# Calculate the response variable
datFA_sim <- trials.sim %>% 
  mutate(Mx_Ind_s = beta_0 + T_0s + (beta_1 + T_1s) * X_i + e_si) %>%
  select(subj_id, obs_n, Tegaderm, X_i, Mx_Ind_s)

# Visualizing and checking simulated data, comparing it to collected data

# Graph 1 - simulated data
ggplot(datFA_sim, aes(Tegaderm, Mx_Ind_s, color = Tegaderm)) +
  geom_hline(yintercept = beta_0) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9))

# Graph 1 - collected data
ggplot(datFA, aes(Tegaderm, Mx_Ind, color = Tegaderm)) +
  geom_hline(yintercept = 0.8924) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9))

# Graph 2 - simulated data
ggplot(datFA_sim, aes(Tegaderm, Mx_Ind_s, color=factor(subj_id))) +
  geom_jitter(width=0, size=2) + 
  ggtitle('Forearm') +
  guides(color=guide_legend('Subject', 
                            title.theme = (element_text(size = 11)),
                            label.theme = (element_text(size = 11)))) +
  facet_wrap(. ~ subj_id, nrow = 3) +
  theme_bw() + 
  ylab('Indentation depth (mm)') + 
  xlab ('Skin condition')

# Graph 2 - collected data
ggplot(datFA, aes(Tegaderm, Mx_Ind, color=Participant)) +
  geom_jitter(width=0, size=2) + 
  ggtitle('Forearm') +
  guides(color=guide_legend('Subject', 
                            title.theme = (element_text(size = 11)),
                            label.theme = (element_text(size = 11)))) +
  facet_wrap(. ~ Participant, nrow = 3) +
  theme_bw() + 
  ylab('Indentation depth (mm)') + 
  xlab ('Skin condition')

# Creating a mixed effects model for the simulated data 
mod_datFA_sim <- lmer(Mx_Ind_s ~ Tegaderm + (1 + Tegaderm | subj_id), data = datFA_sim)
summary(mod_datFA_sim)
mod_datFA_simsum <- summary(mod_datFA_sim)

# Parameters in the mixed effects model for the simulated data
# Number of subjects
mod_datFA_simsum$ngrps
# Random factors
mod_datFA_simsum$varcor
# Fixed factors
mod_datFA_simsum$coefficients

# Parameters in the mixed effects model for the collected data
# Number of subjects
mod_datFA_sum$ngrps
# Random factors
mod_datFA_sum$varcor
# Fixed factors
mod_datFA_sum$coefficients

# Get a tidy table of results
broom.mixed::tidy(mod_datFA_sim) %>%
  mutate(sim = c(beta_0, beta_1, tau_0,
                 rho, tau_1, sigma)) %>%
  select(1:3, 9, 4:8)

# Set up the custom data simulation function
my_sim_datFA <- function(
  n_subj = 6, # number of subjects
  n_bare = 30, # number of bare stimuli 
  n_film = 30, # number of film stimuli
  obs_n = 5, # number of repetitions, observation number: 1-5
  beta_0 = 0.8924, # grand mean
  beta_1 = -0.2189, # effect of film
  tau_0 = 0.2061, # by-subject random intercept sd
  tau_1 = 0.07218, # by-subject random slope sd
  rho = -0.9361, # correlation between intercept and slope
  sigma = 0.06389 # residual (standard deviation) 
){
  items <- data.frame(obs_n = seq_len(obs_n), 
                      Tegaderm = rep(c("bare","film"), c(n_bare, n_film)),
                      X_i = rep(c(0, 1), c(n_bare, n_film)))
                      
  # Variance-covariance matrix
  cov_mx <- matrix(
    c(tau_0^2, rho * tau_0 * tau_1,
      rho * tau_0 * tau_1, tau_1^2 ),
    nrow = 2, byrow = TRUE)
  
  subjects <- data.frame(subj_id = seq_len(n_subj), 
                         MASS::mvrnorm(n = n_subj, 
                                       mu = c(T_0s = 0, T_1s = 0),
                                       Sigma = cov_mx))
  
  crossing(subjects, items) %>%
    mutate(e_si = rnorm(nrow(.), mean = 0, sd = sigma),
           Mx_Ind_s = beta_0 + T_0s + (beta_1 + T_1s) * X_i + e_si) %>%
    select(subj_id, obs_n,Tegaderm ,X_i, Mx_Ind_s)
  }

# Simulate, analyze, and return a table of parameter estimates
single_run <- function(...) {
  
# . . . is a shortcut that forwards any arguments to 
  # my_sim_data(), the function created above

  dat_sim <- my_sim_datFA(...)
  mod_sim <- lmer(Mx_Ind_s ~ Tegaderm + (1 + Tegaderm | subj_id), data = dat_sim)
  broom.mixed::tidy(mod_sim) 
}

# Run one model with default parameters
single_run()

# Run simulations and save to a file
n_runs <- 1000 # use at least 1000 to get stable estimates
sims_FA <- purrr::map_df(1:n_runs, ~ single_run()) 
write_csv(sims_FA, "sims_FA.csv")

# Read saved simulation data
sims_FA <- read_csv("/Users/maidamuhtar/Documents/R/sims_FA.csv", col_types = 
                   cols( 
                     group = col_factor(ordered = TRUE),
                     term = col_factor(ordered = TRUE)
                     ))
                   sims_FA %>%
  filter(effect == "fixed") %>%
  select(term, estimate, p.value)
                   
# POWER ANALYSIS

# Calculate mean estimates and power for specified alpha
alpha <- 0.05 
sims_FA %>%
filter(effect == "fixed") %>%
group_by(term) %>%
summarize(
  mean_estimate = mean(estimate),
  mean_se = mean(std.error),
  power = mean(p.value < alpha),
  .groups = "drop"
  )

# Visualizing
sims_FA_intercept <- sims_FA %>%
  filter(term == "(Intercept)")

sims_FA_slope <- sims_FA %>%
  filter(term == "Tegadermfilm")

# Intercept - max indentation on bare 
# Blue dashed line = measured mean max indentation bare skin
# Black dashed line = simulated mean max indentation bare skin
ggplot(sims_FA_intercept, aes(estimate)) +
  geom_histogram(color="green", fill="white") +
  geom_vline(aes(xintercept=mean(estimate)),
             color="black", linetype="dashed", size=1) +
  geom_vline(aes(xintercept= 0.89524),
             color="blue", linetype="dashed", size=1) +
  theme_light()

# Slope - effect of film
# Blue dashed line = measured effect of film on max indentation
# Black dashed line = simulated effect of film on max indentation
ggplot(sims_FA_slope, aes(estimate)) +
  geom_histogram(color="red", fill="white") +
  geom_vline(aes(xintercept=mean(estimate)),
             color="black", linetype="dashed", size=1) +
  geom_vline(aes(xintercept= - 0.2189),
             color="blue", linetype="dashed", size=1) +
  theme_light()









