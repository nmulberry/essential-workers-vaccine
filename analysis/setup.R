library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(furrr)
library(cowplot)
library(lubridate)
library(ggpubr)
no_cores <- availableCores() - 1
plan(multicore, workers = no_cores)
theme_set(theme_light())

source('../vaccine-model.R') # main model & simulation funcs
source('../utils.R') 
source('../contact-matrix.R')

# BUILD CONTACT MATRIX
p_ess <- c(0.0,0.0,0.1705,0.2043,0.1675,0.1536,0.1642,0.1069, 0.0) 
age_demo_by_five <- as.matrix(readr::read_csv('~/essential-workers-vaccine/data/Population_Estimates.csv'))
generate_contact_matrix(p_ess, age_demo_by_five) #build

age_demo <- readRDS("~/essential-workers-vaccine/generated-data/age_demographics_essential.rds")
mu_home <- readRDS("~/essential-workers-vaccine/generated-data/mu_home_essential.rds")
mu_work <- readRDS("~/essential-workers-vaccine/generated-data/mu_work_essential.rds")
mu_school <- readRDS("~/essential-workers-vaccine/generated-data/mu_school_essential.rds")
mu_other <-  readRDS("~/essential-workers-vaccine/generated-data/mu_other_essential.rds")

# RATES (10-yr age bins, 0-9,10-19,...,80+)
YLL_vec <- readRDS("~/essential-workers-vaccine/data/yll_vec_CAN.RData")
IFR <- readRDS("~/essential-workers-vaccine/data/ifr_vec_CAN.RData")
IHR <- readRDS("~/essential-workers-vaccine/data/ihr_vec_CAN.RData")

# Symptom duration lognormal. log sigma=0.8. log  mu are 1.9, 2.2, 2.5, 2.8,
# for ages under 30, 30-39, 40-49, over 50.
# means are 9, 12, 17 and 22
SDUR_vec = c(9,9,9, 12, 17, 22,22,22,22)  # symptom dur Sudre et al and Paul

LCR = c(0.037, 0.037, 0.037, 0.079,0.15,0.252, 0.252, 0.252, 0.252) # from Paul from Sudre
longDUR = c(rep(40.3,3), 42.5, 45.6, rep(49.7,4)) 
HOSPDUR = c(4.692308,3.380952,5.260870,4.304348,7.212598,7.874346,9.092857,12.752475,12.964346)


# ADD ESSENTIAL WORKERS TO IFR 
IFR <- c(IFR, IFR[3:8])


# GLOBAL MODEL PARAMETERS
n <- length(age_demo)
pop_total <- age_demo[n] #pop BC
age_demo <- age_demo[1:n-1]
N_i <-  pop_total*age_demo
num_groups <- length(age_demo)


sero_none <- rep(0, num_groups) # no prior immunity
startDate <- lubridate::ymd("2021-01-01")

I_0 <- c(133.0, 241.0, 488.8, 364.0, 315.2, 243.2, 203.2, 111.2, 155.0, 122.2,  91.0,  78.8,
   60.8,  50.8,  27.8) # estimated prevalence by age in BC in Jan 2021

u_var <- c(0.38, 0.4, 0.79, 0.86, 0.8, 0.82, 0.88,1.0, 1.0,
  0.79, 0.86, 0.8, 0.82, 1.0,1.0) # rel. susceptibility (this works better with validation than one in bubar)

percent_vax <- 1.0 # just a limit we can't exceed

H = c(0.0,0.0,0.3,0.2,0.2,0.2,0.15,0.15,0.15,0.3,0.2,0.2,0.2,0.15,0.15)*N_i # Hesitancy
#H = c(0.0,0.0,0.3,0.2,0.2,0.2,0.01,0.01,0.01,0.3,0.2,0.2,0.2,0.01,0.01)*N_i # Hesitancy
##########################
# DEFINE MAIN STRATEGIES
##########################

labels <- c('A: Oldest to Youngest', 
            'B: 80+, 20-79',
            'C: 80+, EW, 70-79,...', 
            'D: 80+, EW, 20-79',
            'E: 80+, EW, 70-79, 20-69',
            'F: 80+, 70-79, EW, 20-69',
            'G: 80+, EW, 70-79, 60-69, 20-59')

strategies <- list(list(9,c(8,15),c(7,14), c(6,13), c(5,12), c(4,11), c(3,10)), 
                   list(9,3:15),
                   list(9, 10:15, 8, 7,6,5,4,3),
                   list(9, 10:15, 3:8),
                   list(9, 10:15, 8, 3:7),
                   list(9,c(8,15), 10:14, 3:7),
                   list(9,10:15, 8,7,3:6))


#######################
# DEFINE MAIN SCENARIO
# This is the main piecewise
# scenario. Start with first phase
# of vaccinating 80+ slowly. Then
# move on to a second phase of
# faster vaccination 
########################
run_over_scen_2 = function(R, ve, vp, scen,alpha=0.0){
   T1 <- 60
   T2 <- 210
   # Initial stage (vax all 80+)
   R_init <- 1.05
   n <- age_demo[9]/T1
   C <- construct_C_from_prem(home=mu_home, work=mu_work, school=mu_school, other=mu_other, u=u_var,
                              target_R0=R_init, in_school=TRUE, alpha_factor=alpha)

   df0 <- run_sim_basic(C, I_0=I_0, percent_vax =1.0, strategy=list(9), num_perday=n,
                        v_e = rep(ve, num_groups), v_p=rep(vp, num_groups),
                        u = u_var, num_days=T1, with_essential=TRUE, H=H) 
   # Final stage
   n <- sum(age_demo[-9])/T2
   C <- construct_C_from_prem(home=mu_home, work=mu_work, school=mu_school, other=mu_other, u=u_var,
                              target_R0=R, in_school=TRUE, alpha_factor=alpha)
   df <- run_sim_restart(C, df_0=tail(df0, n=1), percent_vax =1.0, strategy= strategies[[scen]], num_perday=n,
                         v_e = rep(ve, num_groups), v_p=rep(vp, num_groups),
                         u = u_var, num_days=T2, with_essential=TRUE, H=H)
   # combine 
   df$time <- df$time+T1+1
   df <- combine_age_groups(rbind(df0,df))
   # add pars
   df$R <- R
   df$ve <- ve
   df$vp <- vp
   df$type <- labels[[scen]]
   df$scen <- scen
   df$alpha <- alpha
   return(df)}


########################
# Another scenario
# no piecewise
# constant R and rate (n)
########################

run_over_scen = function(R, ve, vp, scen, n, alpha=0.0){
   T <- 270 

   C <- construct_C_from_prem(home=mu_home, work=mu_work, school=mu_school, other=mu_other, u=u_var,
      target_R0=R, in_school=TRUE, alpha_factor=alpha)
    
   df <- run_sim_basic(C, I_0=I_0, percent_vax =1.0, strategy= alt_strategies[[scen]], num_perday=n,
                           v_e = rep(ve, num_groups),v_p=rep(vp, num_groups),
                           u = u_var, num_days=T, with_essential=TRUE, H=H)

   df <- combine_age_groups(df)
   df$R0 <- R
   df$ve <- ve
   df$vp <- vp
   df$type <- alt_labels[[scen]]
   df$scen <- scen
   df$alpha <- alpha
   return(df)}

