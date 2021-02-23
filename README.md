# Strategies for vaccine allocation with essential workers

This analysis accompanies the text, *Vaccine Rollout Strategies: The Case for Vaccinating Essential Workers Early* (Mulberry et al., 2021). 
The preprint is [here](https://github.com/nmulberry/essential-workers-vaccine/vaccine_rollout_preprint.pdf).

The core model and simulation functions are adapted from [Bubar et al.](https://github.com/kbubar/vaccine_prioritization). 

## Code Layout

* analysis : contains code for reproducing the figures from the paper.
  * analysis/setup.R : defines key rates & the scenarios used in the paper.
* data : contains the input data used in our analysis.
* vaccine-model.R : core model and simulation functions.
* utils.R :  helper functions for plotting & summarizing output.
* contact-matrix.R : functions for building the contact matrix with essential workers. 

## Simple Example
```R
# define key parameters according to BC rates & demographics, 
#   & build the 4 contact matrices with essential workers
source('setup.R') 

#--- Define variables
ve <- 0.75 
vp <- 0.9
T <- 200  # simulation days
n <- sum(age_demo)/T # assume that everyone vaccinated within T days
R <- 1.3 #R value for entire simulation

# now construct contact matrix with target R
C <- construct_C_from_prem(home=mu_home, work=mu_work, school=mu_school, other=mu_other, u=u_var,
                          target_R0=R, in_school=TRUE)

# define a strategy as a list of age groups 
# age groups are:
# 0-9, 10-19,20-29,...,80+, 20-29e,...,70-79e 
# so essential workers corresponding to indices 10:15

S <- list(9, 10:15, 8, 7,6,5,4,3) # 80+, EW, 70-79,...,20-29

# run
df1 <- run_sim_basic(C, I_0=I_0, percent_vax =1.0, strategy=S, num_perday=n,
                    v_e = rep(ve, num_groups), v_p=rep(vp, num_groups),
                    u = u_var, num_days=T, with_essential=TRUE, H=H) 

# try another strategy
S <- list(9,c(8,15),c(7,14), c(6,13), c(5,12), c(4,11), c(3,10)) # oldest to youngest (EWs are vaccinated with their age group)

df2 <- run_sim_basic(C, I_0=I_0, percent_vax =1.0, strategy=S, num_perday=n,
                    v_e = rep(ve, num_groups), v_p=rep(vp, num_groups),
                    u = u_var, num_days=T, with_essential=TRUE, H=H) 

# compare trajectories
trajectories <- compare_sims(sim1 = df1, 
                             sim2 = df2,
                             name1 = '80+, EW, 70-79,...,20-29', 
                             name2 = 'Oldest to Youngest', 
                             startDate=ymd("2021-01-01"), 
                             textsize = 16)

ggarrange(plotlist=trajectories)

```
