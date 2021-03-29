source('setup.R')
tsize <- 16

#--- Define variables
ve <- 0.75 
vp <- 0.9
T <- 90  # simulation days
n <- sum(age_demo)/T # assume that everyone vaccinated within T days

# now construct contact matrix (any R, will be rescaled)
C <- construct_C_from_prem(home=mu_home, work=mu_work, school=mu_school, other=mu_other, u=u_var,
                          target_R0=1.0, in_school=TRUE)

# define a strategy as a list of age groups 
# age groups are:
# 0-9, 10-19,20-29,...,80+, 20-29e,...,70-79e 
# so essential workers corresponding to indices 10:15

S <- list(9, 10:15, 8, 7,6,5,4,3) # 80+, EW, 70-79,...,20-29

# run
R_vec <- get_R_vec(R1=1.5, R2=2.5, start_ramp=30, end_ramp=60, ndays=T)
df1 <- run_sim_ramp_R(C, I_0=I_0, R_vec=R_vec,percent_vax =1.0, strategy=S, num_perday=n,
                    v_e = rep(ve, num_groups), v_p=rep(vp, num_groups),
                    u = u_var, num_days=T, with_essential=TRUE, H=H) 


# run with constant R vec
R_vec <- get_R_vec(R1=3.5, R2=1.5, start_ramp=30, end_ramp=60, ndays=T)
df2 <- run_sim_ramp_R(C, I_0=I_0, R_vec=R_vec, percent_vax =1.0, strategy=S, num_perday=n,
                    v_e = rep(ve, num_groups), v_p=rep(vp, num_groups),
                    u = u_var, num_days=T, with_essential=TRUE, H=H) 

# compare trajectories
trajectories <- compare_sims(sim1 = df1, 
                             sim2 = df2,
                             name1 = '80+, EW, 70-79,...,20-29', 
                             name2 = 'Oldest to Youngest', 
                             startDate=ymd("2021-01-01"), 
                             textsize = tsize)

ggarrange(plotlist=trajectories, align="v")

