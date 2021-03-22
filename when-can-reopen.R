# when can we reopen to R = x ? 

source('setup_ON.R')


phac.accept = c(72,72, 75.18, 78.25, 84.7, 84.7, 84.70)
# 15-19, 20-44, 45-54, 55-64, 65+, 75+, 80+ 
# these will give a hes like this:
H = c(0,0,0.28, 0.28, 0.235, 0.2, 0.153, 0.153, 0.153) 
H=c(H, H[3:8])*N_i # with EW and scaled to ON pop size 

H=0.25*H
# rollout info from PHAC had 3.4M# (population, not doses) in Q1, 12.75M in Q2, vague on Q3 and Q4 

# i think i can continue w 80+ in jan - feb at a low rate
# but then in Q2 they have some EW and 60-79 l probably model as 60-79, then EW, 
# then rest in Q3 to be done by september, much like the basic BC rollout

# so the question is: 
# - pick a strategy, and a baseline R to go with it, using 2-phase sim 
# - simulate until time T
# - at that time, make a much higher R 
# - detect whether cases still rise a lot, or is it now OK to reopen without a huge rise
# or with only a small rise 
# see run_over_scen_3 below - can look manually, get intuition, then expand to 
# code that returns a reopening -ready  time if desired. 

# need to look at what the threshold is for hosps 
# phac has 37 bed / 100K but that is acute and ICU, not sure if it is all hosp
# but for now let's use it. that is a hosp of  37*pop_total/1e5 # of 5500
# but Dec numbers on the website were more like 1600, probably a more reasonable threshold

test1 = run_over_scen_3(1.1, 2.25, trytime = 120,T2=720)
test2 = run_over_scen_3(1.1, 2.25, trytime = 210,T2=720)

tmp = compare_sims(test1, test2,textsize = 10)
ggarrange(plotlist = tmp[1:2], nrow=1)

run_over_scen_3 = function(R1, R2, trytime,T2=210, ve=0.75, vp=0.9,  scen=1,alpha=0.0){
    T1 <- 60
#    T2 <- 360
    # Initial stage (vax all 80+)
    R_init <- 1.05
    n <- age_demo[9]/T1
    C <- construct_C_from_prem(home=mu_home, work=mu_work, school=mu_school, other=mu_other, u=u_var,
                               target_R0=R_init, in_school=TRUE, alpha_factor=alpha)
    
    df0 <- run_sim_basic(C, I_0=I_0, percent_vax =1.0, strategy=list(9), num_perday=n,
                         v_e = rep(ve, num_groups), v_p=rep(vp, num_groups),
                         u = u_var, num_days=T1, with_essential=TRUE, H=H) 

        # next stage : R moves to R1 and rate increases
    n <- sum(age_demo[-9])/210
    C <- construct_C_from_prem(home=mu_home, work=mu_work, school=mu_school, other=mu_other, u=u_var,
                               target_R0=R1, in_school=TRUE, alpha_factor=alpha)
    df1 <- run_sim_restart(C, df_0=tail(df0, n=1), percent_vax =1.0, strategy= strategies[[scen]], num_perday=n,
                          v_e = rep(ve, num_groups), v_p=rep(vp, num_groups),
                          u = u_var, num_days=trytime, with_essential=TRUE, H=H)

        # final stage L R moves to R2, rate stays the same 
    C <- construct_C_from_prem(home=mu_home, work=mu_work, school=mu_school, other=mu_other, u=u_var,
                               target_R0=R2, in_school=TRUE, alpha_factor=alpha)
    df <- run_sim_restart(C, df_0=tail(df1, n=1), percent_vax =1.0, strategy= strategies[[scen]], num_perday=n,
                          v_e = rep(ve, num_groups), v_p=rep(vp, num_groups),
                          u = u_var, num_days=T2-trytime, with_essential=TRUE, H=H)
    # combine 
    df1$time = df1$time+T1+1
    df$time = df$time+T1+trytime+1
    
    df <- combine_age_groups(rbind(df0,df1,df))

    # add pars
    df$R1 <- R1
    df$R2 <- R2
    df$trytime = trytime
    df$ve <- ve
    df$vp <- vp
    df$type <- labels[[scen]]
    df$scen <- scen
    df$alpha <- alpha
    return(df)}


