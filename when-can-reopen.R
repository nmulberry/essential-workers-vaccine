# when can we reopen to R = x ? 

# source('setup_ON.R')
source('setup.R') # bc model 

 # H=0.25*H
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

# ----
test1 = run_over_scen_3(1.1, 2.5, trytime = 120,T2=720)
test2 = run_over_scen_3(1.1, 2.5, trytime = 180,T2=720)
tmp = compare_sims(test1, test2,textsize = 10, 
                   name1 = "July reopening", name2 = "Sept. reopening",scale_y = TRUE)
ggarrange(plotlist = tmp[1:2], nrow=1)
ggsave("reopen-july-sept-slowvax.pdf", height = 4, width = 12)
# ----

phac.accept = c(72,72, 75.18, 78.25, 84.7, 84.7, 84.70)
# 15-19, 20-44, 45-54, 55-64, 65+, 75+, 80+ 
# these will give a hes like this:
H = c(0,0,0.28, 0.28, 0.235, 0.2, 0.153, 0.153, 0.153) 
H=c(H, H[3:8])*N_i # with EW and scaled to ON pop size 

# my best guess for the PHAC rollout 
strategies <- list( list(9, c(8,15), c(7, 10, 11, 12, 13, 14), 6,5,4,3))


# ---- figure 2: compare 2 to 2.5 reopening in september 
test1nk = run_over_scen_3(1.1, 2.5, ve = 0.8, vp=0.75, 
                          trytime = 180,T2=720,speedup = 1)
test2nk = run_over_scen_3(1.1, 2, ve = 0.8, vp=0.75,
                          trytime = 180,T2=720,speedup = 1)
tmp = compare_sims(filter(test1nk, time <=600), filter(test2nk, time<=600),
                   textsize = 10, 
                   name1 = "Reopen to R = 2.5", 
                   name2 = "Reopen to R = 2")
ggarrange(plotlist = tmp[1:3], nrow=1)
ggsave("reopen-sept-rcompare.pdf", height = 4, width = 12)

# ---- figure 3: bar plots illustrating WHO gets infected ---- 
dd=total_cases_origin(filter(test1nk, time >= 180 , time <=600))
yscaler = 1e5/pop_total
p1 = ggplot(dd, aes(x=age_band, y=cases*yscaler, fill=Protection)) + 
    geom_bar(stat="identity",alpha=0.6) +
    ylab("Infections per 100K total pop") +
    theme(axis.title.x = element_blank())
p2 = ggplot(dd, aes(x=age_band, y=hosp*yscaler, fill=Protection)) +
    geom_bar(stat="identity",alpha=0.6)+
    ylab("Hospitalizations per 100K total pop") +
    theme(axis.title.x = element_blank())
ggarrange(p1, p2, nrow=1, common.legend = TRUE, legend = "bottom")
ggsave("end-who-protect1.pdf", width=10, height = 5) 
sum(dd$hosp)*yscaler


# ---- figure 4: again but with vaccinating kids ----
strategies <- list( list(9, c(8,15), c(7, 10, 11, 12, 13, 14), 6,5,4,3, 2))
H[1:2]=c(0,0.7)*N_i[1:2]
test1k = run_over_scen_3(1.1, 2.5, ve = 0.8, vp=0.75,
                         trytime = 180,T2=720,speedup =1)
test2k = run_over_scen_3(1.1, 2, ve = 0.8, vp=0.75, 
                         trytime = 180,T2=720,speedup = 1)
tmp = compare_sims(filter(test1nk, time <=700), filter(test1k, time<=700),
                   textsize = 10, 
                   name1 = "Reopen to R = 2.5", 
                   name2 = "Reopen to R = 2.5, with 10-19")
toprow = ggarrange(plotlist = tmp[1:3], nrow=1)

dd=total_cases_origin(filter(test1k, time >= 180 , time <=700))
yscaler = 1e5/pop_total
p1 = ggplot(dd, aes(x=age_band, y=cases*yscaler, fill=Protection)) + 
    geom_bar(stat="identity",alpha=0.6) +
    ylab("Infections per 100K total pop") +
    theme(axis.title.x = element_blank())
p2 = ggplot(dd, aes(x=age_band, y=hosp*yscaler, fill=Protection)) +
    geom_bar(stat="identity",alpha=0.6)+
    ylab("Hospitalizations per 100K total pop") +
    theme(axis.title.x = element_blank())
botrow = ggarrange(p1, p2, nrow=1, common.legend = TRUE, legend = "bottom")
ggarrange(toprow,botrow,nrow = 2)
sum(dd$hosp)*yscaler
ggsave("reopen-kids-prot.pdf", width = 12,height = 8)
#


# ---- now VOC issues ---- 
strategies <- list( list(9, c(8,15), c(7, 10, 11, 12, 13, 14), 6,5,4,3, 2))
H[1:2]=c(0,0.7)*N_i[1:2]
test1voc = run_over_scen_3(1.1, 2.5,ve = 0.5,vp=0.75, trytime = 180,T2=720,speedup = 2)
test2voc = run_over_scen_3(1.1, 2,ve = 0.5,vp=0.75,  trytime = 180,T2=720,speedup = 2)
tmp = compare_sims(filter(test1k, time <=600), filter(test1voc, time<=600),
                   textsize = 10, 
                   name1 = "Reopen to R = 2.5", 
                   name2 = "Reopen to R = 2.5 w VOC")
toprow = ggarrange(plotlist = tmp[1:3], nrow=1)

dd=total_cases_origin(filter(test1voc, time >= 180 , time <=600))
yscaler = 1e5/pop_total
p1 = ggplot(dd, aes(x=age_band, y=cases*yscaler, fill=Protection)) + 
    geom_bar(stat="identity",alpha=0.6) +
    ylab("Infections per 100K total pop") +
    theme(axis.title.x = element_blank())
p2 = ggplot(dd, aes(x=age_band, y=hosp*yscaler, fill=Protection)) +
    geom_bar(stat="identity",alpha=0.6)+
    ylab("Hospitalizations per 100K total pop") +
    theme(axis.title.x = element_blank())
botrow = ggarrange(p1, p2, nrow=1, common.legend = TRUE, legend = "bottom")
ggarrange(toprow,botrow,nrow = 2)
ggsave("reopen-kids-voc.pdf", width = 12,height = 8)



# try ramp
test = run_over_scen_ramp(1.3, 2.5, trytime = 60)
tmp = compare_sims(test1, test,textsize = 10, 
                   name1 = "july no ramp", name2 = "ramp")
ggarrange(plotlist = tmp[1:2], nrow=1)





# define 3-stage function 
run_over_scen_3 = function(R1, R2, trytime,T2=210,
                           ve=0.75, vp=0.9,  scen=1,speedup=1, alpha=0.0){
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
    n <- speedup*sum(age_demo[-9])/210
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
    df$r1 <- R1
    df$r2 <- R2
    df$trytime = trytime
    df$ve <- ve
    df$vp <- vp
    df$type <- labels[[scen]]
    df$scen <- scen
    df$alpha <- alpha
    return(df)
    }

run_over_scen_ramp = function(R1, R2,ramptime=60, trytime,T2=210,
                           ve=0.75, vp=0.9,  scen=1,speedup=1, alpha=0.0){
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
    n <- speedup*sum(age_demo[-9])/210
    C <- construct_C_from_prem(home=mu_home, work=mu_work, school=mu_school, other=mu_other, u=u_var,
                               target_R0=R1, in_school=TRUE, alpha_factor=alpha)
    df1 <- run_sim_restart(C, df_0=tail(df0, n=1), percent_vax =1.0, strategy= strategies[[scen]], num_perday=n,
                           v_e = rep(ve, num_groups), v_p=rep(vp, num_groups),
                           u = u_var, num_days=trytime, with_essential=TRUE, H=H)
    
    # final stage L R moves to R2, rate stays the same 
    C <- construct_C_from_prem(home=mu_home, work=mu_work, school=mu_school, other=mu_other, u=u_var,
                               target_R0=1.0, in_school=TRUE)
    
#     C <- construct_C_from_prem(home=mu_home, work=mu_work, school=mu_school, other=mu_other, u=u_var,
 #                              target_R0=R2, in_school=TRUE, alpha_factor=alpha)
    R_vec <- get_R_vec(R1=R1, R2=R1, start_ramp=10, end_ramp=10+ramptime, ndays=T2-trytime)
    
      df <- run_sim_restart_ramp_R(C=C, df_0=tail(df1, n=1),  R_vec=R_vec, percent_vax =1.0, strategy= strategies[[scen]], num_perday=n,
                          v_e = rep(ve, num_groups), v_p=rep(vp, num_groups),
                          u = u_var, num_days=T2-trytime, with_essential=TRUE, H=H)
    # combine 
    df1$time = df1$time+T1+1
    df$time = df$time+T1+trytime+1
    
    df <- combine_age_groups(rbind(df0,df1,df))
    
    # add pars
    df$r1 <- R1
    df$r2 <- R2
    df$trytime = trytime
    df$ve <- ve
    df$vp <- vp
    df$type <- labels[[scen]]
    df$scen <- scen
    df$alpha <- alpha
    return(df)
}








