# ontario FSA stuff 
source('setup_ON.R')

# data and files - quick guide
# AgeFSA_forBrief -- FSA names, risk group from 1 to 10, number of 80+ in the FSA
# note M5H is like 3 blocks of downtown TO and has 27 80+s whereas L4C is a 
# large area of Richmond Hill with 4000+ of them 
# justDeciles is a reduced version of this 

# Copy of COVID19 ... 
# FSA, PHU, population, .. lasti s portion EW, 15 cols before that are the age_demo vector
# this one gives us population, and age_demo

# the AgeFSA_forBrief one gives us risk group

# steps will be

# (1) set up a new model for each FSA including its age_demo vector, its population size, 
# its initial condition I0, and its overall vaccination rate 

# generate_contact_matrix has to have a new age_demo_by_fives that has the total population 
# for each FSA and the p_ess for each FSA in its argument 


# (1a) test the setup 

# (2) each FSA i guess gets a strategy, depending on how low we go down the age rollout

# (3) run the model for each FSA

# (4) add the total numbers of infections, deaths etc from each run 

# to set this up I will need a function like run_over_scen which takes in the 
# FSA and runs the whole above setup, returning a big df with FSA in it 

# then use summarize() to add the FSAs up 

# ---- testing ----
fsa1=readxl::read_excel("~/Dropbox/COVID/ON-vaccine/AgeFSA_forBrief.xlsx")
fsa2=readxl::read_excel("~/Dropbox/COVID/ON-vaccine/Copy of COVID19 Modeling - occupation by FSA - 2021-03-04 SHARED (COVID data removed).xlsx")




# change each of these for each fsa 
# p_ess <- c(0.0,0.0, 0.4,0.35,0.35,0.35,0.25,0.1,0.0) # est.
get_p_ess <- function(fsa_name, fsa_data) {
    ii = match(fsa_name, fsa_data$FSA) 
    myvec = as.numeric(fsa_data[ii, 16:30]) # 1:9 portion non-EW. 10:15 portion EW of 20s - 70s .
    numew = c(0,0, myvec[10:15],0)
    return( numew/(numew+myvec[1:9]))
} 

run_over_scen_fsa = function(R,fsa_name,r, scen,alpha=0.0, use_fsa= TRUE){
    T1 <- 60
    T2 <- 60
    T3 <- 140
#    r=0.4 # percent per day
    # set up FSA specific stuff 
    age_demo_by_five <- as.matrix(readr::read_csv('~/essential-workers-vaccine/data/Population_Estimates_ON.csv'))
    this_p_ess = get_p_ess(fsa_name, fsa2)
    this_pop_fsa = fsa2$Pop[match(fsa_name, fsa2$FSA)]
    newmats = generate_contact_matrix_fsa(this_p_ess,this_pop_fsa, age_demo_by_five) #build
    age_demo <- newmats$age_demo
    mu_home = newmats$mu_home
    mu_school = newmats$mu_school
    mu_other = newmats$mu_other
    mu_work = newmats$mu_work
    
    # GLOBAL MODEL PARAMETERS
    n <- length(age_demo)
    pop_total <- age_demo[n] #pop of the FSA, now
    age_demo <- age_demo[1:n-1]
    N_i <-  pop_total*age_demo
    num_groups <- length(age_demo)
    
    # I0 comes from setup - here will have to scale I0 by pop size of the FSA and by risk
    # of the FSA
    risk=fsa1$`Neighbourhood Risk Group`[match(fsa_name, fsa1$FSA)]
   scalefactor =1
   if (risk==1 & use_fsa) {scalefactor =3}
    if (risk==2 & use_fsa) {scalefactor = 2} 
    I_0_fsa = scalefactor*I_0*pop_total/sum(age_demo_by_five)
    
    H = c(0.0,0.0,0.2,0.2,0.2,0.2,0.15,0.15,0.15,0.2,0.2,0.2,0.2,0.15,0.15)*N_i # Hesitancy
    # could update according to PHAC stuff --> moved to 0.2 not 0.3 in 20s 
    
    # Initial stage (vax all 80+)
    R_init <- 1.05
    n <- age_demo[9]/T1
    C <- construct_C_from_prem(home=mu_home, work=mu_work, school=mu_school, other=mu_other, u=u_var,
                               target_R0=R_init, in_school=TRUE, alpha_factor=alpha)
    
    df0 <- run_sim_basic(C, I_0=I_0, percent_vax =1.0, strategy=list(9), num_perday=n,
                         v_e = rep(ve, num_groups), v_p=rep(vp, num_groups),
                         u = u_var, num_days=T1, with_essential=TRUE, H=H) 
    # second stage
# NOTE NEED TO FIX THESE IF I WANT TO USE THIS - some NA where logical needed blah
    if (use_fsa & (risk ==1 | risk ==2)) { n = 3*r } else {n =r}
#    n <- sum(age_demo[-9])/T2
    C <- construct_C_from_prem(home=mu_home, work=mu_work, school=mu_school, other=mu_other, u=u_var,
                               target_R0=R, in_school=TRUE, alpha_factor=alpha)
    df1 <- run_sim_restart(C, df_0=tail(df0, n=1), percent_vax =1.0, strategy= strategies[[scen]], num_perday=n,
                          v_e = rep(ve, num_groups), v_p=rep(vp, num_groups),
                          u = u_var, num_days=T2, with_essential=TRUE, H=H)
    # third stage
    if (use_fsa & (risk ==1 | risk ==2)) {
        n = (pop_total - 3*r*pop_total*T1)/T2 
    } else {
            n =( pop_total-r*pop_total*T1)/T2 }
    #    n <- sum(age_demo[-9])/T2
    C <- construct_C_from_prem(home=mu_home, work=mu_work, school=mu_school, other=mu_other, u=u_var,
                               target_R0=R, in_school=TRUE, alpha_factor=alpha)
    df <- run_sim_restart(C, df_0=tail(df1, n=1), percent_vax =1.0, strategy= strategies[[scen]], num_perday=n,
                          v_e = rep(ve, num_groups), v_p=rep(vp, num_groups),
                          u = u_var, num_days=T3, with_essential=TRUE, H=H)
    # combine 
    df1$time = df1$time+T1+1
    df$time = df$time+T1+T2+1

    df <- combine_age_groups(rbind(df0,df1,df))
    # add pars
    df$R <- R
    df$fsa = fsa_name
    df$r = r
  #  df$ve <- ve
  #  df$vp <- vp
    df$type <- labels[[scen]]
    df$scen <- scen
    df$alpha <- alpha
    return(df)
    }



pars <- crossing(R=1.3,fsa_name = fsa2$FSA,r=0.004,  scen=1) 
res <- pars %>%  future_pmap_dfr(run_over_scen_fsa, .progress=TRUE)






