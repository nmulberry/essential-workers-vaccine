source('setup.R')
#-----------------------------------------------#

paramss <- crossing(R=c(1.15,1.3,1.5), alpha=0.0, ve=c(0.6,0.75,0.8), vp=c(0.9), scen=c(1,2,3,4,5))

# original results
mu_home <- readRDS("mu_home_essential.rds")
mu_work <- readRDS("mu_work_essential.rds")
mu_school <- readRDS("mu_school_essential.rds")
mu_other <-  readRDS("mu_other_essential.rds")
res0 <- params %>%  future_pmap_dfr(run_over_scen_2, .progress=TRUE)

res0 <-  res0 %>%
    group_by(type, R0, alpha, ve,vp) %>%
    nest() %>%
    summarize(cases=map_dbl(data, total_cases),
              hosp=map_dbl(data, total_hosp,hosp_efficacy=vp),
              deaths=map_dbl(data, total_deaths),
              long = map_dbl(data, total_long, hosp_efficacy=vp))%>%
    mutate(prop_e='original')

# Half the number of essential workers
generate_contact_matrix(0.5*p_ess, age_demo_by_five)
mu_home <- readRDS("mu_home_essential.rds")
mu_work <- readRDS("mu_work_essential.rds")
mu_school <- readRDS("mu_school_essential.rds")
mu_other <-  readRDS("mu_other_essential.rds")
res_half <- params %>%  future_pmap_dfr(run_over_scen_2, .progress=TRUE)

res_half <-  res_half %>%
    group_by(type, R0, alpha, ve,vp) %>%
    nest() %>%
    summarize(cases=map_dbl(data, total_cases),
              hosp=map_dbl(data, total_hosp,hosp_efficacy=vp),
              deaths=map_dbl(data, total_deaths),
              long = map_dbl(data, total_long, hosp_efficacy=vp))%>%
    mutate(prop_e='half')

# Double the number of essential workers
generate_contact_matrix(2*p_ess, age_demo_by_five)
mu_home <- readRDS("mu_home_essential.rds")
mu_work <- readRDS("mu_work_essential.rds")
mu_school <- readRDS("mu_school_essential.rds")
mu_other <-  readRDS("mu_other_essential.rds")
res_double <- params %>%  future_pmap_dfr(run_over_scen_2, .progress=TRUE)


res_double<-  res_double %>%
    group_by(type, R0, alpha, ve,vp) %>%
    nest() %>%
    summarize(cases=map_dbl(data, total_cases),
              hosp=map_dbl(data, total_hosp,hosp_efficacy=vp),
              deaths=map_dbl(data, total_deaths),
              long = map_dbl(data, total_long, hosp_efficacy=vp))%>%
    mutate(prop_e='double')


res_tot <- rbind(rbind(res_half, res_double), res0)

res_tot$prop_e <- factor(res_tot$prop_e , levels=c("half", "original", "double"))



g1 <- ggplot(filter(res_tot, R0==1.15), 
  aes(x=prop_e, group=type, col=type))+
  geom_line(aes(y=cases), size=1.5)+
  geom_point(aes(y=cases), size=2)+
  facet_wrap(~ve, labeller=label_both, scales='free')+ 
  scale_fill_brewer(palette = "Dark2")+ 
  scale_color_brewer(palette = "Dark2")+
  theme(text=element_text(size=16))+
  theme(panel.spacing.x=unit(1.5, "lines"))+
  labs(x='', y='Infections', col='Strategy') +
            theme(axis.title.x = element_blank(),text=element_text(size=16))+
            theme(panel.spacing.x=unit(1.5, "lines") ,
                          axis.text.x = element_blank())

g2 <- ggplot(filter(res_tot, R0==1.15), 
  aes(x=prop_e, group=type, col=type))+
  geom_line(aes(y=hosp), size=1.5)+
  geom_point(aes(y=hosp), size=2)+
  facet_wrap(~ve, labeller=label_both, scales='free')+ 
  scale_fill_brewer(palette = "Dark2")+ 
  scale_color_brewer(palette = "Dark2")+
  theme(text=element_text(size=16))+
  labs(x='', y='Hospitalizations', col='Strategy') +
  theme(panel.spacing.x=unit(1.5, "lines"))+
            theme(axis.title.x = element_blank(),text=element_text(size=16))+theme(panel.spacing.x=unit(1.5, "lines") ,
                          axis.text.x = element_blank())

g3 <- ggplot(filter(res_tot, R0==1.15), 
  aes(x=prop_e, group=type, col=type))+
  geom_line(aes(y=deaths), size=1.5)+
  geom_point(aes(y=deaths), size=2)+
  facet_wrap(~ve, labeller=label_both, scales='free')+ 
  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+
  theme(text=element_text(size=16))+
  theme(panel.spacing.x=unit(1.5, "lines"))+
  labs(x='Amount of Essential Workers', y='Deaths', col='Strategy') +
            theme(text=element_text(size=16))+theme(panel.spacing.x=unit(1.5, "lines") ,
                          axis.text.x = element_text(angle = 35,hjust = 1))

ggarrange(g1,g2,g3, ncol=1, common.legend=TRUE, legend='bottom', align='v', heights = c(1,1,1.5))
ggsave("figures/sensitivity-prope.pdf", width=14,height=10)


###############################
# RANDOM RESMAPLING
###############################
perturb_C <- function(C){
  new_C <- C
  # create a symmetric pertrubabtion matrix
  alpha <- matrix(rnorm(nrow(C)*nrow(C), mean=0, sd=0.3), ncol=nrow(C))
  alpha <- 0.5*(alpha+t(alpha))
  alpha[which(alpha < -1)] = -1
  alpha[which(alpha > 1)] = 1
  new_C <- C + alpha*C

  return(new_C)
}




pars <- crossing(R=c(1.15,1.3,1.5), alpha=0.0, ve = 0.75, vp = 0.9, scen=c(1,2,3,4,5), iter=1:50)   

run_over_C_2 = function(R, alpha,ve, vp, scen,iter){
   T1 <- 60
   T2 <- 210
   # Initial stage (vax all 80+)
   R_init <- 1.05
   n <- age_demo[9]/T1
   C <- construct_C_from_prem(home=mu_home, work=mu_work, school=mu_school, other=mu_other, u=u_var,
                              target_R0=10, in_school=TRUE, alpha_factor=alpha)
   C <- perturb_C(C)
   C <- C*R_init/compute_R0(u_var,C)

   df0 <- run_sim_basic(C, I_0=I_0, percent_vax =1.0, strategy=list(9), num_perday=n,
                        v_e_type="aorn", v_e = rep(ve, num_groups),
                        u = u_var, num_days=T1, v_p=rep(vp, num_groups), with_essential=TRUE, H=H) 
   # Final stage
   n <- sum(age_demo[-9])/T2
   C <- C*R/compute_R0(u_var,C)
   df <- run_sim_restart(C, df_0=tail(df0, n=1), percent_vax =1.0, strategy= strategies[[scen]], num_perday=n,
                         v_e_type="aorn", v_e = rep(ve, num_groups),
                         u = u_var, num_days=T2, v_p=rep(vp, num_groups), with_essential=TRUE, H=H)
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
   df$iter <- iter
   return(df)}


res <- pars %>%  
  future_pmap_dfr(run_over_C_2, .progress=TRUE,.options = furrr_options(seed = TRUE))

res <- res %>%
    group_by(type, R, alpha, ve,vp,iter) %>%
    nest() %>%
    summarize(cases=map_dbl(data, total_cases),
              hosp=map_dbl(data, total_hosp,hosp_efficacy=vp),
              deaths=map_dbl(data, total_deaths),
              long = map_dbl(data, total_long, hosp_efficacy=vp))%>%
    group_by(type, R, alpha, ve,vp) %>%
    summarize(cases_50=quantile(cases, probs=0.5), deaths_50=quantile(deaths, probs=0.5), hosp_50=quantile(hosp, probs=0.5),
      cases_5=quantile(cases, probs=0.05), deaths_5=quantile(deaths, probs=0.05), hosp_5=quantile(hosp, probs=0.05),
      cases_95=quantile(cases, probs=0.95), deaths_95=quantile(deaths, probs=0.95), hosp_95=quantile(hosp, probs=0.95))


g1 <- ggplot(res, aes(x=type, fill=type, y=cases_50))+
  geom_col()+
  geom_errorbar(aes(ymin=cases_5, ymax=cases_95),  col='black', width=.2)+ 
  scale_fill_brewer(palette = "Dark2")+ 
  facet_wrap(~ R, ncol=3, scales="free",labeller=label_both)+
  scale_color_brewer(palette = "Dark2")+
  theme(text=element_text(size=16))+
  theme(panel.spacing.x=unit(1.5, "lines"))+
  labs(x='', y='Infections', fill='Strategy') +
            theme(text=element_text(size=16))+theme(panel.spacing.x=unit(1.5, "lines") ,
                          axis.text.x = element_blank())

g2 <- ggplot(res, aes(x=type, fill=type, y=hosp_50))+
  geom_col()+
  geom_errorbar(aes(ymin=hosp_5, ymax=hosp_95),  col='black', width=.2)+ 
  scale_fill_brewer(palette = "Dark2")+ 
  scale_color_brewer(palette = "Dark2")+
  facet_wrap(~ R,  ncol=3, scales="free", labeller=label_both)+
  theme(text=element_text(size=16))+
  theme(panel.spacing.x=unit(1.5, "lines"))+
  labs(x='', y='Hospitalizations', col='Strategy') +
            theme(text=element_text(size=16))+theme(panel.spacing.x=unit(1.5, "lines") ,
                          axis.text.x = element_blank())


g3 <- ggplot(res, aes(x=type, fill=type,  y=deaths_50))+
  geom_col()+
  geom_errorbar(aes(ymin=deaths_5, ymax=deaths_95),  col='black', width=.2)+ 
  scale_fill_brewer(palette = "Dark2")+ 
  facet_wrap(~ R,  ncol=3, scales="free",labeller=label_both)+
  scale_color_brewer(palette = "Dark2")+
  theme(text=element_text(size=16))+
  theme(panel.spacing.x=unit(1.5, "lines"))+
  labs(x='', y='Deaths', fill='Strategy') +
            theme(text=element_text(size=16))+theme(panel.spacing.x=unit(1.5, "lines") ,
                          axis.text.x = element_text(angle = 35,hjust = 1))

ggarrange(g1,g2,g3, ncol=1, align='v', common.legend=TRUE, legend='none', heights=c(1,1,1.5))
ggsave("figures/sensitivity-random-C.pdf")


