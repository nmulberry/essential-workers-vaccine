source('setup_ON.R')

# parameter space
pars <- crossing(R=c(1.15, 1.3), ve = c(0.75,0.9), vp = 0.9, scen=c(1,2,3,4))

# RUN (according to piecewise scenario)
res <- pars %>%  future_pmap_dfr(run_over_scen_2, .progress=TRUE)

############################
# FIGURE 1 (trajectories)
#############################
# Look at trajectories
trajA <- compare_sims(sim1 = filter(res, R==1.3 & scen==1 & ve==0.75), 
                             sim2=filter(res, R==1.3 & scen ==4 & ve==0.75),
                             name1=labels[1], name2=labels[4], startDate=startDate, 
                             textsize = 12)

trajB <- compare_sims(sim1 = filter(res, R==1.3 & scen==2 & ve==0.75), 
                             sim2=filter(res, R==1.3 & scen==3 & ve==0.75),
                             name1=labels[2], name2=labels[3], startDate=startDate, 
                             textsize = 12)

# Look at number vaccinated
gg_vax <- res %>% 
    filter(R==1.15) %>% # note: all R values will give approx. the same plot
    group_by(type)%>% 
    nest()%>%
    summarize(plot=map(data, display_prop_vax, startDate, type, textsize=12))


fig1a = ggarrange(gg_vax$plot[[1]]+ggtitle(labels[1]),
                  gg_vax$plot[[2]]+ggtitle(labels[2]),
                  gg_vax$plot[[3]]+ggtitle(labels[3]),
                  gg_vax$plot[[4]]+ggtitle(labels[4]),
                  ncol=4, nrow=1, common.legend=TRUE, legend="bottom")

# set y lims for second trajectory manually 
trajB[[1]]=trajB[[1]]+ylim(c(0,16000))
trajA[[2]]=trajA[[2]]+ylim(c(0,2600))
# trajB[[3]]=trajB[[3]]+ylim(c(0,16000))
trajB[[4]]=trajB[[4]]+ylim(c(0,45000))

# Arrange
fig1b = ggarrange(ggarrange(plotlist=trajA, nrow=1, ncol=4, widths = c(1,1,1,1)),
           ggarrange(plotlist=trajB, nrow=1, ncol=4, widths = c(1,1,1,1)),
          nrow=2)
 
ggarrange(fig1a, fig1b, nrow=2,heights = c(1, 1.6))
ggsave("../figures/fig-trajectories_on1.pdf", width = 18, height = 12)      


##########################
# FIG 2 (bar plots)
###########################
tsize <- 12

#---SUMMARIZE
res2 <- res %>% 
    group_by(type, R, ve,vp) %>%
    nest() %>%
    summarize(cases=map_dbl(data, total_cases),
              hosp=map_dbl(data, total_hosp,hosp_efficacy=rep(vp, 9)),
              deaths=map_dbl(data, total_deaths),
              long = map_dbl(data, total_long, hosp_efficacy=rep(vp, 9)),
              t_turn = map_dbl(data, time_to_decr),
              cases_turn = map_dbl(data,cases_immunity),
              vax_immunity = map_dbl(data, vax_immunity)
)

R_vec <- c(1.15,1.3) # R vals to plot

g1 <- ggplot(filter(res2, R %in% R_vec), aes(x=ve, y=cases, group=type, fill=type))+
  geom_col(position='dodge', alpha=0.8)+ 
  scale_fill_brewer(palette = "Dark2")+ 
  facet_wrap(~ R,  ncol=4,labeller=label_both)+
  scale_color_brewer(palette = "Dark2")+
  theme(text=element_text(size=tsize))+
  theme(panel.spacing.x=unit(1.5, "lines"))+
  labs(x='Efficacy against Infection', y='Infections', fill='Strategy') +
            theme(text=element_text(size=tsize))+
            theme(panel.spacing.x=unit(1.5, "lines"),
                          axis.text.x = element_text(angle = 35,hjust = 1))+
  scale_x_continuous(breaks=c(0.6,0.75,0.9))

g2 <- ggplot(filter(res2, R %in% R_vec), aes(x=ve, y=hosp, group=type, fill=type))+
  geom_col(position='dodge',alpha=0.8)+ 
  scale_fill_brewer(palette = "Dark2")+ 
  facet_wrap(~ R,  ncol=4,labeller=label_both)+
  scale_color_brewer(palette = "Dark2")+
  theme(text=element_text(size=tsize))+
  theme(panel.spacing.x=unit(1.5, "lines"))+
  labs(x='Efficacy against Infection', y='Hospitalizations', fill='Strategy') +
            theme(text=element_text(size=tsize))+theme(panel.spacing.x=unit(1.5, "lines") ,
                          axis.text.x = element_text(angle = 35,hjust = 1))+
  scale_x_continuous(breaks=c(0.6,0.75,0.9))

g3 <- ggplot(filter(res2, R %in% R_vec), aes(x=ve, y=deaths, group=type, fill=type))+
  geom_col(position='dodge',alpha=0.8)+ 
  scale_fill_brewer(palette = "Dark2")+ 
  facet_wrap(~ R,  ncol=4, ,labeller=label_both)+
  scale_color_brewer(palette = "Dark2")+
  theme(text=element_text(size=tsize))+
  theme(panel.spacing.x=unit(1.5, "lines"))+
  labs(x='Efficacy against Infection', y='Deaths', fill='Strategy') +
            theme(text=element_text(size=tsize))+theme(panel.spacing.x=unit(1.5, "lines") ,
                          axis.text.x = element_text(angle = 35,hjust = 1))+
  scale_x_continuous(breaks=c(0.6,0.75,0.9))

g4 <- ggplot(filter(res2, R %in% R_vec), aes(x=ve, y=long, group=type, fill=type))+
  geom_col(position='dodge',alpha=0.8)+ 
  scale_fill_brewer(palette = "Dark2")+ 
  facet_wrap(~ R,  ncol=4,,labeller=label_both)+
  scale_color_brewer(palette = "Dark2")+
  theme(text=element_text(size=tsize))+
  theme(panel.spacing.x=unit(1.5, "lines"))+
  labs(x='Efficacy against Infection', y='Long COVID', fill='Strategy') +
            theme(text=element_text(size=tsize))+theme(panel.spacing.x=unit(1.5, "lines") ,
                          axis.text.x = element_text(angle = 35,hjust = 1))+
  scale_x_continuous(breaks=c(0.6,0.75,0.9))

bars <- ggarrange(g1,g2,g3,g4, ncol=2, nrow=2, common.legend=TRUE, legend='bottom', align='v')

ggsave('../figures/fig-barplots-ON.pdf', width=14, height=10)

#####################
# QALYs and Cost
#####################
source('qalys_cost.R')





#####################
# explore a higher rate of vaccination 
#####################
# parameter space
pars <- crossing(R=c(1.15, 1.3), ve =0.75, ratefactor=c(1,1.5), scen=c(1,2,3,4))

# RUN (according to piecewise scenario)
res <- pars %>%  future_pmap_dfr(run_over_scen_rate, .progress=TRUE)
res$ratefactor = as.factor(res$ratefactor)
res2 <- res %>% 
  group_by(type, R, ratefactor,vp) %>%
  nest() %>%
  summarize(cases=map_dbl(data, total_cases),
            hosp=map_dbl(data, total_hosp,hosp_efficacy=rep(vp, 9)),
            deaths=map_dbl(data, total_deaths),
            long = map_dbl(data, total_long, hosp_efficacy=rep(vp, 9)),
            t_turn = map_dbl(data, time_to_decr),
            cases_turn = map_dbl(data,cases_immunity),
            vax_immunity = map_dbl(data, vax_immunity)
  )

g1 <- ggplot(filter(res2, R %in% R_vec), aes(x=ratefactor, y=cases, group=type, fill=type))+
  geom_col(position='dodge', alpha=0.8)+ 
  scale_fill_brewer(palette = "Dark2")+ 
  facet_wrap(~ R,  ncol=4,labeller=label_both)+
  scale_color_brewer(palette = "Dark2")+
  theme(text=element_text(size=tsize))+
  theme(panel.spacing.x=unit(1.5, "lines"))+
  labs(x='Speedup', y='Infections', fill='Strategy') +
  theme(text=element_text(size=tsize))   #+
 # theme(panel.spacing.x=unit(1.5, "lines"),
 #       axis.text.x = element_text(angle = 35,hjust = 1)) 

g2 <- ggplot(filter(res2, R %in% R_vec), aes(x=ratefactor, y=hosp, group=type, fill=type))+
  geom_col(position='dodge',alpha=0.8)+ 
  scale_fill_brewer(palette = "Dark2")+ 
  facet_wrap(~ R,  ncol=4,labeller=label_both)+
  scale_color_brewer(palette = "Dark2")+
  theme(text=element_text(size=tsize))+
  theme(panel.spacing.x=unit(1.5, "lines"))+
  labs(x='Speedup', y='Hospitalizations', fill='Strategy') +
  theme(text=element_text(size=tsize))# +theme(panel.spacing.x=unit(1.5, "lines") ,
   #                                          axis.text.x = element_text(angle = 35,hjust = 1))+
 # scale_x_continuous(breaks=c(0.6,0.75,0.9))

g3 <- ggplot(filter(res2, R %in% R_vec), aes(x=ratefactor, y=deaths, group=type, fill=type))+
  geom_col(position='dodge',alpha=0.8)+ 
  scale_fill_brewer(palette = "Dark2")+ 
  facet_wrap(~ R,  ncol=4, ,labeller=label_both)+
  scale_color_brewer(palette = "Dark2")+
  theme(text=element_text(size=tsize))+
  theme(panel.spacing.x=unit(1.5, "lines"))+
  labs(x='Speedup', y='Deaths', fill='Strategy') +
  theme(text=element_text(size=tsize)) #+theme(panel.spacing.x=unit(1.5, "lines") ,
#                                             axis.text.x = element_text(angle = 35,hjust = 1))+
#  scale_x_continuous(breaks=c(0.6,0.75,0.9))

g4 <- ggplot(filter(res2, R %in% R_vec), aes(x=ratefactor, y=long, group=type, fill=type))+
  geom_col(position='dodge',alpha=0.8)+ 
  scale_fill_brewer(palette = "Dark2")+ 
  facet_wrap(~ R,  ncol=4,,labeller=label_both)+
  scale_color_brewer(palette = "Dark2")+
  theme(text=element_text(size=tsize))+
  theme(panel.spacing.x=unit(1.5, "lines"))+
  labs(x='Speedup', y='Long COVID', fill='Strategy') +
  theme(text=element_text(size=tsize))# +theme(panel.spacing.x=unit(1.5, "lines") ,
          #                                   axis.text.x = element_text(angle = 35,hjust = 1))+
 # scale_x_continuous(breaks=c(0.6,0.75,0.9))

bars <- ggarrange(g1,g2,g3,g4, ncol=2, nrow=2, common.legend=TRUE, legend='bottom', align='v')

ggsave('../figures/fig-barplots-rate-ON.pdf', width=14, height=10)

