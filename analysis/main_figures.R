source('setup.R')

# parameter space
pars <- crossing(R=c(1.15, 1.3, 1.5), ve = c(0.6,0.75,0.9), vp = 0.9, scen=c(1,2,3,4,5))

# RUN (according to piecewise scenario defined in run_over_scen_2)
res <- pars %>%  future_pmap_dfr(run_over_scen_2, .progress=TRUE)

############################
# FIGURE 1 (trajectories)
#############################
# Look at trajectories
trajA <- compare_sims(sim1 = filter(res, R==1.3 & scen==1 & ve==0.75), 
                             sim2=filter(res, R==1.3 & scen ==2 & ve==0.75),
                             name1=labels[1], name2=labels[2], startDate=startDate, 
                             textsize = 16)

trajB <- compare_sims(sim1 = filter(res, R==1.3 & scen==3 & ve==0.75), 
                             sim2=filter(res, R==1.3 & scen==4 & ve==0.75),
                             name1=labels[3], name2=labels[4], startDate=startDate, 
                             textsize = 16)

# Look at number vaccinated
gg_vax <- res %>% 
    filter(R==1.15) %>% # note: all R values will give approx. the same plot
    group_by(type)%>% 
    nest()%>%
    summarize(plot=map(data, display_prop_vax, startDate, type, textsize=16))


fig1a = ggarrange(gg_vax$plot[[1]]+ggtitle("A: oldest first"),
                  gg_vax$plot[[2]]+ggtitle("B: 80+, 20-79"),
                  gg_vax$plot[[3]]+ggtitle("C: 80+, EW, 70-79,..."),
                  gg_vax$plot[[4]]+ggtitle("D: 80+, EW, 20-79"),
                  ncol=4, nrow=1, common.legend=TRUE, legend="bottom")


# Arrange
fig1b = ggarrange(ggarrange(plotlist=trajA, nrow=1, ncol=4, widths = c(1,1,1,1)),
           ggarrange(plotlist=trajB, nrow=1, ncol=4, widths = c(1,1,1,1)),
          nrow=2)
 
ggarrange(fig1a, fig1b, nrow=2,heights = c(1, 1.6))
ggsave("figures/fig-trajectories.pdf", width = 15, height = 10)      


##########################
# FIG 2 (bar plots)
###########################
tsize <- 16

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
  theme(text=element_text(size=16))+
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
            theme(text=element_text(size=16))+theme(panel.spacing.x=unit(1.5, "lines") ,
                          axis.text.x = element_text(angle = 35,hjust = 1))+
  scale_x_continuous(breaks=c(0.6,0.75,0.9))

g3 <- ggplot(filter(res2, R %in% R_vec), aes(x=ve, y=deaths, group=type, fill=type))+
  geom_col(position='dodge',alpha=0.8)+ 
  scale_fill_brewer(palette = "Dark2")+ 
  facet_wrap(~ R,  ncol=4, ,labeller=label_both)+
  scale_color_brewer(palette = "Dark2")+
  theme(text=element_text(size=16))+
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
  theme(text=element_text(size=16))+
  theme(panel.spacing.x=unit(1.5, "lines"))+
  labs(x='Efficacy against Infection', y='Long COVID', fill='Strategy') +
            theme(text=element_text(size=tsize))+theme(panel.spacing.x=unit(1.5, "lines") ,
                          axis.text.x = element_text(angle = 35,hjust = 1))+
  scale_x_continuous(breaks=c(0.6,0.75,0.9))

bars <- ggarrange(g1,g2,g3,g4, ncol=2, nrow=2, common.legend=TRUE, legend='bottom', align='v')

ggsave('figures/fig-barplots.pdf', width=14, height=10)

#####################
# QALYs and Cost
#####################

data = readxl::read_xlsx("~/essential-workers-vaccine/data/tupper_QALY_dollars.xlsx")

new_yll <- data$`QALYs lost to death` # goes into q1
caseqalyfactor <- data$`QALD lost per community case` # into q2
hospqfactor <- data$`burden of disease hospital stay` # into q3
chronicq <- data$`QALYs lost from chronic (max 25 years)` # into q4; supp use other col

qalys_res <- res %>% 
    group_by(type, R, ve) %>%
    nest() %>%
    summarize(Deaths=map_dbl(data, qaly1, new_yll),
            Cases=map_dbl(data, qaly2,SDUR_vec[1:9],caseqalyfactor),
            Hospitalization=map_dbl(data, qaly3,IHR[1:9],HOSPDUR[1:9],hospqfactor),
            Chronic = map_dbl(data,qaly4,0.015, YLL_vec[1:9],IHR[1:9],chronicq)) %>%
    pivot_longer(cols=4:7, names_to="Source",
                            values_to ="qaly")

qalys_res$Source <- factor(qalys_res$Source, 
    levels=c("Deaths", "Hospitalization","Chronic","Cases"))

ggplot(data=qalys_res, aes(x=type,y=qaly,fill=Source))+
    geom_bar(stat="identity",alpha = 0.8)+
    facet_grid(R~ve,scales= "free_y", labeller=label_both)+
    theme(panel.spacing.x=unit(1.5, "lines"),
          axis.text.x = element_text(angle = 35,hjust = 1),
          axis.title.x = element_blank(),text=element_text(size=tsize))+
    ylab("QALY lost")+
    scale_fill_brewer(palette = "Paired") 

ggsave("figures/qalybars.pdf",height=7, width=10)

# COST
chroniccost <- data$`cost of chronic condition max 25 years`
vp <- rep(unique(res$vp), 9)
cost_res <- res %>% 
    group_by(type, R, ve) %>%
    nest() %>%
    summarize(qq=map_dbl(data, cost1, new_yll),
              Hospitalization=map_dbl(data, cost2,IHR[1:9],vp),
              Chronic = map_dbl(data,cost3, IHR[1:9],chroniccost)) %>% 
    rename(QALY = qq) %>%
    pivot_longer(cols=4:6, names_to="Source",
                           values_to ="Cost") %>% 
    group_by(R, ve, Source) %>%
    mutate(Savings = max(Cost)-Cost)%>% #strategy A is always the most expensive here
    filter(type != labels[1]) # get rid of Strategy A

cost_res$Source <- factor(cost_res$Source, 
    levels=c("QALY", "Hospitalization","Chronic"))

ggplot(data=cost_res, aes(x=type,y=Cost/1e6,fill=Source))+
            geom_bar(stat="identity",alpha = 0.8)+
            facet_grid(R~ve,scales= "free_y", labeller=label_both)+
            theme(panel.spacing.x=unit(1.5, "lines") ,
                  axis.text.x = element_text(angle = 35,hjust = 1),
                  axis.title.x = element_blank(),text=element_text(size=tsize))+
    ylab("Cost (millions)")+
    scale_fill_brewer(palette = "Paired") 

ggsave("figures/costbars.pdf",height=7, width=10)

ggplot(data=cost_res, aes(x=type,y=Savings/1e6,fill=Source))+
            geom_bar(stat="identity",alpha = 0.8)+
            facet_grid(R~ve,scales= "free_y", labeller=label_both)+
            theme(panel.spacing.x=unit(1.5, "lines") ,
                  axis.text.x = element_text(angle = 35,hjust = 1),
                  axis.title.x = element_blank(),text=element_text(size=tsize))+
    ylab("Savings (millions)")+
    scale_fill_brewer(palette = "Paired") 


ggsave("savingsbars.pdf",height=7, width=10)



