source('setup.R')
# Define strategies
labels <- c('A: Oldest to Youngest', 
            'B: 80+, 20-79',
            'C: 80+, EW, 70-79,...', 
            'D: 80+, EW, 20-79',
            'E: 80+, EW, 70-79,20-69')

strategies <- list(list(9,c(8,15),c(7,14), c(6,13), c(5,12), c(4,11), c(3,10)), 
                   list(9,3:15),
                   list(9, 10:15, 8, 7,6,5,4,3),
                   list(9, 10:15, 3:8),
                   list(9, 10:15, 8, 3:7))

# parameter space
pars <- crossing(R=c(1.05,1.15,1.5,1.8,2.0,2.2), ve = c(0.6,0.75,0.9), vp = 0.9, scen=c(1,2,3,4,5))

# RUN (according to piecewise scenario defined in run_over_scen_2)
res <- pars %>%  future_pmap_dfr(run_over_scen_2, .progress=TRUE)

#---SUMMARIZE
res2 <- res %>% 
    group_by(type, R, ve,vp) %>%
    nest() %>%
    summarize(t_turn = map_dbl(data, time_to_decr),
              cases_turn = map_dbl(data,cases_immunity),
              vax_turn = map_dbl(data, vax_immunity))

#---PLOT
tsize <- 16

g1 <- ggplot(res2, aes(x=R, y=cases_turn/pop_total, group=type, col=type))+
  geom_line(size=1.5)+
  geom_point(size=2)+
  facet_wrap(~ve, labeller=label_both,ncol=3)+ 
  scale_color_brewer(palette = "Dark2")+
  theme(text=element_text(size=tsize))+
  theme(panel.spacing.x=unit(1.5, "lines"))+
  labs(x='R',y='Proportion of Population',col='Strategy')
  
g2 <- ggplot(res2, aes(x=R, y=(cases_turn+vax_turn)/pop_total, group=type, col=type))+
  geom_line(size=1.5)+
  geom_point(size=2)+
  facet_wrap(~ve, labeller=label_both,ncol=3)+ 
  scale_color_brewer(palette = "Dark2")+
  theme(text=element_text(size=tsize))+
  theme(panel.spacing.x=unit(1.5, "lines"))+
  labs(x='R',y='Proportion of Population',col='Strategy')

ggarrange(g1,g2, ncol=1, nrow=2, common.legend=TRUE, legend='bottom', labels=c("A","B"))
ggsave("figures/immunity.pdf", width=16, height=10)
  
