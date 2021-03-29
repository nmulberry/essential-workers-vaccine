source('setup.R')


###############
# VE vs VP
################
pars <- crossing(R=c(1.15, 1.3, 1.5), ve = seq(0.2,0.9,by=0.1), vp = seq(0.2,0.9,by=0.1), scen=c(1,2,3), n=0.003)
res <- pars %>%  future_pmap_dfr(run_over_scen, .progress=TRUE)

#---SUMMARIZE & GET OPT STRATS
res2 <- res %>% 
    group_by(type, R, ve,vp) %>%
    nest() %>%
    summarize(deaths=map_dbl(data, total_deaths)) %>%
    group_by(R,ve,vp) %>%
    filter(deaths == min(deaths))

res2$type <- factor(res2$type, levels=labels[1:3])
    
#--- PLOT
ggplot(res2, aes(x=ve,y=vp, fill=type))+
  geom_tile()+
  facet_wrap(~ R, labeller=label_both, ncol=3)+
  labs(x='Efficacy against infection', y='Efficacy against disease', fill='Strategy')+
  scale_fill_brewer(palette = "Dark2", drop=FALSE)

ggsave(paste0(PATH,"figures/heatmap_ve_vp.pdf"))
  



###############
# VE vs n
################
pars <- crossing(R=c(1.15, 1.3, 1.5), ve = seq(0.2,0.9,by=0.1), ve = 0.9, scen=c(1,2,3), n=seq(0.003,0.01, by=0.001))
res <- pars %>%  future_pmap_dfr(run_over_scen, .progress=TRUE)

#---SUMMARIZE & GET OPT STRATS
res2 <- res %>% 
    group_by(type, R, ve,vp,n) %>%
    nest() %>%
    summarize(deaths=map_dbl(data, total_deaths)) %>%
    group_by(R,ve,vp,n) %>%
    filter(deaths == min(deaths))

res2$type <- factor(res2$type, levels=labels[1:3])
    
#--- PLOT
ggplot(res2, aes(x=ve,y=n, fill=type))+
  geom_tile()+
  facet_wrap(~ R, labeller=label_both, ncol=3)+
  labs(x='Efficacy against infection', y='Vaccination rate', fill='Strategy')+
  scale_fill_brewer(palette = "Dark2", drop=FALSE)

ggsave(paste0(PATH,"figures/heatmap_ve_n.pdf"))
  
