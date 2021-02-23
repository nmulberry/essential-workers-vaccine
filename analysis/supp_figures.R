source('setup.R')


###########################
# SENSITIVITY W.R.T. VE
############################
pars <- crossing(R=c(1.15, 1.3, 1.5), ve = seq(0.2,0.9,by=0.1), vp = 0.9, scen=c(1,2,3,4,5))
res <- pars %>%  future_pmap_dfr(run_over_scen_2, .progress=TRUE)
#---SUMMARIZE
res2 <- res %>% 
    group_by(type, R, ve,vp,alpha) %>%
    nest() %>%
    summarize(cases=map_dbl(data, total_cases),
              hosp=map_dbl(data, total_hosp,hosp_efficacy=rep(vp, 9)),
              deaths=map_dbl(data, total_deaths),
              long = map_dbl(data, total_long, hosp_efficacy=rep(vp, 9))
)

g1 <- ggplot(res2, 
  aes(x=ve, group=type, col=type))+
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

g2 <- ggplot(res2, 
  aes(x=ve, group=type, col=type))+
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

g3 <- ggplot(res2, 
  aes(x=ve, group=type, col=type))+
  geom_line(aes(y=deaths), size=1.5)+
  geom_point(aes(y=deaths), size=2)+
  facet_wrap(~ve, labeller=label_both, scales='free')+ 
  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+
  theme(text=element_text(size=16))+
  theme(panel.spacing.x=unit(1.5, "lines"))+
  labs(x='Vaccine Efficacy Against Infection', y='Deaths', col='Strategy') +
            theme(text=element_text(size=16))+theme(panel.spacing.x=unit(1.5, "lines") ,
                          axis.text.x = element_text(angle = 35,hjust = 1))

ggarrange(g1,g2,g3, ncol=1, common.legend=TRUE, legend='bottom', align='v', heights = c(1,1,1.5))
ggsave("figures/sensitivity-ve.pdf", width=14,height=10)






###########################
# SENSITIVITY W.R.T. ALPHA
############################
pars <- crossing(R=c(1.15, 1.3, 1.5), ve = 0.75, vp = 0.9, scen=c(1,2,3,4,5), alpha=seq(0.0,1.0, by=0.25))
res <- pars %>%  future_pmap_dfr(run_over_scen_2, .progress=TRUE)


#---SUMMARIZE
res2 <- res %>% 
    group_by(type, R, ve,vp,alpha) %>%
    nest() %>%
    summarize(cases=map_dbl(data, total_cases),
              hosp=map_dbl(data, total_hosp,hosp_efficacy=rep(vp, 9)),
              deaths=map_dbl(data, total_deaths),
              long = map_dbl(data, total_long, hosp_efficacy=rep(vp, 9))
)

g1 <- ggplot(res2, 
  aes(x=alpha, group=type, col=type))+
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

g2 <- ggplot(res2, 
  aes(x=alpha, group=type, col=type))+
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

g3 <- ggplot(res2, 
  aes(x=alpha, group=type, col=type))+
  geom_line(aes(y=deaths), size=1.5)+
  geom_point(aes(y=deaths), size=2)+
  facet_wrap(~ve, labeller=label_both, scales='free')+ 
  scale_fill_brewer(palette = "Dark2")+ scale_color_brewer(palette = "Dark2")+
  theme(text=element_text(size=16))+
  theme(panel.spacing.x=unit(1.5, "lines"))+
  labs(x='Proportion Nonessential Work Contact', y='Deaths', col='Strategy') +
            theme(text=element_text(size=16))+theme(panel.spacing.x=unit(1.5, "lines") ,
                          axis.text.x = element_text(angle = 35,hjust = 1))

ggarrange(g1,g2,g3, ncol=1, common.legend=TRUE, legend='bottom', align='v', heights = c(1,1,1.5))
ggsave("figures/sensitivity-alpha.pdf", width=14,height=10)




