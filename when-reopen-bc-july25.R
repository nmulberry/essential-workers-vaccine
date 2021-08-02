# when can we reopen to R = x ? 

# source('setup_ON.R')
setwd("analysis/")
source('setup.R') # bc model 
setwd("../")
library("viridis") 
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

# make something up that I think will kind of get BC's actual vaccination 
# as of may 12, sort of right 
labels <- c('A: Oldest to Youngest', 
            'B: BC approximate') 
# 80, 70, (EW, 50, 60), (20-40) 

strategies <- list(list(9,c(8,15),c(7,14), c(6,13), c(5,12), c(4,11), c(3,10)), 
                   list(9,c(8,15), c(7,14), c(6,13, 5,12,11,10), c(3,4),2))





# ----
test1 = run_over_scen_3(1.1, 2.5, trytime = 120,T2=720)
test2 = run_over_scen_3(1.1, 2.5, trytime = 180,T2=720)
tmp = compare_sims(test1, test2,textsize = 10, 
                   name1 = "July reopening", name2 = "Sept. reopening",scale_y = FALSE)
ggarrange(plotlist = tmp[1:2], nrow=1)
ggsave("reopen-july-sept-slowvax.pdf", height = 4, width = 12)


# ----
ascFrac=0.75
textsize=12
phac.accept = c(72,72, 75.18, 78.25, 84.7, 84.7, 84.70)
# 15-19, 20-44, 45-54, 55-64, 65+, 75+, 80+ 
# these will give a hes like this:
H = c(0,0.3,0.23, 0.23, 0.233, 0.2, 0.153, 0.153, 0.153) 
# july 25, using numbesr from https://covid19tracker.ca/provincevac.html?p=BC, 
# let's just try 
H = c(0, 0.35, 0.3, 0.25, 0.2 ,0.2, 0.13, 0.07, 0.04) 
H=c(H, H[3:8])*N_i # with EW and scaled to ON pop size 
# ve=0.8
# vp=0.75


########## sims 1: best current guesses, 70% ve, overall 0.88
# and CURRENT VACCINATION 
# dat = readRDS("~/BC-dat.rds") 
ve=0.7
vp=0.6

dat = readRDS("~/BC-dat-july.rds") 
H = c(0, 0.35, 0.3, 0.25, 0.2 ,0.2, 0.13, 0.07, 0.04) 
H=c(H, H[3:8])*N_i # with EW and scaled to ON pop size 

test1 = run_over_scen_4(R1 = 1.23, R2=1.0, R3=2.8,ve=0.7, 
                        vp=0.6,trytime = 82, Tfinal=365,scen = 2 )
# o1 = extract_cases_deaths(test1,LCFAC = 1);
p1a = plot_incid_data(test1,headertext="Short term: 80% vacc",textsize = 9) + 
 scale_x_date( date_breaks = "months",date_labels = "%b-%d",
                limits=c(ymd("2021-02-15"), ymd("2021-08-30")))
p1a
p1b=plot_incid_data(test1, headertext="Medium term: 80% vacc",textsize = 9) +
  scale_x_date( date_breaks = "months",date_labels = "%b-%d", 
                limits = c(ymd("2021-01-15"), ymd("2021-11-30")))#  +  xlim(c(ymd("2021-01-15"), ymd("2021-09-15")))

ggarrange(p1a,p1b, common.legend = TRUE, nrow=1,legend = "bottom")

# ggsave("~/incid-reopen-0.6.pdf", width = 10, height = 5)

 display_prop_vax(test1, label = "BC")

 
 ########## sims 2: best current guesses, 70% ve, overall 0.88
 # and OPTIMISTIC  VACCINATION but it takes longer 
 H = c(0,0.20,0.14, 0.14, 0.13, 0.12, 0.1, 0.07, 0.04) 
 H=c(H, H[3:8])*N_i # with EW and scaled to ON pop size 
 
 
  test2 = run_over_scen_4(R1 = 1.23, R2=1.0, R3=2.8,ve=0.7, 
                         vp=0.6,trytime = 82,speedup=0.925, Tfinal=365,scen = 2 )
 # o1 = extract_cases_deaths(test1,LCFAC = 1);
 
 p2a = plot_incid_data(test2,headertext="Short term: 90% vacc",textsize = 9) + 
   scale_x_date( date_breaks = "months",date_labels = "%b-%d",
                 limits=c(ymd("2021-01-15"), ymd("2021-08-30")))
 p2a
 p2b=plot_incid_data(test2, headertext="Medium term: 90% vacc",textsize = 9) +
   scale_x_date( date_breaks = "months",date_labels = "%b-%d", 
                 limits = c(ymd("2021-01-15"), ymd("2021-11-30")))#  +  xlim(c(ymd("2021-01-15"), ymd("2021-09-15")))
 
 ggarrange(p1b,p2b+ylim(c(0,6000)), common.legend = TRUE, nrow=1,legend = "bottom")
 
 # ggsave("~/incid-reopen-0.6.pdf", width = 10, height = 5)
 
 
 ### explore a few comparison plots 
 
 ggarrange(p1a, # + xlim(c(ymd("2021-01-15"), ymd("2021-09-31"))), 
           p2a ,# + xlim(c(ymd("2021-01-15"), ymd("2021-09-31")))+ylim(c(0,6000)),
           common.legend = TRUE, nrow=1, legend = "bottom")
 
 ggarrange(p1b + xlim(c(ymd("2021-01-15"), ymd("2021-10-31"))), 
           p2b + xlim(c(ymd("2021-01-15"), ymd("2021-10-31")))+ylim(c(0,18000)),
                       common.legend = TRUE, nrow=1, legend = "bottom")
 
 
 
 tmp = compare_sims(test1,test2,
                    textsize = 10, 
                    name1 = "80%", 
                    name2 = "90%",scale_y = FALSE)
 
 ggarrange(plotlist = tmp[1:2], nrow=1,common.legend = TRUE, legend="bottom")
 
 
 
# left here from before - all below here is old 
 
######### BC REPORT 2 : 0.7 eff against infection , 3 reopenings May 25 
# same thing with lower efficacy against infection 0.7, given we all have a single dose 
vp=0.5 # NOTE now in as default scaling in extracting the incidence
test1 = run_over_scen_4(R1 = 1.23, R2=0.95, R3=1.7,ve=0.7, vp=0.5, trytime = 45, Tfinal=365,scen = 2 )
test2 = run_over_scen_4(R1 = 1.23, R2=0.95, R3=2,ve=0.7, vp=0.5,trytime = 45, Tfinal=365,scen = 2 )
test3 = run_over_scen_4(R1 = 1.23, R2=0.95, R3=2.3,ve=0.7, vp=0.5,trytime = 45, Tfinal=365,scen = 2 )
# o1 = extract_cases_deaths(test1,LCFAC = 1);
p1 = plot_incid_data(test1,headertext="Acceptance 0.75-0.85, Reopen R=1.7",textsize = 9)
p2 = plot_incid_data(test2,headertext="Acceptance 0.75-0.85, Reopen R=2",textsize = 9)
p3 = plot_incid_data(test3,headertext="Acceptance 0.75-0.85, Reopen R=2.3",textsize = 9)
ggarrange(p1,p2,p3, common.legend = TRUE, nrow=1,legend = "bottom")
ggsave("~/incid-reopen-0.7.pdf", width = 10, height = 5) # USED THIS ONE 



########## BC REPORT 3 : 0.7 eff against infection , 3 reopenings June 15 

test1 = run_over_scen_4(R1 = 1.23, R2=0.95, R3=1.7,ve=0.7, vp=0.5,trytime = 65, Tfinal=365,scen = 2 )
test2 = run_over_scen_4(R1 = 1.23, R2=0.95, R3=2,ve=0.7, vp=0.5,trytime = 65, Tfinal=365,scen = 2 )
test3 = run_over_scen_4(R1 = 1.23, R2=0.95, R3=2.3,ve=0.7, vp=0.5,trytime = 65, Tfinal=365,scen = 2 )
# o1 = extract_cases_deaths(test1,LCFAC = 1);
p1 = plot_incid_data(test1,headertext="Acceptance 0.75-0.85, Reopen R=1.7",textsize = 9)
p2 = plot_incid_data(test2,headertext="Acceptance 0.75-0.85, Reopen R=2",textsize = 9)
p3 = plot_incid_data(test3,headertext="Acceptance 0.75-0.85, Reopen R=2.3",textsize = 9)
ggarrange(p1,p2,p3, common.legend = TRUE, nrow=1,legend = "bottom")
ggsave("~/incid-reopen-june-0.7.pdf", width = 10, height = 5)

display_prop_vax(test1, label = "BC")


######### BC REPORT 4 : 0.7 eff against infection , 3 reopenings July 25 

# same thing with lower efficacy against infection 0.7, given we all have a single dose 

test1 = run_over_scen_4(R1 = 1.23, R2=0.9, R3=1.7,ve=0.7, vp=0.55, trytime = 85, Tfinal=365,scen = 2 )
test2 = run_over_scen_4(R1 = 1.23, R2=0.9, R3=2,ve=0.7, vp=0.55,trytime = 85, Tfinal=365,scen = 2 )
test3 = run_over_scen_4(R1 = 1.23, R2=0.9, R3=2.3,ve=0.7, vp=0.55,trytime = 85, Tfinal=365,scen = 2 )
# o1 = extract_cases_deaths(test1,LCFAC = 1);
p1 = plot_incid_data(test1,headertext="Acceptance 0.75-0.85, Reopen R=1.7",textsize = 9)
p2 = plot_incid_data(test2,headertext="Acceptance 0.75-0.85, Reopen R=2",textsize = 9)
p3 = plot_incid_data(test3,headertext="Acceptance 0.75-0.85, Reopen R=2.3",textsize = 9)
ggarrange(p1,p2,p3, common.legend = TRUE, nrow=1,legend = "bottom")
ggsave("~/incid-reopen-june-0.7.pdf", width = 10, height = 5)



# ---- stuff RG asked for ---- >> begin here but see reopening-script-rg.R 





# ---- figure 2: compare 2 to 2.5 reopening in september 
test1nk = run_over_scen_3(1.1, 2.5, ve = 0.8, vp=0.75, 
                          trytime = 180,T2=720,speedup = 1)
test2nk = run_over_scen_3(1.1, 2, ve = 0.8, vp=0.75,
                          trytime = 180,T2=720,speedup = 1)
tmp = compare_sims(filter(test1nk, time <=600), filter(test2nk, time<=600),
                   textsize = 10, 
                   name1 = "Reopen to R = 2.5", 
                   name2 = "Reopen to R = 2")
ggarrange(plotlist = tmp[1:3], nrow=1,common.legend = TRUE, legend="bottom")
ggsave("reopen-sept-rcompare.pdf", height = 4, width = 12)

singpans  = one_sim_plot(test1nk, name1 = "Reopen in September", scale_y = FALSE) 
ggarrange(plotlist = singpans[1:3], nrow=1)

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
sum(dd$hosp)*yscaler # 504


# ---- figure 4: again but with vaccinating kids ----
strategies <- list( list(9, c(8,15), c(7, 10, 11, 12, 13, 14), 6,5,4,3, 2))
H[1:2]=c(0,0.3)*N_i[1:2]
test1k = run_over_scen_3(1.1, 2.5, ve = 0.8, vp=0.75,
                         trytime = 180,T2=720,speedup =1)
test2k = run_over_scen_3(1.1, 2, ve = 0.8, vp=0.75, 
                         trytime = 180,T2=720,speedup = 1)
tmp = compare_sims(filter(test1nk, time <=700), filter(test1k, time<=700),
                   textsize = 10, 
                   name1 = "Reopen R = 2.5", 
                   name2 = "Reopen R = 2.5, w. 10-19")
toprow = ggarrange(plotlist = tmp[1:3], nrow=1,
                   common.legend = TRUE, legend="bottom")
toprow
ggsave("reopen-kids-protworks.pdf", width = 12,height = 4)

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
H = c(0,0,0.28, 0.28, 0.235, 0.2, 0.153, 0.153, 0.153) 
H=c(H, H[3:8])*N_i # with EW and scaled to ON pop size 

H[1:2]=c(0,0.3)*N_i[1:2]
test1voc = run_over_scen_3(1.1, 2.5,ve = 0.5,vp=0.75, trytime = 180,T2=720,speedup = 2)
test2voc = run_over_scen_3(1.1, 2,ve = 0.5,vp=0.75,  trytime = 180,T2=720,speedup = 2)
tmp = compare_sims(filter(test1k, time <=600), filter(test1voc, time<=600),
                   textsize = 10, 
                   name1 = "Reopen R = 2.5", 
                   name2 = "Reopen R = 2.5 w VOC")
toprow = ggarrange(plotlist = tmp[1:3], nrow=1,
                   common.legend = TRUE, legend="top")
toprow 


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
ggarrange(toprow,botrow,nrow = 2,heights = c(1.3,1))
ggsave("reopen-kids-voc.pdf", width = 12,height = 8)


# we decided to do a version of figure 5 that has 
# better H comparing no VOC to VOC


# keep Paul's current simple calculation but adjust to not include 0-9 
# and mirror what we do in the text. 



# ---- now VOC with better vaccine acceptance, still w kids ---- 
strategies <- list( list(9, c(8,15), c(7, 10, 11, 12, 13, 14), 6,5,4,3, 2))
H[1:2]=c(0,0.3)*N_i[1:2]
H = 0.5*H # hesitancy down by a factor of 2 
test1vocLH = run_over_scen_3(1.1, 2.5,ve = 0.5,vp=0.75, trytime = 180,T2=720,speedup = 2)
test2vocLH = run_over_scen_3(1.1, 2,ve = 0.5,vp=0.75,  trytime = 180,T2=720,speedup = 2)
tmp = compare_sims(filter(test1voc, time <=600), filter(test1vocLH, time<=600),
                   textsize = 10, 
                   name1 = "R = 2.5 w VOC", 
                   name2 = "R = 2.5 w VOC, low hes")
toprow = ggarrange(plotlist = tmp[1:3], nrow=1,
                   common.legend = TRUE, legend="top")
toprow

dd=total_cases_origin(filter(test1vocLH, time >= 180 , time <=600))
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
ggarrange(toprow,botrow,nrow = 2,heights = c(1.3,1))
ggsave("reopen-kids-voc-lowhes.pdf", width = 12,height = 8)


# ---- herd immunity calculations figure ---- 

# reset H, ve, vp, strategy 
phac.accept = c(72,72, 75.18, 78.25, 84.7, 84.7, 84.70)
# 15-19, 20-44, 45-54, 55-64, 65+, 75+, 80+ 
# these will give a hes like this:
H = c(0,0,0.28, 0.28, 0.235, 0.2, 0.153, 0.153, 0.153) 
H=c(H, H[3:8])*N_i # with EW and scaled to ON pop size 
ve=0.8
vp=0.75
# my best guess for the PHAC rollout 
strategies <- list( list(9, c(8,15), c(7, 10, 11, 12, 13, 14), 6,5,4,3))

labels = "PHAC strategy" 

pars <- crossing(R=c(1.05,1.15,1.3, 1.5,1.8,2.0,2.2, 2.5, 3, 3.5), ve = 0.8, vp = 0.75, scen=1)
# NOTE I  have set T2 in run_.. 2 to be 110 not 210 !! 

# RUN (according to piecewise scenario defined in run_over_scen_2)
res <- pars %>%  future_pmap_dfr(run_over_scen_2, .progress=TRUE)

#---SUMMARIZE
mycases_imm <- function(df){
  t <- time_to_decr(df)
  n <- rowSums(df[t, grep("R|D|I|E", names(df))]) # all infections
  return(n)
}
myvax_imm <- function(df){
  t <- time_to_decr(df)
  n <- rowSums(df[t, grep("v", names(df))]) # NOT the x bc they weren't protected
  return(n)
}
double_counted <- function(df){ 
  t <- time_to_decr(df)
  n <- rowSums(df[t, grep("Rv", names(df))]) # subtract these dudes 
  return(n)
} 

res2 <- res %>% 
  group_by(type, R, ve,vp) %>%
  nest() %>%
  summarize(t_turn = map_dbl(data, time_to_decr),
            cases_turn = map_dbl(data,mycases_imm) - map_dbl(data, double_counted),
            vax_turn = map_dbl(data, myvax_imm))
res2$type = "Age and contact model"
theorydf = data.frame(R=res2$R, ve = 0.8, vp =0.75, type="Theory", t_turn = NA, 
                      cases_turn=0, vax_turn=pop_total*(1-1/res2$R))

resboth = rbind(res2, theorydf)

tsize <- 12


ggplot(resboth, aes(x=R, y=(cases_turn+vax_turn)/pop_total, group=type, col=type))+
  geom_line(size=1.5)+
  geom_point(size=2)+
 # facet_wrap(~ve, labeller=label_both,ncol=3)+ 
 # scale_color_brewer(palette = "Dark2")+
  theme(text=element_text(size=tsize))+
  theme(panel.spacing.x=unit(1.5, "lines"), legend.position = "bottom")+
  labs(x='R',y='Proportion of Population',col='')

ggsave(file = "herd-immunity-compare-simple.pdf", height = 5, width = 7)
# ---- end herd imm figure 








# try ramp
test = run_over_scen_ramp(1.3, 2.5, trytime = 60)
tmp = compare_sims(test1, test,textsize = 10, 
                   name1 = "july no ramp", name2 = "ramp")
ggarrange(plotlist = tmp[1:2], nrow=1)










