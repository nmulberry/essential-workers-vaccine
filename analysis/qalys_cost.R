tsize <- 16
R_vec <- c(1.15,1.3) # R vals to plot
data = readxl::read_xlsx(paste0(PATH,"data/tupper_QALY_dollars.xlsx"))

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
                            values_to ="qaly") %>%
    mutate(Source = case_when(Source == "Cases" ~ "Community Infections",
                              TRUE ~ Source))

qalys_res$Source <- factor(qalys_res$Source, 
    levels=c("Deaths", "Hospitalization","Chronic","Community Infections"))

ggplot(data=filter(qalys_res, R %in% R_vec), aes(x=type,y=qaly,fill=Source))+
    geom_bar(stat="identity",alpha = 0.8)+
    facet_grid(R~ve,scales= "free_y", labeller=label_both)+
    theme(panel.spacing.x=unit(1.5, "lines"),
          axis.text.x = element_text(angle = 35,hjust = 1),
          axis.title.x = element_blank(),text=element_text(size=tsize))+
    ylab("Health loss (QALYs)")+
    scale_fill_brewer(palette = "Paired") 

ggsave(paste0(PATH,"figures/qalybars.pdf"),height=7, width=10)

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
    mutate(Savings = max(Cost)-Cost)

cost_res$Source <- factor(cost_res$Source, 
    levels=c("QALY", "Hospitalization","Chronic"))

ggplot(data=filter(cost_res, R %in% R_vec), aes(x=type,y=Cost/1e6,fill=Source))+
            geom_bar(stat="identity",alpha = 0.8)+
            facet_grid(R~ve,scales= "free_y", labeller=label_both)+
            theme(panel.spacing.x=unit(1.5, "lines") ,
                  axis.text.x = element_text(angle = 35,hjust = 1),
                  axis.title.x = element_blank(),text=element_text(size=tsize))+
    ylab("NMB loss (millions)")+
    scale_fill_brewer(palette = "Paired") 

ggsave(paste0(PATH,"figures/costbars.pdf"),height=7, width=10)


ggplot(data=filter(cost_res, R %in% R_vec, type != labels[1]), aes(x=type,y=Savings/1e6,fill=Source))+
            geom_bar(stat="identity",alpha = 0.8)+
            facet_grid(R~ve,scales= "free_y", labeller=label_both)+
            theme(panel.spacing.x=unit(1.5, "lines") ,
                  axis.text.x = element_text(angle = 35,hjust = 1),
                  axis.title.x = element_blank(),text=element_text(size=tsize))+
    ylab("NMB savings (millions)")+
    scale_fill_brewer(palette = "Paired") 


ggsave(paste0(PATH,"figures/savingsbars.pdf"),height=7, width=10)


 
