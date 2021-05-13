
extract_totals_all <- function(output) {
    return( extract_cases_deaths(output) %>% group_by(time) %>%
                summarise(Prevalence=sum(cases), 
                          Deaths=sum(newdeaths),
                          Hospitalized = sum(hosp),
                          Incidence = sum(incid)))
}
one_sim_plot <- function(sim1, name1="1", 
                         startDate = ymd("2021-01-10"), LCFAC=1,
                         stages=NULL,textsize=16, scale_y=TRUE) {
if (length(LCFAC)==1) {LCFAC1 = LCFAC; LCFAC2=LCFAC} 
if (length(LCFAC)==2) {LCFAC1 = LCFAC[1]; LCFAC2=LCFAC[2]} 

# oo = extract_cases_deaths(sim1,LCFAC = LCFAC1); 
oo = extract_totals_all(sim1)
oo$scen = name1

oo$date = startDate + oo$time 
yscaler = ifelse(scale_y == TRUE, 1e5/pop_total, 1)
ystr = ifelse(scale_y, "per 100K", "")
p1=  ggplot(data = oo, aes(x=date, y=Incidence*yscaler))+theme_light()+
    geom_area(position="stack",alpha=0.7, fill = "blue")+ #guides(fill=FALSE)+
    theme(axis.title.x = element_blank(), 
          text=element_text(size=textsize),
          strip.text.x = element_text(size = textsize+1))+
    ylab(paste("Incicence", ystr))+scale_x_date( date_labels = "%m-%y")

p2=  ggplot(data = oo, aes(x=date, y=Hospitalized*yscaler))+theme_light()+
    geom_area(position="stack",alpha=0.7, fill = "red")+# guides(fill=FALSE)+
    theme(axis.title.x = element_blank(), 
          text=element_text(size=textsize),
          strip.text.x = element_text(size = textsize+1))+
    ylab(paste("Hospitalizations", ystr))+scale_x_date( date_labels = "%m-%y")

p3 = ggplot(data = oo, aes(x=date, y=Deaths*yscaler))+
    theme_light()+
    geom_area(position="stack",alpha=0.7, fill = "grey")+ #guides(fill=FALSE)+
    theme(axis.title.x = element_blank(), 
          text=element_text(size=textsize),
          strip.text.x = element_text(size = textsize+1))+
    ylab(paste("Daily deaths",ystr))+scale_x_date( date_labels = "%m-%y")



return(list(p1,p2,p3))
}


