library(rvest)
url <- "https://www.canada.ca/en/public-health/services/diseases/2019-novel-coronavirus-infection/prevention-risks/covid-19-vaccine-treatment/vaccine-rollout.html"
tables <- read_html(url) %>% html_table()
delivery_data <- bind_rows(tables[[2]] %>% 
                               pivot_longer(-`Distribution location`) %>%
                               mutate(type="Pfizer/BioNTech") ,
                           tables[[3]] %>% 
                               pivot_longer(-`Distribution location`) %>%
                               mutate(type="Moderna")) %>%
    mutate(Date=as.Date(paste0(gsub("^.+- *","",name)," 2021"),format="%d %b %Y"),
           value=gsub(",","",value) %>% as.integer(),
           `Distribution location`=gsub(" +"," ",`Distribution location`))


vaccine_age <- read_csv("https://health-infobase.canada.ca/src/data/covidLive/vaccination-coverage-byAgeAndSex.csv",col_types = cols(.default="c")) %>%
    mutate_at(vars(matches("num")),as.numeric)
vaccine_age = vaccine_age %>% mutate(week_end = ymd(week_end))
vbc = filter(vaccine_age, prename=="British Columbia") %>% 
    group_by(prename,week_end) %>% 
    summarise(total = sum(numtotal_atleast1dose))
vbc
ggplot(vbc,  aes(x=week_end, y=total)) + geom_point() # + ylim(c(0, 300000))

# vaccination rate 
mean(diff(vbc$total)[1:8])/(7*5e6) # 0.0043 in period 1 but that's per *week*?! what the heck

mean(diff(vbc$total)[9:16])/(7*5e6)#  % half a percent per day 


