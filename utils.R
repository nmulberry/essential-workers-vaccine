
#######################################
# Helper functions for summarizing output
#######################################

total_vaccinated = function(df){
  ind <- grep("v|x", names(df))
  n <- sum(tail(df, n=1)[,ind])
  return(n)
}

total_deaths = function(df) {
  # D compartments
 ind <- grep("D", names(df))
 n <- sum(tail(df, n=1)[,ind])
 return(n)
}

total_cases = function(df) {
  # R + D compartments
 ind <- grep("D|R", names(df))
 n <- sum(tail(df, n=1)[,ind])
 return(n)
}


total_cases_origin = function(df, ifr=IFR, ihr = IHR, hosp_efficacy=vp) {
# note -- the D group doesn't track whether the individual was vaccinated or not
# this makes counting the cases and splitting by whether they were vaccinated 
  # challenging
  
# however, the deaths are just like the R but R has (1-IFR) and D has IFR
# so if I take the Rs and multiply by 1/(1-IFR) i will get the totals 

  ind = grep("R", names(df))
  ind2=grep("Rv", names(df)) # Rv 
  ind3=grep("Rx", names(df)) # Rx 
  ind1=setdiff(ind, c(ind2, ind3)) # just R 
  ifr = ifr[1:length(ind1)] # 9 age groups, df is collected to not have ew, IFR has all 15
  ihr = ihr[1:length(ind1)]
  lastrow =tail(df,n=1)
age_names = c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+")
dd1 =data.frame(age_band =age_names , 
                cases= as.numeric(lastrow[ind1]*(1/(1-ifr))), 
                hosp = as.numeric(lastrow[ind1]*ihr), 
                 Protection = "Unvaccinated") # unvax, by age 
dd2 = data.frame(age_band =age_names , # these are the Rvs - vax after exposure
                 cases= as.numeric(lastrow[ind2]*(1/(1-ifr))),
                 #hosp = as.numeric(lastrow[ind2]*ihr*(1-hosp_efficacy)), 
                 hosp = as.numeric(lastrow[ind2]*ihr),  # note remove eff here
                Protection="Vaccinated after exposure") # vax, by age 
dd3 = data.frame(age_band =age_names, 
                 cases= as.numeric(lastrow[ind3]*(1/(1-ifr))),
                 hosp = as.numeric(lastrow[ind3]*ihr*(1-hosp_efficacy)), 
                 Protection="Vaccinated but not protected") # vax unprot, by age 
return(rbind(dd1,dd2,dd3))
}

total_hosp = function(df, ihr = IHR, hosp_efficacy=vp) {
  ind = grep("R", names(df))
#   scalevec=c(IHR, (1-hosp_efficacy)*IHR, (1-hosp_efficacy)*IHR) 
    scalevec=c(IHR, IHR, (1-hosp_efficacy)*IHR) # Rv are not protected by efficacy
  htot  <- scalevec*tail(df, n=1)[,ind]
  return(sum(htot))
}


total_long = function(df, lcr = LCR, hosp_efficacy=vp,ltfac=1) {
  ind = grep("R", names(df))
#   scalevec=c(LCR, (1-ltfac*hosp_efficacy)*LCR, (1-ltfac*hosp_efficacy)*LCR) 
    scalevec=c(LCR, (1-ltfac)*LCR, (1-ltfac*hosp_efficacy)*LCR)  # Rv are not protected by efficacy
  ltot  <- scalevec*tail(df, n=1)[,ind]
  return(sum(ltot))
}

total_YLL = function(df,YLL_vec){
  # get total YLL from output
 ind <- grep("D", names(df))
 df <- tail(df, n=1)[,ind] # deaths
 n <- sum(df*YLL_vec)
 return(n)
}

total_YL_infect = function(df, YLL_vec){
  # get total "years of life after infection" from output
  ind <- grep("R", names(df))
   R_tot <- tail(df, n=1)[,ind] # R+Rx+Rv
   n <- sum(R_tot*YLL_vec)
 return(n)
}

# qaly lost due to death 
qaly1 = function(df, new_chris_yll) {
  ind <- grep("D", names(df))
  df <- tail(df, n=1)[,ind] # deaths
  ds=df[1:9]
  return(sum(ds*new_chris_yll))
 }

 # qaly lost due to cases:
qaly2 = function(df, SDUR_vec,caseqalyfactor) { 
  ind <- grep("D|R", names(df))[1:36] # AAA HARD CODED... BEWARE
  # now I need total cases by age
  cases=tail(df, n=1)[,ind]
  cases=cases[1:9]+cases[10:18]+cases[19:27]+cases[28:36]
#   return( sum(caseqalyfactor*cases*SDUR_vec/365) )# community cases
     return( sum(caseqalyfactor*cases/365) )# community cases GUESSING
 }

qaly3 = function(df,IHR, HOSPDUR,hospq ) { 
  ind <- grep("D|R", names(df))[1:36] # AAA HARD CODED... BEWARE
  # now I need total cases by age
  cases=tail(df, n=1)[,ind]
  cases=cases[1:9]+cases[10:18]+cases[19:27]+cases[28:36]
  q3 = sum(hospq*cases*IHR*as.numeric(HOSPDUR)/365)+ sum(0.1*cases*IHR)
}

qaly4 = function(df,discount_rate=0.015, YLL_vec, IHR, chronicq ) {
  ind <- grep("D|R", names(df))[1:36] # AAA HARD CODED... BEWARE
  # now I need total cases by age
  cases=tail(df, n=1)[,ind]
  cases=cases[1:9]+cases[10:18]+cases[19:27]+cases[28:36]
  r=discount_rate # paul says there is a canadian doc, cadth, with this number in it
#  q4= sum(chronicq*cases*(1/r)*(1-exp(-r*pmin(10, YLL_vec)))*IHR*0.1 )
return( sum(chronicq*cases*IHR*0.2 )) # GUESSING 
}


# cost per qaly
cost1 = function(df, yll ) {
 totqaly = qaly1(df,yll)+
   qaly2(df,SDUR_vec[1:9],caseqalyfactor)+
   qaly3(df,IHR, HOSPDUR , hospq= hospqfactor)+
   qaly4(df,discount_rate = 0.015, YLL_vec, IHR, chronicq)
 return(totqaly*30000)
}
  
# cost per hosp
cost2=function(df, ihr = IHR, vp){
  tothosp = total_hosp(df, IHR,hosp_efficacy = vp)
  return(tothosp*1392*8)
}

# cost per chronic
cost3 = function(df, IHR, chroniccost = chroniccost) {
  ind <- grep("D|R", names(df))[1:36] # AAA HARD CODED... BEWARE
  # now I need total cases by age
  cases=tail(df, n=1)[,ind]
  cases=cases[1:9]+cases[10:18]+cases[19:27]+cases[28:36]
     return( sum(cases*IHR*0.2*chroniccost )) # GUESSING 
}



combine_age_groups = function(sim){
  # take a dataframe with essential worker compartments
  # and add together the corresponding age groups
  #-------------------------------------------#
  e_ind <- grep("e", names(sim)) 
  e_ind <- e_ind[e_ind >1] #time/date will always pop up
  if (length(e_ind)==0) {
    return(sim)}
  else {
    newsim <- sim %>%
      mutate(S3=S3+Se3,S4=S4+Se4,S5=S5+Se5, S6=S6+Se6, S7=S7+Se7,S8=S8+Se8,
        Sv3=Sv3+Sve3,Sv4=Sv4+Sve4,Sv5=Sv5+Sve5, Sv6=Sv6+Sve6, Sv7=Sv7+Sve7,Sv8=Sv8+Sve8,
        Sx3=Sx3+Sxe3,Sx4=Sx4+Sxe4,Sx5=Sx5+Sxe5, Sx6=Sx6+Sxe6, Sx7=Sx7+Sxe7,Sx8=Sx8+Sxe8,
        E3=E3+Ee3,E4=E4+Ee4,E5=E5+Ee5, E6=E6+Ee6, E7=E7+Ee7,E8=E8+Ee8,
        Ev3=Ev3+Eve3,Ev4=Ev4+Eve4,Ev5=Ev5+Eve5, Ev6=Ev6+Eve6, Ev7=Ev7+Eve7,Ev8=Ev8+Eve8,
        Ex3=Ex3+Exe3,Ex4=Ex4+Exe4,Ex5=Ex5+Exe5, Ex6=Ex6+Exe6, Ex7=Ex7+Exe7,Ex8=Ex8+Exe8,
        I3=I3+Ie3,I4=I4+Ie4,I5=I5+Ie5, I6=I6+Ie6, I7=I7+Ie7,I8=I8+Ie8,
        Iv3=Iv3+Ive3,Iv4=Iv4+Ive4,Iv5=Iv5+Ive5, Iv6=Iv6+Ive6, Iv7=Iv7+Ive7,Iv8=Iv8+Ive8,
        Ix3=Ix3+Ixe3,Ix4=Ix4+Ixe4,Ix5=Ix5+Ixe5, Ix6=Ix6+Ixe6, Ix7=Ix7+Ixe7,Ix8=Ix8+Ixe8,
        R3=R3+Re3,R4=R4+Re4,R5=R5+Re5, R6=R6+Re6, R7=R7+Re7,R8=R8+Re8,
        Rv3=Rv3+Rve3,Rv4=Rv4+Rve4,Rv5=Rv5+Rve5, Rv6=Rv6+Rve6, Rv7=Rv7+Rve7,Rv8=Rv8+Rve8,
        Rx3=Rx3+Rxe3,Rx4=Rx4+Rxe4,Rx5=Rx5+Rxe5, Rx6=Rx6+Rxe6, Rx7=Rx7+Rxe7,Rx8=Rx8+Rxe8,
        D3=D3+De3,D4=D4+De4,D5=D5+De5, D6=D6+De6, D7=D7+De7,D8=D8+De8)%>%
    select(-all_of(e_ind))
  }


  return(newsim)
}

remove_essential = function(x){
  # helper func to combine age demo back to 9 grps
  if (length(x)==0){
    return(x)
  }else{
  new_x <- rep(0, 9)
  new_x[1:2] <- x[1:2]
  new_x[3:8] <- x[3:8]+x[10:15]
  new_x[9] <- x[9]  
  return(new_x)}
}



get_num_vax <- function(sim){
  # total number vaccinated over time by age
  res <- data.frame(date=startDate+sim$time)
  ind <- grep("v|x", names(sim))
  sim <- sim[,ind]
  res$num <- rowSums(sim)
  return(res)

}

add_hosp_lc <-  function(output, ihr = IHR,lcr = LCR,
           hosp_efficacy = vp, lcrfac=1, 
           dI=1/5, dur_hosp = as.numeric(HOSPDUR), dur_ltc = longDUR, num_groups = 9) {

    Is    = output[,grep("I", colnames(output))] # this is 27 columns of Is

  # Recoveries + hosp + long covid ~=  those coming into R from I
  rec =dI*Is[,1:num_groups] # NEW recoveries
  recv = dI*Is[,(num_groups+1):(2*num_groups)]
  recx = dI*Is[,(2*num_groups+1):(3*num_groups)]

H = 0*rec; Hv = 0*recv; Hx = 0*recx # init
colnames(H) <- colnames(Hv) <- colnames(Hx) <- paste("H", 1:num_groups, sep="")
if (length(hosp_efficacy)==1) hosp_efficacy = rep(hosp_efficacy, length(H))
L = 0*rec; Lv = 0*recv; Lx = 0*recx # init
colnames(L) <- colnames(Lv) <- colnames(Lx) <- paste("L", 1:num_groups, sep="")

  # and in the h matrix , each column is a function of the corresponding column of rec
  # and it needs two other inputs: the hosp rate (age dep) and the hosp duration (age dep)
  for (n in 1:num_groups) {
    H[,n] = get_hosp(rec[,n], thishr = ihr[n], thisdh = dur_hosp[n])
  #  Hv[,n] = get_hosp(recv[,n], thishr = (1-hosp_efficacy[n])*ihr[n], thisdh = dur_hosp[n])
      Hv[,n] = get_hosp(recv[,n], thishr = ihr[n], thisdh = dur_hosp[n])
    Hx[,n] = get_hosp(recx[,n], thishr = (1-hosp_efficacy[n])*ihr[n], thisdh = dur_hosp[n])

    L[,n] = get_hosp(rec[,n], thishr = lcr[n], thisdh = dur_ltc[n])
 #    Lv[,n] = get_hosp(recv[,n], thishr = (1-lcrfac*hosp_efficacy[n])*lcr[n], thisdh = dur_ltc[n])
    Lv[,n] = get_hosp(recv[,n], thishr = lcr[n], thisdh = dur_ltc[n])
    Lx[,n] = get_hosp(recx[,n], thishr = (1-lcrfac*hosp_efficacy[n])*lcr[n], thisdh = dur_ltc[n])
  }
H = H+Hv+Hx; L = L+Lv+Lx
return(cbind(output, H, L))
}


extract_incidence <- function(output,d_E= 1/3) {
  Es    = output[,grep("E", colnames(output))] # this is 27 columns of Es
  L=ncol(Es)/3 # number of cols in each of E, Ev and Ex
  # add them up. The incidence into I is those coming into I from E
  incid =d_E*(Es[,1:L] +  Es[,(L+1):(2*L)] + Es[,(2*L+1):(3*L)])
  mydatC = cbind(select(output, time), incid) %>%
    rename( `0-9`=E1, `10-19`=E2, `20-29`=E3, `30-39`=E4,
            `40-49`=E5, `50-59`=E6, `60-69`=E7, `70-79`=E8, `80+`=E9)   %>%
    pivot_longer(cols = 2:10, names_to = "age_band", values_to = "incid") %>%
    mutate(age_band = as.factor(age_band))
  return(mydatC)
}


extract_totals <- function(output) {
  return( extract_cases_deaths(output) %>% group_by(time) %>%
      summarise(totcases=sum(cases), totdeaths=sum(deaths)))
}

extract_cases_deaths <- function(output, LCFAC=1) {

  allcases = select(output,  I1, I2, I3, I4, I5, I6, I7, I8, I9) +
    select(output,  Iv1, Iv2, Iv3, Iv4, Iv5, Iv6, Iv7, Iv8, Iv9)+
    select(output,  Ix1, Ix2, Ix3, Ix4, Ix5, Ix6, Ix7, Ix8, Ix9)

    mydatC = cbind(select(output, time), allcases) %>%
        rename( `0-9`=I1, `10-19`=I2, `20-29`=I3, `30-39`=I4,
                `40-49`=I5, `50-59`=I6, `60-69`=I7, `70-79`=I8, `80+`=I9)   %>%
        pivot_longer(cols = 2:10, names_to = "age_band", values_to = "cases") %>%
        mutate(age_band = as.factor(age_band))

    mydatD = output %>%  # mutate(D1=0*time, D2=0*time, D3=0*time) %>%
        select( time,D1, D2, D3,  D4, D5, D6, D7, D8, D9) %>%
        rename( `0-9`=D1, `10-19`=D2, `20-29`=D3, `30-39`=D4,
                `40-49`=D5, `50-59`=D6, `60-69`=D7, `70-79`=D8, `80+`=D9) %>%
        pivot_longer(cols = 2:10, names_to = "age_band", values_to = "deaths") %>%
        mutate(age_band = as.factor(age_band))

        mydatDI = output %>%  # mutate(D1=0*time, D2=0*time, D3=0*time) %>%
      select( time,D1, D2, D3,  D4, D5, D6, D7, D8, D9) %>%
      mutate(D1=c(0, diff(D1)), D2=c(0, diff(D2)),
             D3=c(0, diff(D3)), D4=c(0, diff(D4)),
             D5=c(0, diff(D5)), D6=c(0, diff(D6)),
             D7=c(0, diff(D7)), D8=c(0, diff(D8)),
             D9=c(0, diff(D9))) %>%
      rename( `0-9`=D1, `10-19`=D2, `20-29`=D3, `30-39`=D4,
              `40-49`=D5, `50-59`=D6, `60-69`=D7, `70-79`=D8, `80+`=D9) %>%
      pivot_longer(cols = 2:10, names_to = "age_band", values_to = "deaths") %>%
      mutate(age_band = as.factor(age_band))

    mydatE = add_hosp_lc(output, hosp_efficacy=unique(output$vp)) %>% 
      select(time, H1, H2, H3, H4, H5, H6, H7, H8, H9) %>%
      rename( `0-9`=H1, `10-19`=H2, `20-29`=H3, `30-39`=H4,
              `40-49`=H5, `50-59`=H6, `60-69`=H7, `70-79`=H8, `80+`=H9) %>%
      pivot_longer(cols = 2:10, names_to = "age_band", values_to = "hosp") %>%
      mutate(age_band = as.factor(age_band))

    mydatF = add_hosp_lc(output, hosp_efficacy=unique(output$vp), lcrfac=LCFAC) %>%
      select(time, L1, L2, L3, L4, L5, L6, L7, L8, L9) %>%
      rename( `0-9`=L1, `10-19`=L2, `20-29`=L3, `30-39`=L4,
              `40-49`=L5, `50-59`=L6, `60-69`=L7, `70-79`=L8, `80+`=L9) %>%
      pivot_longer(cols = 2:10, names_to = "age_band", values_to = "long") %>%
      mutate(age_band = as.factor(age_band))
    

    if (all(mydatC$age_band==mydatD$age_band) &
        all(mydatC$time==mydatD$time)) { mydatC$deaths = mydatD$deaths}
  mydatC$incid = extract_incidence(output)$incid
  mydatC$hosp= mydatE$hosp
  mydatC$newdeaths = mydatDI$deaths
  mydatC$long = mydatF$long
    return(mydatC)
}

# like extract_cases_deaths, but keeping track of whether they're vax'd or not



get_hosp <- function(r_col, thishr, thisdi= 1/5,thisdh = 10) {
  hosp = 0*r_col # init
  for (k in 2:length(hosp)) {
    hosp[k] = max(0,hosp[k-1]+thishr*r_col[k-1] - (1/thisdh)*hosp[k-1])
  }
  return(hosp)
}

time_to_decr = function(df){
  ind <- grep("I", names(df))
  turn_point <- which(diff(rowSums(df[,ind]))<0)
  turn_point <- turn_point[turn_point > 10] # ignore transients...
  turn_point <- min(turn_point)
  return(df$time[turn_point])
}

herd_immunity = function(df){
  t <- time_to_decr(df)
  if (t==0){
    n=0}
  else {
    ind <- grep("v|x", grep("S", names(df), value=T), invert=T)
    n <- pop_total-rowSums(df[t,ind])
  }
  
  return(n)
}

cases_immunity = function(df){
  t <- time_to_decr(df)
  n <- rowSums(df[t, grep("R|D", names(df))])
  return(n)
}
vax_immunity = function(df){
  t <- time_to_decr(df)
  n <- rowSums(df[t, grep("v|x", names(df))])
  return(n)
}

######################
# PLOTS
#######################
compare_sims <- function(sim1, sim2, name1="1", name2="2",
                             startDate = ymd("2021-01-10"), LCFAC=1,
                             stages=NULL,textsize=16, scale_y=TRUE) {
  if (length(LCFAC)==1) {LCFAC1 = LCFAC; LCFAC2=LCFAC} 
  if (length(LCFAC)==2) {LCFAC1 = LCFAC[1]; LCFAC2=LCFAC[2]} 

  o1 = extract_cases_deaths(sim1,LCFAC = LCFAC1); o1$scen = name1
  o3 = extract_cases_deaths(sim2,LCFAC = LCFAC2); o3$scen= name2
  oo = rbind(o1, o3); oo$scen=as.factor(oo$scen)
  oo$date = startDate + oo$time 
yscaler = ifelse(scale_y == TRUE, 1e5/pop_total, 1)
ystr = ifelse(scale_y, "per 100K", "")
  p1=  ggplot(data = oo, aes(x=date, y=incid*yscaler, fill=age_band))+theme_light()+
    facet_wrap(~scen,nrow = 1) +
    geom_area(position="stack",alpha=0.7)+guides(fill=FALSE)+
    theme(axis.title.x = element_blank(),text=element_text(size=textsize))+
    ylab(paste("Incidence", ystr))
  
  p2=  ggplot(data = oo, aes(x=date, y=hosp*yscaler, fill=age_band))+theme_light()+
    facet_wrap(~scen,nrow = 1) +
    geom_area(position="stack",alpha=0.7)+guides(fill=FALSE)+
    theme(axis.title.x = element_blank(),text=element_text(size=textsize))+
    ylab(paste("Hospitalizations", ystr))
  
  p3 = ggplot(data = oo, aes(x=date, y=newdeaths*yscaler, fill=age_band))+
    facet_wrap(~scen,nrow = 1) +
    theme_light()+
    geom_area(position="stack",alpha=0.7)+guides(fill=FALSE)+
    theme(axis.title.x = element_blank(),text=element_text(size=textsize)) +
    ylab(paste("Daily deaths",ystr))
  p4 = ggplot(data = oo, aes(x=date, y=long*yscaler, fill=age_band))+
    facet_wrap(~scen,nrow = 1) +
    theme_light()+
    geom_area(position="stack",alpha=0.7)+guides(fill=FALSE)+
    theme(axis.title.x = element_blank(),text=element_text(size=textsize)) +
    ylab(paste("Long covid",ystr))
  
  
  if (!is.null(stages)){
    shade_phases <- data.frame(xstart=stages[c(TRUE,FALSE)], xend=stages[c(FALSE,TRUE)])
    p1 <- p1 +   geom_rect(data = shade_phases, aes(xmin = xstart, xmax = xend, 
                                                    ymin = -Inf, ymax = Inf), alpha = 0.15, inherit.aes=FALSE)
    p2 <- p2 +    geom_rect(data = shade_phases, aes(xmin = xstart, xmax = xend, 
                                                     ymin = -Inf, ymax = Inf), alpha = 0.15, inherit.aes=FALSE)
    p3 <- p3 +     geom_rect(data = shade_phases, aes(xmin = xstart, xmax = xend, 
                                                      ymin = -Inf, ymax = Inf), alpha = 0.15, inherit.aes=FALSE)
    p4 <- p4 +   geom_rect(data = shade_phases, aes(xmin = xstart, xmax = xend, 
                                                    ymin = -Inf, ymax = Inf), alpha = 0.15, inherit.aes=FALSE)
    
  }
  
  return(list(p1,p2,p3,p4))
}


display_prop_vax <- function(sim,startDate=lubridate::ymd("2021-01-01"),
                             label, stages=NULL, textsize=16){
# plot total number vaccinated over time per age group
# mostly for debugging/double checking sims

  N_i <- remove_essential(N_i)
  t <- sim$time
  ind <- grep("v|x", names(sim))
  nvax <- data.frame(time = sim$time)
  sim <- sim[,ind]

  tot_by_age <- function(g,dat, t,N_i){
    ind <- grep(g, names(dat))
    tot <- rowSums(dat[,ind])
    prop <- 1-(N_i[g] - tot)/(N_i[g])
    age_groups <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
    return(data.frame(time=t, Total=tot, Proportion=prop, Age=age_groups[g]))
  }

  res <- map_df(1:length(N_i), tot_by_age,sim,t, N_i) 
  res$date <- startDate+res$time

  # as a proportion of ppl in that age group
  g2 <- ggplot(res, aes(x=date, y=Proportion, group=Age, col=Age)) +
    geom_line(size=2)+
    labs(x='',y='', col='') +
            theme(legend.position = "right",
                 axis.title.x = element_blank(),text=element_text(size=textsize))+
      ggtitle(label)

    if (!is.null(stages)){
      shade_phases <- data.frame(xstart=stages[c(TRUE,FALSE)], xend=stages[c(FALSE,TRUE)])
      g2 <- g2 +   geom_rect(data = shade_phases, aes(xmin = xstart, xmax = xend, 
                 ymin = -Inf, ymax = Inf), alpha = 0.15, inherit.aes=FALSE)

  }

  return(g2)

}


plot_grid_barplots <- function(manysims, numcol = NULL,tsize = 16, includeLong=FALSE){
  # another plotting function
  # takes multiple simulation output (defined by type, R0)
  # creates barplot grid
  manysims <- manysims %>% 
    group_by(type, R0) %>%
    nest() %>%
    summarize(cases=map_dbl(data, total_cases),
              hosp=map_dbl(data, total_hosp,hosp_efficacy=ve_p),
              deaths=map_dbl(data, total_deaths),
              long = map_dbl(data, total_long, hosp_efficacy=ve_p))
     

  p1 <- ggplot(manysims, aes(x=type,y=cases, fill=veff))+
    geom_bar(stat="identity",alpha=0.75) + guides(fill=FALSE) +
    facet_wrap('R0', labeller = label_both, scales="free", ncol=numcol)+labs(x='', y='Cases')+
       scale_fill_brewer(palette = "Dark2")+
            theme(axis.title.x = element_blank(),text=element_text(size=tsize))+theme(panel.spacing.x=unit(1.5, "lines") ,
                          axis.text.x = element_text(angle = 35,hjust = 1),
                          axis.title.x = element_blank(),legend.position = "none")

  p2 <- ggplot(manysims, aes(x=type,y=hosp,group=type,fill=type))+
    geom_bar(stat="identity",alpha=0.75) + guides(fill=FALSE) +
    scale_fill_brewer(palette = "Dark2")+
    facet_wrap('R0', labeller = label_both, scales="free", ncol=numcol)+labs(x='', y='Hospitalizations')+
            theme(axis.title.x = element_blank(),text=element_text(size=tsize))+theme(panel.spacing.x=unit(1.5, "lines") ,
                          axis.text.x = element_text(angle = 35,hjust = 1),
                          axis.title.x = element_blank(),legend.position = "none") +theme(axis.text.x = element_blank())
  p3 <- ggplot(manysims, aes(x=type,y=deaths,fill=type, group=type))+
    geom_bar(stat="identity",alpha=0.75) + guides(fill=FALSE) +
    scale_fill_brewer(palette = "Dark2")+
    facet_wrap('R0', labeller = label_both, scales="free", ncol=numcol)+labs(x='', y='Deaths')+
            theme(axis.title.x = element_blank(),text=element_text(size=tsize))+theme(panel.spacing.x=unit(1.5, "lines") ,
                          axis.text.x = element_text(angle = 35,hjust = 1),
                          axis.title.x = element_blank(),legend.position = "none")
if (includeLong ==TRUE) { 
  p4 <-  ggplot(manysims, aes(x=type,y=long,fill=type, group=type))+
    geom_bar(stat="identity",alpha=0.75) + guides(fill=FALSE) +
    scale_fill_brewer(palette = "Dark2")+
    facet_wrap('R0', labeller = label_both, scales="free", ncol=numcol)+labs(x='', y='Long COVID')+
    theme(axis.title.x = element_blank(),text=element_text(size=tsize))+theme(panel.spacing.x=unit(1.5, "lines") ,
                               axis.text.x = element_text(angle = 35,hjust = 1),
                                axis.title.x = element_blank(),legend.position = "none")
return(list(p1,p2,p3,p4))  } else { 
  return(list(p1,p2,p3))} 
}

