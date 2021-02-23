
run_sim_basic = function(C, I_0, percent_vax, strategy, num_perday, v_e, v_p,
                       u = u_var, sero = sero_none, syn_sero_compartments = NA, 
                        num_days=365,H=rep(0,num_groups), 
                        with_essential=FALSE,set_vax_prop=NULL,p_grp=NULL){
  # BASED ON run_sim from Bubar et al.
  # Vaccine rollout is continuous until all vaccines are distributed
  # strategy: list describing the target groups in order
  # e.g. strategy=list(9,8,7,6,3:5,2,1)
  #      strategy=list(9,8,3,4,5,7,6,2,1)
  if (!is.list(strategy)){
    warning("strategy should be a list of target groups.")
  }
  # Check dimensions
  if (!length(v_p)==num_groups){
    warning("wrong dimensions in v_p")
  }
  if (!length(v_e)==num_groups){
    warning("wrong dimensions in v_e")
  }
  if (!length(H)==num_groups){
    warning("wrong dimensions in H")
  }
  if (!nrow(C)==num_groups){
    warning("wrong dimensions in C")
  }
  if (!length(sero)==num_groups){
    warning("wrong dimensions in sero")
  }
  if (!length(I_0)==num_groups){
    warning("wrong dimensions in I_0")
  }


  # Disease Tranmission
  d_E <- 1/3 # incubation period (E -> I), ref: Davies
  d_I <- 1/5 # recovery period (I -> R), ref: Davies

  # Initialize
  E_0 <- Ev_0 <- Ex_0 <- Sv_0 <- Sx_0 <- Iv_0 <- Ix_0 <- Rv_0 <- Rx_0 <- D_0 <- rep(0, num_groups)
  R_0 <- N_i * sero
  E_0 = (3/5)*(I_0) # was 0 but this makes incidence start at 0 and prevalence start by falling
  S_0 <- N_i - I_0 - R_0 - E_0
  num_stages <- length(strategy)

  # starting group to vaccinate
  stage <- 1
  groups <- strategy[[stage]]

  compartments_initial <- c(S_0,Sv_0,Sx_0,E_0,Ev_0,Ex_0,I_0,Iv_0,Ix_0,R_0,Rv_0,Rx_0,D_0)

  if (length(syn_sero_compartments) > 1){
    compartments_initial <- syn_sero_compartments
  }
  vax_supply <- percent_vax*pop_total

  out <- move_vaccinated_BC(compartments_initial, strategy=strategy, stage=1, num_perday, vax_supply,
                                           v_e, H,set_vax_prop, p_grp)
  stage <- out$stage
  compartments_aftervax <- out$states

  parameters = list(u=u, C=C, d_E=d_E, d_I=d_I, v_e=v_e,  v_p=v_p)

  running = TRUE
  t <- c(0:1)
  df <- as.data.frame(deSolve::lsoda(compartments_aftervax, t, calculate_derivatives_BC, parameters,
                                     hmax = 0.01))
  df[1,] <- c(0, compartments_initial)
  vax_supply <- vax_supply - num_perday*pop_total
  vax_supply[vax_supply < 0] <- 0
  t <- t + 1

  while(running == TRUE){
    compartments_initial <- as.numeric(df[t[2], -(1)])
    # update
    out <- move_vaccinated_BC(compartments_initial, strategy=strategy, stage=stage, num_perday, vax_supply,
                                             v_e, H,set_vax_prop,p_grp)
    stage <- out$stage
    compartments_aftervax <- out$states

    parameters = list(u=u, C=C, d_E=d_E, d_I=d_I, v_e=v_e, v_p=v_p)
    temp <- as.data.frame(deSolve::lsoda(compartments_aftervax, t, calculate_derivatives_BC, parameters,
                                         hmax = 0.01))
    row.names(temp) <- t+1
    temp <- temp[-(1),]

    df <- rbind(df, temp)

    vax_supply <- vax_supply - num_perday*pop_total

    vax_supply[vax_supply < 0] <- 0

    if (vax_supply == 0){
      remaining_t = c(t[2]:num_days)
      inits <- as.numeric(df[t[2]+1, -(1)])
      parameters = list(u=u, C=C, d_E=d_E, d_I=d_I, v_e=v_e, v_p=v_p)
      temp <- as.data.frame(deSolve::lsoda(inits, remaining_t, calculate_derivatives_BC, parameters,
                                           hmax = 0.01))
      row.names(temp) <- remaining_t+1
      temp <- temp[-(1),]

      df <- rbind(df, temp)
      running = FALSE
    } else if (t[2] == num_days){
      running = FALSE
    } else {
      t <- t + 1
    }
  }

  df <- add_names(df, with_essential)

  return(df)
}


add_names = function(df, with_essential){
  if (with_essential){

    names(df) <- c("time", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9",
                  "Se3", "Se4", "Se5", "Se6", "Se7", "Se8",
                 "Sv1", "Sv2", "Sv3", "Sv4", "Sv5", "Sv6", "Sv7", "Sv8", "Sv9",
                  "Sve3", "Sve4", "Sve5", "Sve6", "Sve7", "Sve8",
                 "Sx1", "Sx2", "Sx3", "Sx4", "Sx5", "Sx6", "Sx7", "Sx8", "Sx9",
                  "Sxe3", "Sxe4", "Sxe5", "Sxe6", "Sxe7", "Sxe8", 
                 "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9",
                  "Ee3", "Ee4", "Ee5", "Ee6", "Ee7", "Ee8",
                 "Ev1", "Ev2", "Ev3", "Ev4", "Ev5", "Ev6", "Ev7", "Ev8", "Ev9",
                   "Eve3", "Eve4", "Eve5", "Eve6", "Eve7", "Eve8",
                 "Ex1", "Ex2", "Ex3", "Ex4", "Ex5", "Ex6", "Ex7", "Ex8", "Ex9",
                  "Exe3", "Exe4", "Exe5", "Exe6", "Exe7", "Exe8",
                 "I1", "I2", "I3", "I4", "I5", "I6", "I7", "I8", "I9",
                  "Ie3", "Ie4", "Ie5", "Ie6", "Ie7", "Ie8",
                 "Iv1", "Iv2", "Iv3", "Iv4", "Iv5", "Iv6", "Iv7", "Iv8", "Iv9",
                  "Ive3", "Ive4", "Ive5", "Ive6", "Ive7", "Ive8", 
                 "Ix1", "Ix2", "Ix3", "Ix4", "Ix5", "Ix6", "Ix7", "Ix8", "Ix9",
                  "Ixe3", "Ixe4", "Ixe5", "Ixe6", "Ixe7", "Ixe8",
                 "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9",
                  "Re3", "Re4", "Re5", "Re6", "Re7", "Re8",
                 "Rv1", "Rv2", "Rv3", "Rv4", "Rv5", "Rv6", "Rv7", "Rv8", "Rv9",
                   "Rve3", "Rve4", "Rve5", "Rve6", "Rve7", "Rve8", 
                 "Rx1", "Rx2", "Rx3", "Rx4", "Rx5", "Rx6", "Rx7", "Rx8", "Rx9",
                    "Rxe3", "Rxe4", "Rxe5", "Rxe6", "Rxe7", "Rxe8", 
                 "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9",
                  "De3", "De4", "De5", "De6", "De7", "De8")

  } else {
    names(df) <- c("time", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9",
                 "Sv1", "Sv2", "Sv3", "Sv4", "Sv5", "Sv6", "Sv7", "Sv8", "Sv9",
                 "Sx1", "Sx2", "Sx3", "Sx4", "Sx5", "Sx6", "Sx7", "Sx8", "Sx9",
                 "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9",
                 "Ev1", "Ev2", "Ev3", "Ev4", "Ev5", "Ev6", "Ev7", "Ev8", "Ev9",
                 "Ex1", "Ex2", "Ex3", "Ex4", "Ex5", "Ex6", "Ex7", "Ex8", "Ex9",
                 "I1", "I2", "I3", "I4", "I5", "I6", "I7", "I8", "I9",
                 "Iv1", "Iv2", "Iv3", "Iv4", "Iv5", "Iv6", "Iv7", "Iv8", "Iv9",
                 "Ix1", "Ix2", "Ix3", "Ix4", "Ix5", "Ix6", "Ix7", "Ix8", "Ix9",
                 "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9",
                 "Rv1", "Rv2", "Rv3", "Rv4", "Rv5", "Rv6", "Rv7", "Rv8", "Rv9",
                 "Rx1", "Rx2", "Rx3", "Rx4", "Rx5", "Rx6", "Rx7", "Rx8", "Rx9",
                 "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9")
  }

  return (df)
}



move_vaccinated_BC = function(x, strategy, stage, num_perday, vax_supply, v_e, H, set_vax_prop, p_grp) {


  groups <- strategy[[stage]]
  if (is.null(p_grp)){
    priority <- groups
  } else {
    priority <- p_grp[[stage]]
  }
  # move those who are vaccinated in a given day
  num_compartment <- 13
  num_groups <- length(x)/num_compartment
  S    <- as.matrix(x[1:num_groups])
  Sv   <- as.matrix(x[(1*num_groups+1):(2*num_groups)])
  Sx   <- as.matrix(x[(2*num_groups+1):(3*num_groups)])
  E    <- as.matrix(x[(3*num_groups+1):(4*num_groups)])
  Ev   <- as.matrix(x[(4*num_groups+1):(5*num_groups)])
  Ex   <- as.matrix(x[(5*num_groups+1):(6*num_groups)])
  I    <- as.matrix(x[(6*num_groups+1):(7*num_groups)])
  Iv   <- as.matrix(x[(7*num_groups+1):(8*num_groups)])
  Ix   <- as.matrix(x[(8*num_groups+1):(9*num_groups)])
  R    <- as.matrix(x[(9*num_groups+1):(10*num_groups)])
  Rv   <- as.matrix(x[(10*num_groups+1):(11*num_groups)])
  Rx   <- as.matrix(x[(11*num_groups+1):(12*num_groups)])
  D    <- as.matrix(x[(12*num_groups+1):(13*num_groups)])

  people_to_vax <- sum(S[groups] + E[groups] + R[groups])

  if (people_to_vax > 0) {
    vax_proportion <- rep(0, num_groups)

    if (is.null(set_vax_prop)){
       vax_proportion[groups] <- (S[groups] + E[groups] + R[groups])/people_to_vax
     } else{
       vax_proportion[groups] <- set_vax_prop[[stage]]    
    }

    if (vax_supply >= num_perday*pop_total){
      nvax <- num_perday*pop_total
    } else {
      nvax <- vax_supply
    }

    vax_distribution <- nvax*vax_proportion
    vax_eligible <- pmax((S+E+R)-H, 0)

    if (any(vax_distribution > vax_eligible)){
      # make sure everyone in the PRIORITY age groups are vaccinated    
      if (!all(vax_distribution[priority] > vax_eligible[priority])){
        temp <- vax_distribution
        temp[vax_distribution > vax_eligible] <- vax_eligible[vax_distribution > vax_eligible]
        leftover_vax <- sum(vax_distribution - temp)
        
        full_groups <- which(vax_distribution > vax_eligible)
        leftover_groups <- groups[!groups %in% full_groups]
        people_to_vax <- sum(vax_eligible[leftover_groups])
        
        if (people_to_vax > 0){
          # split up leftovers evenly
          vax_proportion <- rep(0, num_groups)
          vax_proportion[leftover_groups] <- vax_eligible[leftover_groups]/people_to_vax
          
          vax_leftover_dist <- leftover_vax*vax_proportion
          vax_distribution <- temp + vax_leftover_dist
        }
      }

      if (any(vax_distribution > vax_eligible)){
        temp <- vax_distribution
        temp[vax_distribution > vax_eligible] <- vax_eligible[vax_distribution > vax_eligible]
        leftover_vax <- sum(vax_distribution - temp)
        
        full_groups <- which(vax_distribution > vax_eligible)
        leftover_groups <- groups[!groups %in% full_groups]
        people_to_vax <- sum(vax_eligible[leftover_groups])
      
        if (people_to_vax > 0){
          vax_proportion <- rep(0, num_groups)
          vax_proportion[leftover_groups] <- vax_eligible[leftover_groups]/people_to_vax
          
          vax_leftover_dist <- leftover_vax*vax_proportion
          vax_distribution <- temp + vax_leftover_dist
        }
      }
      # don't go over the number eligible
      if (any(vax_distribution > vax_eligible)){ 
        vax_distribution[vax_distribution > vax_eligible] <-  vax_eligible[vax_distribution > vax_eligible]
        }

      # Now check if we need to update the stage
      if (all(vax_distribution[priority] >= vax_eligible[priority])){
         if (!(stage == length(strategy))) {
          stage <- stage+1
      }
    } 

    }

    alpha <- vax_distribution/(S+E+R)
    alpha[alpha == Inf] <- 0 # no people in S,E,R
    alpha[is.nan(alpha)] <- 0 # no vax left and no people in S,E,R
  } else {alpha <- 0}


  if(any(alpha > 1)){print("ERROR: alpha > 1 in move_vaccinated_BC")}

  # all-or-nothing
  dS  <- -as.matrix(S*alpha)
  dSv <- as.matrix(S*alpha*v_e)
  dSx <- as.matrix(S*alpha*(1-v_e))

  dE  <- -as.matrix(E*alpha) 
  dEv <- as.matrix(E*alpha*v_e)
  dEx <- as.matrix(E*alpha*(1-v_e))


  dR <- -as.matrix(R*alpha) 
  dRv <- as.matrix(R*alpha*v_e)
  dRx <- as.matrix(R*alpha*(1-v_e))

  # update compartments
  S    <- S + dS
  Sv   <- Sv + dSv
  Sx   <- Sx + dSx
  E    <- E + dE
  Ev   <- Ev + dEv
  Ex   <- Ex + dEx
  R    <- R + dR
  Rv   <- Rv + dRv
  Rx   <- Rx + dRx

  # output updated compartments
  out <- list(stage=stage, states=c(S,Sv,Sx,E,Ev,Ex,I,Iv,Ix,R,Rv,Rx,D))

}


run_sim_restart = function(C, df_0, percent_vax, strategy, num_perday,  v_e, v_p,
                       u = u_var, sero = sero_none, syn_sero_compartments = NA, num_days=365, H=rep(0,num_groups), 
                      with_essential=FALSE,set_vax_prop=NULL, p_grp=NULL){
  # EXACT SAME as run_sim_basic except initialize with a dataframe (from previous output)
  # df_0 should just be the last time stamp from previous sim
  if (!is.list(strategy)){
    warning("strategy should be a list of target groups.")
  }
  # Check dimensions
  if (!length(v_p)==num_groups){
    warning("wrong dimensions in v_p")
  }
  if (!length(v_e)==num_groups){
    warning("wrong dimensions in v_e")
  }
  if (!length(H)==num_groups){
    warning("wrong dimensions in H")
  }
  if (!nrow(C)==num_groups){
    warning("wrong dimensions in C")
  }
  if (!length(sero)==num_groups){
    warning("wrong dimensions in sero")
  }

  # Disease Tranmission
  d_E <- 1/3 # incubation period (E -> I), ref: Davies
  d_I <- 1/5 # recovery period (I -> R), ref: Davies

  num_stages <- length(strategy)
  compartments_initial <- as.numeric(df_0[, -(1)])
  df_0 <- df_0[, -(1)]
  # find out what stage we are on
  # needs to account for hesitancy
  check <- function(i){
    return(sum(df_0[1:num_groups][strategy[[i]]]) > sum(H[strategy[[i]]]))
  }
  stage <- min(which(
    sapply(1:num_stages, check),
    arr.ind=T))

  # starting group to vaccinate
  groups <- strategy[[stage]]
 
  # calculate new vax_supply (take out those already vaccinated)
  inds <- grep("v|x", names(df_0))
  vax_supply <- percent_vax*pop_total -
    sum(df_0[1,inds])

  out <- move_vaccinated_BC(compartments_initial, strategy=strategy, stage=stage, num_perday, vax_supply,
                                           v_e, H,set_vax_prop, p_grp)
  stage <- out$stage
  compartments_aftervax <- out$states

  parameters = list(u=u, C=C, d_E=d_E, d_I=d_I, v_e=v_e, v_p = v_p)

  running = TRUE
  t <- c(0:1)
  df <- as.data.frame(deSolve::lsoda(compartments_aftervax, t, calculate_derivatives_BC, parameters,
                                     hmax = 0.01))
  df[1,] <- c(0, compartments_initial)
  vax_supply <- vax_supply - num_perday*pop_total
  vax_supply[vax_supply < 0] <- 0
  t <- t + 1

  while(running == TRUE){
    compartments_initial <- as.numeric(df[t[2],-(1)]) # remove time column
    # update
    out <- move_vaccinated_BC(compartments_initial, strategy=strategy, stage=stage, num_perday, vax_supply,
                                             v_e, H,set_vax_prop,p_grp)
    stage <- out$stage
    compartments_aftervax <- out$states


    parameters = list(u=u, C=C, d_E=d_E, d_I=d_I, v_e=v_e, v_p=v_p)

    temp <- as.data.frame(deSolve::lsoda(compartments_aftervax, t, calculate_derivatives_BC, parameters,
                                         hmax = 0.01))

    row.names(temp) <- t+1
    temp <- temp[-(1),]

    df <- rbind(df, temp)

    vax_supply <- vax_supply - num_perday*pop_total
    vax_supply[vax_supply < 0] <- 0

    if (vax_supply == 0){
      remaining_t = c(t[2]:num_days)
      inits <- as.numeric(df[t[2]+1, -(1)])
      parameters = list(u=u, C=C, d_E=d_E, d_I=d_I, v_e=v_e, v_p=v_p)
      temp <- as.data.frame(deSolve::lsoda(inits, remaining_t, calculate_derivatives_BC, parameters,
                                           hmax = 0.01)) #and people still end up in the vaccinated categories
      row.names(temp) <- remaining_t+1
      temp <- temp[-(1),]

      df <- rbind(df, temp)
      running = FALSE
    } else if (t[2] == num_days){
      running = FALSE
    } else {
      t <- t + 1
    }
  }

  df <- add_names(df, with_essential)

  return(df)
}


calculate_derivatives_BC=function(t, x, parameters){
  # x is a vector of length (# model compartment types)*(# age groups)
  # S, E, I, R etc. are vectors of length num_groups
  # v_e : transmission Efficacy
  # v_p : IFR protection
  num_compartment <- 13
  num_groups <- length(x)/num_compartment
  S    <- as.matrix(x[1:num_groups])
  Sv   <- as.matrix(x[(1*num_groups+1):(2*num_groups)])
  Sx   <- as.matrix(x[(2*num_groups+1):(3*num_groups)])
  E    <- as.matrix(x[(3*num_groups+1):(4*num_groups)])
  Ev   <- as.matrix(x[(4*num_groups+1):(5*num_groups)])
  Ex   <- as.matrix(x[(5*num_groups+1):(6*num_groups)])
  I    <- as.matrix(x[(6*num_groups+1):(7*num_groups)])
  Iv   <- as.matrix(x[(7*num_groups+1):(8*num_groups)])
  Ix   <- as.matrix(x[(8*num_groups+1):(9*num_groups)])
  R    <- as.matrix(x[(9*num_groups+1):(10*num_groups)])
  Rv   <- as.matrix(x[(10*num_groups+1):(11*num_groups)])
  Rx   <- as.matrix(x[(11*num_groups+1):(12*num_groups)])
  D    <- as.matrix(x[(12*num_groups+1):(13*num_groups)])

  S[S < .Machine$double.eps] <- 0
  Sv[Sv < .Machine$double.eps] <- 0
  Sx[Sx < .Machine$double.eps] <- 0
  E[E < .Machine$double.eps] <- 0
  Ev[Ev < .Machine$double.eps] <- 0
  Ex[Ex < .Machine$double.eps] <- 0
  I[I < .Machine$double.eps] <- 0
  Iv[Iv < .Machine$double.eps] <- 0
  Ix[Ix < .Machine$double.eps] <- 0
  R[R < .Machine$double.eps] <- 0
  Rv[Rv < .Machine$double.eps] <- 0
  Rx[Rx < .Machine$double.eps] <- 0
  D[D < .Machine$double.eps] <- 0

  u <- parameters$u
  C <- parameters$C
  d_E <- parameters$d_E
  d_I <- parameters$d_I
  v_e <- parameters$v_e
  v_p <- parameters$v_p

  lambda <- as.matrix((C)%*%as.matrix((I+Iv+Ix)/N_i) * as.matrix(u) )

  # all-or-nothing
  dSv <- rep(0, num_groups)
  dEv <- -d_E*as.matrix(Ev)
  

  dS  <- -as.matrix(S*lambda)
  dSx <- -as.matrix(Sx*lambda)

  dE  <- as.matrix(S*lambda) - d_E*as.matrix(E)
  dEx <- as.matrix(Sx*lambda) - d_E*as.matrix(Ex)

  dI  <- as.matrix(E*d_E) - as.matrix(I*d_I)
  dIv <- as.matrix(Ev*d_E) - as.matrix(Iv*d_I)
  dIx <- as.matrix(Ex*d_E) - as.matrix(Ix*d_I)

  dR  <- as.matrix(I*d_I*(1-IFR))
  dRv <- as.matrix(Iv*d_I*(1-IFR)*v_p)
  dRx <- as.matrix(Ix*d_I*(1-IFR)*v_p)

  dD  <- as.matrix(I*d_I*IFR + Iv*d_I*IFR*(1-v_p) + Ix*d_I*IFR*(1-v_p))

  out <- c(dS,dSv,dSx,dE,dEv,dEx,dI,dIv,dIx,dR,dRv,dRx,dD)
  list(out)
}



compute_R0 = function(u, C, n=num_groups){
  gamma <- 1/5 # recovery period (I -> R)
  # Davies NGM
  Du <- diag(u, n)
  Dy <- diag(1/gamma, n)
  
  NGM <- Du %*% C %*% Dy
  R0  <- abs(eigen(NGM)$values[1])
  return(R0)
}


construct_C_from_prem = function(home, work, school, other, u,target_R0, in_school=TRUE, alpha_factor=0.0){
  # Construct C matrix from the four PREM matrices
  # achieve a target R0 by playing with the weights of each matrix
  # if !in_school, then no school contacts
  # essential work contacts are assumed to be as "normal"
  # non-essential work contacts occur at rate alpha_factor of normal (fixed)
  #------------------------

  # is school in or not
  school <- in_school*school 

  # eliminate non-essential work
  work[1:9,1:9] <- alpha_factor*work[1:9,1:9]

  C <- school+home+work+other
  return(C*target_R0/compute_R0(u,C))
}




