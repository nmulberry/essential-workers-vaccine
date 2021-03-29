
# Generate Contact Matrix with Essential Workers
# p_ess: The proportion of essential workers by age group
# age_demo_by_fives: Age demographics by 5-yr age bins
# Reads in the contact matrices :
#     "mu_home.csv", "mu_work.csv", "mu_school.csv", "mu_other.csv"
# Saves the contact matrices:
#      "mu_home_essential.csv", "mu_work_essential.csv", "mu_school_essential.csv", "mu_other_essential.csv"
generate_contact_matrix = function(p_ess, age_demo_by_fives){
  age_demo <- age_toTens(age_demo_by_five)
  pop_total <- sum(age_demo)

  # Re-organize the demographics
  num_per_age <- age_demo[1:9]
  num_essential <- num_per_age*p_ess
  new_age_demo <-c(num_per_age - num_essential, #num wfh
    num_essential[3:8], #num essential (20-79)
    pop_total*pop_total)/pop_total 

  # proportion of each grp over its age grp
  prop_by_grp <- c(1-num_essential/num_per_age, num_essential[3:8]/num_per_age[3:8])
  num_per_grp <- head(new_age_demo,-1)*pop_total


  # LOAD IN PREM CONTACT MATRICES
  mu_home <- as.matrix(read.csv(paste0(PATH, "data/mu_home.csv"), sep=",", header=FALSE))
  mu_work <- as.matrix(read.csv(paste0(PATH, "data/mu_work.csv"), sep=",", header=FALSE))
  mu_school <- as.matrix(read.csv(paste0(PATH,"data/mu_school.csv"), sep=",", header=FALSE))
  mu_other <- as.matrix(read.csv(paste0(PATH,"data/mu_other.csv"), sep=",", header=FALSE))


  # put each matrix into 10 yr age bins using BC demographics
  # add 80+ compartment according to Bubar model
  mu_home <- add_80bin(matrix_toTens(mu_home, age_demo_by_five))
  mu_work <- add_80bin(matrix_toTens(mu_work, age_demo_by_five))
  mu_school <- add_80bin(matrix_toTens(mu_school, age_demo_by_five))
  mu_other <- add_80bin(matrix_toTens(mu_other, age_demo_by_five))


  # Look at est. total num contacts & enforce symmetry here
  mu_home_total <- 0.5*(mu_home*num_per_age + t(mu_home*num_per_age))
  mu_work_total <- 0.5*(mu_work*num_per_age + t(mu_work*num_per_age))
  mu_school_total <- 0.5*(mu_school*num_per_age + t(mu_school*num_per_age))
  mu_other_total <- 0.5*(mu_other*num_per_age + t(mu_other*num_per_age))


  # Add essential worker comparments
  mu_home_essential <- matrix_addEssential(mu_home_total, prop_by_grp)/num_per_grp
  mu_school_essential <- matrix_addEssential(mu_school_total, prop_by_grp)/num_per_grp
  mu_other_essential <- matrix_addEssential(mu_other_total, prop_by_grp)/num_per_grp
  mu_work_essential <- matrix_addEssential(mu_work_total, prop_by_grp)

  # For 70-79, and 80+..move all work compartments to essential 
  # also fudge
  mu_work_essential['70-79', 10:14] <- 1*(mu_work_essential['70-79', 10:14] + mu_work_essential['70-79', 3:7])
  mu_work_essential[10:14,'70-79'] <- 1*(mu_work_essential[10:14,'70-79'] + mu_work_essential[3:7,'70-79'])
  mu_work_essential['70-79', 3:7] <- 0.0
  mu_work_essential[3:7,'70-79'] <- 0.0
  mu_work_essential['80+', 10:15] <- 1*(mu_work_essential['80+', 10:15] + mu_work_essential['80+', 3:8])
  mu_work_essential[10:15,'80+'] <- 1*(mu_work_essential[10:15,'80+'] + mu_work_essential[3:8,'80+'])
  mu_work_essential['80+', 3:8] <- 0.0
  mu_work_essential[3:8,'80+'] <- 0.0

  mu_work_essential <- mu_work_essential/num_per_grp
  # SAVE
  saveRDS(mu_home_essential, paste0(PATH, "generated-data/mu_home_essential.rds"))
  saveRDS(mu_work_essential, paste0(PATH, "generated-data/mu_work_essential.rds"))
  saveRDS(mu_school_essential, paste0(PATH,"generated-data/mu_school_essential.rds"))
  saveRDS(mu_other_essential, paste0(PATH, "generated-data/mu_other_essential.rds"))
  saveRDS(new_age_demo, paste0(PATH,"generated-data/age_demographics_essential.rds"))

}

generate_contact_matrix_fsa = function(p_ess_fsa,pop_total_fsa, age_demo_by_fives){
p_ess= p_ess_fsa
    age_demo <- age_toTens(age_demo_by_five)
  pop_total_on <- sum(age_demo)
 age_demo = round(age_demo*pop_total_fsa/pop_total_on)
 pop_total = sum(age_demo)
  
  # Re-organize the demographics
  num_per_age <- age_demo[1:9]
  num_essential <- num_per_age*p_ess
  new_age_demo <-c(num_per_age - num_essential, #num wfh
                   num_essential[3:8], #num essential (20-79)
                   pop_total*pop_total)/pop_total 
  
  # proportion of each grp over its age grp
  prop_by_grp <- c(1-num_essential/num_per_age, num_essential[3:8]/num_per_age[3:8])
  num_per_grp <- head(new_age_demo,-1)*pop_total
  
  
  # LOAD IN PREM CONTACT MATRICES
  mu_home <- as.matrix(read.csv(paste0(PATH,"data/mu_home.csv"), sep=",", header=FALSE))
  mu_work <- as.matrix(read.csv(paste0(PATH,"data/mu_work.csv"), sep=",", header=FALSE))
  mu_school <- as.matrix(read.csv(paste0(PATH,"data/mu_school.csv"), sep=",", header=FALSE))
  mu_other <- as.matrix(read.csv(paste0(PATH,"data/mu_other.csv"), sep=",", header=FALSE))
  
  
  # put each matrix into 10 yr age bins using BC demographics
  # add 80+ compartment according to Bubar model
  mu_home <- add_80bin(matrix_toTens(mu_home, age_demo_by_five))
  mu_work <- add_80bin(matrix_toTens(mu_work, age_demo_by_five))
  mu_school <- add_80bin(matrix_toTens(mu_school, age_demo_by_five))
  mu_other <- add_80bin(matrix_toTens(mu_other, age_demo_by_five))
  
  
  # Look at est. total num contacts & enforce symmetry here
  mu_home_total <- 0.5*(mu_home*num_per_age + t(mu_home*num_per_age))
  mu_work_total <- 0.5*(mu_work*num_per_age + t(mu_work*num_per_age))
  mu_school_total <- 0.5*(mu_school*num_per_age + t(mu_school*num_per_age))
  mu_other_total <- 0.5*(mu_other*num_per_age + t(mu_other*num_per_age))
  
  
  # Add essential worker comparments
  mu_home_essential <- matrix_addEssential(mu_home_total, prop_by_grp)/num_per_grp
  mu_school_essential <- matrix_addEssential(mu_school_total, prop_by_grp)/num_per_grp
  mu_other_essential <- matrix_addEssential(mu_other_total, prop_by_grp)/num_per_grp
  mu_work_essential <- matrix_addEssential(mu_work_total, prop_by_grp)
  
  # For 70-79, and 80+..move all work compartments to essential 
  # also fudge
  mu_work_essential['70-79', 10:14] <- 1*(mu_work_essential['70-79', 10:14] + mu_work_essential['70-79', 3:7])
  mu_work_essential[10:14,'70-79'] <- 1*(mu_work_essential[10:14,'70-79'] + mu_work_essential[3:7,'70-79'])
  mu_work_essential['70-79', 3:7] <- 0.0
  mu_work_essential[3:7,'70-79'] <- 0.0
  mu_work_essential['80+', 10:15] <- 1*(mu_work_essential['80+', 10:15] + mu_work_essential['80+', 3:8])
  mu_work_essential[10:15,'80+'] <- 1*(mu_work_essential[10:15,'80+'] + mu_work_essential[3:8,'80+'])
  mu_work_essential['80+', 3:8] <- 0.0
  mu_work_essential[3:8,'80+'] <- 0.0
  
  mu_work_essential <- mu_work_essential/num_per_grp

  return(list(mu_home = mu_home_essential,
              mu_work = mu_work_essential,
              mu_school = mu_school_essential, 
              mu_other = mu_other_essential,
              age_demo = new_age_demo))
}





age_toTens = function(x){
  # helper func to convert to tens
  res <- rep(0, 9)
  j <- 1
  for (i in 1:9*2){
    res[j] <- x[i]+x[i-1]
    j <- j+1
  }
  # get 80+ bin
  res[9] <- res[9] + x[19]

  return(res)
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

matrix_addEssential = function(C, p_grp){
  # C should be by tens, with 80+ (9x9) 
  # C(i,j) should be total number of contacts between age classes i & j
  # evenly distribute contacts according to demog data
  # p_grp is the prop. of the total age grp over all grps
  C_e <- matrix(nrow=15, ncol=15,0)
  colnames(C_e) <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", 
    "60-69", "70-79", "80+", "20-29e", "30-39e", "40-49e", "50-59e", "60-69e", "70-79e")
  rownames(C_e) <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", 
    "60-69", "70-79","80+", "20-29e", "30-39e", "40-49e", "50-59e", "60-69e", "70-79e")
  # Now copy & extend C
  C_e[1:9,1:9] <- C
  C_e[10:15,1:9] <- C[3:8,]
  C_e[,10:15] <- C_e[,3:8]
  # Now scale by p_grp to allocate contacts evenly
  C_e <- t(t(C_e)*p_grp)*p_grp
  
  return(C_e)
}


matrix_toTens <- function(C_byfives, pop_demo) {
  # Modified from Bubar et al.
  # https://github.com/kbubar/vaccine_prioritization/tree/v1 
  l <- dim(C_byfives)[1]
  C_bytens <- matrix(nrow = l/2, ncol = l/2)

  col_count <- 1
  for (col in seq(1,l,by = 2)){
    row_count <- 1
    for (row in seq(1,l, by = 2)){
      p1 <- pop_demo[row]
      p2 <- pop_demo[row + 1]
      C_bytens[row_count, col_count] <- ((C_byfives[row, col] + C_byfives[row, col + 1])*p1 +
        (C_byfives[row + 1, col] + C_byfives[row + 1, col + 1])*p2)/(p1+p2)
      row_count <- row_count + 1
    }
    col_count <- col_count + 1
  }
  colnames(C_bytens) <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79")
  rownames(C_bytens) <- colnames(C_bytens)

  return(C_bytens)
}

add_80bin <- function(C_bytens){
  # Modified from Bubar et al.
  # https://github.com/kbubar/vaccine_prioritization/tree/v1 
  # INPUT:
  #  C_bytens: Contact matrix with 10 year age bins, with last bin = 70-79
  #
  # OUTPUT:
  # C_new: Same C matrix with new row & col for 80+ age bin

  bin_70 <- dim(C_bytens)[1]
  bin_80 <- bin_70 + 1
  bin_60 <- bin_70 - 1

  C_new <- matrix(nrow = bin_80, ncol = bin_80)
  colnames(C_new) <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
  rownames(C_new) <- colnames(C_new)

  # fill rows for 80+ with same data from 70-79
  C_new[1:bin_70, 1:bin_70] <- C_bytens
  C_new[2:bin_80, bin_80] <- C_bytens[1:bin_70, bin_70]
  C_new[bin_80, 2:bin_80] <- C_bytens[bin_70, 1:bin_70]

  C_new[1, bin_80] <- C_new[1, bin_70]
  C_new[bin_80, 1] <- C_new[bin_70, 1]

  # Assumption: 80+ have similar (but less) contact rates as 70-79, with increased contact with 70-80+ (long term living facilities)
  # Implementation: Decrease contact for bins 0 - 69 for 80+  by 10% then split these contacts between 70-79 & 80+

  # col 80+
  to_decrease <- 0.1 * C_new[1:bin_60, bin_80]
  C_new[1:bin_60, bin_80] <- 0.9 * C_new[1:bin_60, bin_80]

  to_increase <- sum(to_decrease)/2
  C_new[bin_70, bin_80] <- C_new[bin_70, bin_80] + to_increase
  C_new[bin_80, bin_80] <- C_new[bin_80, bin_80] + to_increase

  # row 80+
  to_decrease <- 0.1 * C_new[bin_80, 1:bin_60]
  C_new[bin_80, 1:bin_60] <- 0.9 * C_new[bin_80, 1:bin_60]

  to_increase <- sum(to_decrease)/2
  C_new[bin_80, bin_70] <- C_new[bin_80, bin_70] + to_increase
  C_new[bin_80, bin_80] <- C_new[bin_80, bin_80] + to_increase

  # fudge a little bit
  C_new[bin_80, bin_80] <- 3*C_new[bin_80, bin_80]
  C_new[bin_70, bin_70] <- 2*C_new[bin_70, bin_70]
  C_new[bin_80, bin_70] <- 2*C_new[bin_80, bin_70]
  C_new[bin_70, bin_80] <- 2*C_new[bin_70, bin_80]
  return(C_new)
}



