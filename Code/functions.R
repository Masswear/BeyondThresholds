###############################################################################################################
#
# This script contains the custom functions to generate simulated refill histories with 6 distinct patterns.
# The distinct patterns are achieved by variations in the delays of refill events. The parameters for each
# group are defined in the offset function, which is encapsulated in the med.events.sample function.
#
# In addition, it contains the clustering-function, which provides error handling for kml to make sure that the 
# simulation doesn't break if cluster analysis fails and the CMA.dichot-function, which is a helper function to
# dichotomize the LCMA estimates
#
# Copyright (C) 2018 Samuel Allemann
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License v3
# as published by the Free Software Foundation.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, you can find it here:
# https://www.gnu.org/licenses/gpl-3.0.en.html
#
#################################################################################################################


## function to create a vector of numbers following a logistic function - used to create groups 3 and 4
logistic <- function(x, L=0, S, D, h=1, B = expression(x-D)) {
  y <- lapply(x, FUN = function(x) (h-(5*L))/(1+exp(1)^(S*eval(B))) + L)
  unlist(y)
}

## functions to create refill patterns for different groups
med.events.sample <- function(ntot,
                              start.date = "01.01.2022",
                              tot.duration = 2*365,
                              disp.durations = c(30, 60, 90),
                              dist.durations = c(0.3, 0.5, 0.2),
                              dist = c(0.1, 0.2, 0.2, 0.2, 0.2, 0.1)){
  #function to create offset for different groups
  offset <- function(group){
    
    # Group 1: “High adherence” with an average CMA9 of around 95%.
    # Individuals in this group continuously refill with a normally distributed delay 
    if(group == 1){
      L <- -0.1
      U <- 0.2
      
      m <- 0.05
      s <- 0.1
      
      pL <- pnorm(L, mean=m, sd=s)
      pU <- pnorm(U, mean=m, sd=s)
      
      offset <- qnorm(runif(n, pL, pU), mean=m, sd=s)
    }
    
    # Group 2: “Erratic adherence” with a median CMA9 between 50% and 90%.
    # Individuals in this group continuously refill with a normally distributed delay
    if(group == 2){
      L <- -0.2
      U <- 1.2
      
      m <- 0
      s <- 1
      
      pL <- pnorm(L, mean=m, sd=s)
      pU <- pnorm(U, mean=m, sd=s)
      
      offset <- qnorm(runif(n, pL, pU), mean=m, sd=s)
    }
    
    # Group 3: “Gradual decline” with increasingly delayed refills.
    # Delays increase linearly for each refill by a factor of 2/(total number of refills) with added noise 
    if(group == 3){
      L <- 0.5
      U <- 1.5
      
      m <- 1
      s <- 1
      
      pL <- pnorm(L, mean=m, sd=s)
      pU <- pnorm(U, mean=m, sd=s)
      
      offset1 <- qnorm(runif(n, pL, pU), mean=m, sd=s)
      
      t <- 2/n * 1:n
      
      offset <- t * offset1
    }
    
    # Group 4: “Intermittent adherence” with a change between high and low adherence at regular intervals.
    # Delays follow a logistic function with a sinus term over time with added noise
    if(group == 4){
      L <- 0.8
      U <- 1.2
      
      m <- 1
      s <- 0.1
      
      pL <- pnorm(L, mean=m, sd=s)
      pU <- pnorm(U, mean=m, sd=s)
      
      offset1 <- qnorm(runif(n, pL, pU), mean=m, sd=s)
      
      # D <- sample((n/6):(n/4+1), 1)
      
      offset <- logistic(x=1:n, L = 0.05, S = 10, D = n, B = expression(sin(2*x-D))) * offset1
    }
    
    # Group 5: “Partial drop-off” with high adherence initially and partial drop after some time.
    # Delays follow a logistic function over time with added noise identical to Group 3.
    if(group == 5){
      L <- 0.5
      U <- 1.5
      
      m <- 1
      s <- 1
      
      pL <- pnorm(L, mean=m, sd=s)
      pU <- pnorm(U, mean=m, sd=s)
      
      offset1 <- qnorm(runif(n, pL, pU), mean=m, sd=s)
      
      offset <- logistic(x=1:n, L = 0.05, S = -15, D = n/3) * offset1
    }
    
    # Group 6: “Non-persistence” with one or two refills after the initial fill and no refills afterwards.
    # Delays are normally distributed.
    if(group == 6){
      n <- sample(2:3, 1)
      L <- -0.2
      U <- 0.8
      
      m <- 0.3
      s <- 1
      
      pL <- pnorm(L, mean=m, sd=s)
      pU <- pnorm(U, mean=m, sd=s)
      
      offset <- qnorm(runif(n, pL, pU), mean=m, sd=s)
    }
    
    offset
  }
  
  #function to create refill pattern for a number of patients per group
  refills <- function(x, group){
    initial_fill <- 30
    
    offsets <- offset(group = group)
    
    durations <- sample(disp.durations, size = length(offsets)-1, replace = TRUE, prob = dist.durations)
    durations <- c(initial_fill, durations)
    
    date <- as.Date(start.date, format = "%d.%m.%Y")
    
    refill_dates <- c(date, date + cumsum(durations + round(offsets*durations, 0)))
    
    dt <- data.table(GROUP = group, PATIENT_ID = x, DATE = head(refill_dates, -1), DURATION = durations)
    
    dt
  }
  
  #initialize
  ID_last <- 0
  sample <- NULL
  
  #calculate mean duration of dispense
  mean.duration <- sum(disp.durations*dist.durations)
  
  #calculate number of required dispenses
  n <- ceiling((tot.duration/mean.duration) * 1.5)
  
  #for each group, create a sample based on the distribution given in dist
  for(i in 1:5){
    num_pat <- round(dist[i] * ntot, 0)
    ID_first <- ID_last + 1
    ID_last <- ID_first + num_pat - 1
    group <- rbindlist(lapply(ID_first:ID_last, refills, group = i))
    
    sample <- rbind(group, sample)
  }
  num_pat <- ntot-ID_last
  ID_first <- ID_last + 1
  ID_last <- ID_first + num_pat - 1
  group <- rbindlist(lapply(ID_first:ID_last, refills, group = 6))
  
  sample <- rbind(group, sample)
  sample
}

##function to compute clusters
clustering <- function(sw, CMA, size){

  #extract assigned groups
  traj <- tryCatch(
    { sw_wide <- dcast(sw, GROUP + PATIENT_ID ~ window_ID, value.var = CMA)
    conv_clust <- clusterLongData(idAll = sw_wide$PATIENT_ID,
                                  time = 1:(nrow(sw)/size),
                                  traj = sw_wide[,3:ncol(sw_wide)])
    ptm <- proc.time()
    kml(conv_clust, nbClusters = 6)
    t <- proc.time() - ptm
    
    #extract assigned groups
    matrix_traj <- cbind.data.frame(GROUP = as.numeric(sw_wide$GROUP),
                                    traj = as.numeric(getClusters(conv_clust, 6))
    )
    
    
    cm <- as.matrix(table(Predicted = matrix_traj$traj,
                          Actual = matrix_traj$GROUP))
    
    matrix_traj$traj <- mapvalues(matrix_traj$traj,
                                  from = 1:6,
                                  to = max.col(cm))
    
    list(traj = matrix_traj$traj,
         time = t[3],
         aRand = adjustedRand(matrix_traj$GROUP, matrix_traj$traj, randMethod = "HA"),
         Rand = adjustedRand(matrix_traj$GROUP, matrix_traj$traj, randMethod = "Rand")
         )
    },
    error = function(e){
      #message(e)
      list(traj = NA,
           time = NA,
           aRand = NA,
           Rand = NA)
      
    }
  )
  
  traj
  
}

##function to dichotomize LCMA
  CMA.dichot <- function(.duration, .event.interval) {
    if(.event.interval > .duration){
      c(rep(1, .duration), rep(0, .event.interval-.duration))
    } else rep(1, .event.interval)
  }
