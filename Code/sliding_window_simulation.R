###############################################################################################################
#
# This script simulates refill histories and performs longitudinal k-means clustering for different scenarios.
# Parameters to specify are the sample sizes, number of sets per sample size, window sizes and step length.
# The script uses parallel processing with the foreach package.
# In the default version, it will use a SNOW cluster (with the package doSNOW) for compatibility with Windows,
# but it can easily be customized to run with other parallel backends as well.
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

library(plyr)
library(AdhereR)
library(data.table)
library(clues)
library(kml)
library(evobiR)
library(foreach)
library(doSNOW)
#library(doParallel)

source("functions.R")

#set parameters for foreach
sample.sizes <- c(1000) #the sample sizes to run the simulation on
n.sets <- 100 #the number of iterations to run per sample size
window.sizes <- c(7, 14, seq(30, duration.total, 30) #the window sizes to analyze per LCMA method
overlaps <- seq(0.1,1,0.1) #the step lengths to analyze, as fraction of the window size per LCME method

#set initial parameters for simulation
sample.size <- 0
start.date = as.Date("01.01.2022", format = "%d.%m.%Y")
duration.total <- 2*365 #the (minimum) duration of the simulated refills, this will be used as observation window

#initialize clusters; otherwise kml does not work within loops
conv_clust <- NULL

#create parallel backend, here for a SNOW cluster
clus <- makeCluster(max(1, detectCores() -2))
registerDoSNOW(clus)

#export necessary objects and packages
clusterExport(clus,c("sample.sizes",
                     "n.sets",
					 "window.sizes",
					 "overlaps",
					 "sample.size",
                     "start.date",
                     "duration.total",
                     "conv_clust",
                     "logistic",
                     "med.events.sample",
                     "clustering",
					 "CMA.dichot"))
clusterEvalQ(clus, { library(data.table); library(kml); library(evobiR); library(clues); library(AdhereR); library(plyr); library(foreach); library(doParallel) })

#run simulation --> running this for loop will run the simulation for the specified number of samples, which could take a lot of time!

results <- foreach(n = 1:(length(sample.sizes)*n.sets),
                   .verbose = TRUE,
                   .errorhandling = "pass") %dopar% {
                     
  set.seed(n)
                     
  #set/update parameters
  sample.size <- sample.sizes[ceiling(n/n.sets)]
  
  #create random proportions of dispensing durations summing to 1
  num <- sample(0:10, 3, replace = TRUE)
  if(sum(num) == 0) num <- sample(c(1,0,0))
  disp.durations <- c(30, 60, 90)
  dist.durations <- num/sum(num)
  mean.duration <- sum(disp.durations*dist.durations)
  
  #create random proportions of group sizes summing to 0.8  
  num <- sample(1:3, 4, replace = TRUE)
  dist_groups <- num/sum(num)-0.05
  dist <- c(0.1, dist_groups, 0.1)
  
  #create sample
  sample <- med.events.sample(ntot=sample.size,
                              start.date = start.date,
                              tot.duration = duration.total,
                              disp.durations = disp.durations,
                              dist.durations = dist.durations,
                              dist = dist)
  

  sample$GROUP <- as.factor(sample$GROUP)

  setkey(sample, PATIENT_ID)

  unique <- unique(sample[,.(PATIENT_ID, GROUP)])

  #calculate gaps, with or without carryover
  observation_window <- duration.total

  df_sample <- as.data.frame(sample)

  gap_days <- compute.event.int.gaps(data = df_sample,
                                     ID.colname = "PATIENT_ID",
                                     event.date.colname = "DATE",
                                     event.duration.colname = "DURATION",
                                     observation.window.duration = observation_window,
                                     followup.window.duration = observation_window*2,
                                     remove.events.outside.followup.window = FALSE,
                                     keep.event.interval.for.all.events = TRUE,
                                     keep.window.start.end.dates = TRUE,
                                     return.data.table = TRUE,
                                     carryover.within.obs.window = FALSE,
                                     carryover.into.obs.window = FALSE)

  #recalculate interval and gap days for last event
  gap_days[!is.na(.TDIFF2), event.interval.2:= .TDIFF2]
  gap_days <- gap_days[DATE < start.date+duration.total]
  gap_days[,gap.days.2:= event.interval.2 - (DURATION + .CARRY.OVER.FROM.BEFORE)]
  gap_days[gap.days.2 < 0, gap.days.2 := 0]

  #calculate availability for each event interval
  gap_days[,CMA_carryover := DURATION/(DURATION + gap.days.2)]
  gap_days[is.na(CMA_carryover),CMA_carryover := DURATION/(DURATION + gap.days)]

  gap_days <- gap_days[,c(1:6, 19:21), with = FALSE]
  gap_days[,CMA_no_carryover := DURATION/event.interval.2]
  gap_days[is.na(CMA_no_carryover),CMA_no_carryover := DURATION/event.interval]

## create matrix with probabilities of availability for each day

  # create a vector with all dates for the total duration
  Dates <- seq.Date(start.date, by = 1, length.out = duration.total)

  # Option 1: fill with CMA for interval between two refills
  DT <- gap_days[rep(1:.N,event.interval)][,Indx:=1:.N,by=.(PATIENT_ID, DATE)]

  # Option 2: fill with 1 for the first consecutive period per interval that is covered by the duration

  CMA_dichot <- unlist(mapply(CMA.dichot,
                              .duration = gap_days[,DURATION],
                              .event.interval = gap_days[,event.interval]))

  DT$CMA_dichot <- CMA_dichot

  DT$DATE <- rep(Dates, times = nrow(unique))

  #calculate default aRand value for simple k-means clustering based on average adherence (CMA9)
  CMA9 <- CMA9(data = sample,
               ID.colname = "PATIENT_ID",
               event.date.colname = "DATE",
               event.duration.colname = "DURATION",
               observation.window.duration = observation_window,
               followup.window.duration = observation_window*2)

  CMA9_sample <- getCMA(CMA9)

  CMA9_clusters <- kmeans(CMA9_sample$CMA, 6)
  CMA9_groups <- data.table(sample_id = n,
                            GROUP = unique$GROUP,
                            clus = CMA9_clusters$cluster)
  CMA9_aRand <- adjustedRand(CMA9_clusters$cluster, as.numeric(unique$GROUP), randMethod = "HA")
  
  #create correlation matrix to relabel groups
  cm <- as.matrix(table(Predicted = CMA9_groups$clus,
                        Actual = CMA9_groups$GROUP))

  CMA9_groups$clus <- mapvalues(CMA9_groups$clus,
                                from = 1:6,
                                to = max.col(cm))

  ### the actual partitioning with kml

  # run kml on all combinations of window sizes and lag times.
  # This is the part that takes most of the time
  test <- foreach(window.size = window.sizes
  ),
  .combine = rbind,
  .verbose = TRUE,
  .errorhandling = "pass") %do% {
    
    do.call(rbind, lapply(overlaps, FUN = function(i) {
      
      test <- tryCatch(
        {
          #create sliding windows using LCMA2 without carryover
          sw_1 <- DT[,SlidingWindow(mean, CMA_no_carryover,
                                    window.size, floor(window.size*i)), by = .(PATIENT_ID, GROUP)]
          #dichotomized version of LCMA2
          sw_1[V1 >= 0.8, LCMA2_threshold := 1][V1 < 0.8, LCMA2_threshold := 0]
          
          #create sliding windows using LCMA1
          sw_2 <- DT[,SlidingWindow(mean, CMA_dichot,
                                    window.size, floor(window.size*i)), by = .(PATIENT_ID, GROUP)]
          #dichotomized version of LCMA1
          sw_2[V1 >= 0.8, LCMA1_threshold := 1][V1 < 0.8, LCMA1_threshold := 0]
          
          sw <- cbind(sw_1, sw_2[,3:4])
          setnames(sw, c("PATIENT_ID",
                         "GROUP",
                         "LCMA2",
                         "LCMA2_threshold",
                         "LCMA1",
                         "LCMA1_threshold"))
          
          #create a table with as many rows as days per patient
          sw[,window_ID:=rep(1:(nrow(sw)/nrow(unique)), times = nrow(unique))]
          
          #use kml
          
          LCMA2 <- clustering(sw = sw, CMA = "LCMA2", size = nrow(unique))
          LCMA2_threshold <- clustering(sw = sw, CMA = "LCMA2_threshold", size = nrow(unique))
          LCMA1 <- clustering(sw = sw, CMA = "LCMA1", size = nrow(unique))
          LCMA1_threshold <- clustering(sw = sw, CMA = "LCMA1_threshold", size = nrow(unique))
          
          matrix_traj <- cbind.data.frame(sample_id = n,
                                          win_size = window.size,
                                          mult_step_leng = i,
                                          PATIENT_ID = unique$PATIENT_ID,
                                          GROUP = unique$GROUP,
                                          clus_CMA9 = CMA9_groups$clus,
                                          traj_LCMA2 = LCMA2$traj,
                                          traj_LCMA2_threshold = LCMA2_threshold$traj,
                                          traj_LCMA1 = LCMA1$traj,
                                          traj_LCMA1_threshold = LCMA1_threshold$traj)
          
          res <- c(sample_id = n,
                   win_size = window.size,
                   mean_duration = mean.duration,
                   mean_gap = mean(DT$gap.days.2, na.rm = TRUE),
                   step_leng = floor(window.size*i),
                   mult_step_leng = i,
                   aRand_LCMA1 = LCMA1$aRand,
                   time_LCMA1 = LCMA1$time,
                   aRand_LCMA1_threshold = LCMA1_threshold$aRand,
                   time_LCMA1_threshold = LCMA1_threshold$time,
                   aRand_LCMA2 = LCMA2$aRand,
                   time_LCMA2 = LCMA2$time,
                   aRand_LCMA2_threshold = LCMA2_threshold$aRand,
                   time_LCMA2_threshold = LCMA2_threshold$time,
                   CMA9_aRand = CMA9_aRand,
                   Duration_Var = var(sample$DURATION),
                   Duration_SD = sd(sample$DURATION),
                   sample_size = nrow(unique))
          
          
          list(results = res, trajs = matrix_traj)
        },
        error = function(e){
          matrix_traj <- cbind.data.frame(sample_id = n,
                                          win_size = window.size,
                                          mult_step_leng = i,
                                          PATIENT_ID = unique$PATIENT_ID,
                                          GROUP = unique$GROUP,
                                          clus_CMA9 = CMA9_groups$clus,
                                          traj_LCMA2 = NA,
                                          traj_LCMA2_threshold = NA,
                                          traj_LCMA1 = NA,
                                          traj_LCMA1_threshold = NA)
          
          res <- c(sample_id = n,
                   win_size = window.size,
                   mean_duration = mean.duration,
                   mean_gap = mean(DT$gap.days.2, na.rm = TRUE),
                   step_leng = floor(window.size*i),
                   mult_step_leng = i,
                   aRand_LCMA2 = NA,
                   time_LCMA2 = NA,
                   aRand_LCMA2_threshold = NA,
                   time_LCMA2_threshold = NA,
                   aRand_LCMA1 = NA,
                   time_LCMA1 = NA,
                   aRand_LCMA1_threshold = NA,
                   time_LCMA1_threshold = NA,
                   CMA9_aRand = CMA9_aRand,
                   Duration_Var = var(sample$DURATION),
                   Duration_SD = sd(sample$DURATION),
                   sample_size = nrow(unique))
          
          list(results = res, trajs = matrix_traj)
          
        }
      )
      
      
      test
    }))
    
  }
  test
  
                   }

snow::stopCluster(clus)

