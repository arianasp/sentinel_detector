#First pass at a sentinel detector
#
#The detector uses a combination of ACC-based behavioral states with audio labels to (attempt to) detect bouts of sentinel behavior
#
#For each time point, we look in a time window t_window (currently 31 seconds) centered on that time point
#and count up the fraction of the seconds where that individual was labeled as vigilant.
#We also look within the same time window and count up how many (focal) sn's and cc's the individual produced, 
#then normalize this to calculate the call rate for each call type.
#
#We then define an individual as being in the "sentinel" state at time t if:
#1. The fraction of time spent vigilant within the time window exceeds 0.95 (vig_frac_thresh)
#2. The individual's call rate for sn's in the time window exceeds min_sn_thresh (here set to 0 - so the individual has to give at least 1 short note)
#3. The individual's call rate for cc's in the time window falls below max_cc_thresh (here set to 0.001 meaning no cc's are allowed)
#
#Once we have cateogrized each time point at being in the "sentinel" state or not, we 
#find contiguous bouts of sentinel states for each individual. 
#The resulting bouts are stored in the data frame sentinel_bouts, with onset and offset times (t0 and tf)
#correspoding to the beginning of the first time window and the end of the last time window 
#
#Finally, we connect adjacent bouts to one another if they are less than max_time_connect_bouts time apart

#-----PARAMS-----
#input data file
infile <- '/mnt/EAS_shared/meerkat/working/processed/combined/ZU_2021_allstreams.RData'

#where to save the output table
outfile <- '/mnt/EAS_ind/astrandburg/sentinel_detector/first_pass_sentinel_bouts_2025-03-18.csv'

#detecting whether an individual is a sentinel at any given time point
vig_frac_thresh <- 0.95 #minimum fraction of the time window that has to be designated as vigilant for the individual to be able to be sentinel
max_cc_thresh <- 0.001 #maximum fraction of cc's within a time window to allow the individual to be designated as sentinel
min_sn_thresh <- 0 #minimum fraction of sn's within a time window to consider an individual as sentinel

#connecting nearby bouts of sentinel into one
max_time_connect_bouts <- 120 #adjacent bouts of sentinel will be combined if they are within this time period of one another

#-----LOAD DATA-----
load(infile)

#-----BASIC INFO-----
n_inds <- nrow(xs)
n_times <- ncol(xs)
n_days <- length(breaks)
breaks <- c(breaks, n_times+1)

#-----PROCESS-----

#define each individual at each time point as sentinel-like or not, according to the rules described above
sentinel <- behavs_array_smoothed[,,'Vigilance'] > vig_frac_thresh &
  calls_array_smoothed[,,'cc'] <= max_cc_thresh &
  calls_array_smoothed[,,'sn'] > min_sn_thresh

#for each individual on each day, get bouts of contiguous sentinel-like behavior
sentinel_bouts <- data.frame()
dt <- (params$time_window - 1) / 2 #smoothing window
for(i in 1:n_inds){
  for(j in 1:n_days){
    t_idxs <- seq(breaks[j],(breaks[j+1]-1))
    if(sum(!is.na(sentinel[i,t_idxs]))>0){
      sent_runs <- rle(sentinel[i,t_idxs])
      curr_t_idx <- t_idxs[1]
      for(k in 1:length(sent_runs$lengths)){
        if(!is.na(sent_runs$values[k])){
          if(sent_runs$values[k]==T){
            sentinel_bouts <- rbind(sentinel_bouts, data.frame(ind = i, t0 = curr_t_idx - dt, tf = curr_t_idx + sent_runs$lengths[k] + dt - 1))
          }
        }
        curr_t_idx <- curr_t_idx + sent_runs$lengths[k]
      }
    }
  }
}

#calculate some information about the bouts

#get number of calls given during the bouts
sentinel_bouts[,params$call_types] <- NA
for(i in 1:nrow(sentinel_bouts)){
  t_idxs <- seq(sentinel_bouts$t0[i], sentinel_bouts$tf[i])
  ind <- sentinel_bouts$ind[i]
  for(c in params$call_types){
    sentinel_bouts[i,c] <- sum(calls_array[ind, t_idxs, c],na.rm=T)
  }
}

#some additional basic information
sentinel_bouts$id_code <- ids$id_code[sentinel_bouts$ind]
sentinel_bouts$t0_UTC <- timestamps[sentinel_bouts$t0]
sentinel_bouts$tf_UTC <- timestamps[sentinel_bouts$tf]
sentinel_bouts$date <- lubridate::date(sentinel_bouts$t0_UTC)

#match back to the original audio file to enable manual checking
sentinel_bouts$ind_date <- paste0(sentinel_bouts$id_code,'_',sentinel_bouts$date)
calls$ind_date <- paste0(calls$ind, '_', calls$date)
sentinel_bouts$wav <- calls$wavFileName[match(sentinel_bouts$ind_date, calls$ind_date)]

#find the time of the first call during the bout, in the audio file
#then use this to calculate the start time of the bout of sentinel-like behavior
calls$t0_UTC <- as.POSIXct(calls$t0GPS_UTC, tz = 'UTC')
sentinel_bouts$tf_file <- sentinel_bouts$t0_file <- seconds_to_period(0)
for(i in 1:nrow(sentinel_bouts)){
  wav <- sentinel_bouts$wav[i]
  t0 <- sentinel_bouts$t0_UTC[i]
  first_call <- calls[which(calls$t0_UTC >= t0 & calls$wavFileName == wav)[1],]
  first_call_time <- first_call$t0File
  first_call_time_sec <- cocomo::parse_audition_time(first_call_time)
  time_diff_to_bout_start <- as.numeric(difftime(first_call$t0_UTC, t0, units = 'secs'))
  bout_start_time_in_wav_file <- first_call_time_sec - time_diff_to_bout_start
  sentinel_bouts$t0_file[i] <- lubridate::seconds_to_period(bout_start_time_in_wav_file)
  sentinel_bouts$tf_file[i] <- lubridate::seconds_to_period(bout_start_time_in_wav_file + sentinel_bouts$tf[i] - sentinel_bouts$t0[i])
}

#connect nearby bouts for each individual
inds <- unique(sentinel_bouts$id_code)
sentinel_bouts_connected <- data.frame()
for(i in 1:length(inds)){
  sentinel_bouts_ind <- sentinel_bouts[which(sentinel_bouts$id_code==inds[i]),]
  sentinel_bouts_ind <- sentinel_bouts_ind[order(sentinel_bouts_ind$t0),]
  idx <- 1
  curr_time <- sentinel_bouts_ind$tf[1]
  for(j in 1:nrow(sentinel_bouts_ind)){
    if((sentinel_bouts_ind$t0[j] - curr_time) < max_time_connect_bouts){
      curr_time <- sentinel_bouts_ind$tf[j]
    } else{
      row <- sentinel_bouts_ind[idx,]
      row$tf <- sentinel_bouts_ind$tf[j-1]
      row$tf_UTC <- sentinel_bouts_ind$tf_UTC[j-1]
      row$tf_file <- sentinel_bouts_ind$tf_file[j-1]
      idx <- j
      curr_time <- sentinel_bouts_ind$tf[j]
      sentinel_bouts_connected <- rbind(sentinel_bouts_connected, row)
    }
    
  }
  
}

#save output
write.csv(sentinel_bouts_connected, file = outfile, quote = F, row.names = F)
