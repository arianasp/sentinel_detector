#Second pass at a sentinel detector - audio only
#
#The detector uses audio labels and gps to (attempt to) detect bouts of sentinel behavior
#
#The detector first detects strings of sn-containing seconds that are separated by a maximum
#gap of max_sn_gap_len (=60 sec). It then filters to only sequences with at least min_sns_in_seq (=10)
#Next, it checks the gps data. It filters out any periods where the displacement over a minute goes over a threshold (uses a double threshold method to get contiguous bouts)
#and breaks any sequences that were interrupted by these high speed periods into multiple sequences.
#From the resulting sequences, it filters out those shorter than min_dur (=120 s)
#Finally, it filters out sequences that contain more than max_cc_rate (=0.01) cc-containing seconds

library(lubridate)

#-----PARAMS-----
#input data file
infile <- '~/Desktop/Meerkat_data_for_hackathon/processed_data/ZU_2021_allstreams.RData'

#where to save the output table
outfile <- '~/Desktop/Meerkat_data_for_hackathon/results/audio_only_sentinel_bouts_2025-03-21_gpssplit.csv'

#detecting whether an individual is a sentinel at any given time point
max_sn_gap_len <- 60
min_sns_in_seq <- 10
min_dur <- 120 #minimum (time) duration
max_cc_rate <- .01 #maximum rate of ccs 
speed_dt <- 60 #time window for computing speed (rolling window centered on timestep) - must be even
speed_thresh_upper <- 10 #threshold for determining if a period of movement is to be excluded (in m / min)
speed_thresh_lower <- 5 #threshold for finding the edges of periods of movement (in m / min)

#-----LOAD DATA-----
load(infile)

#-----BASIC INFO-----
n_inds <- nrow(xs)
n_times <- ncol(xs)
n_days <- length(breaks)
breaks <- c(breaks, n_times+1)

#-----PROCESS-----

#Start by finding sequences of short notes with a maximum gap between short notes
#of max_sn_gap_len
sn_seqs <- data.frame()
for(i in 1:n_inds){
  for(d in 1:n_days){
    t0 <- breaks[d]
    tf <- breaks[d + 1]-1
    #get time series of short notes
    sn_ts <- calls_array[i,t0:tf,'sn']
    
    #get indexes to times with short notes
    sn_idxs <- which(sn_ts > 0) + t0 - 1
    n_sns <- length(sn_idxs)
    if(n_sns>1){
      start_idx <- 1
      sn_count <- 1
      for(end_idx in 2:n_sns){
        if((sn_idxs[end_idx] - sn_idxs[end_idx-1]) > max_sn_gap_len){
          row <- data.frame(ind = i, t0 = sn_idxs[start_idx], tf = sn_idxs[end_idx-1], n = sn_count)
          start_idx <- end_idx
          sn_seqs <- rbind(sn_seqs, row)
          sn_count <- 1
        } else{
          sn_count <- sn_count + 1
        }
      }
    }
  }
}

#filter to only sequences with at least min_sns_in_seq short notes
sn_seqs <- sn_seqs[which(sn_seqs$n >= min_sns_in_seq),]

#If there is any movement in these sequences, remove periods above a threshold speed (double threhsold method used to identify contiguous high-speed chunks of time)
sentinel_bouts <- data.frame()
for(i in 1:nrow(sn_seqs)){
  ind <- sn_seqs$ind[i]
  t0 <- sn_seqs$t0[i]
  tf <- sn_seqs$tf[i]
  
  #if too near start, can't compute speed - skip the speed check
  if(t0 - speed_dt/2 <= 0){
    sentinel_bouts <- rbind(sentinel_bouts, data.frame(ind = ind, t0 = t0, tf = tf))
    next
  }
  
  x_past <- xs[ind,(t0-speed_dt/2):(tf-speed_dt/2)]
  y_past <- ys[ind,(t0-speed_dt/2):(tf-speed_dt/2)]
  x_fut <- xs[ind, (t0+speed_dt/2):(tf + speed_dt / 2)]
  y_fut <- ys[ind, (t0+speed_dt/2):(tf + speed_dt / 2)]
  speed <- sqrt((x_fut - x_past)^2 + (y_fut - y_past)^2) / speed_dt * 60 #speed in m / min
  
  #if there are any NAs skip the speed check
  if(sum(is.na(speed)>0)){
    sentinel_bouts <- rbind(sentinel_bouts, data.frame(ind = ind, t0 = t0, tf = tf))
    next
  }
  above_upper <- speed > speed_thresh_upper
  above_lower <- speed > speed_thresh_lower
  too_fast <- cocomo::get_together_sticky(above_upper, above_lower) #use double threshold to pull out times that are too fast
  
  #if there are no periods of fast movement found, just keep the bout whole
  if(sum(too_fast)==0){
    sentinel_bouts <- rbind(sentinel_bouts, data.frame(ind = ind, t0 = t0, tf = tf))
    next
  }
  subbouts_rle <- rle(too_fast)
  idxs_ends <- cumsum(subbouts_rle$lengths)
  idxs_starts <- c(1,idxs_ends+1)[1:(length(idxs_ends))] 
  subbouts <- data.frame(start = idxs_starts, end = idxs_ends, val = subbouts_rle$values)
  subbouts <- subbouts[which(subbouts$val==FALSE),] #subset to bouts below the speed threshold
  subbouts$start <- subbouts$start + t0 - 1
  subbouts$end <- subbouts$end + t0 - 1
  sentinel_bouts_new <- data.frame(ind = rep(ind,nrow(subbouts)), t0 = subbouts$start, tf = subbouts$end)
  sentinel_bouts <- rbind(sentinel_bouts, sentinel_bouts_new)
      
}

#count up the number of sns, rate of sns, and rate of close calls during the sequence time
sentinel_bouts$dur <- sentinel_bouts$tf - sentinel_bouts$t0
sentinel_bouts$sns <- NA
sentinel_bouts$ccs <- NA
for(i in 1:nrow(sentinel_bouts)){
  sentinel_bouts$sns[i] <- sum(calls_array[sentinel_bouts$ind[i], sentinel_bouts$t0[i]:sentinel_bouts$tf[i],'sn'],na.rm=T)
  sentinel_bouts$ccs[i] <- sum(calls_array[sentinel_bouts$ind[i], sentinel_bouts$t0[i]:sentinel_bouts$tf[i],'cc'],na.rm=T)
}
sentinel_bouts$sn_rate <- sentinel_bouts$sns / sentinel_bouts$dur
sentinel_bouts$cc_rate <- sentinel_bouts$ccs / sentinel_bouts$dur

#filter to minimum duration and max cc rate
sentinel_bouts <- sentinel_bouts[which(sentinel_bouts$dur > min_dur & sentinel_bouts$cc_rate < max_cc_rate),]

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
sentinel_bouts$tf_file <- sentinel_bouts$t0_file <- lubridate::seconds_to_period(0)
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

#compute rates of other call types
sentinel_bouts$al_rate <- sentinel_bouts$soc_rate <- sentinel_bouts$agg_rate <- NA
for(i in 1:nrow(sentinel_bouts)){
  ind <- sentinel_bouts$ind[i]
  t0 <- sentinel_bouts$t0[i]
  tf <- sentinel_bouts$tf[i]
  dur <- sentinel_bouts$dur[i]
  sentinel_bouts$al_rate[i] <- sum(calls_array[ind,t0:tf,'al'], na.rm=T) / dur
  sentinel_bouts$soc_rate[i] <- sum(calls_array[ind,t0:tf,'soc'], na.rm=T) / dur
  sentinel_bouts$agg_rate[i] <- sum(calls_array[ind,t0:tf,'agg'], na.rm=T) / dur
}

#i <- 11
# cocomo::plot_behav_and_calls(behavs, calls_array, behavs_key,
#                              focal_ind = sentinel_bouts$ind[i],
#                              t0 = sentinel_bouts$t0[i],
#                              tf = sentinel_bouts$tf[i],
#                              nonfocal_calls_to_plot <- c('cc','sn','al'),
#                              nonfocal_behavs_to_plot = c('Vigilance'),
#                              smooth_window = params$time_window)

#DISTANCE MOVED
sentinel_bouts$dist_spread <- NA
for(i in 1:nrow(sentinel_bouts)){
  x <- xs[cbind(sentinel_bouts$ind[i], sentinel_bouts$t0[i]:sentinel_bouts$tf[i])]
  y <- ys[cbind(sentinel_bouts$ind[i], sentinel_bouts$t0[i]:sentinel_bouts$tf[i])]
  centr_x <- mean(x, na.rm=T)
  centr_y <- mean(y, na.rm=T)
  sentinel_bouts$dist_spread[i] <- mean(sqrt((x - centr_x)^2 + (y - centr_y)^2),na.rm=T)
  
}

#BREAKING INTO BOUTS BASED ON DISTANCE MOVED

#plot(t0:tf,speed)
#abline(v=subbouts$start,col='green')
#abline(v=subbouts$end,col='red')

#row <- 50
#plot(xs[sentinel_bouts$ind[row],sentinel_bouts$t0[row]:sentinel_bouts$tf[row]],
#        ys[sentinel_bouts$ind[row],sentinel_bouts$t0[row]:sentinel_bouts$tf[row]],asp=1)

#save output
write.csv(sentinel_bouts, file = outfile, quote = F, row.names = F)

