#Second pass at a sentinel detector - audio only
#
#The detector uses ONLY audio labels to (attempt to) detect bouts of sentinel behavior
#
#As a start, we simply remove the ACC-based requirements of the first pass sentinel detector
#
#For each time point, we look in a time window t_window (currently 31 seconds) centered on that time point
#and count up how many (focal) sn's and cc's the individual produced, then normalize this to calculate the call rate for each call type.
#
#We then define an individual as being in the "sentinel" state at time t if:
#1. The individual's call rate for sn's in the time window exceeds min_sn_thresh (here set to 0 - so the individual has to give at least 1 short note)
#2. The individual's call rate for cc's in the time window falls below max_cc_thresh (here set to 0.001 meaning no cc's are allowed)
#
#Once we have cateogrized each time point at being in the "sentinel" state or not, we 
#find contiguous bouts of sentinel states for each individual. 
#The resulting bouts are stored in the data frame sentinel_bouts, with onset and offset times (t0 and tf)
#correspoding to the beginning of the first time window and the end of the last time window 
#
#Finally, we connect adjacent bouts to one another if they are less than max_time_connect_bouts time apart

library(lubridate)

#-----PARAMS-----
#input data file
infile <- '/mnt/EAS_shared/meerkat/working/processed/combined/ZU_2021_allstreams.RData'

#where to save the output table
outfile <- '/mnt/EAS_ind/astrandburg/sentinel_detector/audio_only_sentinel_bouts_2025-03-19.csv'

#detecting whether an individual is a sentinel at any given time point
max_sn_gap_len <- 60
min_sns_in_seq <- 10
min_dur <- 120 #minimum (time) duration
max_cc_rate <- .01 #maximum rate of ccs 

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

#count up the rate of close calls during the sequence time
sn_seqs$dur <- sn_seqs$tf - sn_seqs$t0
for(i in 1:nrow(sn_seqs)){
  sn_seqs$ccs[i] <- sum(calls_array[sn_seqs$ind[i], sn_seqs$t0[i]:sn_seqs$tf[i],'cc'],na.rm=T)
}
sn_seqs$cc_rate <- sn_seqs$ccs / sn_seqs$dur

#filter to minimum duration and max cc rate
sentinel_bouts <- sn_seqs[which(sn_seqs$dur > min_dur & sn_seqs$cc_rate < max_cc_rate),]

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

#save output
write.csv(sentinel_bouts, file = outfile, quote = F, row.names = F)

i <- 11
cocomo::plot_behav_and_calls(behavs, calls_array, behavs_key,
                             focal_ind = sentinel_bouts$ind[i],
                             t0 = sentinel_bouts$t0[i],
                             tf = sentinel_bouts$tf[i],
                             nonfocal_calls_to_plot <- c('cc','sn','al'),
                             nonfocal_behavs_to_plot = c('Vigilance'),
                             smooth_window = params$time_window)
