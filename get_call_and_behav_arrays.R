#Get call and behavior arrays of structure:
#calls_array[i,t,c] = number of calls of type c given by individual i at time step t (will be NA if no labeled audio during this window)
#calls_array_smoothed[i,t,c] = call rate of type c given by individual i at time step t (in a window of time_window timesteps around the current time)
#audio_labelled[i,t] = T or F whether audio was labelled (by machine) for individual i at time step t
#behavs_array[i,t,b] = T or F whether individual i was doing behavior b at time t
#behavs_array_smoothed[i,t,b] = fraction of time spent in behavior b for ind i at time t, over the time window time_window (centered on t)
#behavs_key = data frame indicating which behavior is associated with which index in the 3rd dimension of the behavior arrays

#Save this as well as the following items in one big file:
#xs, ys, ids, timestamps
#behavs

library(lubridate)
library(zoo)

#------PARAMETERS------
#Call types to include in arrays
call_types <- c('cc','sn','mo','ld','agg','al','soc')

#Time window for smoothing of call rates (and behavioral state rates)
time_window <- 31

#Input file paths
gps_file <- '/mnt/EAS_shared/meerkat/working/processed/movement/ZU2021_COORDINATES_all_sessions.RData'
behav_file <- '/mnt/EAS_shared/meerkat/working/processed/behavior/ZU2021_ACC_behavs.RData'
audio_file <- '/mnt/EAS_shared/meerkat/working/processed/acoustic/total_synched_call_tables/ML_synch/ZU2021_ALL_ML_CALLS_SYNCHED.csv'

#Output directory where everything will be stored together in 1 combined file
outdir <- '/mnt/EAS_shared/meerkat/working/processed/combined/'

#store parameters
params <- list()
params$gps_file <- gps_file
params$behav_file <- behav_file
params$audio_file <- audio_file
params$time_window <- time_window
params$call_types <- call_types

#------LOAD DATA------

load(gps_file)
load(behav_file)
calls <- read.csv(audio_file)

#------PREPROCESS------

#convert to standard naming conventions
xs <- ZU2021_allX
ys <- ZU2021_allY
timestamps <- ZU2021_timeLine
ids <- ZU2021_indInfo
ids$id_code <- ids$code
calls$date <- lubridate::date(calls$t0GPS_UTC)
timestamps <- as.POSIXct(timestamps, tz = 'UTC')

#get basic info
n_inds <- nrow(xs)
n_times <- ncol(xs)
n_calltypes <- length(call_types)
n_behavs <- nrow(behavs_key)

#get indexes to recording period starts
#time diffs between timestamps (in sec)
diffs <- as.numeric(difftime(timestamps[2:n_times], timestamps[1:(n_times)-1], units = 'sec'))
breaks <- which(diffs > 1) + 1
breaks <- c(1, breaks)
n_breaks <- length(breaks)

#add time index for calls
calls$ind_idx <- match(calls$ind, ids$id_code)
calls$time_idx <- match(lubridate::round_date(as.POSIXct(calls$t0GPS_UTC, tz='UTC'), unit = 'second'), timestamps)

#get onset and offset of labels for each individual on each date
uniq_dates <- unique(calls$date)
label_efforts <- data.frame(id_code = rep(ids$id_code, each = length(uniq_dates)),
                            date = rep(uniq_dates, length(ids$id_code)))
label_efforts$end_UTC <- label_efforts$start_UTC <- NA
label_efforts$end_idx <- label_efforts$start_idx <- NA
for(i in 1:nrow(label_efforts)){
  labs_ind <- calls[which(calls$ind == label_efforts$id_code[i] & calls$date == label_efforts$date[i]),]
  if(nrow(labs_ind)>0){
    label_efforts$start_UTC[i] <- min(labs_ind$t0GPS_UTC, na.rm=T)
    label_efforts$end_UTC[i] <- max(labs_ind$tendGPS_UTC, na.rm=T)
    label_efforts$start_idx[i] <- min(labs_ind$time_idx, na.rm=T)
    label_efforts$end_idx[i] <- max(labs_ind$time_idx, na.rm=T)
  }
}

#make a matrix that shows whether labeled audio data is available for that individual
audio_labelled <- matrix(F, nrow = n_inds, ncol = n_times)
for(i in 1:nrow(label_efforts)){
  date <- label_efforts$date[i]
  id_code <- label_efforts$id_code[i]
  ind_idx <- which(ids$id_code == id_code)
  t0_idx <- label_efforts$start_idx[i]
  tf_idx <- label_efforts$end_idx[i]
  if(!is.na(t0_idx) & !is.na(tf_idx) & !is.infinite(t0_idx) & !is.infinite(tf_idx)){
    audio_labelled[ind_idx, t0_idx:tf_idx] <- T
  }
}

#get number of calls in each second for each call type, in array format
calls_array <- array(data = 0, dim = c(n_inds, n_times, n_calltypes))
dimnames(calls_array) <- list(ids$id_code, timestamps, call_types)
for(i in 1:nrow(calls)){
  ind_idx <- calls$ind_idx[i]
  call_type <- calls$callType[i]
  time_idx <- calls$time_idx[i]
  foc <- calls$focalType[i]
  if(!is.na(ind_idx) & !is.na(call_type) & call_type %in% call_types & foc == 'F'){
    calls_array[ind_idx, time_idx, call_type] <- calls_array[ind_idx, time_idx, call_type] + 1
  }
}

#replace times when audio was not labelled with NAs
for(i in 1:n_calltypes){
  calls_curr <- calls_array[,,i]
  calls_curr[!audio_labelled] <- NA
  calls_array[,,i] <- calls_curr
}

#get smoothed call rates over a time window
calls_array_smoothed <- array(NA, dim = dim(calls_array), dimnames = dimnames(calls_array))
breaks_with_endpoint <- c(breaks, n_times + 1)
print('smoothing call rates')
for(i in 1:n_breaks){
  print(i)
  t0 <- breaks_with_endpoint[i]
  tf <- breaks_with_endpoint[i+1]-1
  for(ind in 1:n_inds){
    for(c in 1:n_calltypes){
      calls_array_smoothed[ind, t0:tf, c] <- zoo::rollmean(calls_array[ind, t0:tf, c], k = time_window, fill = NA, align = 'center')
    }
  }
}

#convert behavs matrix into behavs array (one hot encoding)
behavs_array <- array(NA, dim = c(n_inds, n_times, nrow(behavs_key)))
for(b in 1:n_behavs){
  for(i in 1:n_inds){
    behavs_array[i,,b] <- behavs[i,] == b
  }
}
dimnames(behavs_array) <- list(ids$id_code, timestamps, behavs_key$behav)

#get smoothed behavioral state fractions over a time window
behavs_array_smoothed <- array(NA, dim = dim(behavs_array), dimnames = dimnames(behavs_array))
breaks_with_endpoint <- c(breaks, n_times + 1)
print('smoothing behavior rates')
for(i in 1:n_breaks){
  print(i)
  t0 <- breaks_with_endpoint[i]
  tf <- breaks_with_endpoint[i+1]-1
  for(ind in 1:n_inds){
    for(b in 1:n_behavs){
      behavs_array_smoothed[ind, t0:tf, b] <- zoo::rollmean(behavs_array[ind, t0:tf, b], k = time_window, fill = NA, align = 'center')
    }
  }
}

setwd(outdir)
save(list = c('xs','ys','ids','timestamps','behavs','behavs_array','behavs_array_smoothed','behavs_key','calls','calls_array','calls_array_smoothed', 'audio_labelled','params','breaks'), file = 'ZU_2021_allstreams.RData')

#PLAYING WITH UMAP
#flatten array
# flat_mat <- matrix(NA, nrow = n_inds*n_times, ncol = n_calltypes)
# for(c in 1:n_calltypes){
#   for(i in 1:n_inds){
#     idx_start <- ((i-1)*n_times+1)
#     idx_end <- idx_start + n_times - 1
#     flat_mat[idx_start:idx_end,c] <- calls_array_smoothed[i,,c]
#   }
# }
# 
# non_nas <- which(!is.na(flat_mat[,1]))
# samps <- sample(non_nas, 10000)
# test <- umap(flat_mat[samps,])
# plot(test$layout[,1], test$layout[,2])
# 
# histo <- hist(calls_array_smoothed[,,'sn']*60, breaks=50)
# plot(histo$mids, log(histo$density))
