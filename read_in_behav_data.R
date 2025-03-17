#Read in and convert to matrix form the behavioral sequence data generated from Pranav's classifier
#Do this for now only for ZU2021 (both deployments)
#save to the output file

library(cocomo)

#PATHS

#path to GPS data
gps_data_path <- '/mnt/EAS_shared/meerkat/working/processed/movement/ZU2021_COORDINATES_all_sessions.RData'
load(gps_data_path)

#get files in directory
input_dirs <- c('/mnt/EAS_ind/pminasandra/data/PredictionsByIndividualDeploymentwise/ZU_2021_1/',
                '/mnt/EAS_ind/pminasandra/data/PredictionsByIndividualDeploymentwise/ZU_2021_2/')

#output file name
outfile <- '/mnt/EAS_shared/meerkat/working/processed/behavior/ZU2021_ACC_behavs.RData'

#PREPROCESS

reality_check_plot <- F

#put GPS data and ids table into standard format
ids <- ZU2021_indInfo
ids$id_code <- ids$code
timestamps2 <- ZU2021_timeLine
timestamps <- as.POSIXct(timestamps2, tz = 'UTC')

out <- import_behavioral_seqs_to_matrix(input_dirs = input_dirs, ids, timestamps)

behavs <- out$behavs
behavs_key <- out$behavs_key

save(list = c('behavs','behavs_key'), file = outfile)


#reality check - compare with GPS speed for each behavior
if(reality_check_plot){
  xs <- ZU2021_allX
  ys <- ZU2021_allY
  n_inds <- nrow(xs)
  n_times <- ncol(xs)
  
  speeds <- matrix(NA, nrow = n_inds, ncol = n_times)
  for(i in 1:n_inds){
    out <- cocomo::get_heading_and_speed_temporal(xs[i,], ys[i,], t_window = 5, seconds_per_time_step = 1)
    speeds[i,] <- out$speeds
  }
  
  bins <- seq(-4,2,.4)
  mids <- bins[2:length(bins)] - diff(bins)/2
  
  speeds_hist <- matrix(nrow = 3, ncol = length(mids))
  for(i in 1:3){
    speeds_hist[i,] <- hist(log(speeds[which(behavs==i)]),breaks=bins, plot=F)$counts
    speeds_hist[i,] <- speeds_hist[i,] / sum(speeds_hist[i,])
  }
  
  
  plot(mids, speeds_hist[1,], xlim = c(-2.5,2), type = 'l', col = 'darkblue', xlab = 'Log speed (m/s)', ylab = 'Probability')
  cols <- c('darkblue','darkorange','darkred')
  for(i in 1:3){
    lines(mids, speeds_hist[i,], col = cols[i],lwd=2)
  }
  
  legend('topright',legend=c('Vigilance','Foraging','Running'),col=cols, lwd=2)
}
