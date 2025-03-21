#

#Paths

sentinel_bouts_csv <- '~/Desktop/Meerkat_data_for_hackathon/results/audio_only_sentinel_bouts_2025-03-20.csv'
rawdata_csv <- '~/Desktop/Meerkat_data_for_hackathon/processed_data/ZU_2021_allstreams.RData'
noise_csv <- '~/Desktop/Meerkat_data_for_hackathon/noise_detector/ZU_VZUF052_RRRT_R07_20210719-20210725_file_11_(2021_07_24-08_44_59)_ASWMUX221002.csv'
#noise_csv <- '~/Desktop/Meerkat_data_for_hackathon/noise_detector/ZU_VZUM053_RRLT_R31_20210714-20210725_file_11_(2021_07_24-08_44_59)_ASWMUX221102.csv'
wavfile <- 'ZU_VZUF052_RRRT_R07_20210719-20210725_file_11_(2021_07_24-08_44_59)_ASWMUX221002.wav'

#load data
sentinel_bouts <- read.csv(sentinel_bouts_csv)
load(rawdata_csv)
noise_bouts <- read.csv(noise_csv)

#get time in seconds into file
sentinel_bouts$t0_file_sec <- lubridate::period_to_seconds(as.period(sentinel_bouts$t0_file))
sentinel_bouts$tf_file_sec <- lubridate::period_to_seconds(as.period(sentinel_bouts$tf_file))

curr_bouts <- sentinel_bouts[which(sentinel_bouts$wav == wavfile),]
curr_bouts$noise_cov <- NA
for(i in 1:nrow(curr_bouts)){
  start <- curr_bouts$t0_file_sec[i]
  end <- curr_bouts$tf_file_sec[i]
  noise_bouts_curr <- noise_bouts[which(noise_bouts$TimeInSeconds >= start & noise_bouts$TimeInSeconds < end),]
  curr_bouts$noise_cov[i] <- sum(noise_bouts_curr$WidthInSeconds) / curr_bouts$dur[i]
}



plot(NULL, xlim = c(start, end), ylim = c(0,2))
arrows(x0 = noise_bouts_curr$TimeInSeconds, x1 = noise_bouts_curr$TimeInSeconds,
       y0 = rep(0, nrow(noise_bouts_curr)), y1 = rep(noise_bouts_curr$Amplitude),length=0)



noise_coverage <- sum(noise_bouts_curr$WidthInSeconds)
