#playing with different exclusion criteria

#EXCLUSION CRITERIA
max_dist <- 5
max_soc_agg_rate <- .05

infile <- '~/Desktop/Meerkat_data_for_hackathon/results/sentinel_detector_audio_only_2025-03-19_validated.csv'

sentinel_bouts <- read.csv(infile)
sentinel_bouts <- sentinel_bouts[which(sentinel_bouts$validation...presence!=''),]

sentinel_bouts$aggsoc_rate <- sentinel_bouts$agg_rate + sentinel_bouts$soc_rate


yes_idxs <- which(sentinel_bouts$validation...presence=='yes')
no_idxs <- which(sentinel_bouts$validation...presence=='no')
sun_idxs <- which(sentinel_bouts$validation...presence=='sunning')


