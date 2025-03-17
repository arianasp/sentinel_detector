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

#put GPS data and ids table into standard format
ids <- ZU2021_indInfo
ids$id_code <- ids$code
timestamps2 <- ZU2021_timeLine
timestamps <- as.POSIXct(timestamps2, tz = 'UTC')

out <- import_behavioral_seqs_to_matrix(input_dirs = input_dirs, ids, timestamps)

behavs <- out$behavs
behavs_key <- out$behavs_key

save(list = c('behavs','behavs_key'), file = outfile)
