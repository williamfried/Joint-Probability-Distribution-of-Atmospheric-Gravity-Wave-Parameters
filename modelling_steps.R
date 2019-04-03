source('station_info.R')

# run this script to build model

# import data
source('import.R')
data = import(dataset)

# fit marginal distributions
source('fitting.R')

# construct and optimize copula
source('copula.R')

# determine optimal number of intervals of zonal/meridional wavelength
source('cross_validation.R')

# evaluate 3 models: independent, copula, complete
source('model_evaluation.R')

# sample from joint probability distribution
source('final_joint_distribution.R')



