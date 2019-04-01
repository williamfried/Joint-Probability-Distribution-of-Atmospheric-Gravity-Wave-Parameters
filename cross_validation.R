source('station_info.R')
source('functions.R')
source('conditional.R')
library('parallel')
setwd(path_to_data)

# determine AIC given testing set, training set and the number of intervals
cross_validation_AIC = function(test_data, train_data, bins)
{
  # extract training data
  U_wavelength_train = as.numeric(as.vector(train_data['U_wavelength',]))
  V_wavelength_train = as.numeric(as.vector(train_data['V_wavelength',]))
  freq_train = as.numeric(as.vector(train_data['freq',]))
  freq_shifted_train = freq_train - (min(freq) - 0.025)

  # train conditional probability distributions
  both_conditional_list_train = tryCatch(cond_dist_both(U_wavelength_train, V_wavelength_train, freq_shifted_train, render_bins(bins, 'U'), render_bins(bins, 'V'), 'cross-validation'), error=function(cond){return(NA)})
  if (is.na(both_conditional_list_train))
  {
    return(NA)
  }
  #if (is(try(cond_dist_both(U_wavelength_train, V_wavelength_train, freq_shifted_train, render_bins(bins, 'U'), render_bins(bins, 'V'), 'cross-validation')), 'try-error'))
  #{
  #  return(NA)
  #}
  #else
  #{
  #  both_conditional_list_train = cond_dist_both(U_wavelength_train, V_wavelength_train, freq_shifted_train, render_bins(bins, 'U'), render_bins(bins, 'V'), 'cross-validation')
  #}
  
  # extract testing data
  vert_wavelength_test = as.numeric(as.vector(test_data['vertical_wavelength',]))
  U_amp_test = as.numeric(as.vector(test_data['U_amp',]))
  V_amp_test = as.numeric(as.vector(test_data['V_amp',]))
  T_amp_test = as.numeric(as.vector(test_data['T_amp',]))
  U_wavelength_test = as.numeric(as.vector(test_data['U_wavelength',]))
  V_wavelength_test = as.numeric(as.vector(test_data['V_wavelength',]))
  freq_test = as.numeric(as.vector(test_data['freq',]))
  freq_shifted_test = freq_test - (min(freq) - 0.025)
  
  # calculate AIC of joint distribution of testing data
  tot = 0
  for (i in 1:dim(test_data)[2])
  {
    copula = dMvdc(c(U_amp_test[i],V_amp_test[i],T_amp_test[i],U_wavelength_test[i],V_wavelength_test[i],vert_wavelength_test[i]), cop_dist)
    if(copula == 0)
    {
      stop('zero copula')
    }
    f = conditional_freq(freq_shifted_test[i], paste(category(U_wavelength_test[i], render_bins(bins, 'U')), category(V_wavelength_test[i], render_bins(bins, 'V'))), both_conditional_list_train, 'density') 
    
    tot = tot + log(copula*f)
  }
  conditional_AIC = AIC(tot, 0)
  return(conditional_AIC)
}

# apply 10-fold cross-validation
cross_validate = function(bins)
{
  # randomly shuffle data
  data_shuffled = data[,sample(ncol(data))]
  
  # create 10 equally size folds 
  folds = cut(seq(1,ncol(data)), breaks=10, labels=F)
  
  AIC_vals = c()
  
  for (i in 1:10)
  {
    test_indices = which(folds==i, arr.ind=T)
    test_data = data_shuffled[,test_indices]
    train_data = data_shuffled[,-test_indices]
    AIC_vals[i] = cross_validation_AIC(test_data, train_data, bins)
    if (is.na(AIC_vals[i]))
    {
      return(NA)
    }
  }
  print(bins)
  return(sum(AIC_vals) + 2*(marginal_parameters_num + length(copula_list) + 2*(length(bins)^2)))
}

# List of number of intervals to try  
bin_trials = 5:20

# parallelize process of testing different number of intervals
# note: this function takes ~15 minutes to run
clusters = makeCluster((detectCores()-1), type='FORK', outfile='cross_validation_progress.txt')
results_vector = unlist(parLapply(clusters, bin_trials, cross_validate))
stopCluster(clusters)
file.remove('cross_validation_progress.txt')

bin_trials = bin_trials[!is.na(results_vector) & !is.infinite(results_vector)]
results_vector = results_vector[!is.na(results_vector) & !is.infinite(results_vector)]
results_vector = results_vector / 1000

# best fit curve of AIC vs. number of intervals
bin_fit = lm(results_vector~poly(bin_trials, floor(length(bin_trials)/1.5), raw=T))
AIC_predicted = predict(bin_fit)

# plot raw data points and overlay best fit curve
par(mfrow=c(1,1))
plot(bin_trials, results_vector, xlab='number of intervals', ylab='AIC (thousands)', main='Optimal number of intervals using cross-validation', cex.lab=1.5, cex.main=1.5)
par(new=T)
lines(bin_trials, AIC_predicted, col='red')

# determine number of intervals that minimizes AIC
optimal_bin_num = bin_trials[which.min(AIC_predicted)]
print(paste('Optimal number of intervals:', optimal_bin_num))
