setwd(path_to_data)
source('functions.R')
library('mixR')
library('stringr')
library('dplyr')

MLE_list = list(normal_MLE, lognormal_MLE, gamma_MLE, weibull_MLE, rayleigh_MLE, gumbel_MLE, burr_MLE, dagum_MLE, invgamma_MLE, genray_MLE, generalized_normal_MLE, polynomial_MLE)        
names(MLE_list) = c('norm', 'lnorm', 'gamma', 'weibull', 'rayleigh', 'gumbel', 'burr', 'dagum', 'invgamma', 'genray', 'generalized_normal', 'polynomial')

cond_dist_both = function(U_parameter, V_parameter, frequency_train, U_intervals, V_intervals, task)
{
  if (task == 'cross-validation')
  {
    distributions = c('gamma', 'weibull', 'rayleigh', 'lnorm')
    #distributions = c('polynomial')
  }
  else if (task == 'modeling')
  {
    #distributions = c('gamma', 'weibull', 'rayleigh', 'lnorm', 'gumbel', 'burr', 'dagum', 'invgamma', 'mixture model')#, 'polynomial')
    distributions = c('polynomial')
  }
  else
  {
    stop('invalid')
  }
  
  conditional_list = list()
  
  for (i in 1:(length(U_intervals)-1))
  {
    #print(i)
    for (j in 1:(length(V_intervals)-1))
    {
      #print(j)
      segment = frequency_train[between(U_parameter, U_intervals[i], U_intervals[i+1]) & between(V_parameter, V_intervals[j], V_intervals[j+1])]
      if (length(distributions) == 1)
      {
        dist = distributions
      }
      else
      {
        fit = fit_dist(segment, F, F, F, distributions, 'freq')
        dist = names(fit)
      }
      
      if (grepl('mixture model', dist))
      {
        conditional_list[[paste(toString(U_intervals[i]), toString(V_intervals[j]))]] = list('dist' = 'mixture model', 'parameters' = mixture_model_MLE(segment, word(dist, 1, 1)))
      }
      else
      {
        conditional_list[[paste(toString(U_intervals[i]), toString(V_intervals[j]))]] = list('dist' = dist, 'parameters' = MLE_list[[dist]](segment))
      }
    }
  }
  return(conditional_list)
}

# split the data into a specified number of intervals that all contain the same number of data points 
render_bins = function(bin_num, type)
{
  if (type == 'U')
  {
    U_wavelength_intervals = unname(quantile(U_wavelength, seq(0,1,length=bin_num+1)))
    U_wavelength_intervals[bin_num+1] = U_wavelength_intervals[bin_num] + 10
    return(U_wavelength_intervals)
  }
  else if (type == 'V')
  {
    V_wavelength_intervals = unname(quantile(V_wavelength, seq(0,1,length=bin_num+1)))
    V_wavelength_intervals[bin_num+1] = V_wavelength_intervals[bin_num] + 10
    return(V_wavelength_intervals)
  }
  else
  {
    stop('invalid')
  }
}

# determine which interval a given zonal/meridional wavelength value falls into 
category = function(x, intervals)
{
  for (i in 1:(length(intervals)-1))
  {
    if (x >= intervals[i] & x < intervals[i+1])
    {
      return(toString(intervals[i]))
    }
  }
}

# total number of parameters associated with maximum likelihood estimates of the conditional frequency distribution for each interval 
conditional_list_param_num = function(cond_list)
{
  tot = 0
  for (i in 1:length(cond_list))
  {
    if (cond_list[[i]]$dist == 'polynomial')
    {
      add = 2
    }
    else if (cond_list[[i]]$dist == 'mixture model')
    {
      add = 3*length(cond_list[[i]]$parameters$fit[[1]])-1
    }
    else
    {
      add = length(cond_list[[i]]$parameters)
    }
    tot = tot + add
  }
  return(tot)
}
