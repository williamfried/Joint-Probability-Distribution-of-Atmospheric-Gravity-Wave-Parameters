setwd(path_to_data)
source('functions.R')

all_dists = c('gamma', 'weibull', 'lnorm', 'rayleigh', 'gumbel', 'burr', 'dagum', 'invgamma', 'genray', 'generalized_normal')

# determine distribution with lowest AIC for all 7 gravity wave parameters
U_amp_fit = fit_dist(U_amp, T, F, F, all_dists, 'U_amp')
V_amp_fit = fit_dist(V_amp, T, F, F, all_dists, 'V_amp')
stop(1)
T_amp_fit = fit_dist(T_amp, F, F, F, all_dists, 'T_amp')
freq_fit = fit_dist(freq_shifted, F, F, F, all_dists, 'freq')
U_wavelength_fit = fit_dist(U_wavelength, F, F, F, all_dists, 'U_wavelength')
V_wavelength_fit = fit_dist(V_wavelength, F, F, F, all_dists, 'V_wavelength')
vert_wavelength_fit = fit_dist(vert_wavelength, F, F, F, all_dists, 'vertical_wavelength')

# caculuate maximum likelihood estimate corresponding to distribution with lowest AIC for all 7 gravity wave parameters
calc_MLE = function(dist, data)
{
  MLE_list = list(normal_MLE, generalized_normal_MLE, lognormal_MLE, gamma_MLE, weibull_MLE, pareto_MLE, rayleigh_MLE, gumbel_MLE, burr_MLE, dagum_MLE, invgamma_MLE, genray_MLE)
  names(MLE_list) = c('norm', 'generalized_normal', 'lnorm', 'gamma', 'weibull', 'pareto', 'rayleigh', 'gumbel', 'burr', 'dagum', 'invgamma', 'genray')
  return(MLE_list[[dist]](data))
}
  
U_amp_MLE = calc_MLE(names(U_amp_fit)[1], U_amp)
V_amp_MLE = calc_MLE(names(V_amp_fit)[1], V_amp)
T_amp_MLE = calc_MLE(names(T_amp_fit)[1], T_amp)
freq_MLE = calc_MLE(names(freq_fit)[1], freq_shifted)
vert_wavelength_MLE = calc_MLE(names(vert_wavelength_fit)[1], vert_wavelength)
U_wavelength_MLE = calc_MLE(names(U_wavelength_fit)[1], U_wavelength)
V_wavelength_MLE = calc_MLE(names(V_wavelength_fit)[1], V_wavelength)

# compute probability density given maximum likelihood estimate
MLE_density = function(x, dist, MLE)
{
  if (dist == 'gamma')
  {
    return(dgamma(x, MLE['shape'][[1]], rate=MLE['rate'][[1]]))
  }
  if (dist == 'lnorm')
  {
    return(dlnorm(x, meanlog=MLE['meanlog'][[1]], sdlog=MLE['sdlog'][[1]]))
  }
  if (dist =='generalized_normal')
  {
    return(dgnorm(x, MLE['mu'][[1]], MLE['alpha'][[1]], MLE['beta'][[1]]))
  }
  if (dist == 'norm')
  {
    return(dnorm(x, MLE['mean'][[1]], MLE['sd'][[1]]))
  }
  if (dist == 'weibull')
  {
    return(dweibull(x, MLE['shape'][[1]], scale=MLE['scale'][[1]]))
  }
  if (dist == 'pareto')
  {
    return(dpareto(x, MLE['scale'][[1]], shape=MLE['shape'][[1]]))
  }
  if (dist == 'rayleigh')
  {
    return(drayleigh(x, MLE['scale'][[1]]))
  }
  if (dist == 'gumbel')
  {
    return(dgumbel(x, MLE[['a']], MLE[['b']]))
  }
  if (dist == 'burr')
  {
    return(dburr(x, MLE[['shape1']], MLE[['shape2']], MLE[['rate']]))
  }
  if (dist == 'dagum')
  {
    return(ddagum(x, MLE[['scale']], MLE[['shape1.a']], MLE[['shape2.p']]))
  }
  if (dist == 'invgamma')
  {
    return(dinvgamma(x, MLE[['shape']], MLE[['rate']]))
  }
  if (dist == 'genray')
  {
    return(dgenray(x, MLE[['scale']], MLE[['shape']]))
  }
  stop('invalid distribution')
}

U_amp_PDF = function(x)
{
  return(MLE_density(x, names(U_amp_fit)[1], U_amp_MLE))
}

V_amp_PDF = function(x)
{
  return(MLE_density(x, names(V_amp_fit)[1], V_amp_MLE))
}

T_amp_PDF = function(x)
{
  return(MLE_density(x, names(T_amp_fit)[1], T_amp_MLE))
}

U_wavelength_PDF = function(x)
{
  return(MLE_density(x, names(U_wavelength_fit)[1], U_wavelength_MLE))
}

V_wavelength_PDF = function(x)
{
  return(MLE_density(x, names(V_wavelength_fit)[1], V_wavelength_MLE))
}

freq_PDF = function(x)
{
  return(MLE_density(x, names(freq_fit)[1], freq_MLE))
}

vert_wavelength_PDF = function(x)
{
  return(MLE_density(x, names(vert_wavelength_fit)[1], vert_wavelength_MLE))
}

# total number of parameters associated with maximum likelihood estimates of all parameters except for frequency
marginal_parameters_num = length(c(U_amp_MLE, V_amp_MLE, T_amp_MLE, U_wavelength_MLE, V_wavelength_MLE, vert_wavelength_MLE))
