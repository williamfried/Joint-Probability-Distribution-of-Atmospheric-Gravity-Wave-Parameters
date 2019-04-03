setwd(path_to_data)
source('explore_functions.R')
library('copula')

# Spearman's rho correlation coefficient between all parameters except for frequency 
U_amp_V_amp = cor(U_amp, V_amp, method='spearman') 
U_amp_T_amp = cor(U_amp, T_amp, method='spearman') 
U_amp_U_wave = cor(U_amp, U_wavelength, method='spearman') 
U_amp_V_wave = cor(U_amp, V_wavelength, method='spearman') 
U_amp_vert = cor(U_amp, vert_wavelength, method='spearman') 
V_amp_T_amp = cor(V_amp, T_amp, method='spearman') 
V_amp_U_wave = cor(V_amp, U_wavelength, method='spearman') 
V_amp_V_wave = cor(V_amp, V_wavelength, method='spearman') 
V_amp_vert = cor(V_amp, vert_wavelength, method='spearman')  
T_amp_U_wave = cor(T_amp, U_wavelength, method='spearman')
T_amp_V_wave = cor(T_amp, V_wavelength, method='spearman') 
T_amp_vert = cor(T_amp, vert_wavelength, method='spearman') 
U_wave_V_wave = cor(U_wavelength, V_wavelength, method='spearman')
U_wave_vert = cor(U_wavelength, vert_wavelength, method='spearman')
V_wave_vert = cor(V_wavelength, vert_wavelength, method='spearman')

# offsets from Spearman's rho correlation coefficients (values are initialized to zero)
copula_list = c(
  'U_amp_V_amp' = 0,
  'U_amp_T_amp' = 0,
  'U_amp_U_wave' = 0,
  'U_amp_V_wave' = 0,
  'V_amp_T_amp' = 0,
  'V_amp_U_wave' = 0,
  'V_amp_V_wave' = 0,
  'T_amp_U_wave' = 0,
  'T_amp_V_wave' = 0,
  'U_wave_V_wave' = 0,
  'U_amp_vert' = 0,
  'V_amp_vert' = 0,
  'T_amp_vert' = 0,
  'U_wave_vert' = 0,
  'V_wave_vert' = 0)

# generate copula given correlation matrix 
generate_copula = function(copula_list)
{
  correlation_matrix = matrix(c(1, U_amp_V_amp+copula_list[['U_amp_V_amp']], U_amp_T_amp+copula_list[['U_amp_T_amp']], U_amp_U_wave+copula_list[['U_amp_U_wave']], U_amp_V_wave+copula_list[['U_amp_V_wave']], U_amp_vert+copula_list[['U_amp_vert']], U_amp_V_amp+copula_list[['U_amp_V_amp']], 1, V_amp_T_amp+copula_list[['V_amp_T_amp']], 
                                V_amp_U_wave+copula_list[['V_amp_U_wave']], V_amp_V_wave+copula_list[['V_amp_V_wave']], V_amp_vert+U_amp_vert+copula_list[['V_amp_vert']], U_amp_T_amp+copula_list[['U_amp_T_amp']], V_amp_T_amp+copula_list[['V_amp_T_amp']], 1, T_amp_U_wave+copula_list[['T_amp_U_wave']], T_amp_V_wave+copula_list[['T_amp_V_wave']], T_amp_vert+copula_list[['T_amp_vert']], 
                                U_amp_U_wave+copula_list[['U_amp_U_wave']], V_amp_U_wave+copula_list[['V_amp_U_wave']], T_amp_U_wave+copula_list[['T_amp_U_wave']], 1, U_wave_V_wave+copula_list[['U_wave_V_wave']], U_wave_vert+copula_list[['U_wave_vert']], U_amp_V_wave+copula_list[['U_amp_V_wave']], V_amp_V_wave+copula_list[['V_amp_V_wave']], T_amp_V_wave+copula_list[['T_amp_V_wave']], 
                                U_wave_V_wave+copula_list[['U_wave_V_wave']], 1, V_wave_vert+copula_list[['V_wave_vert']], U_amp_vert+copula_list[['U_amp_vert']], V_amp_vert+copula_list[['V_amp_vert']], T_amp_vert+copula_list[['T_amp_vert']], U_wave_vert+copula_list[['U_wave_vert']], V_wave_vert+copula_list[['V_wave_vert']], 1), ncol=6) 
  
  normal_copula = normalCopula(param=P2p(correlation_matrix), dim=ncol(correlation_matrix), dispstr='un')
  cop_dist = mvdc(copula=normal_copula, margins=c(names(U_amp_fit), names(V_amp_fit), names(T_amp_fit), names(U_wavelength_fit), names(V_wavelength_fit), 'gnorm'),
                  paramMargins=list(U_amp_MLE, V_amp_MLE, T_amp_MLE, U_wavelength_MLE, V_wavelength_MLE, vert_wavelength_MLE))
}

# density of joint distribution of all parameters except for frequency using copula method and given correlation matrix between the parameters
copula_density = function(copula_list)
{    
  cop_dist = generate_copula(copula_list)
  
  tot = 0
  for (i in 1:dim(data)[2])
  {
    dens = dMvdc(c(U_amp[i],V_amp[i],T_amp[i],U_wavelength[i],V_wavelength[i],vert_wavelength[i]), cop_dist)
    tot = tot + log(dens)
  }
  return(AIC(tot, marginal_parameters_num+length(copula_list)))
}

# optimize values of correlation matrix to minimize AIC of joint distribution of all parameters except for frequency (note: function takes 1-2 hours to run)
corr_optimize_copula = function(increment)
{
  while (T)
  {
    cnt = 0
    
    for (i in 1:length(copula_list))
    {
      print(names(copula_list[i]))
      baseline_AIC = copula_density(copula_list)
      initial = copula_list[[i]]
      copula_list[[i]] = initial + increment
      inc_AIC = copula_density(copula_list)
      copula_list[[i]] = initial - increment
      dec_AIC = copula_density(copula_list)
      
      if (inc_AIC < baseline_AIC & dec_AIC > baseline_AIC)
      {
        sign = 1
      }
      else if (inc_AIC > baseline_AIC & dec_AIC < baseline_AIC)
      {
        sign = -1
      }
      else if(inc_AIC > baseline_AIC & dec_AIC > baseline_AIC)
      {
        copula_list[[i]] = initial
        sign = 'optimized'
      }
      else
      {
        return('error')
      }
      
      if (sign != 'optimized')
      {
        cnt = cnt + 1
        copula_list[[i]] = initial
        prev_AIC = baseline_AIC
        
        while (T)
        {
          copula_list[[i]] = copula_list[[i]] + sign*increment
          AIC_val = copula_density(copula_list)
          if (AIC_val > prev_AIC)
          {
            copula_list[[i]] = copula_list[[i]] - sign*increment
            break
          }
          print(copula_list[[i]])
          print(AIC_val)
          prev_AIC = AIC_val
        }
      }
    }
    if (cnt == 0)
    {
      break
    }
    print(paste('Proportion of correlation coefficients changed: ', toString(cnt), '/', toString(length(copula_list)), sep=''))
    print(copula_list)
  }
  print('copula optimization complete')
  return(copula_list)
}

# optimize correlation matrix between all parameters except for frequency
optimized_copula_list = corr_optimize_copula(0.01)

# create copula object using optimized correlation matrix
cop_dist = generate_copula(optimized_copula_list)
