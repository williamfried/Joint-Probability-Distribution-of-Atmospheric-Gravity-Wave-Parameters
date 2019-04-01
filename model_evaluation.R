source('station_info.R')
source('conditional.R')

# AIC corresponding to models of increasing complexity
# 'independent' models all 7 gravity wave parameters independently of each other
# 'copula' models the joint distribution of all gravity wave parameters except for frequency using a copula and independently models frequency  
# 'complete' models the joint distribution of all gravity wave parameters except for frequency using a copula and models the frequency by conditioning on the zonal and meridional wavelengths  
model_AIC = function(model)
{
  tot = 0
  
  if (model == 'independent')
  {
    for (i in 1:dim(data)[2])
    {
      u = U_amp_PDF(U_amp[i])
      v = V_amp_PDF(V_amp[i])
      t = T_amp_PDF(T_amp[i])
      u_wave = U_wavelength_PDF(U_wavelength[i])
      v_wave = V_wavelength_PDF(V_wavelength[i])
      f = freq_PDF(freq_shifted[i])
      vert = vert_wavelength_PDF(vert_wavelength[i])
      tot = tot + log(u*v*t*u_wave*v_wave*f*vert)
    }
    return(AIC(tot, marginal_parameters_num+length(freq_MLE)))
  }
  
  if (model == 'copula')
  {
    for (i in 1:dim(data)[2])
    {
      copula = dMvdc(c(U_amp[i],V_amp[i],T_amp[i],U_wavelength[i],V_wavelength[i],vert_wavelength[i]), cop_dist)
      f = freq_PDF(freq_shifted[i])
      tot = tot + log(copula*f)
    }
    return(AIC(tot, marginal_parameters_num+length(optimized_copula_list)+length(freq_MLE)))
  }
  
  if (model == 'complete')
  {
    both_conditional_list <<- cond_dist_both(U_wavelength, V_wavelength, freq_shifted, render_bins(optimal_bin_num, 'U'), render_bins(optimal_bin_num, 'V'), 'modeling')
    
    for (i in 1:dim(data)[2])
    {
      copula = dMvdc(c(U_amp[i],V_amp[i],T_amp[i],U_wavelength[i],V_wavelength[i],vert_wavelength[i]), cop_dist)
      f = conditional_freq(freq_shifted[i], paste(category(U_wavelength[i], render_bins(optimal_bin_num, 'U')), category(V_wavelength[i], render_bins(optimal_bin_num, 'V'))), both_conditional_list, 'density')   
      tot = tot + log(copula*f)
    }
    return(AIC(tot, marginal_parameters_num+length(optimized_copula_list)+conditional_list_param_num(both_conditional_list)))
  }
}

independent_AIC = model_AIC('independent')
print(paste('Independent model -- AIC:', toString(as.integer(independent_AIC))))
copula_AIC = model_AIC('copula')
print(paste('Copula model -- AIC:', toString(as.integer(copula_AIC))))
complete_AIC = model_AIC('complete')
print(paste('Complete model -- AIC:', toString(as.integer(complete_AIC))))
