source('station_info.R')
source('functions.R')
library('gplots')

quantities = c('vertical_wavelength','U_wavelength','V_wavelength', 'freq', 'U_amp','V_amp','T_amp')
data_list = list(vert_wavelength, U_wavelength, V_wavelength, freq_shifted, U_amp, V_amp, T_amp)
names(data_list) = quantities
axis_list = axis_limits(data_list)

# scale and shift quantities back to original values
rescale = function(data, quantity)
{
  if (quantity == 'freq')
  {
    data = (data*10^4 - (min(freq) - 0.025)) 
  }
  if (quantity %in% c('U_wavelength','V_wavelength'))
  {
    data = data / 1000
  }
  return(data)
}

# draw random sample of size n from joint probability distribution
sample_joint_distribution = function(n)
{
  # random sample from copula of size n
  six_parameters = rMvdc(n, cop_dist)
  
  U_wavelength_rand <<- six_parameters[,4]
  V_wavelength_rand <<- six_parameters[,5]
  
  # random sample of shifted frequency conditioned on zonal and meridional wavelength values
  freq_shifted_rand = rep(NA, n)
  
  for (i in 1:n)
  {
    print(i)
    f = conditional_freq(1, paste(category(U_wavelength_rand[i], render_bins(optimal_bin_num, 'U')), category(V_wavelength_rand[i], render_bins(optimal_bin_num, 'V'))), both_conditional_list, 'random') 
    freq_shifted_rand[i] = f
  }

  # unshift frequency values
  freq_rand = (freq_shifted_rand + (min(freq) - 0.025)) / 10^4
  
  # unscale zonal and meridional wavelength values  
  six_parameters = cbind(six_parameters[,1:3], apply(six_parameters[,4:5], 2, function(x) x*1000), six_parameters[,6])
  
  seven_parameters = cbind(six_parameters, freq_rand)
  colnames(seven_parameters) = c('U_amp','V_amp', 'T_amp', 'U_wavelength', 'V_wavelength', 'vertical_wavelength', 'freq')
  return(seven_parameters)
}

# compare model KDE to raw data KDE
marginals = function(sample_dataset)
{
  par(mfrow=c(1,1))
  for (i in 1:length(quantities))
  {
    param = sample_dataset[,quantities[i]]
    param = rescale(param, quantities[i])
    y_lim_max = max(adjust_density(param)$y) + 0.05
    plot(adjust_density(param), lwd=lwd_val, xlab=param_list[[quantities[i]]], ylab='Density', cex.lab=1.5, cex.main=1.5, xlim=c(0,axis_list[[quantities[i]]]), ylim=c(0, y_lim_max), col='red', main='')
    par(new=T)
    plot(adjust_density(data_list[[quantities[i]]]), lwd=lwd_val, xlab=param_list[[quantities[i]]], ylab='Density', cex.lab=1.5, cex.main=1.5, xlim=c(0,axis_list[[quantities[i]]]), ylim=c(0, y_lim_max), col='black', main='') 
    if (quantities[i] == 'vertical_wavelength')
    {
      position = 'topleft'
    }
    else
    {
      position = 'topright'
    }
    legend(position, c('model', 'raw data'), col=c('red', 'black'), seg.len=seg.len_val, lwd=lwd_val, y.intersp=0.4, cex = 1, text.width=0.5)
  }
}

# compare model 2D histogram to raw data 2D histogram
scatterplots = function(dataset)
{
  for (i in 1:(length(quantities)-1))
  {
    x = dataset[,quantities[i]]
    x = rescale(x, quantities[i])
    x2 = data_list[[quantities[i]]]
    for (j in (i+1):length(quantities))
    {
      y = dataset[,quantities[j]]
      y = rescale(y, quantities[j])
      y2 = data_list[[quantities[j]]]
      
      par(mfrow=c(1,2))
      if (quantities[i] == 'vertical_wavelength')
      {
        xmin = min(min(x), min(x2))
      }
      else
      {
        xmin = 0
      }
      par(mar=c(4.5,5,4,1))
      hist2d(x, y, xlab=param_list[[quantities[i]]], ylab=param_list[[quantities[j]]], xlim=c(xmin, axis_list[[quantities[i]]]), cex.lab=1.5, cex.main=1.5, ylim=c(0, axis_list[[quantities[j]]]), main=paste('Model: Spearman\'s p:', toString(round(cor(x,y,method='spearman'),2))))
      hist2d(x2, y2, xlab=param_list[[quantities[i]]], ylab=param_list[[quantities[j]]], xlim=c(xmin, axis_list[[quantities[i]]]), cex.lab=1.5, cex.main=1.5, ylim=c(0, axis_list[[quantities[j]]]), main=paste('Raw data: Spearman\'s p:', toString(round(cor(x2,y2,method='spearman'),2))))   
    }
  }
}

sample = sample_joint_distribution(dim(data)[2])
marginals(sample)
scatterplots(sample)

# write sample from joint probability distribution to csv file
write.table(sample, 'joint_distribution_sample.csv', sep=',', append=F, row.names=F)








