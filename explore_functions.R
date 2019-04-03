setwd(path_to_data)
source('station_info.R')
source('functions.R')

years = paste(1998:2008)
months = c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

quantities = c('vertical_wavelength','freq','U_wavelength','V_wavelength','U_amp','V_amp','T_amp')
quantities_full_names = c('vertical wavelength','frequency','zonal wavelength','meridional wavelength','zonal wind amplitude','meridional wind amplitude','temperature amplitude')

param_list = list('zonal wind amplitude (m/s)', 'meridional wind amplitude (m/s)', 'temperature wind amplitude (K)', expression('frequency (' ~ 10^{-4} ~ 's'^{-1} ~ ')'), expression('zonal wavelength (' ~ 10^{3} ~ 'km )'), expression('meridional wavelength (' ~ 10^{3} ~ 'km )'), 'vertical wavelength (km)')
names(param_list) = c('U_amp','V_amp','T_amp','freq','U_wavelength','V_wavelength','vertical_wavelength') 

# marginal distribution of 7 parameters by month, year or station
marginal = function(parameter)
{
  if (parameter == 'month')
  {
    param = months
    color_list = list('blue', 'blue', 'blue', 'blue', 'purple', 'red', 'red', 'red', 'red', 'green', 'blue', 'blue')   
    names(color_list) = param
    csv_name = paste('_month.csv', sep='')
  }
  if (parameter == 'year')
  {
    param = years
    color_list = list('blue', 'orange', 'cyan', 'salmon', 'red', 'purple', 'pink', 'green', 'black', 'brown', 'magenta')   
    names(color_list) = param
    csv_name = paste('_year.csv', sep='')
  }
  if (parameter == 'station')
  {
    color_list = c('blue', 'red', 'green', 'purple', 'orange', 'cyan', 'salmon', 'pink', 'black', 'brown', 'magenta', 'firebrick') 
    
    param = station_names
    color_list = color_list[1:length(station_names)]
    names(color_list) = param
    
    csv_name = paste('_station.csv', sep='')
  }
  
  for (i in 1:length(quantities))
  {
    quantity = quantities[i]
    parameter_data = read.csv(paste(quantity, csv_name, sep=''), header=F, col.names=rep(0,10000), row.names=param)
    max_dens = 0 
    max_val = 0
    for (j in 1:length(param))
    {
      para = param[j]
      param_vec = as.numeric(as.vector(parameter_data[para,]))
      param_vec = param_vec[!is.na(param_vec)]
      
      max_dens_param = max(adjust_density(param_vec)$y)
      quantile_val = quantile(param_vec, 0.99)[[1]]
      
      if (max_dens_param > max_dens)
      {
        max_dens = max_dens_param
      }
      
      if (quantile_val > max_val)
      {
        max_val = quantile_val
      }
    }
    
    max_dens = max_dens + 0.05
    
    for (j in 1:length(param))
    {
      para = param[j]
      param_vec = as.numeric(as.vector(parameter_data[para,]))
      param_vec = param_vec[!is.na(param_vec)]
      
      par(mar=c(4.5,4.5,4,1))
      plot(adjust_density(param_vec), col=color_list[[para]], xlim=c(0,max_val), ylim=c(0,max_dens), xlab=param_list[[quantity]], ylab='density', main=paste('Differences in marginal distribution by', parameter), cex.lab=1.5, cex.main=1.5)
      par(new=T)
    }
    
    if (quantity == 'vertical_wavelength')
    {
      position = 'topleft'
    }
    else
    {
      position = 'topright'
    }
    
    if (parameter == 'month')
    {
      legend(position, c('Nov-Apr', 'May', 'Jun-Sep', 'Oct'), col=c('blue', 'purple', 'red', 'green'), lty='solid', seg.len=0.75, y.intersp=0.4, cex = 1, text.width=max_val/15)
    }
    if (parameter == 'station')
    {
      legend(position, station_names, col=unlist(color_list),  lty='solid', seg.len=0.75, y.intersp=0.4, cex = 1, text.width=max_val/8)
    }
    if (parameter == 'year')
    {
      legend(position, years, col=unlist(color_list),  lty='solid', seg.len=0.75, y.intersp=0.4, cex = 1, text.width=max_val/20)
    }
    par(new=F)
  }
}

# marginal distribution of 7 parameters across core 4, extra 7 or all 11 stations
marginal_all = function()
{
  all_data = read.csv('all.csv', header=F, row.names=quantities)
  for (i in 1:length(quantities))
  {
    quantity = quantities[i]
    param_vec = as.numeric(as.vector(all_data[quantity,]))
    param_vec = param_vec[!is.na(param_vec)]
    plot(adjust_density(param_vec), col='black', xlim=c(quantile(param_vec, 0)[[1]],quantile(param_vec, 0.99)[[1]]), ylim=c(0,(max(adjust_density(param_vec)$y)+0.05)), xlab=param_list[[quantity]], ylab='density', main=paste('Kernel density estimate of ', quantities_full_names[i], ' across all soundings', sep=''), cex.lab=1.4, cex.main=1.5)    
  }
}

# scatterplots for all pairs of 7 parameters across core 4, extra 7 or all 11 stations
correlation = function()
{
  data_all = read.csv('all.csv', header=F, row.names=quantities)
  for (i in 1:(length(quantities)-1))
  {
    x = as.numeric(as.vector(data_all[quantities[i],]))
    for (j in (i+1):length(quantities))
    {
      y = as.numeric(as.vector(data_all[quantities[j],]))
      par(mar=c(5,6,4,1)+.1)
      plot(x, y, xlab=param_list[[quantities[i]]], ylab=param_list[[quantities[j]]], xlim=c(quantile(x, 0.001)[[1]], quantile(x, 0.999)[[1]]), ylim=c(quantile(y, 0.001)[[1]], quantile(y, 0.999)[[1]]), main=paste('Spearman\'s rank correlation coefficient: ', toString(round(cor(x,y,method='spearman'),2))), cex.lab=1.5, cex.main=1.5)
    }
  }
}

# scatterplots for all pairs of 7 parameters by month, year or station
correlation_difference = function(parameter)
{
  winter_months = c(seq(1,4), seq(11,12))
  summer_months = seq(6,9)
  between_months = c(5,10)
  
  if (parameter == 'month')
  {
    main_font = 1.5
    color_list = c('blue', 'blue', 'blue', 'blue', 'purple', 'red', 'red', 'red', 'red', 'green', 'blue', 'blue')  
  }
  else
  {
    main_font = 1.4
    color_list = c('blue', 'red', 'green', 'purple', 'orange', 'cyan', 'salmon', 'pink', 'black', 'brown', 'magenta', 'firebrick') 
  }
  
  for (i in 1:(length(quantities)-1))
  {
    x_data = read.csv(paste(quantities[i], '_', parameter, '.csv', sep=''), header=F, col.names=rep(0,4000))
    for (j in (i+1):length(quantities))
    {
      y_data = read.csv(paste(quantities[j], '_', parameter, '.csv', sep=''), header=F, col.names=rep(0,4000))
      
      x_min_overall = Inf
      x_max_overall = -Inf
      y_min_overall = Inf
      y_max_overall = -Inf
      for (k in 1:dim(x_data)[1])
      {
        x = as.numeric(as.vector(x_data[k,]))
        x = x[!is.na(x)]
        y = as.numeric(as.vector(y_data[k,]))
        y = y[!is.na(y)]
        
        x_min = quantile(x, 0.005)[[1]] 
        x_max = quantile(x, 0.995)[[1]] 
        y_min = quantile(y, 0.005)[[1]] 
        y_max = quantile(y, 0.995)[[1]] 
        
        if (x_min < x_min_overall)
        {
          x_min_overall = x_min
        }
        if (x_max > x_max_overall)
        {
          x_max_overall = x_max
        }
        if (y_min < y_min_overall)
        {
          y_min_overall = y_min
        }
        if (y_max > y_max_overall)
        {
          y_max_overall = y_max
        }
      }
      
      if (parameter == 'month')
      {
        corr_vec_winter = c()
        corr_vec_summer = c()
        corr_vec_between = c()
      }
      else
      {
        corr_vec = c()
      }
      
      for (k in 1:dim(x_data)[1])
      {
        x = as.numeric(as.vector(x_data[k,]))
        x = x[!is.na(x)]
        y = as.numeric(as.vector(y_data[k,]))
        y = y[!is.na(y)]
        
        if (parameter == 'month')
        {
          if (k %in% winter_months)
          {
            corr_vec_winter = c(corr_vec_winter, round(cor(x,y, method='spearman'),2))
          }
          else if (k %in% summer_months)
          {
            corr_vec_summer = c(corr_vec_summer, round(cor(x,y, method='spearman'),2))
          }
          else
          {
            corr_vec_between = c(corr_vec_between, round(cor(x,y, method='spearman'),2))
          }
        }
        else
        {
          corr_vec[k] = round(cor(x,y, method='spearman'),2)
        }
        
        if (k != dim(x_data)[1])
        {
          title = ''
        }
        else
        {
          if (parameter == 'month')
          {
            title_winter = paste(corr_vec_winter, collapse=' ') 
            title_summer = paste(corr_vec_summer, collapse=' ') 
            title_between = paste(corr_vec_between, collapse=' ') 
          }
          else
          {
            title = paste('Spearman\'s rank correlation coefficient by ', parameter, ': ', (paste(corr_vec, collapse=' ')), sep='') 
          }
        }
        par(mar=c(4.5,5,4,1))
        plot(x, y, xlab=param_list[[quantities[i]]], ylab=param_list[[quantities[j]]], xlim=c(x_min_overall, x_max_overall), ylim=c(y_min_overall, y_max_overall), col=color_list[k], cex.lab=1.5)
        
        if (parameter != 'month' | k != dim(x_data)[1])
        {
          title(main=title, cex.main=main_font)
        }
        else
        {
          line_val = 0
          title(main=list(title_winter, col='blue'), adj=0, cex.main=1.5)
          title(main=list(title_summer, col='red'), adj=0.5, cex.main=1.5)
          title(main=list(title_between, col='purple'), adj=0.9, cex.main=1.5)
          title(main=list(paste('mean:', toString(round(mean(corr_vec_winter),2))), col='blue'), adj=0.1, line=line_val, cex.main=1.5)
          title(main=list(paste('mean:', toString(round(mean(corr_vec_summer),2))), col='red'), adj=0.5, line=line_val, cex.main=1.5)
          title(main=list(paste('mean:', toString(round(mean(corr_vec_between),2))), col='purple'), adj=0.89, line=line_val, cex.main=1.5)
        }
        par(new=T)
      }
      if (parameter == 'month')
      {
        legend('topright', c('Nov-Apr', 'May', 'Jun-Sep', 'Oct'), col=c('blue', 'purple', 'red', 'green'), pch=1, y.intersp=0.4, cex = 1, text.width=x_max_overall/20)
      }
      par(new=F)
    }
  }
}
