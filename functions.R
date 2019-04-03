setwd(path_to_data)

library('actuar')
library('VGAM')
library('rmutil')
library('invgamma')
library('extraDistr')
library('fitdistrplus')
library('MASS')
source('gnorm_library.R')
library('mixR')
library('rootSolve')

dburr = actuar::dburr
pburr = actuar::pburr
qburr = actuar::qburr
rburr = actuar::rburr
dinvgamma = invgamma::dinvgamma
pinvgamma = invgamma::pinvgamma
qinvgamma = invgamma::qinvgamma
rinvgamma = invgamma::rinvgamma

# order of polynomial fit
order = 100

param_list = list('zonal wind amplitude (m/s)', 'meridional wind amplitude (m/s)', 'temperature wind amplitude (K)', expression('frequency (' ~ 10^{-4} ~ 's'^{-1} ~ ')'), expression('zonal wavelength (' ~ 10^{3} ~ 'km )'), expression('meridional wavelength (' ~ 10^{3} ~ 'km )'), 'vertical wavelength (km)')
names(param_list) = c('U_amp','V_amp','T_amp','freq','U_wavelength','V_wavelength','vertical_wavelength')  
lwd_val = 2
seg.len_val = 1

# compute AIC given log likelihood and number of model parameters
AIC = function(log_likelihood, param_num)
{
  return(2*param_num - 2*log_likelihood)
}

# calculate the maximum likelihood estimate of a gamma distribution (https://tminka.github.io/papers/minka-gamma.pdf)
optimize_gamma = function(x)
{
  mean_logs = sum(log(x)) / length(x)
  u = mean(x)
  a = 0.5 / (log(u) - mean_logs)
  while(T)
  {
    expr = (1/a) + (mean_logs - log(u) + log(a) - digamma(a)) / ((a^2)*((1/a) - trigamma(a)))
    a_new = 1/expr
    if ((a_new - a) < 1e-7)
    {
      shape = a_new
      rate = a_new / u
      
      tot = 0
      for (i in 1:length(x))
      {
        tot = tot + log(dgamma(x[i], shape, rate=rate))
      }
      
      param_list = list(shape, rate, tot)
      names(param_list) = c('shape', 'rate', 'loglik')
      return(param_list)
    }
    a = a_new
  }
}

# calculate the maximum likelihood estimate of a weibull distribution (https://stats.stackexchange.com/questions/230109/comparing-approaches-of-mle-estimates-of-a-weibull-distribution)
optimize_weibull = function(x)
{
  n = length(x)
  weib1 = function(c) 
  { 
    return(1/c - sum(x^c*log(x))/sum(x^c) + (1/n)*sum(log(x)))
  }
  shape = uniroot(weib1, c(0,10), tol=1e-12)$root  
  scale = ((1/n)*sum(x^shape))^(1/shape)
  
  tot = 0
  for (i in 1:length(x))
  {
    tot = tot + log(dweibull(x[i], shape, scale=scale))
  }
  
  param_list = list(shape, scale, tot)
  names(param_list) = c('shape', 'scale', 'loglik')
  return(param_list)
}

# adjust kernel density estimate such that density is zero for all values less than 0
adjust_density = function(data)
{
  h = density(data, kernel='gaussian')$bw
  w = 1 / pnorm(0, mean=data, sd=h, lower.tail=F)
  d = suppressWarnings(density(data, bw=h, kernel="gaussian", weights=w/length(data)))
  d$y[d$x < 0] = 0
  d$y = d$y[d$x > -0.02]
  d$x = d$x[d$x > -0.02]
  
  return(d)
}

# evalulate inverse CDF given standard uniform random variable
solve_for_zero = function(random, CDF, upper_zero, lower_zero)
{
  upper_bound = upper_zero
  lower_bound = lower_zero
  x_val = mean(c(upper_zero, lower_zero))
  while (T)
  {
    y_val = CDF(x_val) - CDF(lower_zero)
    if (abs(y_val - random) < 1e-4)
    {
      return(x_val)
    }
    if (y_val > random)
    {
      upper_bound = x_val
      x_val = mean(c(x_val, lower_bound))
    }
    else
    {
      lower_bound = x_val
      x_val = mean(c(x_val, upper_bound))
    }
  }
}

# find the x value at which a given PDF drops below zero on the upper end of the distribution
find_upper_zero = function(data, PDF)
{
  i = round(0.85*length(data))
  while (T)
  {
    if (i == length(data))
    {
      return(data[i])
    }
    if (PDF(data[i]) > 0 && PDF(data[i+1]) < 0)
    {
      return(mean(c(data[i], data[i+1])))
    }
    i = i + 1
  }
}

# find the x value at which a given PDF drops below zero on the lower end of the distribution
find_lower_zero = function(data, PDF)
{
  i = round(0.1*length(data))
  while (T)
  {
    if (i == 1)
    {
      if (data[i] < 0.05)
      {
        return(0)
      }
      return(data[i])
    }
    if (PDF(data[i]) > 0 && PDF(data[i-1]) < 0)
    {
      return(mean(c(data[i], data[i-1])))
    }
    i = i - 1
  }
}

# construct CDF or PDF given regression coefficients, the order of the polynomial and the area under the PDF curve 
make_DF_string = function(lm_coeffs, order, type, area)
{
  string = '('
  for (i in 1:(order+1))
  {
    if (!is.na(lm_coeffs[i]))
    {
      if (lm_coeffs[i] > 0 && i != 1)
      {
        sign = ' + '
      }
      else if (i == 1)
      {
        sign = ''
      }
      else
      {
        sign = ' '
      }
      if (type == 'PDF')
      {
        string = paste(string, sign, toString(lm_coeffs[i]), '*x^', toString(i-1), sep='')
      }
      else if (type == 'CDF')
      {
        string = paste(string, sign, toString(lm_coeffs[i]), '*x^', toString(i), '/', toString(i), sep='')
      }
      else
      {
        stop('invalid type')
      }
    }
  }
  string = paste(string, ')/', toString(area), sep='')
  return(string)
}

# convert CDF or PDF to function of x
make_DF = function(string)
{
  func = function(num)
  {
    x = num
    return(eval(parse(text=string)))
  }
  return(func)
}

# identify distribution among all candidate distributions whose corresponding maximum likelihood estimate has the lowest AIC 
fit_dist = function(quantity, plots, all, best, dists, quantity_name)
{
  # fraction of bandwidth to use for kernel density estimation
  adj = 1
  
  main_font = 2
  axis_font = 1.5
  y_offset = 0
  quantile_val = 0.995
  increase = 0.05 * quantile(quantity, quantile_val)
  
  x = seq(0.01, quantile(quantity, quantile_val)+increase, length=1000)
  quantity_limited = quantity[quantity < quantile(quantity, quantile_val)]
  
  AICs = list()
  
  tw = 0.8
  max_max = 0
  
  # fit normal distribution
  if ('norm' %in% dists)
  {
    norm = fitdistr(quantity, 'normal')
    norm_likelihood = unname(norm['loglik'])
    avg = unname(norm['estimate'][[1]]['mean'])
    std_dev = unname(norm['estimate'][[1]]['sd'])
    
    AICs['norm'] = AIC(norm_likelihood[[1]],2)
    
    # plot normal distribution
    if (plots)
    {
      y = dnorm(x, mean=avg, sd=std_dev)
      plot(adjust_density(quantity_limited), lwd=lwd_val, ylim=c(0, max(y, adjust_density(quantity)$y) + y_offset), xlab = unname(param_list[quantity_name]), ylab='density', main=paste('MLE: normal -- corresponding AIC:', toString(round(AIC(norm_likelihood[[1]],2),0))), cex.lab=axis_font, cex.main=main_font)
      points(x, y, type="l", lwd=lwd_val, col='red')
      legend('topright', c('kernel density estimate', 'fitted distribution'), col=c('black', 'red'), seg.len=seg.len_val, lwd=lwd_val, y.intersp=0.4, cex = 1, text.width=1)
    }
  }
  
  # fit generalized normal distribution
  if ('generalized_normal' %in% dists)
  {
    if (is(try(fitdist(quantity, 'gnorm', start=list(mu=mean(quantity), alpha=1.5, beta=3.5), lower=c(0,0.1,0.1), method='mle'), T), 'try-error'))
    {
      AICs['generalized_normal'] = Inf
    }
    else
    {
      fitgennorm = fitdist(quantity, 'gnorm', start=list(mu=mean(quantity), alpha=1.5, beta=3.5), lower=c(0,0.1,0.1), method='mle')
      AICs['generalized_normal'] = fitgennorm$aic
    
      # plot generalized normal distribution
      if (plots)
      {
        mu_est = fitgennorm$estimate[['mu']]
        alpha_est = fitgennorm$estimate[['alpha']]
        beta_est = fitgennorm$estimate[['beta']]
        
        y = dgnorm(x, mu_est, alpha_est, beta_est)
        plot(adjust_density(quantity_limited), lwd=lwd_val, ylim=c(0, max(y, adjust_density(quantity)$y, max_max) + y_offset), xlab = unname(param_list[quantity_name]), ylab='density', main=paste('MLE: generalized normal -- corresponding AIC:', toString(round(fitgennorm$aic,0))), cex.lab=axis_font, cex.main=main_font)
        points(x, y, type="l", lwd=lwd_val, col='red')
        legend('topleft', c('kernel density estimate', 'fitted distribution'), col=c('black', 'red'), seg.len=seg.len_val, lwd=lwd_val, y.intersp=0.4, cex = 1, text.width=0.6)
      }
    }
  }
  
  # fit lognormal distribution
  if ('lnorm' %in% dists)
  {
    lognorm = fitdistr(quantity, 'lognormal')
    lognorm_likelihood = unname(lognorm['loglik'])
    avglog = unname(lognorm['estimate'][[1]]['meanlog'])
    std_dev_log = unname(lognorm['estimate'][[1]]['sdlog'])
    
    AICs['lnorm'] = AIC(lognorm_likelihood[[1]],2)
    
    # plot lognormal distribution
    if (plots)
    {
      y = dlnorm(x, meanlog=avglog, sdlog=std_dev_log)
      plot(adjust_density(quantity_limited), lwd=lwd_val, ylim=c(0, max(y, adjust_density(quantity)$y, max_max) + y_offset), xlab = unname(param_list[quantity_name]), ylab='density', cex.lab=axis_font, cex.main=main_font, main=paste('MLE: lognormal -- corresponding AIC:', toString(round(AIC(lognorm_likelihood[[1]],2),0))))
      points(x, y, type="l", lwd=lwd_val, col='red')
      legend('topright', c('kernel density estimate', 'fitted distribution'), col=c('black', 'red'), seg.len=seg.len_val, lwd=lwd_val, y.intersp=0.4, cex = 1, text.width=tw)
    }
  }
  
  # fit gamma distribution
  if ('gamma' %in% dists)
  {
    if (is(try(fitdistr(quantity, 'gamma', lower=c(0.01,0.01)), T), 'try-error'))
    {
      gamma = optimize_gamma(quantity)
      AICs['gamma'] = AIC(gamma[['loglik']],2)
      gamma_shape = gamma[['shape']]
      gamma_rate = gamma[['rate']]
    }
    else
    {
      gamma = fitdistr(quantity, 'gamma', lower=c(0.01,0.01)) 
      gamma_likelihood = unname(gamma['loglik'])
      gamma_shape = unname(gamma['estimate'][[1]]['shape'])
      gamma_rate = unname(gamma['estimate'][[1]]['rate'])
      
      AICs['gamma'] = AIC(gamma_likelihood[[1]],2)
    }
      
      # plot gamma distribution
      if (plots)
      {
        y = dgamma(x, gamma_shape, rate=gamma_rate)
        plot(adjust_density(quantity_limited), lwd=lwd_val, ylim=c(0, max(y, adjust_density(quantity)$y, max_max) + y_offset), xlab = unname(param_list[quantity_name]), ylab='density', cex.lab=axis_font, cex.main=main_font, main=paste('MLE: gamma -- corresponding AIC:', toString(round(AIC(gamma_likelihood[[1]],2),0))))
        points(x, y, type="l", lwd=lwd_val, col='red')
        legend('topright', c('kernel density estimate', 'fitted distribution'), col=c('black', 'red'), seg.len=seg.len_val, lwd=lwd_val, y.intersp=0.4, cex = 1, text.width=tw)
      }
  }
  
  # fit Weibull distribution
  if ('weibull' %in% dists)
  {
    if (is(try(fitdistr(quantity, 'weibull', lower=c(0.01,0.01)), T), 'try-error'))
    {
      weibull = optimize_weibull(quantity)
      AICs['weibull'] = AIC(weibull[['loglik']],2)
      weibull_shape = weibull[['shape']]
      weibull_scale = weibull[['scale']]
    }
    else
    {
      weibull = fitdistr(quantity, 'weibull', lower=c(0.01,0.01))
      weibull_likelihood = unname(weibull['loglik'])
      weibull_shape = unname(weibull['estimate'][[1]]['shape'])
      weibull_scale = unname(weibull['estimate'][[1]]['scale'])
      
      AICs['weibull'] = AIC(weibull_likelihood[[1]],2)
    }
    
      # plot Weibull distribution
      if (plots)
      {
        y = dweibull(x, weibull_shape, scale=weibull_scale)
        plot(adjust_density(quantity_limited), lwd=lwd_val, ylim=c(0, max(y, adjust_density(quantity)$y, max_max) + y_offset), xlab = unname(param_list[quantity_name]), ylab='density', cex.lab=axis_font, cex.main=main_font, main=paste('MLE: Weibull distribution -- corresponding AIC:', toString(round(AIC(weibull_likelihood[[1]],2),0))))
        points(x, y, type="l", lwd=lwd_val, col='red')
        legend('topright', c('kernel density estimate', 'fitted distribution'), col=c('black', 'red'), seg.len=seg.len_val, lwd=lwd_val, y.intersp=0.4, cex = 1, text.width=tw)
      }
  }
  
  # fit Pareto distribution
  if ('pareto' %in% dists)
  {
    
    pareto_ml = pareto_MLE(quantity)
    pareto_loc = pareto_ml[1]
    pareto_shape = pareto_ml[2]
    
    pareto_likelihood = sum(log(dpareto(quantity, pareto_loc, shape=pareto_shape)))
    
    AICs['pareto'] = AIC(pareto_likelihood[[1]],2)
    
    # plot Pareto distribution
    if (plots)
    {
      y = dpareto(x, pareto_loc, shape=pareto_shape)
      plot(adjust_density(quantity_limited), lwd=lwd_val, ylim=c(0, max(y, adjust_density(quantity)$y) + y_offset), xlab = unname(param_list[quantity_name]), ylab='density', main=paste('MLE: Pareto -- corresponding AIC:', toString(round(AIC(pareto_likelihood[[1]],2),0))), cex.lab=axis_font, cex.main=main_font)
      points(x, y, type="l", lwd=lwd_val, col='red')
      legend('topright', c('kernel density estimate', 'fitted distribution'), col=c('black', 'red'), seg.len=seg.len_val, lwd=lwd_val, y.intersp=0.4, cex = 1, text.width=tw)
    }
    
  }
  
  # fit Rayleigh distribution
  if ('rayleigh' %in% dists)
  {
    rayleigh_ml = rayleigh_MLE(quantity)$scale
    
    rayleigh_likelihood = sum(log(drayleigh(quantity, rayleigh_ml)))
    
    AICs['rayleigh'] = AIC(rayleigh_likelihood[[1]],1)
    
    # plot Rayleigh distribution
    if (plots)
    {
      y = drayleigh(x, rayleigh_ml)
      plot(adjust_density(quantity_limited), lwd=lwd_val, ylim=c(0, max(y, adjust_density(quantity)$y, max_max) + y_offset), xlab = unname(param_list[quantity_name]), ylab='density', cex.lab=axis_font, cex.main=main_font, main=paste('MLE: Rayleigh -- corresponding AIC:', toString(round(AIC(rayleigh_likelihood[[1]],1),0))))
      points(x, y, type="l", lwd=lwd_val, col='red')
      legend('topright', c('kernel density estimate', 'fitted distribution'), col=c('black', 'red'), seg.len=seg.len_val, lwd=lwd_val, y.intersp=0.4, cex = 1, text.width=tw)
    }
    
  }
  
  # fit Gumbel distribution
  if ('gumbel' %in% dists)
  {
    fitgumbel = tryCatch(fitdist(quantity, 'gumbel', start=list(a=2, b=2), method='mle'), error=function(cond){return(NA)})
    if (is.na(fitgumbel)[[1]])
    {
      AICs['gumbel'] = Inf
    }
    else
    {
      fitgumbel = fitdist(quantity, 'gumbel', start=list(a=2, b=2), method='mle')
      AICs['gumbel'] = fitgumbel$aic
      
      # plot Gumbel distribution
      if (plots)
      {
        a_est = fitgumbel$estimate[['a']]
        b_est = fitgumbel$estimate[['b']]
        
        y = dgumbel(x, a_est, b_est)
        plot(adjust_density(quantity_limited), lwd=lwd_val, ylim=c(0, max(y, adjust_density(quantity)$y) + y_offset), xlab = unname(param_list[quantity_name]), ylab='density', main=paste('MLE: Gumble -- corresponding AIC:', toString(round(fitgumbel$aic,0))), cex.lab=axis_font, cex.main=main_font)
        points(x, y, type="l", lwd=lwd_val, col='red')
        legend('topright', c('kernel density estimate', 'fitted distribution'), col=c('black', 'red'), seg.len=seg.len_val, lwd=lwd_val, y.intersp=0.4, cex = 1, text.width=tw)
      }
    }
  }
  
  # fit Burr distribution
  if ('burr' %in% dists)
  {
    fitburr = tryCatch(suppressWarnings(fitdist(quantity, 'burr', start = list(shape1 = 0.3, shape2 = 1, rate = 1), method='mle')), error=function(cond){return(NA)})
    if (is.na(fitburr)[[1]])
    {
      AICs['burr'] = Inf
    }
    else
    {
      AICs['burr'] = fitburr$aic
    
      # fit Burr distribution
      if (plots)
      {
        shape1_est = fitburr$estimate[['shape1']]
        shape2_est = fitburr$estimate[['shape2']]
        rate_est = fitburr$estimate[['rate']]
        
        y = dburr(x, shape1_est, shape2_est, rate=rate_est)
        plot(adjust_density(quantity_limited), lwd=lwd_val, ylim=c(0, max(y, adjust_density(quantity)$y, max_max) + y_offset), xlab=unname(param_list[quantity_name]), ylab='density', cex.lab=axis_font, cex.main=main_font, main=paste('MLE: Burr -- corresponding AIC:', toString(round(fitburr$aic,0))))
        points(x, y, type="l", lwd=lwd_val, col='red')
        legend('topright', c('kernel density estimate', 'fitted distribution'), col=c('black', 'red'), seg.len=seg.len_val, lwd=lwd_val, y.intersp=0.4, cex = 1, text.width=tw)
      }
    }
  }
  
  # fit Dagum distribution
  if ('dagum' %in% dists)
  {
    fitdagum = tryCatch(suppressWarnings(fitdist(quantity, 'dagum', start = list(scale = 2, shape1.a = 5, shape2.p = 0.5), method='mle')), error=function(cond){return(NA)})
    if (is.na(fitdagum)[[1]])
    {
      AICs['dagum'] = Inf
    }
    else
    {
      AICs['dagum'] = fitdagum$aic
  
      # plot Dagum distribution
      if (plots)
      {
        scale_est = fitdagum$estimate[['scale']]
        shape1.a_est = fitdagum$estimate[['shape1.a']]
        shape2.p_est = fitdagum$estimate[['shape2.p']]
        
        y = ddagum(x, scale_est, shape1.a_est, shape2.p_est)
        plot(adjust_density(quantity_limited), lwd=lwd_val, ylim=c(0, max(y, adjust_density(quantity)$y, max_max) + y_offset), xlab = unname(param_list[quantity_name]), ylab='density', cex.lab=axis_font, cex.main=main_font, main=paste('MLE: Dagum distribution -- corresponding AIC:', toString(round(fitdagum$aic,0))))
        points(x, y, type="l", lwd=lwd_val, col='red')
        legend('topright', c('kernel density estimate', 'fitted distribution'), col=c('black', 'red'), seg.len=seg.len_val, lwd=lwd_val, y.intersp=0.4, cex = 1, text.width=tw)
      }
    }
  }
  
  # fit inverse gamma distribution
  if ('invgamma' %in% dists)
  {
    fitinvgamma = tryCatch(fitdist(quantity, 'invgamma', start=list(shape=3, rate=5), method='mle'), error=function(cond){return(NA)})
    if (is.na(fitinvgamma)[[1]])
    {
      AICs['invgamma'] = Inf
    }
    else
    {
      AICs['invgamma'] = fitinvgamma$aic
      
      # plot inverse gamma distribution
      if (plots)
      {
        shape_est = fitinvgamma$estimate[['shape']]
        rate_est = fitinvgamma$estimate[['rate']]
        y = dinvgamma(x, shape_est, rate=rate_est)
        plot(adjust_density(quantity_limited), lwd=lwd_val, ylim=c(0, max(y, adjust_density(quantity)$y, max_max) + y_offset), xlab = unname(param_list[quantity_name]), ylab='density', cex.lab=axis_font, cex.main=main_font, main=paste('MLE: Inverse gamma -- corresponding AIC:', toString(round(fitinvgamma$aic,0))))
        points(x, y, type="l", lwd=lwd_val, col='red')
        legend('topright', c('kernel density estimate', 'fitted distribution'), col=c('black', 'red'), seg.len=seg.len_val, lwd=lwd_val, y.intersp=0.4, cex = 1, text.width=tw)
      }
    }
  }
  
  # fit generalized Rayleigh distribution
  if ('genray' %in% dists)
  {
    fitgenray = tryCatch(suppressWarnings(fitdist(quantity, 'genray', start=list(scale=1, shape=1), method='mle')), error=function(cond){return(NA)})
    if (is.na(fitgenray)[[1]])
    {
      AICs['genray'] = Inf
    }
    else
    {
      AICs['genray'] = fitgenray$aic
      
      # plot generalized Rayleigh distribution
      if (plots)
      {
        scale_est = fitgenray$estimate[['scale']]
        shape_est = fitgenray$estimate[['shape']]
        
        y = dgenray(x, scale=scale_est, shape_est)
        plot(adjust_density(quantity_limited), lwd=lwd_val, ylim=c(0, max(y, adjust_density(quantity)$y, max_max) + y_offset), xlab = unname(param_list[quantity_name]), ylab='density', cex.lab=axis_font, cex.main=main_font, main=paste('MLE: generalized Rayleigh -- corresponding AIC:', toString(round(fitgenray$aic,0))))
        points(x, y, type="l", lwd=lwd_val, col='red')
        legend('topright', c('kernel density estimate', 'fitted distribution'), col=c('black', 'red'), seg.len=seg.len_val, lwd=lwd_val, y.intersp=0.4, cex = 1, text.width=tw)
      }
    }
  }
  
  # fit gamma, Weibull and lognormal mixture model 
  if ('mixture model' %in% dists)
  {
    families = c('gamma', 'weibull', 'lnorm')
    for (i in 1:length(families))
    {
      previous_AIC = Inf
      components = 1
    
      # determine number of mixture model components that minimizes AIC
      while(T)
      {
        fit = eval(substitute(mixfit(.quantity, ncomp=.components, max_iter=250, family=families[i]), list(.quantity=quantity, .components=components)))
        if (is.na(fit$aic) || is.nan(fit$aic) || fit$aic > previous_AIC)
        {
          AICs[paste(families[i], 'mixture model')] = previous_AIC
          
          if (plots)
          {
            plot(adjust_density(quantity), lwd=lwd_val, ylim=c(0, 0.25), xlab = unname(param_list[quantity_name]), ylab='density', main=paste('MLE:', families[i], 'mixture model -- corresponding AIC:', toString(round(previous_AIC,0))), cex.lab=axis_font, cex.main=main_font)
            lines(eval(substitute(density(.previous_fit), list(.previous_fit = previous_fit))), lwd=lwd_val, col = 'red')
            legend('topright', c('kernel density estimate', 'fitted mixture model'), col=c('black', 'red'), seg.len=seg.len_val, lwd=lwd_val, y.intersp=0.4, cex = 1, text.width=1.4)
          }
          break
        }
        previous_fit = fit
        previous_AIC = fit$aic
        components = components + 1
      }
    }
  }
  
  # fit polynomial to kernel density estimate
  if ('polynomial' %in% dists)
  {
    d = adjust_density(quantity)
    
    x_data = d$x[d$x >= 0]
    y_data = d$y[d$x >= 0]
    
    # least squares polynomial fit to kernel density estimate
    density_lm = lm(y_data~poly(x_data, order, raw=T))
    
    # extract coefficients corresponding to polynomial fit 
    lm_coeffs = unname(density_lm$coefficients)
    
    # calculate density values based on polynomial fit
    density_fit = predict(density_lm)
    
    # construct initial PDF and CDF from polynomial fit
    PDF = make_DF(make_DF_string(lm_coeffs, order, 'PDF', 1))
    CDF = make_DF(make_DF_string(lm_coeffs, order, 'CDF', 1))
    
    # calculate area under PDF (should be very close to 1)
    upper_zero = find_upper_zero(x_data, PDF)
    lower_zero = find_lower_zero(x_data, PDF)
    area = CDF(upper_zero) - CDF(lower_zero)
    
    # adjust PDF such that area under curve exactly equals 1
    PDF = make_DF(make_DF_string(lm_coeffs, order, 'PDF', area))
    
    # compute AIC corresponding to polynomial fit
    tot = 0 
    for (i in 1:length(quantity))
    {
      tot = tot + log(PDF(quantity[i]))
    }
    AICs['polynomial'] = AIC(tot, 0)
    
    
    # plot PDF and CDF corresponding to polynomial fit
    if (plots)
    {
      plot(adjust_density(quantity), lwd=lwd_val, col='black', ylim=c(0,max(density_fit)+0.1), xlab=unname(param_list[quantity_name]), ylab='density', cex.lab=axis_font, cex.main=main_font, main='Polynomial fit')#, main=paste('Polynomial -- corresponding AIC: ', toString(round(AIC(tot, 2),0))), cex.lab=axis_font, cex.main=main_font)
      par(new=T)
      lines(x_data, density_fit, lwd=lwd_val, col='red')
      legend('topright', c('kernel density estimate', 'fitted polynomial'), col=c('black', 'red'), seg.len=seg.len_val, lwd=lwd_val, y.intersp=0.4, cex = 1, text.width=1.5)
      
      CDF = make_DF(make_DF_string(lm_coeffs, order, 'CDF', area))
      x_data_positive = x_data[x_data < upper_zero]
      predicted_CDF = CDF(x_data_positive)
      plot(x_data_positive, predicted_CDF, xlab=unname(param_list[quantity_name]), lwd=lwd_val, ylab='CDF', cex.lab=axis_font, type='l')
    }
  }
  
  # plot mixture model 
  if (grepl('mixture model', dists) && sapply(strsplit(dists, ' '), length) == 3)
  {
    previous_AIC = Inf
    components = 1
    while(T)
    {
      fit = eval(substitute(mixfit(.quantity, ncomp=.components, max_iter=250, family=word(dists, 1, 1)), list(.quantity=quantity, .components=components)))
      if (is.na(fit$aic) || is.nan(fit$aic) || fit$aic > previous_AIC)
      {
        AICs[dists[1]] = previous_AIC
        main=paste('MLE:', word(dists, 1, 1), ' mixture model -- corresponding AIC: ', toString(round(previous_AIC,0)))
        plot(adjust_density(quantity), xlab = unname(param_list[quantity_name]), ylab='density', main=paste('MLE:', word(dists, 1, 1), 'mixture model -- corresponding AIC:', toString(round(previous_AIC,0))), cex.lab=axis_font, cex.main=main_font)
        lines(eval(substitute(density(.previous_fit), list(.previous_fit = previous_fit))), col = 'red')
        break
      }
      previous_fit = fit
      previous_AIC = fit$aic
      components = components + 1
    }
  }
  
  # plot distribution that minimizes AIC
  if (best)
  {
    return(fit_dist(quantity, T, F, F, c(names(AICs[which.min(AICs)])), quantity_name))
  }
  if (all)
  {
    return(AICs) 
  }
  else
  {
    return(AICs[which.min(AICs)])
  }
}

# maximum likelihood estimate of gamma distribution
gamma_MLE = function(x)
{
  if (is(try(fitdistr(x, 'gamma', lower=c(0.01,0.01)), T), 'try-error'))
  {
    return(optimize_gamma(x)[1:2])
  }
  else
  {
    gamma = fitdistr(x, 'gamma', lower=c(0.01,0.01)) 
    gamma_shape = unname(gamma['estimate'][[1]]['shape'])
    gamma_rate = unname(gamma['estimate'][[1]]['rate'])
    gamma_list = list(gamma_shape, gamma_rate)
    names(gamma_list) = c('shape','rate')
    return(gamma_list)
  }
}

# maximum likelihood estimate of lognormal distribution
lognormal_MLE = function(x)
{
  lognorm = fitdistr(x, 'lognormal')
  avglog = unname(lognorm['estimate'][[1]]['meanlog'])
  std_dev_log = unname(lognorm['estimate'][[1]]['sdlog'])
  lognormal_list = list(avglog, std_dev_log)
  names(lognormal_list) = c('meanlog', 'sdlog')
  return(lognormal_list)
}

# maximum likelihood estimate of generalized normal distribution
generalized_normal_MLE = function(x)
{
  fitgennorm = fitdist(x, 'gnorm', start=list(mu=mean(x), alpha=1.5, beta=3.5), lower=c(0,0.1,0.1), method='mle')   
  mu_est = fitgennorm$estimate[['mu']]
  alpha_est = fitgennorm$estimate[['alpha']]
  beta_est = fitgennorm$estimate[['beta']]
  generalized_normal_list = list(mu_est, alpha_est, beta_est)
  names(generalized_normal_list) = c('mu', 'alpha', 'beta')
  return(generalized_normal_list)
}

# maximum likelihood estimate of Weibull distribution
weibull_MLE = function(x)
{
  if (is(try(fitdistr(x, 'weibull', lower=c(0.01,0.01)), T), 'try-error'))
  {
    return(optimize_weibull(x)[1:2])
  }
  else
  {
    weibull = fitdistr(x, 'weibull', lower=c(0.01,0.01))
    weibull_shape = unname(weibull['estimate'][[1]]['shape'])
    weibull_scale = unname(weibull['estimate'][[1]]['scale'])
    weibull_list = list(weibull_shape, weibull_scale)
    names(weibull_list) = c('shape', 'scale')
    return(weibull_list)
  }
}

# maximum likelihood estimate of Pareto distribution
pareto_MLE = function(x)
{
  n = length(x)
  m = min(x)
  a = n/sum(log(x)-log(m))
  pareto_list = list(m,a)
  names(pareto_list) = c('scale', 'shape')
  return(pareto_list) 
}

# maximum likelihood estimate of Rayleigh distribution
rayleigh_MLE = function(x)
{
  rayleigh_list = list(sqrt((sum(x^2 / 2)) / length(x)))
  names(rayleigh_list) = c('scale')
  return(rayleigh_list)
}

# maximum likelihood estimate of normal distribution
normal_MLE = function(x)
{
  u = mean(x)
  sigma = sd(x)
  normal_list = list(u, sigma)
  names(normal_list) = c('mean', 'sd')
  return(normal_list)
}

# maximum likelihood estimate of Gumbel distribution
gumbel_MLE = function(x)
{
  fitgumbel = fitdist(x, 'gumbel', start=list(a=1, b=1), method='mle')
  a_est = fitgumbel$estimate[['a']]
  b_est = fitgumbel$estimate[['b']]
  gumbel_list = list(a_est, b_est)
  names(gumbel_list) = c('a','b')
  return(gumbel_list)
}

# maximum likelihood estimate of Burr distribution
burr_MLE = function(x)
{
  fitburr = suppressWarnings(fitdist(x, 'burr', start = list(shape1 = 0.3, shape2 = 1, rate = 1), method='mle'))
  shape1_est = fitburr$estimate[['shape1']]
  shape2_est = fitburr$estimate[['shape2']]
  rate_est = fitburr$estimate[['rate']]
  burr_list = list(shape1_est, shape2_est, rate_est)
  names(burr_list) = c('shape1', 'shape2', 'rate')
  return(burr_list)
}

# maximum likelihood estimate of Dagum distribution
dagum_MLE = function(x)
{
  fitdagum = suppressWarnings(fitdist(x, 'dagum', start = list(scale = 2, shape1.a = 5, shape2.p = 0.5), method='mle'))
  scale_est = fitdagum$estimate[['scale']]
  shape1.a_est = fitdagum$estimate[['shape1.a']]
  shape2.p_est = fitdagum$estimate[['shape2.p']]
  dagum_list = list(scale_est, shape1.a_est, shape2.p_est)
  names(dagum_list) = c('scale', 'shape1.a', 'shape2.p')
  return(dagum_list)
}

# maximum likelihood estimate of inverse gamma distribution
invgamma_MLE = function(x)
{
  fitinvgamma = fitdist(x, 'invgamma', start=list(shape=3, rate=5), method='mle')
  shape_est = fitinvgamma$estimate[['shape']]
  rate_est = fitinvgamma$estimate[['rate']]
  invgamma_list = list(shape_est, rate_est)
  names(invgamma_list) = c('shape', 'rate')
  return(invgamma_list)
}

# maximum likelihood estimate of generalized Rayleigh distribution
genray_MLE = function(x)
{
  fitgenray = suppressWarnings(fitdist(quantity, 'genray', start=list(scale=1, shape=1), method='mle'))
  scale_est = fitgenray$estimate[['scale']]
  shape_est = fitgenray$estimate[['shape']]
  genray_list = list(scale_est, shape_est)
  names(genray_list) = c('scale', 'shape')
  return(genray_list)
}

# maximum likelihood estimate of mixture model from given family (gamma, Weibull or lognormal)
mixture_model_MLE = function(x, family)
{
  previous_AIC = Inf
  components = 1
  while(T)
  {
    fit = eval(substitute(mixfit(.x, ncomp=.components, family=family), list(.x=x, .components=components)))
    if (is.na(fit$aic) || is.nan(fit$aic) || fit$aic > previous_AIC)
    {
      mixture_model_list = list(previous_fit, family)
      names(mixture_model_list) = c('fit', 'family')
      return(mixture_model_list)
    }
    previous_fit = fit
    previous_AIC = fit$aic
    components = components + 1
  }
}

# polynomial fit 
polynomial_MLE = function(x)
{
  d = adjust_density(x)
  x_data = d$x[d$x >= 0]
  y_data = d$y[d$x >= 0]
  
  density_lm = lm(y_data~poly(x_data, order, raw=T))
  density_fit = predict(density_lm)
  lm_coeffs = unname(density_lm$coefficients)
  
  PDF = make_DF(make_DF_string(lm_coeffs, order, 'PDF', 1))
  CDF = make_DF(make_DF_string(lm_coeffs, order, 'CDF', 1))
  
  upper_zero = find_upper_zero(x_data, PDF)
  lower_zero = find_lower_zero(x_data, PDF)
  area = CDF(upper_zero) - CDF(lower_zero)
  
  PDF = make_DF(make_DF_string(lm_coeffs, order, 'PDF', area))
  CDF = make_DF(make_DF_string(lm_coeffs, order, 'CDF', area))

  upper_zero = find_upper_zero(x_data, PDF)
  lower_zero = find_lower_zero(x_data, PDF)
  
  polynomial_list = list(PDF, CDF, upper_zero, lower_zero)
  names(polynomial_list) = c('PDF', 'CDF', 'upper_zero', 'lower_zero')
  
  return(polynomial_list)
}

# calculate conditional density of frequency or generate conditional random values of frequency given zonal and meridional wavelengths
conditional_freq = function(f, cat, conditional_list, type)
{
  dist = conditional_list[[cat]]$dist
  
  if (type == 'density')
  {
    func_list = list(dgamma, dweibull, dlnorm, drayleigh, dnorm, dgumbel, dburr, ddagum, dinvgamma, dgenray, dgnorm)
    first_parameter = f
  }
  else if (type == 'random')
  {
    func_list = list(rgamma, rweibull, rlnorm, rrayleigh, rnorm, rgumbel, rburr, rdagum, rinvgamma, rgenray, rgnorm)
    first_parameter = 1
  }
  else
  {
    stop('invalid type')
  }
  
  names(func_list) = c('gamma', 'weibull', 'lnorm', 'rayleigh', 'norm', 'gumbel', 'burr', 'dagum', 'invgamma', 'genray', 'generalized_normal')
  
  if (dist == 'gamma')
  {
    return(func_list$gamma(first_parameter, conditional_list[[cat]]$parameters$shape, conditional_list[[cat]]$parameters$rate))
  }
  else if (dist == 'weibull')
  {
    return(func_list$weibull(first_parameter, conditional_list[[cat]]$parameters$shape, conditional_list[[cat]]$parameters$scale))
  }
  else if (dist == 'lnorm')
  {
    return(func_list$lnorm(first_parameter, conditional_list[[cat]]$parameters$meanlog, conditional_list[[cat]]$parameters$sdlog))
  }
  else if (dist == 'rayleigh')
  {
    return(func_list$rayleigh(first_parameter, conditional_list[[cat]]$parameters$scale))
  }
  else if (dist == 'norm')
  {
    return(func_list$norm(first_parameter, conditional_list[[cat]]$parameters$mean, conditional_list[[cat]]$parameters$sd))
  }
  else if (dist == 'gumbel')
  {
    return(func_list$gumbel(first_parameter, conditional_list[[cat]]$parameters$a, conditional_list[[cat]]$parameters$b))
  }
  else if (dist == 'burr')
  {
    return(func_list$burr(first_parameter, conditional_list[[cat]]$parameters$shape1, conditional_list[[cat]]$parameters$shape2, conditional_list[[cat]]$parameters$rate))
  }
  else if (dist == 'dagum')
  {
    return(func_list$dagum(first_parameter, conditional_list[[cat]]$parameters$scale, conditional_list[[cat]]$parameters$shape1.a, conditional_list[[cat]]$parameters$shape2.p))
  }
  else if (dist == 'invgamma')
  {
    return(func_list$invgamma(first_parameter, conditional_list[[cat]]$parameters$shape, conditional_list[[cat]]$parameters$rate))
  }
  else if (dist == 'genray')
  {
    return(func_list$genray(first_parameter, conditional_list[[cat]]$parameters$scale, conditional_list[[cat]]$parameters$shape))
  }
  else if (dist == 'generalized_normal')
  {
    return(func_list$generalized_normal(first_parameter, conditional_list[[cat]]$parameters$mu, conditional_list[[cat]]$parameters$alpha, conditional_list[[cat]]$parameters$beta))
  }
  else if (dist == 'mixture model')
  {
    if (type == 'density')
    {
      return(eval(substitute(density(.conditional_list[[.cat]]$parameters$fit, smoothness=1, from=.f, to=.f)$y, list(.conditional_list=conditional_list, .cat=cat, .f=f))))
    }
    else if (type == 'random')
    {
      family = conditional_list[[cat]]$parameters$family
      fit = conditional_list[[cat]]$parameters$fit
      
      if (family == 'gamma')
      {
        conversion = to_mu_sd_gamma(fit$alpha, fit$lambda)
        return(rmixgamma(1, fit$pi, conversion$mu, conversion$sd))
      }
      if (family == 'weibull')
      {
        conversion = to_mu_sd_weibull(fit$k, fit$lambda)
        return(rmixweibull(1, fit$pi, conversion$mu, conversion$sd))
      }
      if (family == 'lnorm')
      {
        conversion = to_mu_sd_lnorm(fit$mulog, fit$sdlog)
        return(rmixlnorm(1, fit$pi, conversion$mu, conversion$sd))
      }
      stop('invalid distribution family')
    }
  }
  else if (dist == 'polynomial')
  {
    if (type == 'density')
    {
      return(conditional_list[[cat]]$parameters$PDF(f))
    }
    else if (type == 'random')
    {
      return(solve_for_zero(runif(1), conditional_list[[cat]]$parameters$CDF, conditional_list[[cat]]$parameters$upper_zero, conditional_list[[cat]]$parameters$lower_zero))
    }
    stop('invalid argument')
  }
  stop('invalid distribution')
}

axis_limits = function(data_list)
{
  axis_list = list()
  for (i in 1:length(data_list))
  {
    axis_list[names(data_list)[i]] = quantile(data_list[[i]],0.995)[[1]]
  }
  return(axis_list)
}

