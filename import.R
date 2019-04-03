setwd(path_to_data)
source('station_info.R')

import = function(dataset)
{
  # import raw data 
  quantities = c('vertical_wavelength','freq','U_wavelength','V_wavelength','U_amp','V_amp','T_amp') 
  data = read.csv(dataset, header=F, row.names=quantities)
  #data = read.csv('poly2.csv', header=F, row.names=quantities)
  
  # extract 7 inferred gravity wave parameters
  vert_wavelength <<- as.numeric(as.vector(data['vertical_wavelength',]))
  freq <<- as.numeric(as.vector(data['freq',]))
  
  # shift frequency values close to zero to facilitate modeling 
  freq_shifted <<- freq - (min(freq) - 0.025)
  
  # scale zonal and meridional wavelength values down by a factor of 1000 to facilitate modeling 
  U_wavelength <<- as.numeric(as.vector(data['U_wavelength',]))
  #U_wavelength <<- U_wavelength / 1000
  V_wavelength <<- as.numeric(as.vector(data['V_wavelength',]))
  #V_wavelength <<- V_wavelength / 1000
  
  U_amp <<- as.numeric(as.vector(data['U_amp',]))
  V_amp <<- as.numeric(as.vector(data['V_amp',]))
  T_amp <<- as.numeric(as.vector(data['T_amp',]))
  
  return(data)
}

