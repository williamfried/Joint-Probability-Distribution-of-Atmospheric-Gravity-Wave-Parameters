source('station_info.R')
source('explore_functions.R')
setwd(path_to_data)

# options for each argument of function:
# category: 'marginal', 'correlation'
# by: 'month', 'year', 'station', 'none'
explore = function(category, by)
{
  if (category == 'marginal')
  {
    if (by == 'none')
    {
      return(marginal_all())
    }
    else
    {
      return(marginal(by))
    }
    stop('invalid function arguments')
  }
  if (category == 'correlation')
  {
    if (by == 'none')
    {
      return(correlation())
    }
    else
    {
      return(correlation_difference(by))
    }
    stop('invalid function arguments')
  }
  stop('invalid function arguments')
}

explore('marginal', 'station')

