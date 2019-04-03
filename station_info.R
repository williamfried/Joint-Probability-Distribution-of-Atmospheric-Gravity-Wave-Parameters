# define path to folder where CSV files of the inferred gravity wave parameters are stored
path_to_data = '/Users/willfried/Desktop/test/Joint-Probability-Distribution-of-Atmospheric-Gravity-Wave-Parameters/'

# define vector of station names (needs to be in same order as stations are define in station_info.py)
station_names = c('03020', '23050', '23047', '23023')

# specify the CSV file that contains data that will be used to construct the joint probability distribution
# For Nov-Apr, use 'winter.csv'; for Jun-Sept, use 'summer.csv'; for May and Oct, use 'all.csv'
dataset = 'summer.csv'