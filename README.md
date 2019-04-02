Perform the following steps to construct the joint probability distribution of the seven atmospheric gravity wave parameters.
The data retrieval and inference steps are carried out in Python, and the modeling step is carried out in R. 

Python:
1. Determine the geographical region for which you want to characterize the atmospheric gravity wave activity.
2. Use the following two resources to determine which upper air stations are located within the given geographical region:
https://www.weather.gov/upperair/nws_upper 
https://www.sparc-climate.org/data-centre/data-access/us-radiosonde/us-upper-air-station-details/
3. Insert the WBAN numbers corresponding to the chosen upper air stations in the _stations_ list in _station_info.py_.
4. Look up the latitude corresponding to each upper air station and insert these values into the _latitudes_ list in _station_info.py_. Make sure the latitudes line up with the corresponding stations.
5. Specify the path to the folder where you want to store the raw radiosonde data in the _path_to_radiosonde_data_ variable in _station_info.py_. 
6. Specify the path to the folder where you want to store the inferred gravity wave parameters in the _path_to_csv_files_ variable in _station_info.py_. 
7. Run _FTP.py_ for each upper air station. This retrieves and unzips the radiosonde data from 1998-2008 and stores each sounding as a DAT file in the _path_to_radiosonde_data_ folder. 
8. Run _clean.py_ for each upper air station. This iterates through all the DAT files and removes files that don't meet specifications (see _Clearning Data_ section in report for more details about the criteria that are used for removing DAT files). 
9. Run _execute_inference.py_ for each upper air station. This infers the gravity wave parameters for each radiosonde soundings and stores these quantities in a massive dictionary.
10. Run _exploration.py_ to explore how the inferred gravity wave parameters vary on a sounding-to-sounding basis as well as by time of day (see _Sounding-to-sounding Variation_ and _Variation by Time of Day_ sections in report for more details). 
11. Run _export_to_R.py_ to transfer the inferred gravity wave paramters into CSV files that can be imported into R. These CSV files are stored in the _path_to_csv_files_ folder. Four types of CSV files are created: one type contains all the gravity wave parameters from all the stations, while the other three types consist of the gravity wave parameters categorized by year, month and station.

R:
1. Define the _path_to_data_ variable in _station_info.R_ to be the same as _path_to_csv_files_. 
2. Make sure that all the R files are stored in the _path_to_data_ folder. 
3. Define the _station_names_ variable in _station_info.R_ to be the same as the _stations_ list. Make sure that the radiosonde stations are listed in the same order as they are in the _stations_ list!
4. Insert the name of the CSV file that contains the data that will be used to construct the joint probability distribution as the argument of the _import_ function in _modeling steps.R_. The name of this CSV file should be one of: _all.csv_, _summer.csv_ and _winter.csv_. 
5. Run the _explore_ function in _explore.R_ to confirm that the assumptions described in the report in the _Overall Structure_, _Variation by Radiosonde Station_, _Variation by Year_ and _Variation by Month_ sections are accurate. 
6. Run _modeling_steps.R_. This will perform the following steps:
  * Import and organize data from CSV file.
  * Model marginal distributions of seven gravity wave parameters.
  * Construct copula and optimize correlation coefficients using coordinate descent algorithm. (Note: this step can take over an hour to run depending on the amount of data used to build the model.)
  * Perform cross-validation to determine the optimal number of intervals for determining the conditonal frequency distribution.
  * Calculate the AIC corresponding to each of the following three models of increasing complexity:
    1. Models all 7 gravity wave parameters independently of each other.
    2. Models the joint distribution of all gravity wave parameters except for frequency using a copula and independently models frequency. 
    3. Models the joint distribution of all gravity wave parameters except for frequency using a copula and models the frequency by conditioning on the zonal and meridional wavelengths. 
  * Draw sample of given size from the joint probability distribution and compare these samples to the empirical marginal distributions and correlation stuctures of the inferred gravity wave paramters. This sample is then stored in a CSV file where it can be further analyzed.

