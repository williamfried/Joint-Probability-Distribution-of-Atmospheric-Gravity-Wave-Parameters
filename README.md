# Senior-Thesis

Python:
1. Decide on the geographical region where you want to gain information about the atmospheric gravity wave activity.
2. Use the following two resources to determine which upper air stations are located within the given geographical region:
https://www.weather.gov/upperair/nws_upper 
https://www.sparc-climate.org/data-centre/data-access/us-radiosonde/us-upper-air-station-details/
3. Insert the WBAN number corresponding to the chosen upper air stations in the "stations" list in "station_info.py".
4. Look up the latitude corresponding to each upper air station and insert these values into the "latitudes" list in "station_info.py". Make sure the latitudes line up with the corresponding stations.
5. Specify the path to the folder where you want to store the raw radiosonde data in the "path_to_radiosonde_data" variable in "station_info.py". 
6. Specify the path to the folder where you want to store the inferred gravity wave parameters in the "path_to_csv_files" variable in "station_info.py". 
7. Run "FTP.py" for each upper air station. This retrieves and unzip the radiosonde data from 1998-2008 and stores each sounding as a DAT file in the "path_to_radiosonde_data" folder. 
8. Run "clean.py" for each upper air station. This iterates through all the DAT files and removes files that don't meet specifications (see "Clearning Data" section in report for more details on the criteria that are used for removing DAT files). 
9. Run "execute_inference.py" for each upper air station. This infers the gravity wave parameters for each radiosonde soundings and stores these quantities in a massive dictionary.
10. Run "exploration.py" to explore how the inferred gravity wave parameters vary on a sounding-to-sounding basis as well as by time of day (see "Sounding-to-sounding Variation" and "Variation by Time of Day" sections in report for more details). 
11. Run "export_to_R.py" to transfer the inferred gravity wave paramters into CSV files that can be imported into R. These CSV files are stored in the "path_to_csv_files" folder. Four types of CSV files are created: one type contains all the gravity wave parameters from all the stations, while the other three types consist of the gravity wave parameters separated by year, month and station.

R:
1. Define the "path_to_data" variable in "explore.R" to be the same as "path_to_csv_files". 
2. Define the "station_names" variable in "explore_functions.R" to be the same as the "stations" list. Make sure that the radiosonde stations are listed in the same order as they are in the "stations" list!
3. Run the "explore" function in "explore.R" to confirm that the assumptions described in the report in the "Overall Structure", "Variation by Radiosonde Station", "Variation by Year" and "Variation by Month" sections are accurate. 
4. Insert the name of the CSV file that contains the data that will be used to construct the joint probability distribution as the argument of the "import" function in "modeling steps.R". The name of this CSV file should be one of: "all.csv", "summer.csv" and "winter.csv". 
5. Run "modeling_steps.R". This will perform the following steps:
  a. Import and organize data from CSV file.
  b. Model marginal distributions of seven gravity wave parameters.
  c. Create copula and optimize correlation coefficients using coordinate descent algorithm. (Note: this step can take over an hour to run depending on the amount of data in the CSV file.)
  d. Perform cross-validation to determine optimal number of intervals for determining the conditonal frequency distribution.
  e. Calculate the AIC corresponding to each of the following three models of increasing complexity:
    i.  Models all 7 gravity wave parameters independently of each other.
    ii. Models the joint distribution of all gravity wave parameters except for frequency using a copula and independently models frequency. 
    iii. Models the joint distribution of all gravity wave parameters except for frequency using a copula and models the frequency by conditioning on the zonal and meridional wavelengths. 
  f. Draw sample of given size from the joint probability distribution and compare these samples to the empirical marginal distributions and correlation stuctures of the inferred gravity wave paramters. This sample is then stored in a CSV where it can be further analyzed.

