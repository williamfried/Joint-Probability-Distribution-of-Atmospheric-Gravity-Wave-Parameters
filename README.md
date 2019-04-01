# Senior-Thesis

1. Decide on the geographical region where you want to learn about the atmospheric gravity wave activity.
2. Use the following two resources to determine which upper air stations are located within the given geographical region:
https://www.weather.gov/upperair/nws_upper 
https://www.sparc-climate.org/data-centre/data-access/us-radiosonde/us-upper-air-station-details/
3. Insert the WBAN number corresponding to the chosen upper air stations in the "stations" list in "station_info.py".
4. Look up the latitude corresponding to each upper air station and insert these values into the "latitudes" list in "station_info.py". Make sure the latitudes line up with the corresponding stations.
5. Specify the path to the folder where you want to store the raw radiosonde data in the "path_to_radiosonde_data" variable in "station_info.py". 
6. Specify the path to the folder where you want to store the inferred gravity wave parameters in the "path_to_csv_files" variable in "station_info.py". 
7. Run FTP.py for each upper air station. This will retrieve and unzip the radiosonde data from 1998-2008 and store it in the "path_to_radiosonde_data" folder. 

FTP.py to downl
