from ftplib import FTP
import gzip
import shutil
import os 
import sys
from station_info import *

station_code = sys.argv[1]

years = [str(x + 1998) for x in range(11)]

# connect to FTP server
ftp=FTP('sparc-ftp1.ceda.ac.uk')
ftp.login()

# download and unzip radiosonde readings 
path = path_to_radiosonde_data + station_code
if not os.path.exists(path):
	os.mkdir(path)
os.chdir(path)

for year in years:
	ftp.cwd('/sparc/hres/6_second/' + year + '/' + station_code + '/')
	filenames = ftp.nlst()

	for filename in filenames:
		print(filename)
		local_filename = os.path.join(path, filename)
		with open(local_filename ,'wb') as f:
			ftp.retrbinary('RETR ' + filename, f.write)

		with gzip.open(filename, 'rb') as f_in:
			with open(filename[0:-3], 'wb') as f_out:
				shutil.copyfileobj(f_in, f_out)
			
ftp.quit()

# for code in 94703 14684 14607 54762
# for code in 14929 94980 13995
# do
# python ftp.py $code &
# done 
		